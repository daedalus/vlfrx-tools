//
//  Copyright (c) 2010 Paul Nicholson
//  All rights reserved.
//
//  Redistribution and use in source and binary forms, with or without
//  modification, are permitted provided that the following conditions
//  are met:
//  1. Redistributions of source code must retain the above copyright
//     notice, this list of conditions and the following disclaimer.
//  2. Redistributions in binary form must reproduce the above copyright
//     notice, this list of conditions and the following disclaimer in the
//     documentation and/or other materials provided with the distribution.
//
//  THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR
//  IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
//  OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
//  IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT,
//  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
//  NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
//  DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
//  THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
//  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
//  THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//

#include "config.h"
#include "vtport.h"
#include "vtlib.h"

#define EARTH_RAD 6371.0

// XXX Todo: implement ellipsoidal distance function
#define ERAD_a 6378.137     // WGS84 equatorial radius
#define ERAD_b 6356.7523142 // WGS84 polar radius

static int BFLAG = FALSE;        // -b option: calculate bearing
static int DFLAG = FALSE;        // -d option: calculate destination point
static int GFLAG = FALSE;        // -g option: output great circle
static int SFLAG = FALSE;
static int ncpu = 1;             // -U option: number of worker threads to use
static int OMODE_ISO = FALSE;    // Set by -o iso: ISO output format timestamps
static int OMODE_EXT = FALSE;    // Set by -o ext: extended output
static int OMODE_BRC = FALSE;    // Set by -o brc: output .brc files

static double filter_envelop = 0;  // -f envelop= option
static double filter_nearmax = 0;  // -f nearmax= option

#define MAXSITES 200

//
//  If a site or measurement doesn't provide a sigma for its TOGAs or bearings,
//  then these are the values used.
//

#define DEFAULT_AT_SIGMA 20e-6    // Arrival time sigma, seconds
#define DEFAULT_AZ_SIGMA 5        // Azimuth sigma, degrees

static double max_residual = 1; // Set by -r option:  units of sigma

//
//  Group velocity, km/sec.
//  This is a nominal value suitable for both day and night use at around 13kHz.
//
//  Actually, 0.9872 seems to work better for some reason.
//

static double CVLF = 300e3 * 0.9922;
static int cvlf_given = FALSE;

//
//  Some little utilities.
//

static void usage( void)
{
   fprintf( stderr,
       "usage:  vtspot [options] meas1 [meas2 ...] \n"
       "\n"
       "options:\n"
       "  -v            Increase verbosity\n"
       "  -L name       Specify logfile\n"
       "  -U ncpu       Number of worker threads to use (default 1)\n"
       "\n"
       "  -c factor     Velocity factor (default 0.9922)\n"
       "\n"
       "Measurements:\n"
       "\n"
       "    T/location/timestamp[/sigma]        (sigma in seconds)\n"
       "                                (defaults to 50e-6 sigma)\n"
       "    B/location/bearing[/sigma]          (sigma in degrees)\n"
       "                            (defaults to 10 degrees sigma)\n"
       "\n"
       "Matching:\n"
       " -m location=file    Specify TOGA file for the receiver location\n"
       " -n minsites         Minimum number of sites to attempt solution\n"
       "                     (default is all sites)\n"
       " -r residual         Max residual to accept solution (default 1.0)\n"
       " -f envelop=degrees  Apply surrounding filter rule\n"
       " -f nearmax=km       Apply distance limit to nearest receiver\n"
       " -o ext              Extended output records\n"
       " -o iso              Output ISO timestamps (default numeric)\n"
       "\n"
       "Utility functions:\n"
       " -b location1 location2 [timestamp]  Range and bearing of 2 from 1\n"
       " -d location bearing range_km        Destination from location\n"
       " -s location timestamp               Calculate Sun azimuth, elevation\n"
       " -g location1 location2 [step]       Great circle points\n"
     );
   exit( 1);
}

static double constrain( double a, double low, double high)
{
   double r = high - low;
   while (a < low) a += r;
   while (a >= high) a -= r;
   return a;
}

///////////////////////////////////////////////////////////////////////////////
// N-vector Operations                                                       //
///////////////////////////////////////////////////////////////////////////////

typedef double V3[3];
typedef double A3[3][3]; 

static V3 V3north = { 0, 0, 1 };
// static A3 A0 = { { 0, 0, 0}, {0, 0, 0}, {0, 0, 0} };
static A3 AI = { { 1, 0, 0}, {0, 1, 0}, {0, 0, 1} };

static inline int v3_equal( V3 v1, V3 v2)
{
   return v1[0] == v2[0] && v1[1] == v2[1] && v1[2] == v2[2];
}

static inline void v3_add( V3 v1, V3 v2, V3 vr)
{
   vr[0] = v1[0] + v2[0];
   vr[1] = v1[1] + v2[1];
   vr[2] = v1[2] + v2[2];
}

#if 0   // Not used right now
static void v3_sub( V3 v1, V3 v2, V3 vr)
{
   vr[0] = v1[0] - v2[0];
   vr[1] = v1[1] - v2[1];
   vr[2] = v1[2] - v2[2];
}
#endif

// Dot product
static inline double v3_dot( V3 v1, V3 v2)
{
   return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
}

// Cross product, v1 x v2 -> vc
static inline void v3_cross( V3 v1, V3 v2, V3 vc)
{
   vc[0] = v1[1] * v2[2] - v1[2] * v2[1];
   vc[1] = v1[2] * v2[0] - v1[0] * v2[2];
   vc[2] = v1[0] * v2[1] - v1[1] * v2[0];
}

// Magnitude
static inline double v3_mag( V3 v)
{
   return sqrt( v3_dot( v, v));
}

// Multiply by a scalar
static inline void v3_scale( V3 v, double s, V3 r)
{
   int i;
   for (i = 0; i < 3; i++) r[i] = v[i] * s;
}

// Make into unit vector
static inline void v3_normalise( V3 v)
{
   v3_scale( v, 1/v3_mag( v), v);
}

// Random normalised vector
static inline void v3_random( V3 v)
{
   v[0] = rand(); v[1] = rand(); v[2] = rand();
   v3_normalise( v);
}

// Cross product matrix from vector
static inline void a3_cpmat( V3 v, A3 a)
{
   a[0][0] = 0;     a[0][1] = -v[2];      a[0][2] = v[1];
   a[1][0] = v[2];  a[1][1] = 0;          a[1][2] = -v[0];
   a[2][0] = -v[1]; a[2][1] = v[0];       a[2][2] = 0;
}

static void a3_scalar_mul( A3 a, double s, A3 r)
{
   int i, j;
   for (i = 0; i < 3; i++) for (j = 0; j < 3; j++) r[i][j] = a[i][j] * s;
}

#if 0   // Not used right now
static void a3_copy( A3 s, A3 d)
{
   memcpy( d, s, sizeof(A3));
}
#endif

static void v3_copy( V3 s, V3 d)
{
   memcpy( d, s, sizeof( V3));
}

// Matrix sum, a + b -> r and r can be the same matrix as a and/or b
static void a3_sum( A3 a, A3 b, A3 r)
{
   int i, j;
   for (i = 0; i < 3; i++) for (j = 0; j < 3; j++) r[i][j] = a[i][j] + b[i][j];
}

// Matrix multiplication a x b -> r
static void a3_mul( A3 a, A3 b, A3 r)
{
   int i, j, k;

   for (i = 0; i < 3; i++)
      for (j = 0; j < 3; j++)
      {
         r[i][j] = 0;
         for (k = 0; k < 3; k++) r[i][j] += a[i][k] * b[k][j];
      }
}

#if 0   // Not used right now
static double a3_det( A3 X)
{
   return
    X[0][0]*X[1][1]*X[2][2] + X[0][1]*X[1][2]*X[2][0] +
    X[0][2]*X[1][0]*X[2][1] - X[0][2]*X[1][1]*X[2][0] -
    X[0][1]*X[1][0]*X[2][2] - X[0][0]*X[1][2]*X[2][1];
}
#endif

// Compute the unit vector normal to v1 and v2
static inline void v3_unit_normal_to( V3 v1, V3 v2, V3 vn)
{
   v3_cross( v1, v2, vn);
   v3_normalise( vn);
}

#if 0   // Not used right now
static void v3_outer_prod( V3 v1, V3 v2, A3 r)
{
   int i, j;
   for (i = 0; i < 3; i++) for (j = 0; j < 3; j++) r[i][j] = v1[i] * v2[j];
}
#endif

static void v3_transform( V3 v, A3 rot, V3 r)
{
   int i, j;

   for (i = 0; i < 3; i++)
   {
      r[i] = 0;
      for (j = 0; j < 3; j++) r[i] += rot[i][j] * v[j];
   }
}

// Compute rotation matrix R from unit axis vector v and angle phi
static void a3_rot( V3 v, double phi, A3 r)
{
   A3 k;  a3_cpmat( v, k);
   A3 t;  a3_scalar_mul( k, sin(phi), t);  a3_sum( AI, t, r);
   a3_mul( k, k, t);  a3_scalar_mul( t, 1 - cos(phi), t); a3_sum( r, t, r);
}

///////////////////////////////////////////////////////////////////////////////
// N-vector Calculations                                                     //
///////////////////////////////////////////////////////////////////////////////

//
//  Convert N-vector to latitude/longitude and format into a string.
//
//  If 's' is supplied, it forms the destination and return value.
//  With NULL 's', a pointer to a static string is returned, so should be used
//  no more than once in a printf argument list.
//

static char * v3_string( V3 v, char *s)
{
   static char temp[50];

   if (!s) s = temp;

   double lat = atan2( v[2], sqrt(v[0] * v[0] + v[1] * v[1]));
   double lon = atan2( v[1], v[0]);

   if (lat > M_PI/2)
   {
      lat = M_PI - lat;
      lon += M_PI;
   }
   else
   if (lat < -M_PI/2)
   {
      lat = -M_PI - lat;
      lon += M_PI;
   }

   lon = constrain( lon, -M_PI, M_PI);
   sprintf( s, "%.3f,%.3f", lat*180/M_PI, lon*180/M_PI);
   return s;
}

//
//  Convert latitude/longitude to N-vector.
//

static void v3_make( double lat, double lon, V3 v)
{
   v[0] = cos( lat) * cos( lon);
   v[1] = cos( lat) * sin( lon);
   v[2] = sin( lat);
}

//
//  Range in km between two points v1 and v2.
//

static inline double v3_range( V3 v1, V3 v2)
{
   V3 xp; v3_cross( v1, v2, xp);
   return atan2( sqrt( v3_dot( xp, xp)), v3_dot( v1, v2)) * EARTH_RAD;
}

//
//  Bearing of forepoint v2 from standpoint v1.
//

static double v3_bearing( V3 v1, V3 v2)
{
   V3 x1;  v3_unit_normal_to( V3north, v1, x1);
   V3 x2;  v3_unit_normal_to( v2, v1, x2);

   V3 xp;  v3_cross( x1, x2, xp);
   V3 xpn; v3_copy( xp, xpn); v3_normalise( xpn);
   return -v3_dot( v1, xpn) * atan2( v3_mag( xp), v3_dot( x1, x2));
}

//
//  Compute forepoint vf from standpoint vs along initial bearing b radians.
//  Distance 'a' is in units of Earth radii.
//

static void destination_point( V3 vs, double b, double a, V3 vf)
{
   A3 r0;  a3_rot( vs, -b, r0);
   V3 v0;  v3_transform( V3north, r0, v0);
   V3 n0;  v3_unit_normal_to( vs, v0, n0);
   A3 r1;  a3_rot( n0, a, r1);
   v3_transform( vs, r1, vf);
}

//
//  Parse a latitude or longitude string and return the angle in radians.
//

static double parse_coord( char *str)
{
   double a;
   char temp[50], *p;

   strcpy( temp, str);
   int sign = 1;

   if ((p = strchr( temp, 'S')) ||
       (p = strchr( temp, 'W')) ||
       (p = strchr( temp, 's')) ||
       (p = strchr( temp, 'w')))
   {
      *p = 0;  sign = -sign;
   }

   if ((p = strchr( temp, 'N')) ||
       (p = strchr( temp, 'E')) ||
       (p = strchr( temp, 'n')) ||
       (p = strchr( temp, 'e')))
   {
      *p = 0;
   }

   int fd, fm;
   double fs;

   if (sscanf( temp, "%d:%d:%lf", &fd, &fm, &fs) == 3 ||
       sscanf( temp, "%d:%d", &fd, &fm) == 2)
      a = fd + fm/60.0 + fs/3600.0;
   else
   if (sscanf( temp, "%lf", &a) != 1)
      VT_bailout( "bad coordinate [%s]", temp);

   return a * sign * M_PI/180;
}

///////////////////////////////////////////////////////////////////////////////
//  Sun Position & Path Illumination                                         //
///////////////////////////////////////////////////////////////////////////////

static void sunpos( V3 v, timestamp T, double *solaz, double *solel)
{
   //
   //  Time.
   //

   int ut = timestamp_secs( T);
 
   double JD = ut/86400.0 + 2440587.5;

   VT_report( 3, "JD %.7f", JD);

   double TU = JD - 2451545.0;
   double ERA = 2 * M_PI * (0.779057273 + 1.002737812 * TU);

   double GST = ERA; // - EPREC(t)   // Radians

   GST = GST - 2 * M_PI * (int)(GST/(2*M_PI));

   //
   //  RA and DEC of the Sun.
   //

   double mean_longitude = 280.46 + 0.9856474 * TU;   // Degrees

   double mean_anomaly = 357.528 + 0.9856003 * TU;   // Degrees

   double ecliptic_longitude =
     mean_longitude +
     1.915 * sin( mean_anomaly * M_PI/180) +
     0.020 * sin( 2 * mean_anomaly * M_PI/180);     // Degrees

   double e_of_e = 23.43691 * M_PI/180;   // Obliquity of the ecliptic, radians

   // RA and DEC in radians

   double RA = atan2( cos( e_of_e) * sin( ecliptic_longitude * M_PI/180),
                      cos( ecliptic_longitude * M_PI/180));

   double DEC = asin( sin( e_of_e) * sin( ecliptic_longitude * M_PI/180));

   VT_report( 3, "RA %.3f deg %.5f hrs", RA * 180/M_PI, RA * 12/M_PI);
   VT_report( 3, "DEC %.3f deg", DEC * 180/M_PI);
 
   //
   //  Local coordinates.
   //

   double lat = atan2( v[2], sqrt(v[0] * v[0] + v[1] * v[1]));
   double lon = atan2( v[1], v[0]);

   if (lat > M_PI/2)
   {
      lat = M_PI - lat;
      lon += M_PI;
   }
   else
   if (lat < -M_PI/2)
   {
      lat = -M_PI - lat;
      lon += M_PI;
   }

   lon = constrain( lon, -M_PI, M_PI);
  
   //
   //  Local time and hour angle.
   //
 
   double d = (ut - 946728000)/86400.0;
   double h = (ut % 86400)/3600.0;
   VT_report( 3, "d %.4f h %.4f", d, h);

   double LST = 100.46 * M_PI/180 + 0.985647 * d * M_PI/180 + lon + h * M_PI/12;
   VT_report( 3, "LST %.3f deg %.3f hours", LST * 180/M_PI, LST * 12/M_PI);
   double LHA = LST - RA;   // Radians

   // http://www.stargazing.net/kepler/altaz.html

   double alt = asin( sin(DEC)*sin(lat) + cos(DEC)*cos(lat)*cos(LHA));
   double az = acos( (sin(DEC) - sin(alt)*sin(lat))/cos(alt)/cos(lat));
   if (sin(LHA) > 0) az = 2 * M_PI - az;

   if (solaz) *solaz = az;
   if (solel) *solel = alt;
}

static double path_illumination( V3 v1, V3 v2, timestamp T)
{
   int n, k = 0;

   for (n = 0; n <= 19; n++)
   {
      double s = n/19.0;

      V3 t1; v3_scale( v1, 1 - s, t1);
      V3 t2; v3_scale( v2, s, t2);

      V3 t3; v3_add( t1, t2, t3); v3_normalise( t3);
      
      double el;
      sunpos( t3, T, NULL, &el);
      if (el > 0) k++;
   }

   return k/20.0;
}
 
///////////////////////////////////////////////////////////////////////////////
//  Spots File                                                               //
///////////////////////////////////////////////////////////////////////////////

//
//  Load the 'spots' file into memory.  The spots file is a plain text file
//  which maps symbolic site names into geographic coordinates.
//

static struct SPOT {
   char *name;
   V3 v;
}
 *spots = NULL;

static int nspots = 0;

static void alloc_spot( char *name, V3 v)
{
   VT_report( 3, "alloc spot [%s] %s", name, v3_string( v, NULL));

   spots = VT_realloc( spots, (nspots+1) * sizeof( struct SPOT));
   spots[nspots].name = strdup( name);
   v3_copy( v, spots[nspots].v);
   nspots++;
}

static int load_spots_file( char *filename)
{
   FILE *f = fopen( filename, "r");
   if (!f) return FALSE;

   int lino = 0;

   char temp[500], *p, *q, *s;
   while (fgets( temp, 500, f))
   {
      lino++;

      p = strchr( temp, '\r'); if (p) *p = 0;
      p = strchr( temp, '\n'); if (p) *p = 0;
      p = strchr( temp, ';'); if (p) *p = 0;

      p = temp;
      while (isspace( *p)) p++;
      if (!*p) continue;

      if (isalpha( *p))
      {
         VT_report( 0, "error in %s, line %d", filename, lino);
         continue;
      }

      for (q = p; *q && !isspace( *q); ) q++;
      if (*q) *q++ = 0;
      if ((s = strchr( p, ',')) == NULL)
      {
         VT_report( 0, "error in %s, line %d", filename, lino);
         continue;
      }

      *s++ = 0;
      double lat = parse_coord( p);
      double lon = parse_coord( s);
      V3 v;
      v3_make( lat, lon, v);

      p = q;
      int n = 0;
      while (1)
      {
         while (isspace( *p)) p++;
         if (!*p) break;
         
         for (q = p; *q && !isspace( *q); ) q++;
         if (*q) *q++ = 0;

         alloc_spot( p, v);
         n++;
         p = q;
      }

      if (!n)
      {
         VT_report( 0, "error in %s, line %d", filename, lino);
         continue;
      }
   }

   fclose( f);
   VT_report( 2, "loaded spots file [%s], %d entries", filename, nspots);
   return TRUE;
}

//
//  Attempt to load 'spots' file from the current directory.  If that fails,
//  try from the user's home directory.
//

static void load_spots( void)
{
   if (nspots) return;

   if (load_spots_file( "./spots")) return;
 
   char *home = getenv( "HOME");
   if (!home) return;

   char *path;
   if (asprintf( &path, "%s/spots", home) < 0 || !path) return;
   load_spots_file( path);
   free( path);
}

//
//  A case insensitive search for a symbolic site name in the spots table,
//  returning TRUE if found, with 'v' filled in with the location.
//

static int lookup_spot( char *s, V3 v)
{
   int i;

   for (i = 0; i < nspots; i++)
      if (!strcasecmp( spots[i].name, s))
      {
         v3_copy( spots[i].v, v);
         return TRUE;
      }

   return FALSE;
}

///////////////////////////////////////////////////////////////////////////////
//  Sites                                                                    //
///////////////////////////////////////////////////////////////////////////////

//
//  A list of sites referred to in measurement input records.
//

static struct SITE {
   V3 v;                 // Location of the site
   char *name;           // Symbolic name of the site.  If none, then the
                         // string representation of the lat/lon.
   char *file;           // Source file for sferic TOGAs when matching
   double offset;

   int n_sferic;         // Number of sferics supplied for matching
   int n_matched;        // Number of sferics successfully matched

   double tetotal;       // Overall sum of timing residuals, seconds
}
 sites[MAXSITES];

static int nsites = 0;

static int minsites = 0;     // Set by -n option
static int site_width = 0;   // Max length of site names

//
//  Add a site to the list, if not already listed.  'name' is optional.
//
//  XXX: this function assumes site location is fixed.
//

static struct SITE *add_site( V3 v, char *name)
{
   int i;
   for (i = 0; i < nsites; i++) if (v3_equal( sites[i].v, v)) return sites + i;

   if (nsites == MAXSITES) VT_bailout( "too many sites, max %d", MAXSITES);

   v3_copy( v, sites[nsites].v);
   sites[nsites].name = strdup( name);
   if (strlen( name) > site_width) site_width = strlen( name);

   return sites + nsites++;
}

static void clear_sites( void)
{
   int i;
   for (i = 0; i < nsites; i++)
      if (sites[i].name) free( sites[i].name);
   nsites = 0;
}

//
//  Parse a location.  This can be a symbolic place name and the coordinates
//  are looked up in the spots file (if available), or a latitude/longitude
//  pair.
//
//  Returns a pointer into the sites[] array.
//

static struct SITE *parse_latlon( char *s)
{
   char temp[150], *p;

   strcpy( temp, s);

   V3 v;
   if (isalpha( s[0]))
   {
      if (lookup_spot( s, v)) return add_site( v, temp);
      VT_bailout( "no definition in spots file for [%s]", s);
   }

   if ((p = strchr( temp, ',')) == NULL)
       VT_bailout( "bad lat/long [%s]", s);
   *p++ = 0;

   double lat = parse_coord( temp);
   double lon = parse_coord( p);

   v3_make( lat, lon, v);
   return add_site( v, temp);
}

///////////////////////////////////////////////////////////////////////////////
// Input Measurements                                                        //
///////////////////////////////////////////////////////////////////////////////

//
//  Measurement set structure.
//

struct MSET {

   //
   //  Table of locations and arrival times, usually TOGAs but could also be
   //  trigger times.
   //
   
   struct AT {
      struct SITE *site;
      timestamp toga;
      double sigma;
   } ats[MAXSITES];
   
   int nats;    // Number of elements of ats[]
   
   //
   //  Table of locations and bearings.
   //
   
   struct BR {
      struct SITE *site;
      double bearing;
      double sigma;
   } brs[MAXSITES];
   
   int nbrs;    // Number of elements of brs[]
   
   //
   //  Table of arrival time differences.
   //
   
   struct ATD {
      struct SITE *site1, *site2;
      double atd;                     // Arrival time difference, seconds
      double sigma;                   // Combined sigma of the two arrival times
      double range;                   // Baseline distance, km
   }
    atds[MAXSITES];
   
   int natd;    // Number of elements of atd[]
   
   char *ident;     // User-supplied ident for the measurement set
};
 
static void init_measurement_set( struct MSET *m)
{
   memset( m, 0, sizeof( struct MSET));
}

static void reset_measurement_set( struct MSET *m)
{
   m->nats = 0;
   m->nbrs = 0;
   m->natd = 0;
   if (m->ident) { free( m->ident); m->ident = NULL; }
}

#if 0   // Diagnostic use only

static void print_mset( struct MSET *m)
{
   int i;

   printf( "nats=%d\n", m->nats);
   for (i = 0; i < m->nats; i++)
      printf( "  %s %.6Lf %.6f\n",
          m->ats[i].site->name, m->ats[i].toga, m->ats[i].sigma);
   printf( "natd=%d\n", m->natd);
   for (i = 0; i < m->natd; i++)
      printf( " %s->%s %.6f\n",
         m->atds[i].site1->name, m->atds[i].site2->name, 
         m->atds[i].atd);
} 

#endif

static void parse_measurement( struct MSET *m, char *arg)
{
   char *p1 = strdup( arg), *p2, *p3, *p4 = NULL, *p5 = NULL;

   if (p1[1] != '/') VT_bailout( "cannot parse measurement [%s]", arg);
   p2 = p1 + 2;
   p3 = strchr( p2, '/');
   if (p3)
   {
      *p3++ = 0;
      p4 = strchr( p3, '/');
      if (p4)
      {
         *p4++ = 0;
         p5 = strchr( p4, '/');
         if (p5) *p5++ = 0;
      }
   }

   switch (*p1)
   {
      case 'T':  m->ats[m->nats].site = parse_latlon( p2);
                 if (!p3) VT_bailout( "missing timestamp in [%s]", arg);
                 m->ats[m->nats].toga = VT_parse_timestamp( p3);
                 m->ats[m->nats].sigma = p4 ? atof( p4) : DEFAULT_AT_SIGMA;
                 m->nats++;
                 break;

      case 'B':  m->brs[m->nbrs].site = parse_latlon( p2);
                 if (!p3) VT_bailout( "missing bearing in [%s]", arg);
                 m->brs[m->nbrs].bearing = atof( p3) * M_PI/180;
                 m->brs[m->nbrs].sigma = (p4 ? atof( p4) : 5) * M_PI/180;
                 m->nbrs++;
                 break;

      case 'A':  m->atds[m->natd].site1 = parse_latlon( p2);
                 if (!p3) VT_bailout( "missing location in [%s]", arg);
                 m->atds[m->natd].site2 = parse_latlon( p3);
                 if (!p4) VT_bailout( "missing ATD in [%s]", arg);
                 m->atds[m->natd].atd = atof( p4);
                 m->atds[m->natd].sigma = p5 ? atof( p5) : DEFAULT_AT_SIGMA;
                 m->natd++;
                 break;
                 
      case 'I':  if (m->ident) free( m->ident);
                 m->ident = strdup( p2);
                 break;

      default: VT_bailout( "unknown measurement type [%s]", arg);
   }

   free( p1);
}

///////////////////////////////////////////////////////////////////////////////
//  Propagation Model?                                                       //
///////////////////////////////////////////////////////////////////////////////

static inline double prop_delay( double range)
{
   if (cvlf_given) return range/CVLF;

   double vf = 1.0;

   if (range > 2000) vf = 0.985;
   else
   if (range > 1000) vf = 0.975;
   else
   if (range > 750) vf = 0.965;
   else
   if (range > 500) vf = 0.93;
   else vf = 0.90;

   return range/(300e3 * vf);
}

static inline double group_delay( V3 v1, V3 v2)
{
   return prop_delay( v3_range( v1, v2));
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

//
//  Compute a cost for a proposed position 'v' relative to the current
//  measurement set.  The cost is measured in units of standard deviations.
//  Each measurement in the set is compared with the value it should have if
//  the source is at location v and the result expressed in standard deviations
//  of that measurement.
//
//  We want to emphasise the worst-case deviation so that we can eliminate
//  sferics which don't belong to the intended stroke.  But returning just the
//  cost_max can sometimes cause the Nelder-Mead to hang up or converge
//  very slowly.  As a compromise, cost() returns a combination of mean
//  and worst case cost.
//

static double cost( struct MSET *m, V3 v)
{
   double cost_max = 0;
   double cost_sum = 0;

   int i;

   for (i = 0; i < m->natd; i++)
   {
      struct ATD *a = m->atds + i;
      double t1 = group_delay( v, a->site1->v);
      double t2 = group_delay( v, a->site2->v);
      double e = fabs(t1 - t2 - a->atd);
      double s = e/a->sigma;
      if (s > cost_max) cost_max = s;
      cost_sum += s;
   }

   for (i = 0; i < m->nbrs; i++)
   {
      struct BR *b = m->brs + i;
      double e = v3_bearing( b->site->v, v) - b->bearing;
      if (e > M_PI) e -= 2*M_PI;
      if (e < -M_PI) e += 2*M_PI;
      double s = fabs(e) / b->sigma;
      if (s > cost_max) cost_max = s;
      cost_sum += s;
   }

   //
   //  Return mean cost and max cost, weighted to emphasise the worst case.
   //

   double mean_cost = cost_sum / (m->natd + m->nbrs);
   return (mean_cost + 3 * cost_max)/4;
}

//
//  Estimate the timestamp of a source, given a location and a set of
//  arrival times.
//

static void source_time( struct MSET *m,
                         V3 location, timestamp *tp, double *stddev)
{
   *tp = timestamp_ZERO;
   if (stddev) *stddev = 0;

   if (!m->nats) return;

   timestamp T_base = m->ats[0].toga;
   int i;
   double sum = 0;
   double sumsq = 0;
   for (i = 0; i < m->nats; i++)
   {
      double diff = timestamp_diff( m->ats[i].toga, T_base)
                         - group_delay( location, m->ats[i].site->v);
      sum += diff;
      sumsq += diff * diff;
   }

   double mean = sum/m->nats;
   *tp = timestamp_add( T_base, mean);
   if (stddev) *stddev = sqrt( sumsq/m->nats - mean*mean);
}

//
//  Output a solution location, v.  'n' indicates which solution (0, 1, ...)
//  if there is more than one solution for the measurement set.  'rc' is the
//  residual cost for the solution. 'ns' is the number of sites used for
//  the trilateration.
//
//  stddev is the estimated RMS timing error, if known, otherwise zero.
//

static void output( struct MSET *m,
                    V3 location, timestamp T, double stddev,
                    double rc, int idx, int ns)
{
   // Prefix with an 'A' record identifier if using extended output

   if (OMODE_EXT) printf( "A ");

   // If the measurement set contained I/ident then start the output record
   // with the ident token.

   if (m->ident) printf( "%s ", m->ident);

   // Timestamp, either numeric or ISO format string
   char tstring[30];
   if (OMODE_ISO) VT_format_timestamp( tstring, T);
   else  timestamp_string6( T, tstring);

   printf( "%d %s %s %.2f %.2f %d\n",
          idx, v3_string( location, NULL), tstring, rc, stddev * 1e6, ns);
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

//
//  Parametric function of an ATD hyperbola.  'param' is the input parameter,
//  IS1 and IS2 are the two mirror image points associated with param.
//
//  Ra: Radius of circle centered on a->v2;
//  Rb: Radius of circle centered on a->v1;
//  Ra and Rb in radians subtended at Earth center;
//
//  Difference between Ra and Rb equals the ATD angle so that the circles
//  intersect on the ATD hyperbola.
//
//  Compute the two points IS1, IS2, at which the circles intersect.
//

static int hyperbola_point( struct ATD *a, double param, V3 IS1, V3 IS2)
{
   double Ra = param + (a->range - a->atd * CVLF)/2/EARTH_RAD;
   double Rb = param + (a->range + a->atd * CVLF)/2/EARTH_RAD;

   double CA = cos(Ra);
   double CB = cos(Rb);

   V3 Va;  v3_copy( a->site2->v, Va);
   V3 Vb;  v3_copy( a->site1->v, Vb);

   double Ax = Va[0], Ay = Va[1], Az = Va[2];
   double Bx = Vb[0], By = Vb[1], Bz = Vb[2];

   V3 VX; v3_cross( Va, Vb, VX);

   double SQRT = sqrt(
             -v3_dot(Va,Va)*CB*CB -v3_dot(Vb,Vb)*CA*CA +2*v3_dot(Va,Vb)*CA*CB
                   + v3_dot( VX,VX));

   // Sometimes, limited numeric precision can give a slight -ve argument to
   // sqrt(), when actually it should be zero or nearly zero.
   if (isnan( SQRT))
   {
      if (param) return FALSE;
      SQRT = 0;
   }

   A3 Ar = {
              { +Az*Az+Ay*Ay, -Ax*Ay,      -Ax*Az},
              { -Ax*Ay,       Az*Az+Ax*Ax, -Ay*Az} ,
              { -Ax*Az,       -Ay*Az,      Ay*Ay+Ax*Ax}
           };

   A3 Br = {
             { Bz*Bz+By*By, -Bx*By,       -Bx*Bz} ,
             { -Bx*By,      +Bz*Bz+Bx*Bx, -By*Bz} ,
             { -Bx*Bz,        -By*Bz,     By*By+Bx*Bx}
           };

   V3 FB;   v3_transform( Vb, Ar, FB);
   V3 FA;   v3_transform( Va, Br, FA);

   v3_scale( FA, CA, FA);
   v3_scale( FB, CB, FB);
   V3 F;   v3_add( FA, FB, F);

   v3_scale( VX, +SQRT, IS1);
   v3_add( IS1, F, IS1);

   v3_scale( VX, -SQRT, IS2);
   v3_add( IS2, F, IS2);

   v3_normalise( IS1);
   v3_normalise( IS2);
   return TRUE;
}

///////////////////////////////////////////////////////////////////////////////
//  Pattern Search                                                           //
///////////////////////////////////////////////////////////////////////////////

//
//  The pattern search is used in cases where the Nelder-Mead downhill simplex
//  cannot be used.  For example when there are only three arrival times, or
//  other combinations where there are two exact solutions.
//
//  Called with an approximate location 'v' for a source, pattern_search()
//  tries to improves the solution with reference to the cost() function.
//

#define MING (0.001 * M_PI/180)    // Convergence resolution (pattern search)

static double pattern_search( struct MSET *m, V3 v, double g)
{
   double best_cost = cost( m, v);
   VT_report( 2, "begin pattern_search with %.3f at %s",
                  best_cost, v3_string( v, NULL));

   int N = 0;

   //
   //  Set up a pair of orthogonal axes, a priori there is no preferred
   //  orientation so we just pick something at random.
   //

   V3 vr;   v3_random( vr);
   V3 vx;   v3_unit_normal_to( v, vr, vx);
   V3 vy;   v3_unit_normal_to( v, vx, vy);
 
   while (g >= MING)
   { 
      N++;

      A3 rx1;   a3_rot( vx, +g, rx1);
      A3 ry1;   a3_rot( vy, +g, ry1);
      A3 rx2;   a3_rot( vx, -g, rx2);
      A3 ry2;   a3_rot( vy, -g, ry2);

      //  Evaluate points in orthogonal directions away from the current
      //  point, continue until no further improvement at this scale.

      int n;  // Number of improvements made
      do {
         n = 0;

         V3 vt;
         double s;

         v3_transform( v, rx1, vt); 
         s = cost( m, vt);
         if (s < best_cost) { v3_copy( vt, v); best_cost = s; n++; }
      
         v3_transform( v, rx2, vt); 
         s = cost( m, vt);
         if (s < best_cost) { v3_copy( vt, v); best_cost = s; n++; }
      
         v3_transform( v, ry1, vt); 
         s = cost( m, vt);
         if (s < best_cost) { v3_copy( vt, v); best_cost = s; n++; }
      
         v3_transform( v, ry2, vt); 
         s = cost( m, vt);
         if (s < best_cost) { v3_copy( vt, v); best_cost = s; n++; }
      }
       while (n);

      //  No further improvement so reduce the scale factor
      g *= 0.8;
   }

   VT_report( 2, "pattern search iterations: %d", N);
   return best_cost;
}

//
//  Walk along a bearing line in steps of 'ps' radians, looking for the
//  minimum score.  Call pattern_search() to polish the position.
//

static int walk_bearing( struct MSET *m, double ps)
{
   A3 r0; a3_rot( m->brs[0].site->v, -m->brs[0].bearing, r0);
   V3 v0; v3_transform( V3north, r0, v0);

   V3 n0;  v3_unit_normal_to( m->brs[0].site->v, v0, n0);

   double r;
   double sc = 1e99;
   V3 v;
   for (r = 0.01 * M_PI/180; r<2*M_PI - 2*ps; r += ps)
   {
      A3 r1; a3_rot( n0, r, r1);
      V3 vt; v3_transform( m->brs[0].site->v, r1, vt);

      double st = cost( m, vt);
      if (st < sc) { v3_copy( vt, v); sc = st; }
   } 

   double c = pattern_search( m, v, ps);
   if (c < max_residual)
   {
      output( m, v, timestamp_ZERO, 0, c, 0, nsites);
      return TRUE;
   }

   return FALSE;
}

//
//  Walk along an ATD hyperbola in steps of 'ps' radians.  Find the minimum
//  score positions on each side of the ATD baseline and converge each one.
//

static int walk_atd( struct MSET *m, double ps)
{
   double sc1 = 1e99, sc2 = 1e99;
   int n = 2 * (int)(M_PI/ps); // Number of points to examine on the hyperbola

   static V3 *plist = NULL;
   static double *costs = NULL;

   if (!plist)
   {
      plist = VT_malloc( sizeof( V3) * n);
      costs = VT_malloc( sizeof( double) * n);
   }

   int i;
   for (i = 0; i < n/2; i++)
   {
      double r = ps/2 + i * ps;
      if (!hyperbola_point( m->atds+0, r, plist[n/2 - i], plist[n/2 + i]))
         costs[n/2 - i] = costs[n/2 + i] = 1e99; 
      else
         costs[n/2 - i] = cost( m, plist[n/2 - i]), 
         costs[n/2 + i] = cost( m, plist[n/2 + i]);
   }

   int i1 = -1, i2 = -1;
   for (i = 1; i < n-1; i++)
   {
      if (costs[i] > costs[i-1] ||
          costs[i] > costs[i+1]) continue;

      if (sc1 > sc2 && costs[i] < sc1) { sc1 = costs[i]; i1 = i; }     
      else
      if (sc2 >= sc1 && costs[i] < sc2) { sc2 = costs[i]; i2 = i; }     
   }

   timestamp T;
   double sdev;

   int r = FALSE;
   if (i1 >= 0)
   {
      double c1 = pattern_search( m, plist[i1], ps);
      if (c1 < max_residual)
      {
         source_time( m, plist[i1], &T, &sdev);
         output( m, plist[i1], T, sdev, c1, 0, nsites);
         r = TRUE;
      }
   }

   if (i2 >= 0)
   {
      double c2 = pattern_search( m, plist[i2], ps);
      if (c2 < max_residual)
      {
         source_time( m, plist[i2], &T, &sdev);
         output( m, plist[i2], T, sdev, c2, 1, nsites);
         r = TRUE;
      }
   }

   return r;
}

//
//  Scan an ATD hyperbola and to find an approximate location for the
//  minimum cost.
//

static void scan_atd( struct MSET *m, struct ATD *a, V3 best_location)
{
   double ps = 1.0 * M_PI/180;
   double best_cost = 1e99;

   int n = M_PI/ps; // Number of points to examine on the hyperbola

   int i;
   for (i = 0; i < n; i++)
   {
      double r = ps/2 + i * ps;
      V3 is1, is2;
      if (!hyperbola_point( a, r, is1, is2)) continue;

      double c;

      if ((c = cost( m, is1)) < best_cost)
      {
         best_cost = c;
         v3_copy( is1, best_location);
      }

      if ((c = cost( m, is2)) < best_cost)
      {
         best_cost = c;
         v3_copy( is2, best_location);
      }
   }
}

//
//  Scan a bearing line in steps of 'ps' radians, looking for the
//  minimum cost.
//

static void scan_bearing( struct MSET *m, struct BR *br, V3 best_location)
{
   double ps = 1.0 * M_PI/180;
   A3 r0; a3_rot( br->site->v, -br->bearing, r0);
   V3 v0; v3_transform( V3north, r0, v0);

   V3 n0;  v3_unit_normal_to( br->site->v, v0, n0);

   double r;
   double sc = 1e99;
   for (r = 0.01 * M_PI/180; r<2*M_PI - 2*ps; r += ps)
   {
      A3 r1; a3_rot( n0, r, r1);
      V3 vt; v3_transform( br->site->v, r1, vt);

      double st = cost( m, vt);
      if (st < sc) { v3_copy( vt, best_location); sc = st; }
   }
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

static void output_gc( struct SITE *s1, struct SITE *s2, double step)
{
   //
   //  The great circle between the two sites.
   //

   V3 gc; v3_unit_normal_to( s1->v, s2->v, gc);

   //
   //  Adjust step size to give an integer number of points on the GC
   //

   int istep = 0.5 + 2 * M_PI / step;
   step = 2 * M_PI/istep;

   //
   //  Generate points.
   //

   int i;
   for (i = 0; i < istep; i++)
   {
      A3 r1; a3_rot( gc, i * step, r1);
      V3 p;  v3_transform( s1->v, r1, p);
      printf( "%s\n", v3_string( p, NULL));
   }
}


//
//  This method works with arrival time differences and bearings.  Therefore
//  it can perform a 2-D minimising of the cost function over the lightning
//  position without considering the lightning timestamp.  The timestamp is
//  estimated separately once the location has been determined.
//

static int downhill_simplex( struct MSET *m, V3 location, double *costp)
{
   // print_mset( m);   // Diagnostic use only

   int i;

   //
   //  Three point simplex and cost sorting function.
   //

   struct SIMPLEX {
      V3 v;
      double cost;
   } P[3];

   int cmp_cost( const void *a1, const void *a2)
   {
      struct SIMPLEX *p1 = (struct SIMPLEX *) a1;
      struct SIMPLEX *p2 = (struct SIMPLEX *) a2;
      if (p1->cost < p2->cost) return -1;
      if (p1->cost > p2->cost) return +1;
      return 0;
   }

   //
   //  Initialise the simplex by scanning the ATD hyperbolas of the first
   //  three ATDs, finding the minimum cost location of each.  I have tried
   //  several other methods to choose initial simplex corners and this one
   //  is by far the most robust.
   //
   //  If there are not enough arrival time differences, then use the available
   //  bearings in the same way.
   //

   int ni;
   for (ni = 0; ni < 3 && ni < m->natd; ni++)
      scan_atd( m, m->atds + ni, P[ni].v);
   for (i = 0; ni < 3 && i < 3 && i < m->nbrs; i++, ni++)
      scan_bearing( m, m->brs + i, P[ni].v);

   //
   //  A completely standard Nelder-Mead.
   //

   double alpha = 1.0, gamma = 2.0, rho = 0.5, sigma = 0.5;

   for (i = 0; i < 3; i++) P[i].cost = cost( m, P[i].v);

   int NI = 0;    // Number of iterations
   double previous_cost = 0;

   while (1)
   {
      NI++;
      qsort( P, 3, sizeof( struct SIMPLEX), cmp_cost);

      // Termination conditions 
      if (NI > 50 && previous_cost - P[0].cost < 0.001)  break;
      if (NI > 200) break;

      previous_cost = P[0].cost;

      V3 x0;  v3_add( P[0].v, P[1].v, x0);  v3_normalise( x0);

      // 
      //  Reflection.
      //

      V3 xr;
      xr[0] = x0[0] + alpha * (x0[0] - P[2].v[0]);
      xr[1] = x0[1] + alpha * (x0[1] - P[2].v[1]);
      xr[2] = x0[2] + alpha * (x0[2] - P[2].v[2]);

      double cr = cost( m, xr);
      if (cr >= P[0].cost && cr < P[1].cost)
      {
         v3_copy( xr, P[2].v);
         P[2].cost = cr;
         continue;
      }

      //
      //  Expansion.
      //

      if (cr < P[0].cost)
      {
         V3 xe;
         xe[0] = x0[0] + gamma * (xr[0] - x0[0]);
         xe[1] = x0[1] + gamma * (xr[1] - x0[1]);
         xe[2] = x0[2] + gamma * (xr[2] - x0[2]);
         double ce = cost( m, xe);
         if (ce < cr)
         {
            v3_copy( xe, P[2].v);
            P[2].cost = ce;
         }
         else
         {
            v3_copy( xr, P[2].v);
            P[2].cost = cr;
         }
         continue;
      }

      //
      //  Contraction.
      //

      V3 xc;
      xc[0] = x0[0] + rho * (P[2].v[0] - x0[0]);
      xc[1] = x0[1] + rho * (P[2].v[1] - x0[1]);
      xc[2] = x0[2] + rho * (P[2].v[2] - x0[2]);
      double cc = cost( m, xc);

      if (cc < P[2].cost)
      {
         v3_copy( xc, P[2].v);
         P[2].cost = cc;
         continue;
      }

      //
      //  Shrink.
      //

      for (i = 1; i <= 2; i++)
      {
         P[i].v[0] = P[0].v[0] + sigma * (P[i].v[0] - P[0].v[0]);
         P[i].v[1] = P[0].v[1] + sigma * (P[i].v[1] - P[0].v[1]);
         P[i].v[2] = P[0].v[2] + sigma * (P[i].v[2] - P[0].v[2]);
      }
   }

   if (P[0].cost < max_residual)   // Successful?
   {
      if (costp) *costp = P[0].cost;
      v3_copy( P[0].v, location);
      return TRUE;
   }

   return FALSE;
}

static int process_measurement_set( struct MSET *m)
{
   int i;

   //
   //  For diagnostics, report arrival times and bearings.
   //

   for (i = 0; i < m->nats; i++)
   {
      char temp[30];
      timestamp_string6( m->ats[i].toga, temp);
      VT_report( 2, "at %s %s (%.6f)",
                      v3_string( m->ats[i].site->v, NULL),
                      temp, m->ats[i].sigma);
   }

   for (i = 0; i < m->nbrs; i++)
      VT_report( 2, "br %s %.1f (%.1f)",
                      v3_string( m->brs[i].site->v, NULL),
                      m->brs[i].bearing * 180/M_PI,
                      m->brs[i].sigma * 180/M_PI);

   //
   //  Take the arrival times ats[] of the measurement set and build a list of
   //  all the arrival time differences.   If any ATD exceeds the baseline
   //  this function returns without attempting a solution.
   //

   for (i = 0; i < m->nats - 1; i++)
   {
      struct ATD *a = m->atds + m->natd;
      a->site1 = m->ats[i+0].site;
      a->site2 = m->ats[i+1].site;
      a->atd = timestamp_diff( m->ats[i+0].toga, m->ats[i+1].toga);
      a->range = v3_range( a->site1->v, a->site2->v);
      double tbase = prop_delay( a->range);

      if (fabs( a->atd) > tbase)
      {
         VT_report( 1, "ATD out of range %s %s %.6f limit %.6f",
                           a->site1->name, a->site2->name, a->atd, tbase);
         return FALSE;
      }

      a->sigma = m->ats[i].sigma;
      if (m->ats[i+1].sigma > a->sigma) a->sigma = m->ats[i+1].sigma;

      VT_report( 1, "ATD %s %s %.6f", a->site1->name, a->site2->name, a->atd);
      m->natd++;
   }

   if (m->nats == 0 && m->nbrs == 2)
   {
      // Special case: intersection of two bearings.  Calculate the exact
      // solutions without iteration.

      // Make 2nd point on each GC by rotating V3north around each point
      A3 r0; a3_rot(  m->brs[0].site->v, -m->brs[0].bearing, r0);
      V3 v0; v3_transform( V3north, r0, v0);

      A3 r1; a3_rot(  m->brs[1].site->v, -m->brs[1].bearing, r1);
      V3 v1; v3_transform( V3north, r1, v1);

      // Define GCs by their normal unit vector
      V3 n0;  v3_unit_normal_to( m->brs[0].site->v, v0, n0);
      V3 n1;  v3_unit_normal_to( m->brs[1].site->v, v1, n1);

      V3 it1;  v3_unit_normal_to( n0, n1, it1);
      V3 it2;  v3_scale( it1, -1, it2);

      output( m, it1, timestamp_ZERO, 0, 0, 0, 2);
      output( m, it2, timestamp_ZERO, 0, 0, 1, 2);
      return TRUE;
   }

   //
   //  If enough arrival time differences and bearings, use Nelder-Mead to find
   //  the single best location.
   //

   if (m->natd + m->nbrs >= 3)
   {
      V3 v;
      double c;
      if (downhill_simplex( m, v, &c))
      {
         timestamp T = timestamp_ZERO;
         double sdev = 0;
         source_time( m, v, &T, &sdev);

         output( m, v, T, sdev, c, 0, 0);
         return TRUE;
      }
      return FALSE;
   }

   double ps = 1.0 * M_PI/180;

   // If we have any arrival time differences, then pick the first ATD
   // baseline and scan that for approximate solutions.   Otherwise scan
   // the first bearing line of the measurement set.

   if (m->natd)
   {
      return walk_atd( m, ps);
   }
   else
   if (m->nbrs)
   {
      return walk_bearing( m, ps);
   }

   return FALSE;
}

///////////////////////////////////////////////////////////////////////////////
// TOGA Matching                                                             //
///////////////////////////////////////////////////////////////////////////////

static void parse_matching( char *arg)
{
   char *s = arg;

   // -m site=file

   char *p = strchr( s, '=');
   if (!p) VT_bailout( "invalid -m argument [%s]", arg);
   *p++ = 0;

   struct SITE *site = parse_latlon( s);
   site->file = strdup( p);
   VT_report( 2, "matching: %s from %s", site->name, site->file);
}

static void parse_offset( char *arg)
{
   char *s = arg;

   // -t site=offset

   char *p = strchr( s, '=');
   if (!p) VT_bailout( "invalid -t argument [%s]", arg);
   *p++ = 0;

   struct SITE *site = parse_latlon( s);
   site->offset = atof( p);
   VT_report( 1, "offset: %s %.6f", site->name, site->offset);
}

//
//  An aggregate array of all the sferics to be matched and solved.
//

static struct TOGALIST {
   timestamp toga;    // The timestamp from vttoga
   struct SITE *site; // Link to the site 
   float rms;        // RMS amplitude reported by vttoga
   int16_t azimuth;   // Radians
   int16_t flags;     // Bit flags
}
 *togalist = NULL;

#define HAS_AZIMUTH (1 << 0)
#define USED (1 << 1)

static int nmatch = 0;        // Number of entries in togalist[]
static int match_alloc = 0;   // Allocated entries in togalist[]

static void add_togalist( struct SITE *site, timestamp toga,
                          double rms,
                          int has_az, double azimuth)
{
   // Extend togalist[] in blocks of 1000
   if (nmatch == match_alloc)
   {
      match_alloc += 1000;
      togalist = VT_realloc( togalist, sizeof( struct TOGALIST) * match_alloc);
   }

   struct TOGALIST *p = togalist + nmatch++;

   p->site = site;
   p->toga = toga;
   p->rms = rms;
   p->flags = has_az ? HAS_AZIMUTH : 0;
   p->azimuth = azimuth;
}

static int filter_solution( struct TOGALIST **list, int nlist, uint64_t u,
                            V3 stroke)
{
   struct S5 {
      struct SITE *site;
      double bearing;
      double range;
   } s[64];

   int ns = 0;

   int cmp_bearing( const void *s1, const void *s2)
   {
      struct S5 *p1 = (struct S5 *) s1;
      struct S5 *p2 = (struct S5 *) s2;
      if (p1->bearing < p2->bearing) return -1;
      if (p1->bearing > p2->bearing) return +1;
      return 0;
   }

   if (__builtin_popcountl( u) < minsites) return FALSE;

   //
   //  Link the sferics selected by the bit pattern 'u', into s[].
   //  Compute the range and bearing of the receiver, from the standpoint
   //  of the lightning stroke.
   //

   int i;
   for (i = 0; i < nlist; i++)
      if (u & (1L << i))
      {
         s[ns].site = list[i]->site;
         s[ns].range = v3_range( stroke, s[ns].site->v);
         s[ns].bearing = v3_bearing( stroke, s[ns].site->v);
         ns++;
      }

   if (filter_envelop)
   {
      //  Sort into bearing order so that the receivers are ordered clockwise
      //  around the stroke.

      qsort( s, ns, sizeof( struct S5), cmp_bearing);

      //
      //  Look for any interval greater than filter_envelop degrees.
      //

      for (i = 0; i < ns; i++)
      {
         double db = i < ns-1 ? s[i+1].bearing - s[i].bearing
                              : s[0].bearing - s[i].bearing;

         if (db < 0) db += 2*M_PI;
         if (db > filter_envelop * M_PI/180) return FALSE;
      }
   }

   if (filter_nearmax)
   {
      //
      //  Discard the solution if the nearest receiver is more than
      //  filter_nearmax km distant from the stroke.
      //

      for (i = 0; i < ns; i++) if (s[i].range < filter_nearmax) break;
      if (i == ns) return FALSE;
   }

   //
   // XXX Other filtering options go in here.
   //

   return TRUE;
}

static int prepare_measurement_set( struct MSET *m, 
                                    struct TOGALIST **list, int n,
                                    uint64_t mask)
{
   int i, j;

   reset_measurement_set( m);

   for (i = 0; i < n; i++)
   {
      if ((mask & (1L << i)) == 0) continue;  // This sferic not selected?
    
      m->ats[m->nats].site = list[i]->site;
      m->ats[m->nats].toga = list[i]->toga;
      m->ats[m->nats].sigma = DEFAULT_AT_SIGMA;
      m->nats++;

      if (list[i]->flags & HAS_AZIMUTH)
      {
         m->brs[m->nbrs].site = list[i]->site;
         m->brs[m->nbrs].bearing = list[i]->azimuth;
         m->brs[m->nbrs].sigma = DEFAULT_AZ_SIGMA * M_PI/180;
         m->nbrs++;
      }
   }

   //
   //  Check baseline ATD limits and build list of independent ATDs.
   //

   for (i = 0; i < m->nats - 1; i++)
      for (j = i+1; j < m->nats; j++)
      {
         //  All baselines must be checked, to help remove sferics
         //  that come from different strokes.

         double atd = timestamp_diff( m->ats[i].toga, m->ats[j].toga);
   
         double range = v3_range( m->ats[i].site->v, m->ats[j].site->v);
         double tbase = prop_delay( range);
   
         if (fabs(atd) > 1.001 * tbase) return FALSE;

         //  Add independent ATDs to atds[]

         if (j == i+1)
//         if (i == 0)
         {
            struct ATD *a = m->atds + m->natd;
      
            a->site1 = m->ats[i].site;
            a->site2 = m->ats[j].site;

            a->sigma = m->ats[i].sigma;
            if (m->ats[j].sigma > a->sigma) a->sigma = m->ats[j].sigma;
      
            a->range = range;
            a->atd = atd;
            m->natd++;
         }
      }

   return TRUE;
}

//
//  Attempt to find a stroke solution from the array of up to 64 sferic.
//  All sferics earlier than this list have already been dealt with, either
//  used or exhaustively tried.
//
//  The list is in timestamp order.
//
//  Attempt to find a solution which involves the sferic at the head of the
//  list, plus some number of other later sferics from the list.
//
//  The strategy is to find a solution amongst the first eight sferics by
//  trying all permutations involving the first sferic and four or more others.
//  Then add additional later sferics one at a time, discarding those that
//  don't fit.
//

static pthread_mutex_t output_lock = PTHREAD_MUTEX_INITIALIZER;

static int n_strokes = 0;      // Total number of solutions generated. Doubles
                               // as an index number for output records

static void process_matchlist( struct TOGALIST **list, int nlist)
{
   int i;
   struct MSET m;
   init_measurement_set( &m);

   //  Can only consider up to 64 sferics at a time because we use a 64 bit
   //  word to keep track of the selections

   if (nlist > 64) nlist = 64;

   //
   //  Search permutations of sferics at the head of the list, until some
   //  solution is found.  Permutations of up to max_initial sferics are tried.
   //  We want at least min_initial number of sites in the initial solution.
   //

   int max_initial = 8;
   if (max_initial > nlist) max_initial = nlist;

   int min_initial = 6;
   if (min_initial > minsites) min_initial = minsites;

   V3 best_location;      // The stroke location for the best pattern
   uint64_t best_u = 0;   // The bit pattern which selects the best sferics in
                          // the list
   double best_cost = 0;  // The cost() of the best pattern
   int best_weight = 0;   // Number of sferics in the best pattern: the hamming
                          // weight of best_u

   uint64_t u;  // Bit pattern selecting the list entries currently being tried
 
   uint32_t npu = 1 << max_initial;
   uint32_t npm = (1 << min_initial) - 1;

   for (u = npm; u < npu; u += 2) // Only odd patterns are tried because the
                                  // first entry in the toga list must always be
                                  // present
   {
      int weight = __builtin_popcountl( u);
      if (weight < min_initial) continue;  // We want a minimum number of
                                           // sferics in the initial solution

      //
      //  prepare_measurement_set() will return FALSE if any ATDs exceed a
      //  baseline - this also removes permutations where a site has more than
      //  one sferic in the selection.
      //

      if (prepare_measurement_set( &m, list, max_initial, u) &&
          downhill_simplex( &m, best_location, &best_cost))
      {
         best_u = u;
         best_weight = weight;
         break;
      }
   }

   if (!best_weight) // No initial solution found?
   {
      reset_measurement_set( &m);
      return;   
   }

   //
   //  Extend the solution by adding one additional sferic at a time,
   //  keeping those that still solve.
   //

   for (i = 0; i < nlist; i++)
   {
      if (best_u & (1L << i)) continue;  // Sferic already in the list?

      u = best_u | (1L << i);    //  Add sferic 'i' to the pattern

      // XXX should arrange to continue with the previous simplex instead of
      // starting from scratch

      if (prepare_measurement_set( &m, list, nlist, u) &&
          downhill_simplex( &m, best_location, &best_cost))
      {
         best_weight = __builtin_popcountl( u);
         best_u = u;
      }
   }

   //
   //  Apply filtering rules to the solution.  If a rule fails, the solution
   //  is discarded and the sferics in list[] except the first one, remain
   //  available for further attempts.
   //

   if (!filter_solution( list, nlist, best_u, best_location))
   {
      reset_measurement_set( &m);
      return;
   }

   //
   //  Mark the successful sferics as used in togalist[].  Reload ats[] with
   //  the final selection of sferics - needed by source_time().
   //

   prepare_measurement_set( &m, list, nlist, best_u);

   for (i = 0; i < nlist; i++) if (best_u & (1L << i)) list[i]->flags |= USED;

   timestamp T = timestamp_ZERO;
   double sdev = 0;
   source_time( &m, best_location, &T, &sdev);

   pthread_mutex_lock( &output_lock);

   for (i = 0; i < nlist; i++)
      if (best_u & (1L << i))
      {
         double diff = timestamp_diff( list[i]->toga, T)
                         - group_delay( best_location, list[i]->site->v);
         list[i]->site->tetotal += diff;
         list[i]->site->n_matched++;
      }

   output( &m, best_location, T, sdev, best_cost, n_strokes, best_weight);

   //
   //  Extended output requested?
   //

   if (OMODE_EXT)
   { 
      for (i = 0; i < nlist; i++)
         if (best_u & (1L << i))
         {
            double range = v3_range( best_location, list[i]->site->v);
            double terr =
                timestamp_diff( list[i]->toga, T) - prop_delay( range);

            double illum =
                path_illumination( best_location, list[i]->site->v, T);

            // Timestamp, either numeric or ISO format string
            char tstring[30];
            if (OMODE_ISO) VT_format_timestamp( tstring, list[i]->toga);
            else  timestamp_string6( list[i]->toga, tstring);

            printf( "R %-*s %7.1f %6.1f %5.1f %.3e %.2f %s\n",
               site_width, list[i]->site->name,
               range,
               v3_bearing( best_location, list[i]->site->v) * 180/M_PI,
               terr * 1e6,
               list[i]->rms, illum, tstring);
         }

      printf( "E\n");
   }

   //
   //  Output a .brc file for map plotting.
   //

   if (OMODE_BRC)
   {
      char filename[100];
      sprintf( filename, "brc/%d.brc", n_strokes);
      FILE *f = fopen( filename, "w");
      if (f)
      {
         for (i = 0; i < nlist; i++)
            if (u & (1L << i))
               fprintf( f, "spot %s 0.005\n",
                 v3_string( list[i]->site->v, NULL));
         fprintf( f, "spot %s 0.005 S\n", v3_string( best_location, NULL));
         fclose( f);
      }
   }

   n_strokes++;

   pthread_mutex_unlock( &output_lock);
}

//
//  Load all the sferic data from the files list with -m options.
//
//  The data goes into a single array togalist[] which is then sorted
//  into timestamp order.
//

static double batch_duration = 0;  // Overall time spanned by the batch of 
                                   // TOGA files.

static void load_matching_files( void)
{
   int i;

   for (i = 0; i < nsites; i++)
   {
      struct SITE *site = sites + i;

      FILE *f = fopen( site->file, "r");
      if (!f) VT_bailout( "cannot open %s: %s", site->file, strerror( errno));

      char buff[256];
      while (fgets( buff, 255, f))
      {
         if (strncmp( buff, "H ", 2)) continue;

         long double T;
         double rms;
         double dummy;
        
         int n;               // Number of chars used by the first three fields
         if (sscanf( buff+2, "%Lf %lg %lg%n",
                             &T, &rms, &dummy, &n) != 3)
            VT_bailout( "format error in [%s]", site->file);

         // Optional bearing supplied by the site?
         double azimuth;
         int has_az = sscanf( buff + 2 + n, "%lg", &azimuth) == 1;
            
         timestamp TS = timestamp_compose( (int)T, T - (int)T);

         add_togalist( site, timestamp_add( TS, site->offset),
                       rms, has_az, azimuth * M_PI/180);
         site->n_sferic++;
      }

      fclose( f);
   }

   if (!nmatch) VT_bailout( "no sferics found to match");

   //
   //  Sort the sferics into timestamp order.
   //

   int cmp_togalist( const void *a1, const void *a2)
   {
      struct TOGALIST *p1 = (struct TOGALIST *) a1;
      struct TOGALIST *p2 = (struct TOGALIST *) a2;
      if (timestamp_LT( p1->toga, p2->toga)) return -1;
      if (timestamp_GT( p1->toga, p2->toga)) return +1;
      return 0;
   }
   
   qsort( togalist, nmatch, sizeof( struct TOGALIST), cmp_togalist);

   VT_report( 1, "loaded %d sferics from %d sites", nmatch, nsites);

   batch_duration = timestamp_diff( togalist[nmatch-1].toga, togalist[0].toga);
   VT_report( 1, "duration %.1f seconds", batch_duration);
   VT_report( 1, "average rate %.1f sferics/second", nmatch/batch_duration);
}

struct MATCH_LAUNCH {
   struct TOGALIST *togalist;
   int nmatch;
   int done;
};

static void *task_matching( void *arg)
{
   struct MATCH_LAUNCH *L = (struct MATCH_LAUNCH *)arg;
   int used[MAXSITES];

   int i;

   //
   //  All the sferics from all the sites are merged into togalist[] in
   //  timestamp order.   Work sequentially through the list, trying to find
   //  sferics which belong to the same stroke as the togalist[base] sferic.
   //

   double list_window = prop_delay( M_PI * EARTH_RAD);

   int base;
   for (base = 0; base < L->nmatch - 3; base++)
   {
      //
      //  All sferics togalist[0] to togalist[base-1] have been either used or
      //  exhaustively tried.   Try to find a solution involving togalist[base],
      //  with some number of later sferics.
      //

      struct TOGALIST *bs = L->togalist + base;

      if (bs->flags & USED) continue;  // This TOGA was successful in a previous
                                       // matching, so is not to be used again

      //
      //  Build a list of TOGAs over which to attempt matching.  Up to 64
      //  sferics can be considered.
      //

      for (i = 0; i < nsites; i++) used[i] = FALSE;
      
      struct TOGALIST *matchlist[64];    // Empty list of TOGAs
      int nmatchlist = 0;

      matchlist[nmatchlist++] = bs;  // Add our base entry
      int sites_in_matchlist = 1;    // Number of sites in the matching list
      used[bs->site - sites] = TRUE;

      int j;
      for (j = base+1; j < L->nmatch && nmatchlist < 64; j++)
      {
         if (L->togalist[j].flags & USED) // This TOGA was successful in a
            continue;                     // previous measurement set, so is
                                          // not to be used again

         // Don't include further sferics from the base site
         if (L->togalist[j].site == bs->site) continue; 

         double atd = timestamp_diff( L->togalist[j].toga, bs->toga);
         if (atd > list_window) break;

         double baseline_atd = group_delay( bs->site->v, L->togalist[j].site->v);
         if (atd > 1.0 * baseline_atd) continue;

         matchlist[nmatchlist++] = L->togalist + j;   // Accept into the list

         // Keep a count of sites in the list
         if (!used[L->togalist[j].site - sites]) 
         {
            used[L->togalist[j].site - sites] = TRUE;
            sites_in_matchlist++;
         }
      }

      //
      //  Now matchlist[] contains a list of TOGAs which can plausibly belong
      //  to the same stroke.   A site may contribute more than one sferic
      //  into matchlist[].
      // 

      if (sites_in_matchlist >= minsites)
         process_matchlist( matchlist, nmatchlist);
   }

   L->done = TRUE;
   return 0;
}

static void run_matching( void)
{
   pthread_t pids[ncpu];
   memset( pids, 0, sizeof( pids));

   struct MATCH_LAUNCH L[ncpu];
   double mingap = prop_delay( M_PI * EARTH_RAD);

   int partition_base = 0;
   int cpu = 0;

   if (ncpu == 1)
   {
      //
      //  Only one CPU available so no point in partitioning.
      //

      L[0].togalist = togalist;
      L[0].nmatch = nmatch;
      task_matching( (void *)&L[0]);
   }
   else
   {
      pthread_attr_t pa;
      pthread_attr_init( &pa);
      pthread_attr_setdetachstate( &pa, PTHREAD_CREATE_DETACHED);

      //
      //  Partition togalist[] where we find gaps in the time series of TOGAs.
      //  Run each partition in a separate thread.
      //  To avoid nibbling, ensure at least 500 sferics in a partition.
      //

      while (partition_base < nmatch - 1)
      {
         // Search forward through the togalist[] from partition_base to find a
         // gap in the TOGAs exceeding mingap seconds.
   
         int n = partition_base + 500;
         if (n > nmatch) n = nmatch;
         for (; n < nmatch; n++)
         {
            double dt = timestamp_diff( togalist[n].toga, togalist[n-1].toga);
            if (dt > mingap) break;
         }
   
         // No point in running this partition if it has less than minsites
         // sferics
         if (n - partition_base >= minsites)
         {
            // Find a spare CPU
            while (1)
            {
               for (cpu = 0; cpu < ncpu; cpu++)
                  if (!pids[cpu] || L[cpu].done) break;
      
               if (cpu < ncpu) break;   // Found an expired thread?
      
               usleep( 1000);
            }
      
            L[cpu].togalist = togalist + partition_base;
            L[cpu].nmatch = n - partition_base;
            L[cpu].done = FALSE;
            pthread_create( &pids[cpu], &pa, task_matching, (void *)&L[cpu]);
         }
   
         partition_base = n;
      }
   
      for (cpu = 0; cpu < ncpu; cpu++)
         while (pids[cpu] && !L[cpu].done) usleep( 1000);
      pthread_attr_destroy( &pa);
   }


   int n_used = 0;
   int i;
   for (i = 0; i < nmatch; i++) if (togalist[i].flags & USED) n_used++;

   VT_report( 1, "sferics successfully used: %d (%.1f%%)",
                  n_used, n_used * 100.0/nmatch);
   VT_report( 1, "strokes found: %d", n_strokes);

   if (n_strokes)
      VT_report( 1, "average sferics per stroke: %.1f",
                     n_used/(double) n_strokes);

   if (batch_duration)
      VT_report( 1, "average strokes per second: %.2f",
                        n_strokes/batch_duration);

   //
   //  Extended output?  Add the site summary 'G' records and a final 'Q'
   //  record.
   //

   if (OMODE_EXT)
   {
      for (i = 0; i < nsites; i++)
      {
         struct SITE *s = sites + i;

         printf( "G %-*s %6u %6u %6.1f\n",
             site_width, s->name,
             s->n_sferic, s->n_matched,
             s->n_matched ? s->tetotal/s->n_matched * 1e6 : 0);
      }

      printf( "Q %d %d %d %.1f\n", nmatch, n_used, n_strokes, batch_duration);
   }
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

int main( int argc, char *argv[])
{
   VT_init( "vtspot");

   load_spots();

   while (1)
   {
      int c = getopt( argc, argv, "vL:c:m:n:r:f:o:t:U:gsbd?");

      if (c == 'v') VT_up_loglevel();
      else
      if (c == 'L') VT_set_logfile( "%s", optarg);
      else
      if (c == 'b') BFLAG = 1;
      else
      if (c == 'U') ncpu = atoi( optarg);
      else
      if (c == 'c')
      {
         CVLF = 300e3 * atof( optarg);
         cvlf_given = TRUE;
      }
      else
      if (c == 'd') DFLAG = 1;
      else
      if (c == 's') SFLAG = 1;
      else
      if (c == 'g') GFLAG = 1;
      else
      if (c == 'r') max_residual = atof( optarg);
      else
      if (c == 'n') minsites = atoi( optarg);
      else
      if (c == 'f')   // Filtering options
      {
         if (!strncmp( optarg, "envelop=", 8))
            filter_envelop = atof( optarg + 8);
         else
         if (!strncmp( optarg, "nearmax=", 8))
            filter_nearmax = atof( optarg + 8);
         else
            VT_bailout( "unrecognised -f option [%s]", optarg);
      }
      else
      if (c == 'o')   // Output options
      {
         if (!strcmp( optarg, "iso")) OMODE_ISO = TRUE;
         else
         if (!strcmp( optarg, "ext")) OMODE_EXT = TRUE;
         else
         if (!strcmp( optarg, "brc")) OMODE_BRC = TRUE;
         else
            VT_bailout( "unrecognised -o option [%s]", optarg);
      }
      else
      if (c == 't') parse_offset( strdup( optarg));
      else
      if (c == 'm') parse_matching( strdup( optarg));
      else
      if (c == -1) break;
      else
         usage();
   }

   if (ncpu < 1) VT_bailout( "invalid number of worker threads");

   if (BFLAG)
   {
      if (argc == optind)
      {
         // Read standpoint and forepoint from standard input

         char *p, *q0, *q1, *q2, temp[100];
         while (fgets( temp, 99, stdin))
         {
            // Format: standpoint forepoint rest

            // Remove line terminators
            p = strchr( temp, '\n');  if (p) *p = 0;
            p = strchr( temp, '\r');  if (p) *p = 0;

            // Skip leading white space
            q0 = temp;
            while (*q0 && isspace( *q0)) q0++;
            if (!*q0) VT_bailout( "invalid input format");

            // Standpoint
            q1 = strchr( q0, ' ');
            if (!q1) VT_bailout( "invalid input format");
            *q1++ = 0;
            while (*q1 && isspace( *q1)) q1++;
           
            q2 = strchr( q1, ' ');
            if (q2) *q2++ = 0;
            
            struct SITE *s1 = parse_latlon( q0);  // Standpoint
            struct SITE *s2 = parse_latlon( q1);  // Forepoint

            double d = v3_range( s1->v, s2->v);
            double b = v3_bearing( s1->v, s2->v);
            b = constrain( b, 0, 2*M_PI);
            printf( "%.3f %.1f", d, b * 180/M_PI);

            if (q2 && *q2)
            {
               // Maybe an optional timestamp, plus additional fields to
               // pass through
               timestamp T = timestamp_NONE;

               char *q3 = strchr( q2, ' ');
               if (q3) *q3++ = 0;
               T = VT_parse_timestamp( q2);
               printf( " %.2f %s", path_illumination( s1->v, s2->v, T), q2);
               if (q3 && *q3) printf( " %s\n", q3);
            }
            else printf( "\n");

            clear_sites();
         }
      }
      else
      if (argc - optind == 2)
      {
         // Standpoint and forepoint are next on the command line
         struct SITE *s1 = parse_latlon( argv[optind++]);  // Standpoint
         struct SITE *s2 = parse_latlon( argv[optind++]);  // Forepoint

         double d = v3_range( s1->v, s2->v);
         double b = v3_bearing( s1->v, s2->v);
         b = constrain( b, 0, 2*M_PI);
         printf( "%.3f %.1f\n", d, b * 180/M_PI);
      }
      else
      if (argc - optind == 3)
      {
         // Standpoint and forepoint are next on the command line, followed
         // by a timestamp
         struct SITE *s1 = parse_latlon( argv[optind++]);  // Standpoint
         struct SITE *s2 = parse_latlon( argv[optind++]);  // Forepoint
         timestamp T = VT_parse_timestamp( argv[optind]);
         double d = v3_range( s1->v, s2->v);
         double b = v3_bearing( s1->v, s2->v);
         b = constrain( b, 0, 2*M_PI);
         printf( "%.3f %.1f %.2f\n", d, b * 180/M_PI,
                    path_illumination( s1->v, s2->v, T));
      }
      else usage();
      return 0;
   }

   if (GFLAG)
   {
      if (argc - optind >= 2)
      {
         struct SITE *s1 = parse_latlon( argv[optind++]);
         struct SITE *s2 = parse_latlon( argv[optind++]);

         double step = 2 * M_PI/100;
         if (argc - optind == 1) step = atof( argv[optind++]) * M_PI/180;
         output_gc( s1, s2, step); 
      }
      else usage();

      return 0;
   }

   if (DFLAG)
   {
      if (argc - optind != 3) usage();
      V3 vf;
      struct SITE *s = parse_latlon( argv[optind++]);
      double b = atof( argv[optind++]) * M_PI/180;
      double a = atof( argv[optind])/EARTH_RAD;
      destination_point( s->v, b, a, vf);
      printf( "%s\n", v3_string( vf, NULL));
      return 0;
   }

   if (SFLAG)
   {
      if (argc - optind != 2) usage();

      struct SITE *s = parse_latlon( argv[optind++]);
      timestamp T = VT_parse_timestamp( argv[optind]);
      double solaz, solel;
      sunpos( s->v, T, &solaz, &solel);
      printf( "%.1f %.1f\n", solaz * 180/M_PI, solel * 180/M_PI);
      return 0;
   }

   if (optind < argc)
   {
      // Take a measurement set from the command line
      struct MSET m;
      init_measurement_set( &m);
      reset_measurement_set( &m);
      while (optind < argc) parse_measurement( &m, argv[optind++]);
      process_measurement_set( &m);
      return 0;
   }

   if (nsites)  // Some -m options given?
   {
      time_t tstart = time( NULL);   // To measure running time
      if (!minsites) minsites = nsites;
      load_matching_files();
      run_matching();

      VT_report( 1, "elapsed: %d seconds", (int)(time( NULL) - tstart));
      return 0;
   }

   //
   //  Read measurement sets from stdin, one set per line, and process each.
   //

   char *inbuf = malloc( 4096), *p, *q;

   struct MSET m;
   init_measurement_set( &m);
   while (fgets( inbuf, 4096, stdin))
   {
      if ((p = strchr( inbuf, '\r')) != NULL) *p = 0;
      if ((p = strchr( inbuf, '\n')) != NULL) *p = 0;

      reset_measurement_set( &m);
      p = inbuf;
      while (*p)
      {
         if (isspace( *p)) { p++; continue; }
         for (q = p; *q; q++) if (isspace( *q)) break;
         *q++ = 0; 
         parse_measurement( &m, p);
         p = q;
      }

      process_measurement_set( &m);
   }

   return 0;
}

