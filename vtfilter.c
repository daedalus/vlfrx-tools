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

#include <fftw3.h>

static VTFILE *vtinfile, *vtoutfile;
static char *inname = NULL;
static char *outname = NULL;
static int sample_rate = 0;
static int chans = 0;
static char **eqmap_files = NULL; // List of eqmap files, -e option
static int n_eqmap_files = 0;     // Number of eqmap files
static int DFLAG = 0;             // -d diagnostic option

#define DEFAULT_AN_BW 1
#define DEFAULT_AN_THRESH 6.0

static double an_skip = 2.0;

static double gain = 1.0;
static double *gains;

static int BINS = 8192*2;
static int FFTWID = 0;
static double dF = 0;
static int p_load1 = 0, p_load2 = 0;                       // Two load pointers

static int DEFAULT_AN_LIMIT = 6000;   // Upper limit Hz for auto notching

#define AN_SMOOTH 0.95   // Exponential moving average coefficienct for
                         // bin power averaging
#define AN_WIDTH  20     // Half width Hz of mean reference neighbourhood

fftw_complex *X;
static struct CHANNEL
{
   int unit;
   int an_enable;
   int an_bw;
   double an_thresh;
   int an_limit;
   double ms_acc;
   fftw_complex *filterc;                                // Filter coefficients

   fftw_complex *X;                                   // Working buffer for FFT
   fftw_plan ffp_fwd1, ffp_rev1;              
   fftw_plan ffp_fwd2, ffp_rev2;              

   double *buf1, *buf2;                                 // Two circular buffers
   double *acc;                      // Array of smoothed average signal levels
}
 *channels;   // One of these for each output channel

static double *sinsq;                                     // FT window function

static void init_channel( struct CHANNEL *fp, int unit)
{
   int i;

   fp->unit = unit;
   fp->buf1 = VT_malloc_zero( FFTWID * sizeof( double));

   fp->buf2 = VT_malloc_zero( FFTWID * sizeof( double));

   fp->ffp_fwd1 = fftw_plan_dft_r2c_1d( FFTWID, fp->buf1, X, FFTW_ESTIMATE);
   fp->ffp_rev1 = fftw_plan_dft_c2r_1d( FFTWID, X, fp->buf1, FFTW_ESTIMATE);
   fp->ffp_fwd2 = fftw_plan_dft_r2c_1d( FFTWID, fp->buf2, X, FFTW_ESTIMATE);
   fp->ffp_rev2 = fftw_plan_dft_c2r_1d( FFTWID, X, fp->buf2, FFTW_ESTIMATE);

   fp->acc = VT_malloc_zero( FFTWID * sizeof( double));

   fp->filterc = VT_malloc( (BINS+1) * sizeof( fftw_complex));
   for (i = 0; i <= BINS; i++) fp->filterc[i] = 1;
}

static void init_filter( void)
{
   int i;

   X = VT_malloc( sizeof( fftw_complex) * (BINS+1));

   sinsq = VT_malloc( FFTWID * sizeof( double));
   for (i = 0; i < FFTWID; i++)
   {
      double theta = M_PI * i/(double) FFTWID;
      sinsq[i] = sin( theta) * sin( theta);
   }

   gains = VT_malloc( (BINS+1) * sizeof( double));
}

static void inline smooth_accumulate_double( double *a, double val, 
                                             double factor)
{
   *a = *a * factor + val * (1-factor);
}

static void filter_inner( struct CHANNEL *fp, fftw_plan ffp_forward, 
                                              fftw_plan ffp_reverse)
{
   int i;
   double ms_in;

   //
   //  Do the forward FFT and accumulate the moving average power level in
   //  each bin.
   //

   fftw_execute( ffp_forward);

   for (ms_in = i =0; i <= BINS; i++)
   {
      double a = cabs( X[i]);
      ms_in += a * a;  // Total power
   }

   if (ms_in < an_skip * fp->ms_acc)
      for (i = 0; i <= BINS; i++)
      {
         double a = cabs( X[i]);
         smooth_accumulate_double( &fp->acc[i], a * a, AN_SMOOTH);
      }

   smooth_accumulate_double( &fp->ms_acc, ms_in, AN_SMOOTH);

   //
   //  Start with all bins at unity gain.
   //

   for (i = 0; i <= BINS; i++) gains[i] = 1;

   //
   //  If the autonotch filter is enabled, turn down the gain on any bins whose
   //  accumulated power level exceeds the threshold.
   //

   if (fp->an_enable) 
   {
      // Average power in neighbouring spectrum +/- AN_WIDTH Hz 

      int an_width = AN_WIDTH * (FFTWID / (double) sample_rate);
      for (i = an_width; i < FFTWID/2 - an_width && i < fp->an_limit; i++)
      {
         int j;
         double sum;

         for (sum = 0, j = -an_width; j <= an_width; j++)
            if (j) sum += fp->acc[i+j];
         sum /= 2*an_width;
  
         if (fp->acc[i] > sum * fp->an_thresh)
         {
            int n;
            for (n = -fp->an_bw; n <= fp->an_bw; n++) gains[i+n] = 0;
         }
      }
   }

   //
   //  Apply the gain array and filter coefficients to the array of bins.
   //

   if (!fp->an_enable)
      for (i = 0; i <= BINS; i++) X[i] *= fp->filterc[i] * gains[i] / FFTWID;
   else
   {
      int j, k;

      for (i = 1; i < BINS; i++)
      {
         if (gains[i]) continue;
        
         for (j = i + 1; j < BINS; j++) if (gains[j]) break;
         if (j < BINS)
         {
            complex double p1 = X[i-1];
            complex double p2 = X[j];
            int d = j - (i-1);
            for (k = i; k < j; k++)
               X[k] = p1 + (p2 - p1) * (k-(i-1))/(double)d;
         }
         i = j;
      }

      for (i = 0; i <= BINS; i++) X[i] *= fp->filterc[i] / FFTWID;
   }


   //
   //  Do the reverse FFT.
   //
   fftw_execute( ffp_reverse);
}

//
//  The FFT filter.  This runs two FFTs, overlapping by 50%.  Each
//  has its own circular buffer for input.

static void filter_outer( double *inframe, int *map, double *outframe)
{
   int c;

   for (c = 0; c < chans; c++)
   {
      struct CHANNEL *fp = channels + c;

      outframe[c] = gain * (fp->buf1[p_load1] +
                            fp->buf2[p_load2]);

      double f = inframe[map[c]]; 
      fp->buf1[p_load1] = f * sinsq[p_load1];
      fp->buf2[p_load2] = f * sinsq[p_load2];
   }

   if (++p_load1 == FFTWID)
   {
      p_load1 = 0;

      for (c = 0; c < chans; c++)
      {
         struct CHANNEL *fp = channels + c;
         filter_inner( fp, fp->ffp_fwd1, fp->ffp_rev1);
      }
   }

   if (++p_load2 == FFTWID)
   {
      p_load2 = 0;

      for (c = 0; c < chans; c++)
      {
         struct CHANNEL *fp = channels + c;
         filter_inner( fp, fp->ffp_fwd2, fp->ffp_rev2);
      }
   }
}

static inline complex double cpower( complex double c, int n)
{
   complex double r = 1;
   while (n-- > 0) r *= c;
   return r;
}

static void setup_hpf( struct CHANNEL *fp, int npoles, double corner)
{
   int i;
   if (npoles <= 0) VT_bailout( "filter must specify number of poles");

   for (i = 0; i <= BINS; i++)
   {
      double F = i * dF/corner;   // Normalised frequency
      fp->filterc[i] *= cpower( I*F/(1 + I*F), npoles);
   }

   VT_report( 1, "ch %d hpf: poles=%d corner=%.2f",
                     fp->unit, npoles, corner);
}

static void setup_lpf( struct CHANNEL *fp, int npoles, double corner)
{
   int i;
   if (npoles <= 0) VT_bailout( "filter must specify number of poles");

   for (i = 0; i <= BINS; i++)
   {
      double F = i * dF/corner;   // Normalised frequency
      fp->filterc[i] *= cpower( 1 / (1 + I*F), npoles);
   }

   VT_report( 1, "ch %d lpf: poles=%d corner=%.2f",
                     fp->unit, npoles, corner);
}

static void setup_bpf( struct CHANNEL *fp, double center, double width)
{
   int i;

   for (i = 0; i <= BINS; i++)
   {
      double freq = (i+0.5) * (sample_rate / (double) FFTWID);
      if (freq < center - width/2 ||
          freq > center + width/2)
          {
             fp->filterc[i] = 0;
          }
   }

   VT_report( 1, "ch %d bpf: center=%.2f width=%.2f",
                     fp->unit, center, width);
}

static void setup_bsf( struct CHANNEL *fp, double center, double width)
{
   int i;

   for (i = 0; i <= BINS; i++)
   {
      double freq = (i+0.5) * (sample_rate / (double) FFTWID);
      if (freq >= center - width/2 &&
          freq <= center + width/2)
          {
             fp->filterc[i] = 0;
          }
   }

   VT_report( 1, "ch %d bsf: center=%.2f width=%.2f",
                     fp->unit, center, width);
}

static void parse_filter_args( char *args, struct VT_CHANSPEC *spec)
{
   int type = 0;
   double freq = -1;
   double width = 0;
   int poles = 0;
   int an_bw = DEFAULT_AN_BW;
   double an_thresh = DEFAULT_AN_THRESH;
   double an_limit = DEFAULT_AN_LIMIT;

   while (args && *args)
   {
      char *p = strchr( args, ',');
      if (p) p++;

      if (!strncmp( args, "lp", 2)) type = 1;
      else
      if (!strncmp( args, "hp", 2)) type = 2;
      else
      if (!strncmp( args, "bp", 2)) type = 3;
      else
      if (!strncmp( args, "bs", 2)) type = 4;
      else
      if (!strncmp( args, "an", 2)) type = 5;
      else
      if (!strncmp( args, "f=", 2)) freq = atof( args+2);
      else
      if (!strncmp( args, "w=", 2)) width = atof( args+2);
      else
      if (!strncmp( args, "poles=", 6)) poles = atoi( args+6);
      else
      if (!strncmp( args, "p=", 2)) poles = atoi( args+6);
      else
      if (!strncmp( args, "bw=", 3)) an_bw = atoi( args+3);
      else
      if (!strncmp( args, "th=", 3)) an_thresh = atof( args+3);
      else
      if (!strncmp( args, "ul=",3)) an_limit = atof( args+3);
      else
         VT_bailout( "unrecognised filter option: %s", args);

      args = p;
   }

   if (!type) VT_bailout( "filter must specify an, hp, lp, bp or bs");

   if (type != 5)
      if (freq < 0 || freq > sample_rate/2)
         VT_bailout( "invalid or missing frequency for filter");

   //
   //  Apply this filter to the channels given by 'spec'.  If spec is NULL
   //  then apply to all channels of the output stream.
   //

   if (!spec)
   {
      spec = VT_parse_chanspec( args);
      VT_init_chanspec( spec, vtoutfile);
   }

   int i;
   for (i = 0; i < spec->n; i++)
   {
      int c = spec->map[i];
      struct CHANNEL *fp = channels + c;
      if (type == 1) setup_lpf( fp, poles, freq);
      if (type == 2) setup_hpf( fp, poles, freq);
      if (type == 3) setup_bpf( fp, freq, width);
      if (type == 4) setup_bsf( fp, freq, width);
      if (type == 5)
      {
         fp->an_bw = an_bw;
         fp->an_thresh = an_thresh;
         fp->an_enable = TRUE;
         fp->an_limit = an_limit * (FFTWID / (double) sample_rate);
      }
   }
}

static struct EQMAP
{
   double freq;
   complex double *coeffs;
}
 *eqmap = NULL;

static int eqmap_alloc = 0;
static int eqmap_n = 0;
static double eqmap_freq = 0;

static void load_eqmap( char *file)
{
   eqmap_freq = eqmap_n = 0;

   int i, nr = 0;
   FILE *fh = fopen( file, "r");
   if (!fh) VT_bailout( "cannot open eqmap [%s], %s",
                           file, strerror( errno));

   int read_value( double *v)
   {
      int e = fscanf( fh, " %lf", v);
      if (e == EOF) return 0;
      if (e != 1) VT_bailout( "error in eqmap [%s] line %d", file, nr+1);
      return 1;
   }

   int read_string( char *s)
   {
      int e = fscanf( fh, " %s", s);
      if (e == EOF) return 0;
      if (e != 1) VT_bailout( "error in eqmap [%s] line %d", file, nr+1);
      return 1;
   }

   struct EQMAP *allocate( double freq)
   {
      if (freq < eqmap_freq)
         VT_bailout( "non-monotonic frequency in eqmap %s: %.3e", file, freq);

      if (eqmap_alloc <= eqmap_n)
      {
         eqmap = VT_realloc( eqmap, sizeof( struct EQMAP) * ++eqmap_alloc);
         eqmap[eqmap_n].coeffs = VT_malloc( sizeof( complex double) * chans);
      }

      eqmap_freq = eqmap[eqmap_n].freq = freq;
      return eqmap + eqmap_n++;
   }

   double freq;
   struct EQMAP *ep;

   while (1)
   {
      if (!read_value( &freq)) break;   // Read first field - frequency Hz
      
      if (!eqmap_n && freq > 0)
      {
         // Jam in a DC entry if not supplied from the file
         ep = allocate( 0);
         for (i = 0; i < chans; i++) ep->coeffs[i] = 1;
      }

      ep = allocate( freq);

      for (i = 0; i < chans; i++)
      {
         char temp[50];

         if (!read_string( temp))
            VT_bailout( "incomplete data in eqmap %s line %d", file, nr+1);
         if (!VT_parse_complex( temp, ep->coeffs + i))
            VT_bailout( "bad coefficient [%s] in eqmap %s line %d",
                   temp, file, nr+1);
      }

      nr++;
   }

   fclose( fh);

   if (!nr) VT_bailout( "no records parsed in eqmap %s", file);

   VT_report( 1, "loaded %d eqmap records from %s", nr, file);

   // Append a record for the Nyquist frequency if one wasn't supplied
   if (eqmap_freq < sample_rate / 2.0)
   {
      ep = allocate( sample_rate / 2.0);
      for (i = 0; i < chans; i++) ep->coeffs[i] = 1;
   }
}

static void apply_eqmap( void)
{
   int bin, ne, ch;

   for (ne = bin = 0; bin <= BINS;)
   {
      struct EQMAP *base = eqmap + ne;
      struct EQMAP *next = eqmap + ne + 1;
      double f = bin * dF;
      if (f > next->freq) { ne++; continue; }

      double r = ne == 0 ? 0 : 
                 log( f / base->freq) / log( next->freq / base->freq);

      for (ch = 0; ch < chans; ch++)
      {
         // Interpolate magnitudes logarithmically, the 1e-99 prevents NaN
         // if a zero coefficient is given
         double a1 = fmax( cabs( base->coeffs[ch]), 1e-99);
         double a2 = fmax( cabs( next->coeffs[ch]), 1e-99);
         double ma = exp( log( a1) + (log(a2) - log(a1)) * r);

         // Interpolate phase linearly
         double p1 = carg( base->coeffs[ch]);
         double p2 = carg( next->coeffs[ch]);
         double dp = p2 - p1;
         if (dp >= M_PI) dp -= 2 * M_PI;
         if (dp < -M_PI) dp += 2 * M_PI;
         double pa = p1 + dp * r;

         channels[ch].filterc[bin] *= ma * (cos(pa) + I*sin(pa));
      }

      bin++;
   }
}

//
//  Output filter coefficients.  -d1 option outputs in eqmap format, -d2
//  outputs amplitude and phase in separate columns.
//

static void output_filterc( void)
{
   int bin, ch;
   
   for (bin = 0; bin <= BINS; bin++)
   {
      printf( "%.6f", bin * dF);

      for (ch = 0; ch < chans; ch++)
      {
         struct CHANNEL *cp = channels + ch;

         if (DFLAG == 1)
            printf( " %.6e%+.6ej", creal( cp->filterc[bin]),
                                   cimag( cp->filterc[bin]));
         else
         if (DFLAG == 2)
            printf( " %.6e %.6e", cabs( cp->filterc[bin]),
                                  carg( cp->filterc[bin]) * 180/M_PI);
      }

      printf( "\n");
   }
}

static void usage( void)
{
   fprintf( stderr,
       "usage:  vtfilter [options] input output\n"
       "\n"
       "options:\n"
       "  -v        Increase verbosity\n"
       "  -B        Run in background\n"
       "  -L name   Specify logfile\n"
       "  -g gain   Output gain\n"
       "  -e eqmap  Apply EQ table\n"
       "  -n bins   Frequency bins (default 16384)\n"
       "  -a autonotch_options\n"
       "     th=threshold,bw=notchwidth,ul=upper limit\n"
       "  -h filter_options\n"
       "     Butterworth low pass: lp,f=corner,p=poles\n"
       "     Butterworth high pass: hp,f=corner,p=poles\n"
       "     Brick wall bandpass: bp,f=center,w=width\n"
       "     Brick wall bandstop: bs,f=center,w=width\n"
       "     Autonotch: an,th=threshold,bw=notchwidth,ul=upper limit\n"
     );
   exit( 1);
}

int main(int argc, char *argv[])
{
   VT_init( "vtfilter");

   int i;
   int background = 0;

   // Options for -a and -h are just saved for now, scanned
   // later when we know how many channels we're dealing with

   char **filter_args = NULL;
   struct VT_CHANSPEC **filter_specs = NULL;
   int filter_n = 0;
   struct VT_CHANSPEC *fspec = NULL;

   void add_filter_spec( char *args)
   {
      filter_args = VT_realloc( filter_args, (filter_n + 1) * sizeof( char *));
      filter_args[filter_n] = strdup( args);
   
      filter_specs = VT_realloc( filter_specs,
                             (filter_n + 1) * sizeof( struct VT_CHANSPEC *));
      filter_specs[filter_n] = fspec;
      filter_n++;
   }

   while (1)
   {
      int c = getopt( argc, argv, "vBg:a:h:e:d:n:c:L:?");

      if (c == 'v') VT_up_loglevel();
      else
      if (c == 'B') background = 1;
      else
      if (c == 'L') VT_set_logfile( "%s", optarg);
      else
      if (c == 'g') gain = atof( optarg);
      else
      if (c == 'e')
      {
         eqmap_files = VT_realloc( eqmap_files,
                                   sizeof( char *) * (n_eqmap_files + 1));
         eqmap_files[n_eqmap_files++] = strdup( optarg);
      }
      else
      if (c == 'd') DFLAG = atoi( optarg);
      else
      if (c == 'n') BINS = atoi( optarg);
      else
      if (c == 'a')
      {
         // -a ... is equivalent to -h an,...
         char *s;
         asprintf( &s, "an,%s", optarg); add_filter_spec( s);
      }
      else
      if (c == 'h')
      {
         add_filter_spec( optarg);
      }
      else
      if (c == 'c')
      {
         char *s;
         asprintf( &s, ":%s", optarg); fspec = VT_parse_chanspec( s);
      }
      else
      if (c == -1) break;
      else
         usage();
   }

   if (optind + 2 == argc)
   {
      inname = strdup( argv[optind]);
      outname = strdup( argv[optind+1]);
   }
   else
   if (optind + 1 == argc)
   {
      inname = strdup( argv[optind]);
      outname = strdup( "-");
   }
   else
   if (optind == argc)
   {
      inname = strdup( "-");
      outname = strdup( "-");
   }
   else usage();

   if (background)
   {
      int flags = inname[0] == '-' ? KEEP_STDIN : 0;
      if (outname[0] == '-') flags |= KEEP_STDOUT;
      VT_daemonise( flags);
   }

   struct VT_CHANSPEC *chspec = VT_parse_chanspec( inname);
   vtinfile = VT_open_input( inname);
   if (!vtinfile) 
      VT_bailout( "cannot open input %s: %s", inname, VT_error);

   sample_rate = VT_get_sample_rate( vtinfile);
   FFTWID = 2 * BINS;
   dF = sample_rate/(double) FFTWID;

   VT_init_chanspec( chspec, vtinfile);
   chans = chspec->n;
   VT_report( 1, "channels: %d, sample_rate: %d resolution: %.3f",
                     chans, sample_rate, dF);
   VT_report( 2, "fft width %d, bins %d", FFTWID, BINS);

   vtoutfile = VT_open_output( outname, chans, 0, sample_rate);
   if (!vtoutfile) VT_bailout( "cannot open: %s", VT_error);

   // 
   //  Setup filtering.
   //

   init_filter();

   channels = VT_malloc_zero( chans * sizeof( struct CHANNEL));
   for (i = 0; i < chans; i++) init_channel( channels+i, i+1);

   for (i = 0; i < n_eqmap_files; i++)
   {
      load_eqmap( eqmap_files[i]);
      apply_eqmap();
   }

   for (i = 0; i < filter_n; i++)
      parse_filter_args( filter_args[i], filter_specs[i]);

   for (i = 0; i < chans; i++)
   {
      struct CHANNEL *fp = channels + i;

      if (fp->an_enable)
         VT_report( 1, "channel: %d autonotch bw=%d th=%.2f",
                            i + 1, fp->an_bw, fp->an_thresh);
   }
 
   if (DFLAG == 1 ||
       DFLAG == 2)  // -d1 or -d2 option: Dump the filter transform and exit
   {
      output_filterc(); 
      VT_exit( "done -d");
   }

   //
   //  Main loop.
   //

   double *inframe;
   double *outframe = VT_malloc( sizeof( double) * chans);

   double srcal = 1.0;
   int n = 0, e;
   p_load1 = 0;
   p_load2 = FFTWID/2;

   timestamp T = timestamp_ZERO;

   while (1)
   {
      e = VT_is_block( vtinfile);
      if (e < 0)
      {
         VT_report( 1, "end of input");
         break;
      }

      if (e)  // Starting a new input block?
         srcal = VT_get_srcal( vtinfile);

      if (!vtoutfile->nfb)  // Starting a new output block?
      {
         T = VT_get_timestamp( vtinfile);
         double offset = FFTWID/(srcal * sample_rate);
         VT_set_timebase( vtoutfile, timestamp_add( T, -offset), srcal);
      }

      inframe = VT_get_frame( vtinfile);

      filter_outer( inframe, chspec->map, outframe);

      // Discard the first FFTWID frames as these are just leading zeros.
      if (n < FFTWID) n++;
      else
         VT_insert_frame( vtoutfile, outframe);
   }

   // An empty frame, a dummy for use in flushing
   inframe = VT_malloc_zero( VT_get_chans( vtinfile) * sizeof( double));

   //
   //  If necessary, finish discarding the first FFTWID frames. This kicks in
   //  if the input stream was shorter than FFTWID samples.
   //
   //  n is normally equal to FFTWID, but if the input stream was shorter, 
   //  n will give the length.
   //

   int u;
   for (u = n; u < FFTWID; u++)
      filter_outer( inframe, chspec->map, outframe);

   if (n < FFTWID)
   {
      // First timestamp of output stream will have been set incorrectly, on
      // the assumption that more than FFTWID samples will come through.
      // Now set the correct start time because we know the input length
      double offset = (n-1)/(srcal * sample_rate);
      VT_set_timebase( vtoutfile, timestamp_add( T, -offset), srcal);
   }

   // Drain the filter buffer.   Another min(n,FFTWID) frames to output.
   for (u = 0; u < n; u++)
   {
      filter_outer( inframe, chspec->map, outframe);
      VT_insert_frame( vtoutfile, outframe);
   }

   VT_release( vtoutfile);

   return 0;
}

