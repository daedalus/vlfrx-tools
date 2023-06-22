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

// Important times, all in seconds
static double Tbuf = 0.002;                // Buffer length
static double Tpulse = 0.0012;             // Pulse length analysed
static double Tlead = 200e-6;              // 
static double Tpeaks = 700e-6;             // 
static double Tfft = 0.014;                // FFT input width
static double Tpretrig =  500e-6;          // Pre-trigger period

// Set by command line options
static int EFLAG_B = 0;           // Non-zero if -eB given, multi-band data
static int EFLAG_S = 0;           // Non-zero if -eS given, output spectrum
static int EFLAG_T = 0;           // Non-zero if -eT given, output time domain
static int EFLAG_P = 0;           // Non-zero if -eP given, output pulses
static int EFLAG_Q = 0;           // Non-zero if -eQ given, output quality data
static int EFLAG_R = 0;           // Non-zero if -eR given, output range
static int EFLAG_X = 0;           // Non-zero if -eX given, experimental

static int CFLAG = 0;                 // Non-zero if -c given, calibration mode
static int Nmax = 0;             // From -N option: number of pulses to capture
static double throttle = 0;      // From -r option: target rate triggers/second
static double Etime = 0;            // From -E option, stop after Etime seconds
static int gran = 3600;              // From -G option: output file granularity
static char *outdir = NULL;                 // From -d option: output directory
static timestamp Ttrig;             // From -T option: single shot trigger time
static double limit_u = 0;       // From -u, upper frequency limit to S records

// Variables for auto threshold
static time_t last_T = 0;               // Start time for trigger rate counting
static int last_N = 0;                       // Number of triggers since last_T
#define AT_INT 60                   // Seconds between updates of trigger level

// Limits for accepting a measurement
static double limit_ascore = 0;    // Set by option -f ascore=
static double limit_pscore = 0;    // Set by option -f pscore=
static double limit_tscore = 0;    // Set by option -f tscore=

static double trigger_amplitude = 0;      // Amplitude trigger level, set by -a
static double trigger_impulse = 0;        // Impulse ratio trigger, set by -i

static timestamp Tstart = timestamp_ZERO;          // Timestamp of first sample 
static timestamp T_D = timestamp_ZERO;

static int sample_rate = 0;

static int buflen = 0;                       // Circular buffer length, samples
static int btrig = 0;                            // Pre-trigger period, samples
static int bp = 0;                         // Base pointer into circular buffer
static int FFTWID = 0;                                    // FFT width, samples
static int BINS = 0;                                          // Number of bins
static double DF = 0;                                   // Frequency resolution
static double DT = 0;                                        // Time resolution
static timestamp Tin;                        // Timestamp of most recent sample
static int one_shot = 0;                                // Non-zero if -T given
static int N_out = 0;                      // Number of sferics processed
static double *irbuf = 0;

// Variables for polar operation
static int polar_mode = 0;
static double polar1_align = 0;        // Azimuth of 1st H-loop channel
static double polar2_align = 0;        // Azimuth of 2nd H-loop channel

static int ch_EFIELD = -1;
static int ch_HFIELD1 = -1;
static int ch_HFIELD2 = -1;

#define BUFIDX(P) (((P) + bp + buflen) % buflen)
#define BUF(C,P) channels[C].buf[BUFIDX(P)]

#define EIC_CUTOFF 1670
#define EARTH_RAD 6371.0
#define C 299.792e3

// Rate counters
#define RC_INT 100             // Output interval (seconds) of rate statistics
static time_t rc_start = 0;    // Start time of counters
static int rc_raw = 0;         // Count of raw triggers
static int rc_ascore = 0;      // Number of sferics dropped due to low ascore
static int rc_pscore = 0;      // Number of sferics dropped due to low pscore
static int rc_tscore = 0;      // Number of sferics dropped due to low tscore
static int rc_out = 0;         // Count of output records

static inline void rc_review( timestamp T)
{
   if (timestamp_secs( T) < rc_start + RC_INT) return;

   if (rc_start)
      VT_report( 1,
                 "RC raw %d ascore %d pscore %d tscore %d out %d",
                 rc_raw, rc_ascore, rc_pscore, rc_tscore, rc_out);
 
   rc_start = timestamp_secs( T);
   rc_raw = rc_out = 0;
   rc_ascore = rc_pscore = rc_tscore = 0;
}

static void usage( void)
{
   fprintf( stderr,
       "usage:  vttoga [options] [input]\n"
       "\n"
       "options:\n"
       "  -v            Increase verbosity\n"
       "  -B            Run in background\n"
       "  -L name       Specify logfile\n"
       "  -p polarspec  Specify orientation of input channels\n"
       "  -N count      Examine this many sferics then exit\n"
       "  -E seconds    Stop after so many seconds\n"
       "\nFrequency bands\n"
       "  -F start,end  Frequency range (default 4000,17000)\n"
       "  -M start,end  Mask frequency range (no default)\n"
       "\nTrigger options\n"
       "  -T timestamp  One-shot trigger time\n"
       "  -c            Calibration mode\n"
       "  -a thresh     Amplitude trigger threshold (no default)\n"
       "  -i thresh     Impulse ratio trigger threshold (no default)\n"
       "  -r rate       Auto threshold to this rate/second\n"
       "\nQuality options (experimental)\n"
       "  -f ascore=val Value 0 to 1, spectrum quality limit (default 0)\n"
       "  -f pscore=val Value 0 to 1, phase slope quality limit (default 0)\n"
       "  -f tscore=val Value 0 to 1, TOGA quality limit (default 0)\n"
       "\nOutput options\n"
       "  -e opts       Extended output, opts:\n"
       "                   B   Output multi-band measurements\n"
       "                   S   Output spectrum\n"
       "                   T   Output time domain\n"
       "                   P   Output pulse measurements\n"
       "                   Q   Output quality test results (experimental)\n"
       "  -u freq       Upper frequency limit to -eS output data\n"
       "  -d outdir     Output directory (defaults to stdout)\n"
       "  -G seconds    Output file granularity (default 3600)\n"
    );

   exit( 1);
}

static struct CHAN
{
   double *buf;                                        // Circular input buffer
   double *fft;                                             // FFT input buffer
   complex double *X;                                             // FFT output
   double rms; 
   fftw_plan fp;
} *channels = NULL;

static int nchans = 0;

//
//  Bands in which to measure TOGA, from -F options.
//

static struct BAND {

   // Setup values
   double F1, F2;             // Frequency range to analyse
   int bin1, bin2;            // Start and finish bin numbers, from F1 and F2
   double gcf;                // Geometric center frequency
   int noscore;               // TRUE if this band is not used for scoring
   int do_range;              // TRUE if this band is used for range estimate

   // Measured values
   double energy;
   double phresidual;         // RMS phase residual from linear phase, radians
   double toga_offset;        // Seconds between Tf and toga

   double te;

} *bands = NULL;

static int nbands = 0;

//
//  Bands which are masked, from -M options.
//

static struct MASK {

   double F1, F2;             // Frequency range to mask
   int bin1, bin2;            // Start and finish bin numbers, from F1 and F2
   double gcf;                // Geometric center frequency

} *masks = NULL;

static int nmasks = 0;

static char *binmask = NULL;  // Array of bins to mask
static int *unmasked = NULL;  // List of bins which are not masked
static int nunmasked = 0;     // Number of unmasked bins

static void reset_bands( void)
{
   int band;
   for (band = 0; band < nbands; band++)
   {
      struct BAND *b = bands + band;
      b->energy = b->phresidual = b->toga_offset = 0;
   }
}

#define MAXTDPEAKS 6

struct EVENT {
   timestamp Tf;                   // Timestamp of first sample of buffer
   double energy;                  // Total energy
   double bearing;                 // If in polar mode
   double range_day1;               // Range estimate, km, daytime
   double range_night1;             // Range estimate, km, nighttime
   double range_day2;               // Range estimate, km, daytime
   double range_night2;             // Range estimate, km, nighttime
   double *upb;                    // Unwrapped phase, radians
   double *tsp;                    // Phase at the overall TOGA, radians
   complex double *ca;             // Complex analytic time domain
   double ascore;
   double pscore;
   double tscore;
   double xscore;
   double impulse_ratio;

   struct TDPEAK {
      double a;       // Amplitude of peak
      double t;       // Time offset from Tf
      double f;       // Instantaneous frequency at the peak
   } tdpeaks[MAXTDPEAKS];
   int ntdpeaks;

   double peaks_ratio;
};

static complex double *XX;
static complex double *caX;
static double * car;
static fftw_plan ffca;

static void init_event( struct EVENT *ev)
{
   memset( ev, 0, sizeof( struct EVENT));
   ev->upb = VT_malloc_zero( sizeof( double) * BINS);
   ev->tsp = VT_malloc_zero( sizeof( double) * BINS);
   ev->ca = VT_malloc( sizeof( complex double) * FFTWID);
}

static void reset_event( struct EVENT *ev)
{
   ev->energy = ev->bearing = ev->impulse_ratio =
      ev->ascore = ev->pscore = ev->tscore = 0;

   ev->ntdpeaks = ev->peaks_ratio = 0;
}

///////////////////////////////////////////////////////////////////////////////
//  Output Functions                                                         //
///////////////////////////////////////////////////////////////////////////////

//
//  Output a 'H' record - timestamp and total amplitude.
//  The timestamp in the H record is the earliest TOGA from all the bands.
//

static void output_header( FILE *f, struct EVENT *ev)
{
   double minimum_toga_offset = 1e99;
   int band;
   for (band = 0; band < nbands; band++)
   {
      struct BAND *b= bands + band;
      if (b->toga_offset < minimum_toga_offset)
         minimum_toga_offset = b->toga_offset;
   }

   char temp[30];
   timestamp_string6( timestamp_add( ev->Tf, minimum_toga_offset), temp);

   fprintf( f, "H %s %.3e %.2f",
               temp,               // Timestamp, 1uS resolution
               sqrt( ev->energy),
               ev->impulse_ratio);

   if (polar_mode) fprintf( f, " %.1f", ev->bearing * 180/M_PI);
   fprintf( f, "\n");
}

//
// Extended output records if -e is given.  
//    H is the header record, call output_header() for that;
//    B records for multi-band TOGAs
//    S records for spectrum data;
//    T records for time domain;
//    Q records for quality test results (experimental);
//    R records for range estimate;
//    E is an end marking record;
//

static void output_extended_B( FILE *f, struct EVENT *ev)
{
   char temp[30];

   timestamp_string6( ev->Tf, temp);
   fprintf( f, "B %s %.3e",
               temp,               // Timestamp, 1uS resolution
               ev->energy);

   int band;
   for (band = 0; band < nbands; band++)
   {
      struct BAND *b= bands + band;

      fprintf( f, " %.0f %.0f %.6f %.3e %.1f",
               b->F1,
               b->F2,
               b->toga_offset,
               b->energy,
               b->phresidual * 180/M_PI);
   }
   fprintf( f, "\n");
}

static void output_extended_Q( FILE *f, struct EVENT *ev)
{
   fprintf( f, "Q %.2f %.2f %.2f\n",
               ev->ascore,
               ev->pscore,
               ev->tscore);
}

static void output_extended_S( FILE *f, struct EVENT *ev)
{
   int j;

   int maxbin = limit_u ? limit_u / DF : BINS;
   if (maxbin > BINS) maxbin = BINS;

   for (j = 0; j < maxbin; j++)
   {
      double F = j * DF;
   
      fprintf( f, "S %.1f %d %.3e %.2f\n",
           F,                     // Frequency, Hz
           binmask[j],            // 1 if masked
           cabs(XX[j]),           // Bin amplitude
           ev->tsp[j] * 180/M_PI  // Unwrapped phase, degrees
           );
   }
}

static void output_extended_T( FILE *f, struct EVENT *ev)
{
   int ch, j;

   for (j = 0; j < Tbuf * sample_rate; j++)
   {
      fprintf( f, "T %.6e", j * DT);
      for (ch = 0; ch < nchans; ch++) fprintf( f, " %.4e", BUF(ch, j));

      fprintf( f, " %.4e", cabs( ev->ca[j]));

      double dp = 0, F = 0;
      if (j > 0)
      {
         dp = carg( ev->ca[j]) - carg( ev->ca[j-1]);
         while (dp > M_PI) dp -= 2 * M_PI;
         while (dp < -M_PI) dp += 2 * M_PI;
         F = dp/(DT * 2 * M_PI);
      }
      fprintf( f, " %.4e", F);
      fprintf( f, " %.2f", carg( ev->ca[j])/(2*M_PI));
      fprintf( f, "\n");
   }
}

static void output_extended_P( FILE *fo, struct EVENT *ev)
{
   int j;

   fprintf( fo, "P %.3f", ev->peaks_ratio);
   for (j = 0; j < ev->ntdpeaks; j++)
      fprintf( fo, " %.6e %.3e %.0f",
           ev->tdpeaks[j].t, ev->tdpeaks[j].a, ev->tdpeaks[j].f);
   fprintf( fo, "\n");
}

static void output_extended_R( FILE *fo, struct EVENT *ev)
{
   fprintf( fo, "R %.1f %.1f %.1f %.1f\n",
           ev->range_day1, ev->range_night1,
           ev->range_day2, ev->range_night2);
}

static void output_extended_X( FILE *fo, struct EVENT *ev) // XXX
{
   fprintf( fo, "X %.3f\n", ev->xscore);
}

static void output_record( struct EVENT *ev)
{
   static uint64_t count = 0;

   //
   //  Open an output file if -d given, otherwise use stdout.
   //

   FILE *f;

   if (!outdir) f = stdout;
   else
   {
      char *filename;
      time_t secs = timestamp_secs( timestamp_add( Tin, -buflen * DT));
      secs = (secs / gran) * gran;
      struct tm *tm = gmtime( &secs);

      if (asprintf( &filename, "%s/%02d%02d%02d-%02d%02d%02d",
            outdir,
            tm->tm_year % 100, tm->tm_mon + 1, tm->tm_mday,
            tm->tm_hour, tm->tm_min, tm->tm_sec) < 0)
         VT_bailout( "out of memory");
      if ((f = fopen( filename, "a")) == NULL)
         VT_bailout( "cannot open %s: %s", filename, strerror( errno));
      free( filename);
   }

   output_header( f, ev);

   int ext = FALSE;
   if (EFLAG_B && count % EFLAG_B == 0)
   {
      output_extended_B( f, ev);
      ext = TRUE;
   }
   if (EFLAG_S && count % EFLAG_S == 0)
   {
      output_extended_S( f, ev);
      ext = TRUE;
   }
   if (EFLAG_T && count % EFLAG_T == 0)
   {
      output_extended_T( f, ev);
      ext = TRUE;
   }
   if (EFLAG_P && count % EFLAG_P == 0)
   {
      output_extended_P( f, ev);
      ext = TRUE;
   }
   if (EFLAG_Q && count % EFLAG_Q == 0)
   {
      output_extended_Q( f, ev);
      ext = TRUE;
   }
   if (EFLAG_R && count % EFLAG_R == 0)
   {
      output_extended_R( f, ev);
      ext = TRUE;
   }
   if (EFLAG_X && count % EFLAG_X == 0)
   {
      output_extended_X( f, ev);
      ext = TRUE;
   }

   // If any extended outputs, terminate with an E record
   if (ext) fprintf( f, "E\n");

   // Close the output file, or flush stdout if going to a terminal
   if (f != stdout) fclose( f);
   else
   if (isatty( fileno( f))) fflush( f);

   count++;
}

///////////////////////////////////////////////////////////////////////////////
//  Sferic Analyser                                                          //
///////////////////////////////////////////////////////////////////////////////

//
//  Amplitude spectrum score.
//

static void evaluate_ascore( struct EVENT *ev)
{
   // Span the range of a typical sferic spectrum
   int minbin = 5000/DF;
   int maxbin = 15000/DF;

   // Which unmasked bin has the highest magnitude in the sferic band?
   double pkmax = 0;    // Magnitude of peak bin
   int bin;
   for (bin = minbin; bin < maxbin; bin++)
   {
      if (binmask[bin]) continue;

      double a = cabs( XX[bin]);
      if (a > pkmax) pkmax = a;
   }

   double last = 0;
   double h = 0;
   for (bin = minbin; bin < maxbin; bin++)
   {
      if (binmask[bin]) continue;

      double a = cabs( XX[bin]);
      h += fabs( a - last);
      last = a;
   }
   h += last;

   ev->ascore = h > 0 ? pkmax * 2/h : 1e-5;
}

//
//  Phase spectrum score.
//  Combine the phase residuals of each band, weighted by the RMS amplitude
//  of the band, into an overall score.
//

static void evaluate_pscore( struct EVENT *ev)
{
   double phsum = 0;
   double energy = 0;

   int band;
   for (band = 0; band < nbands; band++)
   {
      struct BAND *b = bands + band;
      if (!b->noscore)
      {
         phsum += b->phresidual * b->energy;
         energy += b->energy;
      }
   }

   double a = phsum/energy;   // Weighted average phase residual
   ev->pscore = (2 * M_PI - a)/(2 * M_PI);
}

//
//  Timing score.
//

static void evaluate_tscore( struct EVENT *ev)
{
   int band;

   double energy_sum = 0;
   for (band = 0; band < nbands - 1; band++)
   {
      struct BAND *b = bands + band;
      if (!b->noscore) energy_sum += b->energy;
   }

   ev->tscore = 1.0;

   int band1;
   for (band1 = 0; band1 < nbands - 1; band1++)
   {
      struct BAND *b1 = bands + band1;
      if (b1->noscore) continue;

      int band2;
      for (band2 = band1 + 1; band2 < nbands; band2++) 
         if (!bands[band2].noscore) break;
      if (band2 == nbands) break;

      struct BAND *b2 = bands + band2;

      //  For a good sferic, lower frequency band TOGAs should be later than
      //  higher bands.   So td will be positive or zero on a good sferic.
      double td = b1->toga_offset - b2->toga_offset;
      if (td < -30e-6) //  Bad interval?
      {
         double emin = fmin( b1->energy, b2->energy);

         ev->tscore *= (energy_sum - emin)/energy_sum;
      }
   }
}

//
//  Calculate analytic time domain waveform using a Hilbert transform.
//

static void make_analytic( struct EVENT *ev)
{
   //
   //  Copy XX to caX, multiply by -j
   //

   int b1 = Tlead * sample_rate;
   int b2 = (Tlead + Tpulse) * sample_rate;

   int bin, i;
   int minbin = 0/DF;
   int maxbin = BINS;

   minbin = 3000/DF;
   maxbin = sample_rate/2/DF;
   memset( caX, 0, sizeof( complex double) * BINS);
   for (bin = minbin; bin < maxbin; bin++)
      caX[bin] = binmask[bin] ? 0 : XX[bin] * -I;
   fftw_execute( ffca);  // Reverse FFT

   //
   //  Copy into the imaginary part of ev->ca
   //

   for (i = 0; i < b1; i++) ev->ca[i] = 0;
   for (i = b1; i < b2; i++) ev->ca[i] = I * car[i] / FFTWID;
   for (i = b2; i < FFTWID; i++) ev->ca[i] = 0;

   //
   //  Copy XX again to caX
   //

   memset( caX, 0, sizeof( complex double) * BINS);
   for (bin = minbin; bin < maxbin; bin++)
      caX[bin] = binmask[bin] ? 0 : XX[bin];
   fftw_execute( ffca);  // Reverse FFT

   //
   //  Copy into the real part of ev->ca
   //

   for (i = b1; i < b2; i++) ev->ca[i] += car[i] / FFTWID;
}

//
//  Estimate range for the band.
//

static double vf( double f1, double f2, double fc)
{
   double wc = fc * 2 * M_PI;
   double w1 = f1 * 2 * M_PI;
   double w2 = f2 * 2 * M_PI;

   double w1_2 = w1*w1;
   double w1_3 = w1*w1*w1;

   double w2_2 = w2*w2;
   double w2_3 = w2*w2*w2;

   double wc_2 = wc*wc;

   double F1 = log( (1+sqrt(1-wc_2/w1_2))*(1 - sqrt(1-wc_2/w2_2)) /
                    ((1+sqrt(1-wc_2/w2_2))*(1 - sqrt(1-wc_2/w1_2))) );

   double G1 = sqrt((w1-wc)*(wc+w1));
   double G2 = sqrt((w2-wc)*(wc+w2));

   double H1 = 8*wc_2 + 6*w1*w2  - 2*w1_2;
   double H2 = 8*wc_2 + 6*w1*w2  - 2*w2_2;

   double vf = -2 * (w2_3-3*w1*w2_2+3*w1_2*w2-w1_3) /
              (3*(w2+w1)*wc_2*F1 + G2*H2 - G1*H1);

//   double wz = wc/sqrt(1 - vf*vf);
//   printf( "least squares %.6f %.2f\n", vf, wz/(2 * M_PI));

   return vf;
}

static void compute_range( struct EVENT *ev, struct BAND *b)
{
   int i, j, k;

   //
   //  Rotate the unwrapped phase ev->upb[], time shifting to the TOGA.
   //

   double tsp[BINS];
   for (j = b->bin1; j <= b->bin2; j++)
      tsp[j] = ev->upb[j] + j * 2 * M_PI * DF * b->toga_offset;

   //
   //  amplitude of each bin.
   //
   double xa[BINS];
   for (j = b->bin1; j <= b->bin2; j++) xa[j] = cabs( XX[j]);

#if 1 // XXX
FILE *ft = fopen( "/tmp/ft.dat", "w");
   for (j = b->bin1; j <= b->bin2; j++)
      fprintf( ft, "%.1f %.3e %.3e %.3e\n", j*DF, ev->upb[j], tsp[j], xa[j]);
fclose( ft);
#endif

   //
   //  Make a list of un-masked bins.
   //

   int list[BINS];
   int nlist = 0;
   list[0] = 0;
   for (j = b->bin1; j <= b->bin2; j++)
      if (!binmask[j]) list[nlist++] = j;
   if (!nlist) return;

   //
   //  Phase slope at the TOGA, d(phase)/d(omega).
   //

   double dphi[BINS];
   memset( dphi, 0, sizeof( double) * BINS);
   double dw = 2 * M_PI * DF;
   for (i = 0; i < nlist; i++)
   {
      if (i < nlist-1)
      {
         j = list[i];  k = list[i+1];
         dphi[j] = (tsp[k] - tsp[j])/(dw * (k - j));
      }
      else
      if (i > 0)
      {
         j = list[i];  k = list[i-1];
         dphi[j] = (tsp[j] - tsp[k])/(dw * (j - k));
      }
   }

   //
   //  Slope of dphi[] against 1/sqrt( 1 - wc^2/w^2) by Theil-Sen regression.
   //

   int nb = b->bin2 - b->bin1 + 1;

   struct TSA {
      double a;
      double w;
   } *tsa = VT_malloc_zero( sizeof( struct TSA) * nb * (nb-1)/2);

   double fc2 = EIC_CUTOFF * EIC_CUTOFF;

   int nts = 0;
   for (j = 0; j < nlist; j++)
      for (k = j+1; k < nlist; k++)
      {
         double fj = list[j] * DF;
         double fk = list[k] * DF;

         double xj = 1/sqrt( 1 - fc2/fj/fj);
         double xk = 1/sqrt( 1 - fc2/fk/fk);

         tsa[nts].a = (dphi[list[k]] - dphi[list[j]])/(xk - xj);
         tsa[nts].w = pow( xa[list[k]] * xa[list[j]], 2.0);
//         tsa[nts].w = xa[list[k]] * xa[list[j]];
//         tsa[nts].w = pow( fmin( xa[list[k]], xa[list[j]]), 4);
         nts++;
      }

   int cmp_tsa( const void *p1, const void *p2)
   {
      double v1 = ((struct TSA *)p1)->a;
      double v2 = ((struct TSA *)p2)->a;

      if (v1 < v2) return -1;
      if (v1 > v2) return 1;
      return 0;
   }

   qsort( tsa, nts, sizeof( struct TSA), cmp_tsa);

   double wsum = 0;
   for (j = 0; j < nts; j++) wsum += tsa[j].w;

   double wt = 0;
   double wtlast = 0;
   double ts_slope = 0;
   for (j = 0; j < nts; j++)
   {
      if (wt > wsum/2)
      {
         if (!j)
            ts_slope = tsa[j].a;
         else
         {
            double q = (wsum/2 - wtlast)/(wt - wtlast);
            double da = tsa[j].a - tsa[j-1].a;
            ts_slope = tsa[j-1].a + da * q;
         }
         break;
      }
      wtlast = wt;
      wt += tsa[j].w;
   }
// fprintf( stderr, "j=%d  nts/2=%d\n", j, nts/2);

   free( tsa);

   ev->range_day1 = -ts_slope * C;

   double ts_sum = 0;
   int ts_sumn = 0;

   for (i = 0; i < nlist; i++)
   {
      j = list[i];
      double f = j * DF;

      double x = 1/sqrt( 1 - fc2/f/f);
      ts_sum += dphi[j] - ts_slope * x;
      ts_sumn++;
   }

   double ts_icept = ts_sum / ts_sumn;

#if 1
// XXX
FILE *fr = fopen( "/tmp/fr.dat", "w");
   for (i = 0; i < nlist; i++)
   {
      j = list[i];
      double fj = j * DF;
      double xj = 1/sqrt( 1 - fc2/fj/fj);
      fprintf( fr, "%.1f %.6e %.6e %.6e %.6e\n",
           fj, xj, dphi[j], ts_icept + xj * ts_slope, xa[j]);
   }
fclose( fr);
#endif

   double vg = C * vf( b->F1, b->F2, EIC_CUTOFF);
// fprintf( stderr, "vf %.5f\n", vf( b->F1, b->F2, EIC_CUTOFF)); // XXX
   ev->range_day2 = ts_icept * vg;
#if 0
   if (ev->range <= 0)
   {
      rc_range++;
      return FALSE;
   }
#endif

/// XXX

   double r;
   double WC = 2 * M_PI * EIC_CUTOFF;
   double minsumsq = 0;
//   double rmin = 0;
   for (r = 100; r <= 10000; r += 100)
   {
      double sumsq = 0;
      for (i = 0; i < nlist; i++)
      {
         j = list[i];
         double w = 2 * M_PI * j * DF;
         double ph = r*w/vg - r*w/300e3 * sqrt(1 - WC*WC/w/w);
         double diff = ph - tsp[j]; 
         while (diff > M_PI) diff -= M_PI;
         while (diff < -M_PI) diff += M_PI;
         sumsq += diff * diff;
      }

      if (!minsumsq || minsumsq > sumsq)
      {
         minsumsq = sumsq;
//         rmin = r;
      }  
   }

//   fprintf( stderr, "r=%.0f\n", rmin); // XXX
}

//
//  Calculations per frequency band:
//  Average phase slope, TOGA from the average phase slope, RMS amplitude.
//

static int compute_band( struct EVENT *ev, struct BAND *b)
{
   int j;

   //
   //  Least squares regression to find the average phase slope.
   //

   double fsum = 0;
   double psum = 0;

   double fn = 0;
   for (j = b->bin1; j <= b->bin2; j++)
      if (!binmask[j])
      {
         double a = cabs( XX[j]);
         fsum += j * DF * a;
         psum += ev->upb[j] * a;
         fn += a;
      }

   double fmean = fsum / fn;
   double pmean = psum / fn;

   double num = 0;
   double denom = 0;

   for (j = b->bin1; j <= b->bin2; j++)
      if (!binmask[j])
      {
         double a = cabs( XX[j]);

         double d = j * DF - fmean;
         num += a * d * (ev->upb[j] - pmean);
         denom += a * d * d;
      }

   double alpha = pmean - num/denom * fmean; // Radians
   double beta = num/denom;                  // Radians per Hz

   //
   //  Offset of the TOGA from the start of the Fourier transform input.
   //


   b->toga_offset = -beta/(2 * M_PI);        // Cycles per Hz = seconds

   //
   //  Calculate:
   //  b->phresidual = RMS phase residual between linear and measured phase.
   //  b->energy = Total energy in this band
   //

   double rs = 0;
   for (j = b->bin1; j <= b->bin2; j++)
      if (!binmask[j])
      {
         double d = ev->upb[j] - (alpha + j * DF * beta);
         double a = cabs( XX[j]);
         rs += d * d * a;

         b->energy += a * a;
      }

   b->phresidual = sqrt( rs/fn);  // Radians 

   if (b->do_range) compute_range( ev, b);

   return TRUE;
}

static void compute_bearing( struct EVENT *ev)
{
   int j;

   if (polar_mode)
   {
      //
      //  Matrix to correct for the loop alignments
      //

      double cos1 = cos( polar1_align);
      double cos2 = cos( polar2_align);
      double sin1 = sin( polar1_align);
      double sin2 = sin( polar2_align);
      double det = sin1*cos2 - cos1*sin2;

      double bsin = 0;
      double bcos = 0;

      //
      //  Calculate the bearing for each frequency bin and do an average
      //  weighted by the RMS amplitude.   Runs over all the bins listed
      //  in all the bands, to produce an overall bearing.
      //

      int band;
      for (band = 0; band < nbands; band++)
         for (j = bands[band].bin1; j <= bands[band].bin2; j++)
         {
            complex double *H1 = channels[ch_HFIELD1].X;
            complex double *H2 = channels[ch_HFIELD2].X;
   
            // N/S and E/W signals, correcting for loop azimuths
            complex double ew = (cos2 * H1[j] - cos1 * H2[j]) * det;
            complex double ns = (-sin2 * H1[j] + sin1 * H2[j]) * det;
   
            double mag_ew = cabs( ew);
            double mag_ns = cabs( ns);
            double pow_ew = mag_ew * mag_ew;
            double pow_ns = mag_ns * mag_ns;
   
            // Phase angle between N/S and E/W
            double phsin = cimag( ns) * creal( ew) - creal( ns) * cimag( ew);
            double phcos = creal( ns) * creal( ew) + cimag( ns) * cimag( ew);
            double a = atan2( phsin, phcos);
   
            // Watson-Watt goniometry to produce cos and sine of 2*bearing.
            double bearing2sin = 2 * mag_ew * mag_ns * cos( a);
            double bearing2cos = pow_ns - pow_ew;
            double pwr = pow_ew + pow_ns;
   
            double weight = pwr;
   
            if (ch_EFIELD < 0)
            {
               // No E-field available, so average the sin,cos of 2*bearing
               bsin += bearing2sin * weight;
               bcos += bearing2cos * weight;
               continue;
            }
   
            // E-field available, compare phase of E with H 
            double bearing180 = atan2( bearing2sin, bearing2cos)/2;
            if (bearing180 < 0) bearing180 += M_PI;
            else
            if (bearing180 >= M_PI) bearing180 -= M_PI;
   
            //  H-field signal in plane of incidence
            complex double or = ew * sin( bearing180) +
                                ns * cos( bearing180);
   
            complex double vr = channels[ch_EFIELD].X[j];
   
            // Phase angle between E and H
            double pha =
                 atan2( cimag( or) * creal( vr) - creal( or) * cimag( vr),
                        creal( or) * creal( vr) + cimag( or) * cimag( vr));
   
            // Reflect the mod 180 bearing to the correct quadrant
            double bearing360 = bearing180;
            if (pha < -M_PI/2 || pha > M_PI/2) bearing360 += M_PI;
   
            // Average the sin,cos of the bearing
            bsin += sin( bearing360) * weight;
            bcos += cos( bearing360) * weight;
         }

      //
      //  Decide on an overall bearing.
      //

      if (ch_EFIELD < 0) // Bearing modulo 180
      {
         ev->bearing = atan2( bsin, bcos)/2;
         if (ev->bearing < 0) ev->bearing += M_PI;
      }
      else  // Bearing modulo 360
      {
         ev->bearing = atan2( bsin, bcos);
         if (ev->bearing < 0) ev->bearing += 2 * M_PI;
      }

      //
      //  Construct XX[] from the signal aligned on the overall bearing.
      //

      for (j = 0; j < BINS; j++)
      {
         complex double *H1 = channels[ch_HFIELD1].X;
         complex double *H2 = channels[ch_HFIELD2].X;

         // N/S and E/W signals, correcting for loop azimuths
         complex double ew = (cos2 * H1[j] - cos1 * H2[j]) * det;
         complex double ns = (-sin2 * H1[j] + sin1 * H2[j]) * det;

         //  H-field signal in plane of incidence
         complex double or = ew * sin( ev->bearing) +
                             ns * cos( ev->bearing);

         XX[j] = or + channels[ch_EFIELD].X[j];
      }
   }
   else
   {
      for (j = 0; j < BINS; j++) XX[j] = channels[0].X[j];
   }
}

struct MOVAV {
   int ns;
   int buflen;
   int bp;
   double sum;
   double *buf;
} ma1, ma2;

static void movav_init( struct MOVAV *m, int len)
{
   m->ns = 0;
   m->buflen = len;
   m->sum = 0;
   m->buf = VT_malloc_zero( m->buflen * sizeof( double));
   m->bp = 0;
}

static inline void movav_update( struct MOVAV *m, double new)
{
   double old = m->buf[m->bp];
   m->buf[m->bp] = new;
   m->bp = (m->bp + 1) % m->buflen;

   if (m->ns < m->buflen) m->ns++;
   m->sum += new - old;
}

static inline double movav_value( struct MOVAV *m)
{
   return m->ns ? m->sum / m->ns : 0;
}

//
//  
//  

static inline double impulse_ratio( struct MOVAV *m1, struct MOVAV *m2)
{
   double energy1 = movav_value( m1);
   double energy2 = movav_value( m2);
   return energy1 ? energy2/energy1 : 0;  
}

static inline int impulse_ratio_peak( int b)
{
   double ir_a = irbuf[BUFIDX(b-1)];
   double ir_b = irbuf[BUFIDX(b)];
   double ir_c = irbuf[BUFIDX(b+1)];

// printf( "irpeak %.2f %.2f %.2f\n", ir_a, ir_b, ir_c); // XXX
   return ir_b > ir_a && ir_b > ir_c; 
}

static void td_peaks( struct EVENT *ev)
{
   void add_tdpeak( int j)
   {
      // Quadratic interpolation to locate time offset of the peak
      double y1 = creal( ev->ca[j-1]);
      double y2 = creal( ev->ca[j]);
      double y3 = creal( ev->ca[j+1]);
      double d = (y3 - y1)/(2 * y2 - y1 - y3)/2;
      double a = y2 - (y1 - y3) * d/4;

      // Insertion-sort into the list of peaks
      int i;
      for (i = 0; i < ev->ntdpeaks; i++)
         if (fabs( a) > fabs( ev->tdpeaks[i].a)) break;

      if (i == MAXTDPEAKS) return;

      if (ev->ntdpeaks > 0 &&
          fabs( a) < 0.1 * fabs( ev->tdpeaks[0].a)) return;

      memmove( ev->tdpeaks + i + 1, ev->tdpeaks + i, 
               sizeof( struct TDPEAK) * (MAXTDPEAKS - i - 1));

      ev->tdpeaks[i].a = a;
      ev->tdpeaks[i].t = (j + d) * DT;

      // Instantaneous frequency at the peak
  
      double dp = carg( ev->ca[j+1]) - carg( ev->ca[j-1]);
      while (dp > M_PI) dp -= 2 * M_PI;
      while (dp < -M_PI) dp += 2 * M_PI;
      dp /= 2;
      ev->tdpeaks[i].f = dp/(DT * 2 * M_PI);

      if (ev->ntdpeaks < MAXTDPEAKS) ev->ntdpeaks++;

      for (i = 0; i < ev->ntdpeaks; i++)
         if (fabs( ev->tdpeaks[i].a) < 0.1 * fabs( ev->tdpeaks[0].a))
            ev->ntdpeaks = i;
   }

   int ns = Tpeaks * sample_rate;
   int b1 = Tlead * sample_rate;
   int j;
   for (j = b1 + 2; j < b1 + ns - 2; j++)
   {
      double v = creal( ev->ca[j]);
      if (v > 0 &&
          v > creal( ev->ca[j-2]) &&
          v > creal( ev->ca[j-1]) &&
          v > creal( ev->ca[j+1]) &&
          v > creal( ev->ca[j+2])) add_tdpeak( j);
      else
      if (v < 0 &&
          v < creal( ev->ca[j-2]) &&
          v < creal( ev->ca[j-1]) &&
          v < creal( ev->ca[j+1]) &&
          v < creal( ev->ca[j+2])) add_tdpeak( j);
   }

   double asum1 = 0, asum2 = 0;
   for (j = 0; j < ev->ntdpeaks; j++)
   {
      double a = ev->tdpeaks[j].a;
      asum1 += fabs( a);
      asum2 += a;
   }

   ev->peaks_ratio = asum1 > 0 ? asum2/asum1 : 0;
}

static void compute_xscore( struct EVENT *ev)
{
   int ns = Tpulse * sample_rate;
   int i1 = Tlead * sample_rate;
   int i2 = i1 + ns;

   int i;
   double sum = 0;
   for (i = i1; i < i2; i++)
   {
      double v = BUF(0, i); // XXX
      sum += v;
   }

   double mean = sum / ns;

   double sumpsq = 0;
   double sumnsq = 0;
   for (i = i1; i < i2; i++)
   {
      double v = BUF(0, i); // XXX
      if (v > mean) sumpsq += v*v;
      else          sumnsq += v*v;
   }

   ev->xscore = (sumpsq - sumnsq)/(sumpsq + sumnsq);
   if (ev->xscore < 0) ev->xscore = -ev->xscore;
}

//
//  Examine the time domain buffer to see if it holds a sferic.
//
//  Determine the TOGA of each band and channel, if possible.
//  Produce an amplitude weighted average TOGA from the channels
//  that gave a measurement.
//
//  If polar operation is called for, work out the bearing.
//

static void process_buffer( struct EVENT *ev)
{
   double ir = irbuf[BUFIDX(btrig)];

   if (trigger_impulse && ir < trigger_impulse) return;

   rc_raw++;

   reset_event( ev);
   ev->impulse_ratio = ir;
   reset_bands();

   // Timestamp of the first sample in the input buffer
   ev->Tf = timestamp_add( Tin, (-buflen + 1) * DT);

   //
   //  Fourier transform each channel.
   //

   int ch;
   for (ch = 0; ch < nchans; ch++)
   {
      struct CHAN *cp = channels + ch;

      //
      //  FFT the pulse and measure the RMS amplitude.  A duration of Tpulse
      //  seconds is analysed, the rest of the buffer duration Tfft is padded
      //  with zeros.
      //

      double sumsq = 0;
      int ns = Tpulse * sample_rate;
      int j;
      int b1 = Tlead * sample_rate;
      for (j = 0; j < b1; j++) cp->fft[j] = 0;
      for (j = b1; j < b1 + ns; j++)
      {
         double v = BUF(ch, j);
         cp->fft[j] = v;
         sumsq += v * v;
      }
      for (; j < FFTWID; j++) cp->fft[j] = 0;

      cp->rms = sqrt(sumsq/ns);
      ev->energy += cp->rms * cp->rms;

      fftw_execute( cp->fp);
   }

   compute_bearing( ev);
   make_analytic( ev);

   //
   //  Extract the bin phases into ev->tsp, rotating the phase to time
   //  shift to the trigger time.  This is where the phase is likely to
   //  be the flattest and therefore give the most reliable unwrapping.
   //

   double g = 0;
   ev->tsp[unmasked[0]] = g;
   int j;
   for (j = 0; j < nunmasked - 1; j++)
   {
      int b1 = unmasked[j];
      int b2 = unmasked[j+1];
      double dp = carg( XX[b2]) + b2 * 2 * M_PI * DF * Tpretrig
                - carg( XX[b1]) - b1 * 2 * M_PI * DF * Tpretrig;
      // carg() returns -pi to +pi, so dp can be -2pi to +2pi
      while (dp > M_PI) dp -= 2*M_PI;
      while (dp < -M_PI) dp += 2*M_PI;
      g += dp;
      ev->tsp[b2] = g;
   }

   //
   //  Time shift the unwrapped phase back to the Tf reference.
   //

   for (j = 0; j < nunmasked; j++)
   {
      int b = unmasked[j];
      ev->upb[b] = ev->tsp[b] - b * 2 * M_PI * DF * Tpretrig;
   }

   compute_xscore( ev);

   //
   //  Process each band.
   //

   int band;
   for (band = 0; band < nbands; band++) compute_band( ev, bands + band);

   //
   //  Compute quality scores.
   //

   evaluate_ascore( ev);
   if (ev->ascore < limit_ascore)
   {
      rc_ascore++;
      return;
   }

   evaluate_pscore( ev);
   if (ev->pscore < limit_pscore)
   {
      rc_pscore++;
      return;
   }

   evaluate_tscore( ev);
   if (ev->tscore < limit_tscore)
   {
      rc_tscore++;
      return;
   }

   td_peaks( ev);

   rc_out++;
   output_record( ev);

   N_out++;
   last_N++;
   if (Nmax && N_out == Nmax) VT_exit( "completed %d", Nmax);
}

//
//  The input buffer is a circular buffer of length buflen indexed through a
//  a macro BUF( channel, P).  P=0 is the oldest sample, P=buflen-1 is the
//  newest.  Tf is the timestamp of BUF(*,0) and Tin is the timestamp of the
//  latest sample BUF(*,buflen-1).
//
//   Tf                                                     Tin
//   |<-- Tpretrig -->|                                     |
//   XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//   |                |
//   P=0              P=btrig
//          

//
//  Called after each input frame.  Test to see if anything triggers and if
//  so, call process_buffer().
//

static int Dflag = FALSE;
static FILE *Dfile = NULL;

static void evaluate_trigger( struct EVENT *ev, int pki, double pkval)
{
   //
   //  Trigger on UT second marks if -c option is given
   //

   if (CFLAG)
   {
      timestamp T = timestamp_add( Tin, (-buflen + 1 + btrig) * DT);
      
      double t = timestamp_frac( T) + Tpretrig;
      if (t < 1 && t + DT >= 1) process_buffer( ev);

      return;
   }

   //
   //  A -T option requests a specific trigger time.
   //

   if (one_shot)
   {
      timestamp T = timestamp_add( Tin, (-buflen + 1 + btrig) * DT);
      if (timestamp_LT( T, Ttrig) &&
          timestamp_GE( timestamp_add( T, DT), Ttrig))
      {
         process_buffer( ev);
         VT_exit( "completed trigger");
      }

      return;
   }

   //
   //  Default behaviour is to trigger whenever the sample at P=btrig has the
   //  maximum absolute value in the input buffer.
   //

   if (!Dflag &&
        !timestamp_is_ZERO( T_D) &&
       timestamp_GE( timestamp_add( Tin, (-buflen + 1) * DT), T_D))
   {
      Dflag = 1;
      Dfile = fopen( "/tmp/d.dat", "w");
   }

   if (Dflag == 1 &&
       timestamp_GE( timestamp_add( Tin, (-buflen + 1) * DT),
                     timestamp_add( T_D, 0.010))) Dflag = 2;
   if (Dflag == 1)
   {
      char temp[50];
      timestamp_string6(  timestamp_add( Tin, (-buflen + 1 + btrig) * DT), temp);
      
      fprintf( Dfile, "%s %.6e %d %.5f\n", 
              temp,  
              BUF(0, btrig),
              pki == btrig ? 1 : 0,
              irbuf[BUFIDX(btrig)]);
   }

   if (pki == btrig)
   {
      //  If trigger_amplitude is set (-a option) then proceed if the amplitude
      //  at offset btrig is greater than the threshold.
   
if (pkval > trigger_amplitude) VT_report( 1, "trigger_amplitude reached: %.3e", pkval); 

      if (!trigger_amplitude || pkval > trigger_amplitude)
         process_buffer( ev);
   }
}

//
//  Revise the trigger threshold every AT_INT seconds.  The threshold is
//  changed by +/-10% according to whether the trigger rate over the last
//  AT_INT seconds is above or below the target rate.
//

static void revise_thresholds( void)
{
   double rate = last_N / (double) AT_INT;

   if (rate < throttle)
   {
      trigger_amplitude *= 0.9;
      trigger_impulse *= 0.9;
   }
   else
   {
      trigger_amplitude *= 1.1;
      trigger_impulse *= 1.1;
   }

   VT_report( 1, "rate %.2f triggers %.3f %.3f",
            rate, trigger_amplitude, trigger_impulse);
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

static void parse_foption( char *args)
{
   while (args && *args)
   {
      char *p = strchr( args, ',');
      if (p) p++;

      if (!strncmp( args, "ascore=", 7)) limit_ascore = atof( args+7);
      else
      if (!strncmp( args, "pscore=", 7)) limit_pscore = atof( args+7);
      else
      if (!strncmp( args, "tscore=", 7)) limit_tscore = atof( args+7);
      else
         VT_bailout( "unrecognised -f option: %s", args);

      args = p;
   }
}

static void parse_band( char *args)
{
   bands = VT_realloc( bands, (nbands + 1) * sizeof( struct BAND));
   struct BAND *b = bands + nbands++;

   VT_parse_freqspec( args, &b->F1, &b->F2);

   if (b->F2 <= b->F1) VT_bailout( "invalid frequency range with %s", optarg);
   char *p = strchr( args, ',');
   if (!p) return;
   p++;

   args = strchr( p, ',');
   if (!args) return;
   args++;

   while (args && *args)
   {
      p = strchr( args, ',');
      if (p) p++;

      if (!strncmp( args, "noscore", 7)) b->noscore = TRUE;
      else
      if (!strncmp( args, "range", 5)) b->do_range = TRUE;
      else VT_bailout( "cannot parse [%s]", args);

      args = p;
   }
}

static void setup_bands( void)
{
   int do_range = 0;
   int band;
   for (band = 0; band < nbands; band++)
   {
      struct BAND *b = bands + band;

      b->gcf = sqrt( b->F1 * b->F2);  // Geometric center frequency

      b->bin1 = b->F1/DF;
      b->bin2 = b->F2/DF;
      if (b->bin2 >= BINS) b->bin2 = BINS-1;

      if (b->do_range) do_range++;
   }

   //
   //  Check the options relating to range estimation.
   //

   if (do_range > 1) VT_bailout( "more than one band configured for range");
   if (do_range && !EFLAG_R)
      VT_bailout( "range estimate requested without -eR option");
   if (EFLAG_R && !do_range)
      VT_bailout( "-eR option given but no band specified for range estimate");

   //
   //  Sort bands into ascending order of center frequency.   The tscore
   //  mechanism requires this.
   //

   int cmp_bands( const void *p1, const void *p2)
   {
      struct BAND *b1 = (struct BAND *) p1,
                  *b2 = (struct BAND *) p2;
      if (b1->gcf > b2->gcf) return +1;
      if (b1->gcf < b2->gcf) return -1;
      return 0;
   }

   qsort( bands, nbands, sizeof( struct BAND), cmp_bands);
}

static void parse_mask( char *args)
{
   masks = VT_realloc( masks, (nmasks + 1) * sizeof( struct MASK));
   struct MASK *m = masks + nmasks++;

   VT_parse_freqspec( optarg, &m->F1, &m->F2);

   if (m->F2 <= m->F1) VT_bailout( "invalid frequency range with %s", optarg);
}

static void setup_masks( void)
{
   int mask;
   for (mask = 0; mask < nmasks; mask++)
   {
      struct MASK *m = masks + mask;

      m->gcf = sqrt( m->F1 * m->F2);  // Geometric center frequency

      m->bin1 = m->F1/DF;
      m->bin2 = m->F2/DF;
      if (m->bin2 >= BINS) m->bin2 = BINS-1;
   }

   //
   //  Sort masks into ascending order of center frequency.   The tscore
   //  mechanism requires this.
   //

   int cmp_masks( const void *p1, const void *p2)
   {
      struct MASK *m1 = (struct MASK *) p1,
                  *m2 = (struct MASK *) p2;
      if (m1->gcf > m2->gcf) return +1;
      if (m1->gcf < m2->gcf) return -1;
      return 0;
   }

   qsort( masks, nmasks, sizeof( struct MASK), cmp_masks);

   //
   //  An array to flag all the masked bins
   //

   binmask = VT_malloc_zero( sizeof( char) * BINS);
   for (mask = 0; mask < nmasks; mask++)
   {
      struct MASK *m = masks + mask;
      int bin;
      for (bin = m->bin1; bin <= m->bin2; bin++) binmask[bin] = 1;
   }

   //
   //  A list to contain all the unmasked bins.
   //

   unmasked = VT_malloc_zero( sizeof( int) * BINS);
   nunmasked = 0;
   int bin;
   for (bin = 0; bin < BINS; bin++)
      if (!binmask[bin]) unmasked[nunmasked++] = bin;
}

static void parse_ext( char *arg)
{
   char *t = strchr( arg, '/');
   int interval = t ? atoi( t + 1) : 1;

   for (; *arg; arg++)
      switch (toupper( *arg))
      {
         case 'B': EFLAG_B = interval;  break;
         case 'S': EFLAG_S = interval;  break;
         case 'T': EFLAG_T = interval;  break;
         case 'P': EFLAG_P = interval;  break;
         case 'Q': EFLAG_Q = interval;  break;
         case 'R': EFLAG_R = interval;  break;
         case 'X': EFLAG_X = interval;  break;
      }
}

int main( int argc, char *argv[])
{
   VT_init( "vttoga");

   int background = 0;
   char *polarspec = NULL;

   while (1)
   {
      int c = getopt( argc, argv, "vBL:F:M:E:d:G:a:i:T:p:N:f:e:r:u:cD:?");
      
      if (c == 'v') VT_up_loglevel();
      else
      if (c == 'B') background = 1;
      else
      if (c == 'E') Etime = atof( optarg);
      else
      if (c == 'L') VT_set_logfile( "%s", optarg);
      else
      if (c == 'F') parse_band( optarg);
      else
      if (c == 'M') parse_mask( optarg);
      else
      if (c == 'd') outdir = strdup( optarg);
      else
      if (c == 'G') gran = atoi( optarg);
      else
      if (c == 'p') polarspec = strdup( optarg);
      else
      if (c == 'T') Ttrig = VT_parse_timestamp( optarg);
      else
      if (c == 'e') parse_ext( strdup( optarg));
      else
      if (c == 'c') CFLAG = 1;
      else
      if (c == 'N') Nmax = atoi( optarg);
      else
      if (c == 'a') trigger_amplitude = atof( optarg);
      else
      if (c == 'i') trigger_impulse = atof( optarg);
      else
      if (c == 'f') parse_foption( strdup( optarg));
      else
      if (c == 'D')
      {
         T_D = VT_parse_timestamp( optarg);
      }
      else
      if (c == 'r') throttle = atof( optarg);
      else
      if (c == 'u') limit_u = atof( optarg);
      else
      if (c == -1) break;
      else
      {
         if (c != '?') VT_report( 0, "unrecognised option [%c]", c);
         usage(); 
      }
   }  
  
   if (argc > optind + 1) usage();
   char *bname = strdup( optind < argc ? argv[optind] : "-");

   if (!one_shot && !CFLAG && !trigger_impulse && !trigger_amplitude)
      VT_bailout( "no triggering option given, needs -a or -i or -T or -C");
 
   if (trigger_impulse && trigger_amplitude)
      VT_bailout( "both -a and -i given");
 
   //  If no -F options were given, put in a default band
   if (!nbands) parse_band( "4000,17000");
  
   if (background)
   {
      int flags = bname[0] == '-' ? KEEP_STDIN : 0;
      flags |= KEEP_STDOUT;
      VT_daemonise( flags);
   }

   struct VT_CHANSPEC *chspec = VT_parse_chanspec( bname);

   VTFILE *vtfile = VT_open_input( bname);
   if (!vtfile) VT_bailout( "cannot open: %s", VT_error);

   VT_init_chanspec( chspec, vtfile);
   nchans = chspec->n;
   sample_rate = VT_get_sample_rate( vtfile);
   VT_report( 1, "channels: %d, sample_rate: %d", nchans, sample_rate);

   if (!timestamp_is_ZERO( Ttrig)) one_shot = 1;

   if (polarspec)
   {
      VT_parse_polarspec( nchans, polarspec,
                          &ch_HFIELD1, &polar1_align,
                          &ch_HFIELD2, &polar2_align,
                          &ch_EFIELD);

      if (ch_HFIELD1 >= 0 && ch_HFIELD2 >= 0) polar_mode = 1;
   }

   //
   //  Set up buffers lengths, etc.
   // 

   buflen = Tbuf * sample_rate;
   btrig = Tpretrig * sample_rate;
   FFTWID = Tfft * sample_rate;
   BINS = FFTWID/2 + 1;
   DF = sample_rate/(double) FFTWID;
   DT = 1/(double)sample_rate;

   VT_report( 2, "buffer length: %d samples trig %d DF=%.3f",
                  buflen, btrig, DF);

   setup_bands();
   setup_masks();
   channels = VT_malloc_zero( sizeof( struct CHAN) * nchans);
   int i;
   for (i = 0; i < nchans; i++)
   {
      struct CHAN *cp = channels + i;

      cp->buf = VT_malloc_zero( sizeof( double) * buflen);
      cp->fft = VT_malloc( sizeof( double) * FFTWID);
      cp->X = VT_malloc( sizeof( complex double) * BINS);
      cp->fp = fftw_plan_dft_r2c_1d( FFTWID, cp->fft, cp->X, FFTW_ESTIMATE);
   }

   XX = VT_malloc( sizeof( complex double) * BINS);
   caX = VT_malloc( sizeof( complex double) * BINS);
   car = VT_malloc( sizeof( double) * FFTWID);
   ffca = fftw_plan_dft_c2r_1d( FFTWID, caX, car, FFTW_ESTIMATE);

   irbuf = VT_malloc_zero( sizeof( double) * buflen);

   double *frame;
   int ch;
   int nbuf = 0;

   Tstart = VT_get_timestamp( vtfile);

   struct EVENT ev;
   init_event( &ev);

   int bpp = -1;
   double pkval = 0;

   int pklim = sample_rate * 0.001;

   movav_init( &ma1, 200e-6 * sample_rate * nchans);
   movav_init( &ma2, 200e-6 * sample_rate * nchans);

   while (1)
   {
      //
      //  Read a frame and add to circular input buffers.  frame is NULL
      //  on end of input.
      //

      Tin = VT_get_timestamp( vtfile); 
      if ((frame = VT_get_frame( vtfile)) == NULL)
      {
         VT_report( 1, "end of input");
         break;
      }

      //
      //  Add new sample to buffer.  Keep track of which buffer index
      //  references the impulse ratio sample.
      //

      if (bpp == bp)   // The established peak is about to scroll out of the
      {                // buffer, so forget it
         bpp = -1;
         pkval = 0;
      }

      if (bpp < 0)     // Need to re-scan the buffer to select a new peak
      {
         int p;
         for (p = pklim-1; p > 1; p--)  // Skips oldest sample, about to be
                                        // discarded
         {
            int i = BUFIDX( p);

            if (trigger_impulse)
            {
               double v = irbuf[i];
               if (v > pkval)
               {
                  pkval = v;
                  bpp = i;
               }
            }
            else
               for (ch = 0; ch < nchans; ch++)
               {
                  double v = fabs( channels[ch].buf[i]);
                  if (v > pkval)
                  {
                     pkval = v;
                     bpp = i;
                  }
               }
         }
      }

      for (ch = 0; ch < nchans; ch++)
      {
         double v = frame[chspec->map[ch]];
         channels[ch].buf[bp] = v;
      }

      if (trigger_impulse && irbuf[BUFIDX(pklim)] > pkval)
      {
         pkval = irbuf[BUFIDX(pklim)];
         bpp = BUFIDX(pklim);
      }
      else
      {
         int i = BUFIDX(pklim);
         for (ch = 0; ch < nchans; ch++)
         {
            double v = fabs( channels[ch].buf[i]);
            if (v > pkval)
            {
               pkval = v;
               bpp = i;
            }
         }
      }

      bp = (bp + 1) % buflen;

      double a1 = 0, a2 = 0;
      int b1 = (Tpulse + 200e-6) * sample_rate + 0.5;
      int b2 = (Tpulse + 500e-6) * sample_rate + 0.5;
      for (ch = 0; ch < nchans; ch++)
      {
         double a;
         a = BUF( ch, b1);
         a1 += a * a;

         a = BUF( ch, b2);
         a2 += a * a;
      }

      movav_update( &ma1, a1);
      movav_update( &ma2, a2);
      irbuf[BUFIDX( (int)((Tpulse + 300e-6) * sample_rate + 0.5))] =
            impulse_ratio( &ma1, &ma2);

      //
      //  Once the buffer is full, start looking for triggers
      //

      if (nbuf < buflen) nbuf++;  
      else
      {
         int pki = (bpp - bp + buflen) % buflen;
         evaluate_trigger( &ev, pki, pkval);
      }

      //
      //  Reconsider the trigger threshold to aim for the target rate.
      //

      time_t t = timestamp_secs( Tin);
      if (throttle && t - last_T > AT_INT)
      {
         revise_thresholds();
         last_T = t;
         last_N = 0;
      }

      //
      //  Finish if reached a given end time.
      //

      if (Etime && timestamp_diff(Tin, Tstart) > Etime)
      {
         VT_report( 1, "completed %f seconds", Etime);
         break;
      }

      rc_review( Tin);
   }

   return 0;
}

