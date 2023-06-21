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

static VTFILE *vtinfile, *vtoutfile;
static char *inname = NULL;
static char *outname = NULL;
static int sample_rate = 0;

static int CFLAG = FALSE;       // Set by -c option: Clip at threshold,
                                // otherwise do blanking
static int DFLAG1 = FALSE;      // Set by -D1 option, output blanking waveform
                                // as an extra channel
static int DFLAG2 = FALSE;      // Set by -D2 option, output moving average
                                // waveform as an extra channel
static int ZFLAG = FALSE;       // Set by -z option, extend blanking out to
                                // previous and next zero-crossings

static double dwelltime = 0;    // From -d option, seconds
static double hfactor = 0;      // Set by -h option, fixed threshold level
static double afactor = 0;      // Set by -a option, auto threshold factor
static double ffactor = 0;      // Set by -f option, target blanking factor
static double speriod = 0;      // Set by -s option, seconds
static double matc = 100.0;     // From -t option: moving average time constant

static double mafac1 = 0;       // Moving average exponential coefficient
static double mafac2 = 0;       // Moving average exponential coefficient
static uint64_t nfp = 0;        // Total number of frames processed

static struct VT_CHANSPEC *espec = NULL;   // -e option: which channels to
                                           // examine for blanking

static struct VT_CHANSPEC *bspec = NULL;   // -b option: which channels to
                                           // apply blanking to

static struct CHAN
{
   double ma;              // Moving average noise floor
   int enable;             // TRUE if set to examine this channel for blanking
   int apply;              // TRUE if blanking is to be applied to this channel
   double outsum;          // Number of unblanked/unclipped output samples
}
 *channels = NULL;         // One for each selected input channel

static int nchans = 0;     // Number of input channels
static int ochans = 0;     // Number of output channels
static int dwellcnt;       // Dwell time sample countdown

//  Frame buffer: circular, length fblen frames.
//
//   <-   Tlead      -> <-          Ttail           ->
//   #################################################
//   ^                 ^
//   fbio              fbex
//   |                 |
//   |                  --- Examination point
//    --- Loading and unloading index
//
//  Tlead and Ttail are windows in which blanking can be extended out from
//  the examination point when looking for zero crossings.

#define TLEAD 0.1   // Seconds
#define TTAIL 0.1   // Seconds

static struct FBUF
{
   double *frame;
   int *zcoefs;    // Per channel coefficients when using zero-crossing
   float gcoef;    // Coefficient to apply to all (enabled) chans in this frame
   double ma;  
}
 *fb;

static int fblen = 0;      // Length of frame buffer
static int nfb = 0;        // Number of frames actually in the frame buffer
static int fbio = 0;       // Load/unload index of frame buffer

//  Blanking coefficient usage:
//
//    out[chan] = in[chan] * gcoef * zcoefs[chan]
//
//  gcoef and zcoef[chan] are initialised to 1
//  One or other or both will be set to zero (or a slope[] value) to apply
//  blanking.
//
//  These coefs allow the actual blanking to be deferred until we are ready
//  to output the frame.  Previously blanked samples can then still be looked
//  at, eg for updating the moving average.  Thus we can apply blanking into
//  future samples.

static int ntlead = 0;     // TLEAD, number of samples
static int nttail = 0;     // TTAIL, number of samples

static double *slope;           // Coefficients for smooth blanking edges
static int nslope = 0;          // Length of speriod in samples
static int active = FALSE;      // TRUE if blanking is activated

//
//  Asymmetric exponential moving average to track the noise floor.  Easy to
//  pull the average down, hard to pull the average up.
//

static inline void update_ma( struct CHAN *cp, double av)
{
   if (av > cp->ma)
      cp->ma = cp->ma * mafac1 + av * (1-mafac1);   // Slow response
   else
      cp->ma = cp->ma * mafac2 + av * (1-mafac2);   // Fast response
}

//
//  Apply blanking to a channel, backwards through the frame buffer starting
//  from the examination index, until reaching the previous zero crossing
//

static void activate_from_zc( int ch, int idx)
{
   int i;
   int s = signbit( fb[idx].frame[ch]);

   for (i = 1; i < ntlead; i++)
   {
      int j = (idx - i + fblen) % fblen;
      if (signbit( fb[j].frame[ch]) != s) return;  // Reached zero crossing?
 
      fb[j].zcoefs[ch] = 0;
   }
}

//
//  Apply blanking to a channel, forwards from the examination point, until
//  the next zero crossing is reached.
//

static void deactivate_to_zc( int ch, int idx)
{
   int i;
   int s = signbit( fb[idx].frame[ch]);

   for (i = 1; i < nttail; i++)
   {
      int j = (idx + i) % fblen;
      if (signbit( fb[j].frame[ch]) != s) return;   // Reached zero crossing?
 
      fb[j].zcoefs[ch] = 0;
   }
}

static void eval_blanking( int idx)
{
   double *frame = fb[idx].frame;
   int ch;

   int trigger = FALSE;  // Will be set TRUE if any of the enabled channels
                         // exceed the threshold

   for (ch = 0; ch < nchans; ch++)
   {
      struct CHAN *cp = channels + ch;

      if (!cp->enable) continue;

      double av = fabs( frame[ch]);

      double th;

      if (afactor)      // We were given a -a option?
      {
         // Update the moving average and set the threshold
         update_ma( cp, av);
         th = afactor * cp->ma;

         if (DFLAG2) fb[idx].ma = cp->ma;
      }
      else 
         th = hfactor;   // Otherwise, use fixed threshold from -h option

      if (av > th) trigger = TRUE;    // Input has exceeded blanking threshold?
   }

   if (trigger)   // Blanking to be activated?
   {
      if (!active && speriod)
      {
         int i;
         for (i = 1; i < nslope; i++)
         {
            float *g = &fb[(idx - i + fblen) % fblen].gcoef;
            *g = fmin( *g, slope[i]);
         }
      }

      // Zero crossing used?  Backdate this activation to the previous zero
      // crossing (independently for each channel)

      if (!active && ZFLAG)
         for (ch = 0; ch < nchans; ch++)
            if (channels[ch].apply) activate_from_zc( ch, idx);

      // Blank one sample, plus however many extra required by the dwell time.
      dwellcnt = 1 + round( dwelltime * sample_rate);

      active = TRUE;
   }

   if (active)
   {
      fb[idx].gcoef = 0;

      dwellcnt--;
      if (!dwellcnt)     // End of activation?
      {
         active = FALSE;

         // Begin start of smooth trailing edge, if -s used.
         // If using zero crossing, extend the activation to the next zero
         // crossing, per channel

         if (speriod)
         {
            int i;
            for (i = 1; i < nslope; i++)
            {
               float *g =  &fb[(idx + i) % fblen].gcoef;
               *g = fmin( *g, slope[i]);
            }
         }
         else
         if (ZFLAG)
            for (ch = 0; ch < nchans; ch++)
               if (channels[ch].apply) deactivate_to_zc( ch, idx);
      }
   }
}

static void eval_clipping( int idx)
{
   double *frame = fb[idx].frame;
   int ch;

   for (ch = 0; ch < nchans; ch++)
   {
      struct CHAN *cp = channels + ch;

      if (!cp->enable) continue;

      double v = frame[ch];
      double av = fabs( v);

      double th;
      if (afactor)      // We were given a -a option?
      {
         // Update the moving average and set the threshold
         update_ma( cp, av);
         th = afactor * cp->ma;
      }
      else th = hfactor;

      if (av > th) frame[ch] = v > 0 ? th : -th;
      else cp->outsum++;
   }
}

static void apply( double *frame, int n)
{
   int ch;

   for (ch = 0; ch < nchans; ch++)
   {
      struct CHAN *cp = channels + ch;

      if (cp->apply)
      {
         double g = fb[n].gcoef;
         if (ZFLAG) g *= fb[n].zcoefs[ch];

         frame[ch] *= g;
         cp->outsum += g;
      }
   }
}

static void final_report( void)
{
   int i;

   VT_report( 1, "frames processed %lld", (long long) nfp);

   for (i = 0; i < nchans; i++)
      if (channels[i].apply)
      {
         // Number of samples affected by blanking or clipping
         double dropsum = nfp - channels[i].outsum;

         // Proportion of the total samples processed
         double dropfactor = dropsum / nfp;

         VT_report( 1, "summary %d %.8f", i+1, dropfactor); 
      }
}

static void usage( void)
{
   fprintf( stderr,
       "usage:  vtblank [options] input output\n"
       "\n"
       "options:\n"
       "  -v         Increase verbosity\n"
       "  -B         Run in background\n"
       "  -L name    Set logfile\n"
       "  -c         Clip at the threshold\n"
       "             (default is to blank the output)\n"
       " -d secs     Dwell time, seconds\n"
       "             (default zero)\n"
       " -h thresh   Threshold amplitude\n"
       " -a factor   Automatic threshold factor\n"
       " -f factor   Target blanking factor\n"
       " -t secs     Time constant for moving average\n"
       "             (default 100 seconds)\n"
       " -s stime    Set smoothing time, seconds\n"
       "             (default zero)\n"
       " -e chanspec Examine only these channels\n"
       " -b chanspec Apply blanking to only these channels\n"
       " -z          Blank between zero crossing\n"
       " -D1         Output blanking waveform as extra channel\n"
     );
   exit( 1);
}

int main( int argc, char *argv[])
{
   VT_init( "vtblank");

   int background = 0;

   while (1)
   {
      int c = getopt( argc, argv, "vBcf:d:h:a:s:t:e:b:L:D:z?");
      
      if (c == 'v') VT_up_loglevel();
      else
      if (c == 'B') background = 1;
      else
      if (c == 'L') VT_set_logfile( "%s", optarg);
      else
      if (c == 'c') CFLAG = 1;
      else
      if (c == 'd') dwelltime = atof( optarg);
      else
      if (c == 'h') hfactor = atof( optarg);
      else
      if (c == 'f') ffactor = atof( optarg);
      else
      if (c == 'a') afactor = atof( optarg);
      else
      if (c == 't') matc = atof( optarg);
      else
      if (c == 's') speriod = atof( optarg);
      else
      if (c == 'e')
      {
         char *s; 
         if (asprintf( &s, ":%s", optarg) > 0) espec = VT_parse_chanspec( s);
      }
      else
      if (c == 'b')
      {
         char *s; 
         if (asprintf( &s, ":%s", optarg) > 0) bspec = VT_parse_chanspec( s);
      }
      else
      if (c == 'z') ZFLAG = TRUE;
      else
      if (c == 'D')
         switch (atoi( optarg))
         {
            case 1: DFLAG1 = TRUE; break;
            case 2: DFLAG2 = TRUE; break;
            default: VT_report( 0, "unrecognised -D argument");
                     usage();
         }
      else
      if (c == -1) break;
      else
         usage();
   }

   if (hfactor && afactor)
      VT_bailout( "options contradict: -a and -h both given");

   if (ZFLAG && CFLAG) VT_bailout( "cannot use -z with -c");
   if (ZFLAG && speriod) VT_bailout( "cannot use -z with -s");

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

   //
   //  Set up input stream.
   //

   struct VT_CHANSPEC *chspec = VT_parse_chanspec( inname);
   vtinfile = VT_open_input( inname);
   if (!vtinfile) 
      VT_bailout( "cannot open input %s: %s", inname, VT_error);

   VT_init_chanspec( chspec, vtinfile);
   nchans = chspec->n;
   if (!nchans) VT_bailout( "no input channels");
   sample_rate = VT_get_sample_rate( vtinfile);

   VT_report( 1, "channels: %d, sample_rate: %d", nchans, sample_rate);
   if (hfactor) VT_report( 1, "hfactor: %.3e", hfactor);
   if (afactor) VT_report( 1, "afactor: %.3f", afactor);

   channels = VT_malloc( sizeof( struct CHAN) * nchans);
   memset( channels, 0, sizeof( struct CHAN) * nchans);

   //
   //  Set up output stream.
   //

   ochans = nchans;
   if (DFLAG1) ochans++;     // One extra channel if -D1 given
   if (DFLAG2) ochans++;     // One extra channel if -D2 given

   vtoutfile = VT_open_output( outname, ochans, 0, sample_rate);
   if (!vtoutfile) VT_bailout( "cannot open: %s", VT_error);

   //
   //  Which channels to examine for blanking threshold and which
   //  channels to apply blanking to.
   //

   int i;
   if (espec)       // A -e option was given?
   {
      VT_report( 2, "espec: %d channels", espec->n );
      for (i = 0; i < espec->n; i++)
      {
         if (espec->map[i] >= nchans)
            VT_bailout( "invalid channel given to -e");
         channels[espec->map[i]].enable = 1;
      }
   }
   else            // All input channels enabled for examination
      for (i = 0; i < nchans; i++) channels[i].enable = 1;

   if (bspec)      // A -b option was given?
   {
      VT_report( 2, "bspec: %d channels", bspec->n );
      for (i = 0; i < bspec->n; i++)
      {
         if (bspec->map[i] >= nchans)
            VT_bailout( "invalid channel given to -e");
         channels[bspec->map[i]].apply = 1;
      }
   }
   else            // Blanking applies to all input channels
      for (i = 0; i < nchans; i++) channels[i].apply = 1;

   for (i = 0; i < nchans; i++)
      VT_report( 2, "channel %d: %s %s", i+1, 
         channels[i].enable ? "examine": "ignore", 
         channels[i].apply  ? "apply" : "bypass");

   // Set up the two moving average coefficients
   mafac1 = exp( -1.0/(matc * sample_rate));
   mafac2 = exp( -1.0/(matc/100 * sample_rate));
 
   // Set up array of coefficients for smooth blanking edges 
   if (speriod)    // -s option given?
   {
      nslope = 1 + speriod * sample_rate;
      slope = VT_malloc( sizeof( double) * nslope);
      for (i = 0; i < nslope; i++)
      {
         double a = i/(double)nslope;
         slope[i] = 1/(1 + exp(-(a-0.46)*12));
      }
   }

   //
   //  Set up the frame buffer.
   //

   double Tlead = TLEAD;
   double Ttail = TTAIL;

   if (speriod > Tlead) Tlead = speriod;  // Maybe -s value longer than TLEAD?
   if (speriod > Ttail) Ttail = speriod;  // Maybe -s value longer than TTAIL?

   ntlead = Tlead * sample_rate;
   nttail = Ttail * sample_rate;
   fblen = 1 + ntlead + nttail;    // Total length of frame buffer

   VT_report( 2, "buffer size %d frames", fblen);

   fb = VT_malloc( sizeof( struct FBUF) * fblen);
   for (i = 0; i < fblen; i++)
   {
      fb[i].frame = VT_malloc_zero( sizeof( double) * ochans);
      fb[i].zcoefs = VT_malloc_zero( sizeof( int) * nchans);
   }

   int *izcoefs = VT_malloc( sizeof( int) * nchans);
   for (i = 0; i < nchans; i++) izcoefs[i] = 1;

   //
   //  Which input channel will -D1 use.
   //

   int d1chan = 0;
   for (i = 0; i < nchans; i++)
      if (channels[i].apply)
      {
         d1chan = i;
         break;
      }

   //
   //  Main loop.
   //

   VT_bailout_hook( final_report);

   void (*eval)(int) = CFLAG ? eval_clipping : eval_blanking;

   while (1)
   {
      int e;
      if ((e = VT_is_block( vtinfile)) < 0)
      {
         VT_report( 1, "end of input");
         break;
      }

      double *inframe = VT_get_frame( vtinfile);
      for (i = 0; i < nchans; i++)
         fb[fbio].frame[i] = inframe[chspec->map[i]];

      fb[fbio].gcoef = 1;
      if (ZFLAG) memcpy( fb[fbio].zcoefs, izcoefs, sizeof( int) * nchans);

      // Evaluate the signal, the evaluation point lags behind the load point
      // by ntlead samples

      eval( (fbio + ntlead) % fblen);

      // Step the buffer in/out index, will then point to the oldest frame
      fbio = (fbio + 1) % fblen;
      if (nfb < fblen) nfb++;    // Buffer still filling?
      else
      {
         // Apply blanking to the oldest frame and send to output

         apply( fb[fbio].frame, fbio);

         if (DFLAG1)
            fb[fbio].frame[ochans-1] =
               ZFLAG ? fb[fbio].gcoef * fb[fbio].zcoefs[d1chan]
                     : fb[fbio].gcoef;
         
         if (DFLAG2)
            fb[fbio].frame[ochans-1] = fb[fbio].ma;
         
         if (!vtoutfile->nfb)  // Starting a new output block?
         {
            // Set output timestamp, always fblen samples earlier than input

            timestamp T = VT_get_timestamp( vtinfile);
            double srcal = VT_get_srcal( vtinfile);
            VT_set_timebase( vtoutfile,
                    timestamp_add( T, -fblen/(double) sample_rate), srcal);
         }

         VT_insert_frame( vtoutfile, fb[fbio].frame);
         nfp++;
      }
   }

   //
   //  Drain the buffer.
   //

   while (nfb--)
   {
      fbio = (fbio + 1) % fblen;
      apply( fb[fbio].frame, fbio);

      if (DFLAG1)
         fb[fbio].frame[ochans-1] =
            ZFLAG ? fb[fbio].gcoef * fb[fbio].zcoefs[d1chan]
                  : fb[fbio].gcoef;
      VT_insert_frame( vtoutfile, fb[fbio].frame);
      nfp++;
   }

   VT_release( vtoutfile);

   final_report();
   return 0;
}

