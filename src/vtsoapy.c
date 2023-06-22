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

#include <SoapySDR/Device.h>
#include <SoapySDR/Formats.h>
#include <SoapySDR/Logger.h>

///////////////////////////////////////////////////////////////////////////////
//  Globals                                                                  //
///////////////////////////////////////////////////////////////////////////////

static char *bname;                               // Output buffer or fifo name
static VTFILE *vtfile;                               // Handle to output stream

static double srcal = 1.0;               // Sample rate calibration coefficient
static unsigned int sample_rate = 0;                     // Nominal sample rate

static char *device = NULL;                                   // From -d option
static double frequency = -1;                                // Tuner frequency

static timestamp timebase = timestamp_ZERO;            // Base of our timestamp
static uint64_t ntb = 0;             // Number of frames since timebase was set

static uint64_t nout = 0;                            // Number of output frames
static double tadj = 0;     // Timing adjustment to be applied per output block
static timestamp Tfalse = timestamp_ZERO;    // False timestamp, from -T option
static int UFLAG = FALSE;                            // Set TRUE with -u option
static int QFLAG = FALSE;                         // -q option: invert Q signal
static char *g_option = NULL;                       // -g option: set SDR gains
static long double TOFFS = 0;                     // Timebase offset, -a option

///////////////////////////////////////////////////////////////////////////////
//  Timestamping                                                             //
///////////////////////////////////////////////////////////////////////////////

#define OFFSET_LIMIT 20e-3 // Reset if more than this many seconds timing error
#define PRE_RUN_SECS 20                         // Number of seconds of pre-run

static timestamp reftime = timestamp_ZERO;               // Reference timestamp
static uint32_t nft = 0;         // Number of frames since reftime last updated

// State variable for timestamping:-
#define STATE_RESET  0    // Transient state - reset everything and begin again
#define STATE_INIT   1    // Pre-run to settle down and get rough sample rate
#define STATE_SETUP  2    // Polish the sample rate
#define STATE_RUN    3    // Normal running

static int state = STATE_INIT;

//
//  Called immediately after a buffer of 'q' frames is received from the
//  device.
//
static void timestamping( int q)
{
   static int n = 0;

   #define AJF_RUN 0.05   // Rate error adjustment factor
   #define AJF_SETUP 0.15
   #define MAX_DEFER 3

   nft += q;

   if (state == STATE_RESET)  
   {
      // When reading from stdin, or if -u option is given, don't bother with
      // timing, go straight to the 'run' state.
      if (UFLAG)
      {
         state = STATE_RUN;
         timebase = reftime =
                          !timestamp_is_ZERO( Tfalse) ? Tfalse : VT_rtc_time();
         nft = 0;
         ntb = nout = 0;
         return;
      }

      reftime = VT_rtc_time();
      nft = 0;
      state = STATE_INIT;
   }

   if (state == STATE_INIT) // Initial pre-run to let things settle down
   {
      if (nft < PRE_RUN_SECS * sample_rate) return;

      // End of pre-run: just get a rough estimate of sample rate
      state = STATE_SETUP;
      n = 0;

      srcal =
         (nft / timestamp_diff(VT_rtc_time(), reftime)) / (double) sample_rate;
      timebase = reftime = VT_rtc_time();
      ntb = nout = 0;
      nft = 0;

      VT_report( 1, "pre-run complete %.2f", srcal * sample_rate);
      return;
   }

   // No sample rate or timebase calibration when using -u
   if (UFLAG) return;

   // Only revise the timing roughly every 10 seconds
   if (nft < sample_rate * 10) return;

   // Expected time to read nft samples, based on the current srcal
   double expected_interval = nft/(double)(srcal * sample_rate);

   // Compare with RTC to calculate a raw sample rate error
   timestamp now = VT_rtc_time();
   double actual_interval = timestamp_diff( now, reftime);

   double err = (expected_interval - actual_interval)/actual_interval;
   double r = fabs( err  * vtfile->bsize/100.0); 

   double raw_rate = nft / actual_interval;

   static int defer = 0;

   if (state == STATE_SETUP) // Initial stabilising of sample rate loop
   {
      // Update smoothed estimate of sample rate calibration factor
      srcal *= 1 + err * AJF_SETUP;

      reftime = now; nft = 0;

      VT_report( 1, "setup %+.3e sr %.2f rr %.1f n=%d r=%.2f",
                    err, sample_rate * srcal, raw_rate, n, r);
      if (r > 20)
      {
         VT_report( 0, "rate error too large, resetting");
         state = STATE_RESET;
      }
      else
      if (r < 0.1 && n >= 0) n += 3;
      else
      if (r < 0.5 && n >= 0) n += 2;
      else
      if (r < 1 && n >= 0) n++;
      else
      if (r < 1) n = 0;
      else 
         n--;

      if (n < -10)
      {
         VT_report( 0, "persistent drift, resetting");
         state = STATE_RESET;
      }
      if (n >= 10)
      {
         // Success.  Timebase has settled enough to start running data.
         state = STATE_RUN;
         timebase = reftime;
         ntb = nout = 0;
      }
   }
   else
   if (state == STATE_RUN)   // Normal operation
   {
      // Allow a couple of bad readings to be ignored - they usually come
      // good again.  Occurs when we didn't get scheduled promptly enough.
      if (r > 0.5 && defer++ <= MAX_DEFER) return;

      // Update smoothed estimate of sample rate calibration factor
      srcal *= 1 + err * AJF_RUN;

      timebase = timestamp_add( timebase, expected_interval); ntb += nft;
      reftime = now; nft = 0;
      defer = 0;
      double offset =
                  timestamp_diff( timebase, reftime); // Timing offset, seconds

      // Timing offset is slewed towards zero by adjusting the timebase
      // by tadj seconds per output data block.
      tadj = offset / (10 * sample_rate/vtfile->bsize);

      // Limit the timebase adjustment to 0.25 samples per block
      if (tadj > 0.25/sample_rate) tadj = 0.25/sample_rate;
      if (tadj < -0.25/sample_rate) tadj = -0.25/sample_rate;

      VT_report( 1, "run %+.3e sr %.3f tadj %+.3e rr %.1f offs %+.3fmS r=%.2f",
                 err, sample_rate * srcal, tadj, raw_rate, 1000 * offset, r);
 
      if (fabs( offset) > OFFSET_LIMIT)
      {
         // Either the R2832U sample rate or the system clock has stepped
         // or is drifting too fast for the control loops to correct.
         // Reset everything and start again.
         VT_report( 0, "timebase error %.3f mS, resetting", 1000 * offset);
         state = STATE_RESET;
      }
   }
}

///////////////////////////////////////////////////////////////////////////////
//  Soapy Interface                                                          //
///////////////////////////////////////////////////////////////////////////////

static SoapySDRDevice *dev = NULL;
static SoapySDRStream *stream = NULL;
static const char * format = NULL;

static void soapy_error_log( const SoapySDRLogLevel level, const char *s)
{
   int i = 1;

   if (level == SOAPY_SDR_FATAL || level == SOAPY_SDR_CRITICAL) i = 0;

   VT_report( i, "soapy: %s", s);
}

static void setup_soapy( int nframes)
{
//   SoapySDRsetLogLevel( 99);

   SoapySDR_registerLogHandler( (SoapySDRLogHandler) soapy_error_log);

   //
   //  Open Soapy device.
   //

   if (!device) device = "";
   
   SoapySDRKwargs stream_args = {0};
   if ((dev = SoapySDRDevice_makeStrArgs( device)) == NULL)
      VT_bailout( "cannot find soapy device");

   if ((vtfile->flags & VTFLAG_FMTMASK) == VTFLAG_INT2)
      format = SOAPY_SDR_CS16;
   else
      format = SOAPY_SDR_CF32;
  
   VT_report( 1, "request format [%s]", format); 
   if (SoapySDRDevice_setupStream( dev, &stream, SOAPY_SDR_RX, format,
                                   NULL, 0, &stream_args) != 0)
      VT_bailout( "cannot initialise device");

   VT_report( 1, "soapy device [%s] format [%s]",
            SoapySDRDevice_getHardwareKey( dev), format);

   //
   //  Set up device.
   //

   if (SoapySDRDevice_setSampleRate( dev, SOAPY_SDR_RX, 0, sample_rate))
      VT_bailout( "cannot set sample rate");

   SoapySDRKwargs args = {0};
   if (SoapySDRDevice_setFrequency( dev, SOAPY_SDR_RX, 0, frequency, &args))
      VT_bailout( "cannot set frequency");

   if (SoapySDRDevice_activateStream(dev, stream, 0, 0, 0))
      VT_bailout( "cannot activate device");

   if (g_option)
   {
      char *s = strdup( g_option);
      while (s && *s)
      {
         char *p = strchr( s, ',');
         if (p) *p++ = 0;

         char *q = strchr( s, '=');
         if (!q) VT_bailout( "cannot parse gain argument [%s]", g_option);
         *q++ = 0;
         double v = atof( q);

         if (SoapySDRDevice_setGainElement( dev, SOAPY_SDR_RX, 0, s, v))
            VT_bailout( "cannot set gain option [%s] = [%s]", s, q);
         VT_report( 1, "gain option [%s] = [%s] set OK", s, q);
         s = p;
      }
   }

   VT_report( 1, "soapy device activated");
}

static int read_data( uint8_t *buff, int nframes)
{
   void *buffs[] = { buff};
   int flags = 0;
   long long timeNs = 0;
                        
   int nr = SoapySDRDevice_readStream( dev, stream, buffs,
                                       nframes, &flags, &timeNs, 1000000);
   if (nr <= 0) VT_bailout( "device read error");

   return nr;
}

static void output_block( uint8_t *buff, int nframes)
{
   if (state != STATE_RUN) return; // Not yet in running state - discard output

   timebase  = timestamp_add( timebase, -tadj);
   VT_set_timebase( vtfile,
         timestamp_add( timebase,
                        TOFFS + (nout - ntb)/(srcal * sample_rate)), srcal);

   VT_next_write( vtfile);

   if (!strcmp( format, SOAPY_SDR_CF32) &&
       (vtfile->flags & VTFLAG_FMTMASK) == VTFLAG_FLOAT4)
   {
      float *d = (float *) VT_data_p( vtfile);
      memcpy( d, buff, vtfile->chans * nframes * 4);
   }
   else
   if (!strcmp( format, SOAPY_SDR_CS16) &&
       (vtfile->flags & VTFLAG_FMTMASK) == VTFLAG_INT2)
   {
      int16_t *d = (int16_t *) VT_data_p( vtfile);
      memcpy( d, buff, vtfile->chans * nframes * 2);
   }
   else
   if (!strcmp( format, SOAPY_SDR_CF32) &&
       (vtfile->flags & VTFLAG_FMTMASK) == VTFLAG_FLOAT8)
   {
      double *d = (double *) VT_data_p( vtfile);
      float *s = (float *) buff;
      int i;
      for (i = 0; i < nframes * 2; i++) d[i] = s[i];
   }

   vtfile->nfb += nframes;
   vtfile->nft += nframes;
   VT_release( vtfile);
   nout += nframes;
}

static void run( void)
{
   uint8_t *buff = VT_malloc( vtfile->bsize * 8);
   VT_report( 1, "blocksize %d chans %d", vtfile->bsize, vtfile->chans);

   while (1)
   {
      int q = read_data( buff, vtfile->bsize);
      if (q <= 0) break;

      timestamping( q);
      output_block( buff, q);
   }

   VT_release( vtfile);
}

///////////////////////////////////////////////////////////////////////////////
//  Main                                                                     //
///////////////////////////////////////////////////////////////////////////////

// Set scheduling priority to the minimum SCHED_FIFO value.
static void set_scheduling( void)
{
   #ifndef HAVE_SCHED_SETSCHEDULER
      int pri = -15;
      if (setpriority( PRIO_PROCESS, getpid(), pri))
      {
         VT_report( 0, "cannot set scheduling priority: %s",
                           strerror( errno));
         return;
      }
   #else
      int pri = sched_get_priority_max( SCHED_FIFO);

      struct sched_param pa;
      pa.sched_priority = pri;
      if (sched_setscheduler( 0, SCHED_FIFO, &pa))
      {
         VT_report( 0, "cannot set scheduling priority: %s",
                          strerror( errno));
         return;
      }
   #endif

   VT_report( 1, "using SCHED_FIFO priority %d", pri);

   if (mlockall( MCL_CURRENT | MCL_FUTURE) < 0)
      VT_report( 0, "unable to lock memory: %s", strerror( errno));
}

static void usage( void)
{
   fprintf( stderr,
       "usage:    vtsoapy [options] buffer_name\n"
       "\n"
       "options:\n"
       "  -v        Increase verbosity\n"
       "  -B        Run in background\n"
       "  -L name   Specify log file\n"
       "\n"
       "  -d device Device selector, comma-separated key=value pairs\n"
       "            eg -d driver=rtlsdr\n"
       "  -d ?      List available devices\n"
       "  -r rate   Sample rate Hz, floating point (no default)\n"
       "            eg -r 2.4e6\n"
       "\n"
       "  -F hertz  Tuner frequency, Hertz\n"
       "  -g gain   Gain options, comma-separated key=value pairs\n"
       "            eg -g LNA=12,MIX=12,AMP=10\n"
       "  -q        Invert the Q signal\n"
       "\n"
       "  -u        No sample rate tracking\n"
       "  -T stamp  Start time when using -u\n"
       "            (default is to use system clock)\n"
       "  -a offset Add timestamp offset, seconds\n"
     );
   exit( 1);
}

static void list_devices( void)
{
   size_t length;

   SoapySDRKwargs *results = SoapySDRDevice_enumerate( NULL, &length);
   int i, j;
   for (i = 0; i < length; i++)
   {
      printf( "device %d:\n", i);
      for (j = 0; j < results[i].size; j++)
         printf("   %s: %s\n", results[i].keys[j], results[i].vals[j]);
   }
   SoapySDRKwargsList_clear(results, length);
}

int main( int argc, char *argv[])
{
   VT_init( "vtsoapy");

   int background = 0;

   while (1)
   {
      int c = getopt( argc, argv, "vBd:r:g:F:T:L:a:uq?");

      if (c == 'v') VT_up_loglevel();
      else
      if (c == 'B') background = 1;
      else
      if (c == 'L') VT_set_logfile( "%s", optarg);
      else
      if (c == 'd') device = strdup( optarg);
      else
      if (c == 'r')
      {
         double fr = atof( optarg);
         sample_rate = fr;
         srcal = fr / sample_rate;
      }
      else
      if (c == 'g') g_option = strdup( optarg);
      else
      if (c == 'T') Tfalse = VT_parse_timestamp( optarg);
      else
      if (c == 'u') UFLAG = TRUE;
      else
      if (c == 'F') frequency = atof( optarg);
      else
      if (c == 'q') QFLAG = TRUE;
      else
      if (c == 'a') TOFFS = strtold( optarg, NULL);
      else
      if (c == -1) break;
      else
         usage();
   }

   if (device && !strcmp( device, "?"))
   {
      list_devices();
      exit( 0);
   }

   if (argc > optind + 1) usage();
   bname = strdup( optind < argc ? argv[optind] : "-");

   if (sample_rate <= 0)
      VT_bailout( "invalid or missing sample rate, needs -r");

   if (frequency < 0)
      VT_bailout( "invalid frequency");

   if (frequency <= 0 )
      VT_bailout( "missing or invalid frequency, needs -F");

   if (background)
   {
      int flags = bname[0] == '-' ? KEEP_STDOUT : 0;
      VT_daemonise( flags);
   }

   VT_report( 1, "buffer name: [%s]", bname);

   vtfile = VT_open_output( bname, 2, 1, sample_rate);
   if (!vtfile) VT_bailout( "cannot create buffer: %s", VT_error);

   setup_soapy( vtfile->bsize);
   set_scheduling();

   state = STATE_RESET;
   run();
   return 0;
}

