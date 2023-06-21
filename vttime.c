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
static double in_srcal = 1.0;

static int timing_chan = 1;  // Selected input channel containing timing signal
                             // set by -c option

static int diag = FALSE;        // Set TRUE by -d option to turn on extra stuff

static double calibration_offset = 0;  // Seconds. Set by mode option c=

static double pulse_width = 0;         // Centroid search half-width, seconds
static int pulse_width_auto = FALSE;   // TRUE if pulse_width set automatically
static int noppm = FALSE;         // If TRUE, does not restrict peak/mean ratio
                                  // set by mode option noppm

static int64_t incnt = 0;                                // Input frame counter
static int state = 0;                                 // 0 = setup, 2 = running

static void (*process_timing)( timestamp, double) = NULL;

static double ppsmad = 0;
static int ppsmad_count = 0;
static double last_interval = 0;

static int holdover_limit = 0;  // Set by -h, maximum seconds for holdover
                                // Zero means holdover forever
static int holdover_count = 0;

static void diagnostic_td( double *buf, int N, int offset)
{
   if (!diag) return;

   double tstep = 1.0/sample_rate;

   FILE *f = fopen( "/tmp/vttime.td", "w");
   if (!f) return;

   int i;
   for (i = 0; i < N; i++)
      fprintf( f, "%d %.7e %.3f\n", i - offset, (i - offset) * tstep, buf[i]);
   fclose( f);
}

static int cmp_double( const void *p1, const void *p2)
{
   double v1 = *(double *)p1;
   double v2 = *(double *)p2;

   if (v1 < v2) return -1;
   if (v1 > v2) return 1;
   return 0;
}

///////////////////////////////////////////////////////////////////////////////
//  Timebases                                                                //
///////////////////////////////////////////////////////////////////////////////

// There are three distinct timestamps handled in this program:
//
// 1/ The timestamp supplied with the input stream.  This is used only to place
//    a capture buffer around the likely position of the 'second'.

// 2/ Our determination (hopefully more accurate) of the input stream's timing.

static timestamp our_timebase = timestamp_ZERO;
static int64_t inbase = 0;   // Input frame count to which our_timebase refers
static double our_srcal = 1.0;
static double our_dt = 0.0;

// 3/ The timestamp of the output stream, after resampling to the nominal sample
//    rate (srcal = 1.0).

static timestamp outbase = timestamp_ZERO;

static int64_t incnt_prev = 0;      // incnt_this from previous capture
static double secmark_prev = 0;     // The secmark_index of the previous pulse

//
//  Rate smoothing.
//

static double rs_TC = 10;
static int rs_once = FALSE;

static double rate_smooth( double in)
{
   static double s1 = 0, s2 = 0;
   static double alpha = 0;

   if (!rs_once)
   {
      s1 = s2 = in;
      alpha = 1 - exp( -1/rs_TC);
      rs_once = TRUE;
      return in;
   }

   s1 += alpha * (in - s1);
   s2 += alpha * (s1 - s2);
 
   return s2;
}

//
//  Timebase offset smoothing.
//

static double ts_TC = 10;
static int ts_once = FALSE;

static double time_smooth( double in)
{
   static double s1 = 0;
   static double alpha = 0;

   if (!ts_once)
   {
      s1 = in;
      alpha = 1 - exp( -1/ts_TC);
      ts_once = TRUE;
      return in;
   }

   s1 += alpha * (in - s1);
   double s2 = alpha * s1;

   s1 -= s2;
   return s2;
}

//
//  Timing queue.  This is a circular list of associations between timestamps
//  and input sample counts.  revise_timebase() inserts a new list entry each
//  second.  Used by the interpolator to determine the timestamp of each input
//  sample.
//

#define TQLEN 100
static struct TQ {
   timestamp our_timebase;
   int64_t inbase;
}
 tq[TQLEN];

static int tqlp = 0;   // Load pointer, used by revise_timebase()
static int tqup = 0;   // Unload pointer, used by the interpolator

//
//  Reset everything to start from scratch.
//
 
static void reset_timebase( void)
{
   state = 0;
   incnt = 0;
   inbase = 0;
   our_srcal = 1.0;
   our_dt = 0.0;
   our_timebase = timestamp_ZERO;
   secmark_prev = 0;
   incnt_prev = 0;
   ppsmad = 0;
   ppsmad_count = 0;
   last_interval = 0;
   rs_once = FALSE;
   ts_once = FALSE;
   tqlp = tqup = 0;
}

//
//  Function revise_timebase() is called after each PPS measurement.  It
//  smooths the intervals and updates our timebase.
//
//  secmark_index: offset in samples of the second mark in the capture buffer
//  incnt_this: input sample count of the first frame of the capture buffer
//  T_start: out timestamp of the first frame of the capture buffer
//  pmr: peak to mean ratio of the timing pulse
//

static void revise_timebase( double secmark_index,
                             int64_t incnt_this,
                             timestamp T_start,
                             double pmr)
{
   //
   //  Sanity check on the second mark index.  Can sometimes be called with
   //  a nan if there's a bad pulse.
   //

   if (isnan( secmark_index))
   {
      VT_report( 0, "invalid secmark, ignoring");
      holdover_count++;
      return;
   }

   //
   //  Calculate pulse interval: number of samples between this second mark
   //  and the previous. incnt_this is the input sample count of the start
   //  of the capture buffer.  secmark_index is the offset of the second mark
   //  into the capture buffer.
   //

   double interval = incnt_this - incnt_prev + secmark_index - secmark_prev;
   secmark_prev = secmark_index;
   incnt_prev = incnt_this;

   //
   //  Only use the raw interval if it is reasonable, say within 1% of the
   //  nominal sample rate.   A seriously bad soundcard crystal could
   //  permanently fail this sanity check.
   //

   if (interval > sample_rate * in_srcal * 1.001 ||
       interval < sample_rate * in_srcal * 0.999)
   {
      VT_report( 0, "wild PPS pulse interval %.6f - skipped", interval);
      holdover_count++;
      return;
   }

   //
   //  Track the mean absolute difference between consecutive PPS intervals.
   //

   #define SRTC  40       // Time constant (seconds) for MAD interval smoothing
   double smfac = exp(-1.0/SRTC);               // Sample rate smoothing factor

   double ad = 0;
   double new_ppsmad = 0;

   if (last_interval)
   {
      ad = fabs(interval - last_interval);
      if (!ppsmad) new_ppsmad = ad;
      else
         new_ppsmad = ppsmad * smfac + ad * (1-smfac);
   }

   //
   //  Sometimes when starting up the first or second interval is bad.
   //  Discard these and reset the timebase.  Otherwise it takes a while
   //  to exponentially recover.
   //

   if (!state && last_interval && fabs( interval - last_interval) > 0.1)
   {
      VT_report( 0, "unreliable PPS interval %.6f - skipped", last_interval);
      reset_timebase();
      return;
   }

   //
   //  If an interval is more than a few times the mean absolute difference
   //  away from the current average rate, then discard it.
   //

   if (state == 2 && ppsmad_count > 50 &&
       fabs(interval - sample_rate * our_srcal) > 5 * ppsmad)
   {
      VT_report( 0, "dubious PPS interval %.6f - skipped", interval);
      holdover_count++;
      return;
   }
  
   //
   //  PPS interval looks good, accept the interval and the new mean absolute
   //  difference.
   //

   ppsmad = new_ppsmad;
   ppsmad_count++;
   last_interval = interval;
 
   //
   //  Update our sample rate estimate.  This is just a moving average of
   //  recent intervals.  rate_smooth() implements a low pass filter
   //  which provides the averaging.
   //

   double current_rate = sample_rate * our_srcal;         
   double new_rate = rate_smooth( interval);
   double rate_error = new_rate - current_rate;  

   //
   //  Inform the interpolator about the new sample rate.
   //

   our_srcal = new_rate/sample_rate;
   our_dt = 1/(sample_rate * our_srcal);

   //
   //  Calculate the time correction.  The arrival time of the PPS is compared
   //  with the arrival time predicted by our timebase.
   //

   double time_error = 0;           // Raw timebase offset in samples
   double smoothed_time_error = 0;  // Low pass filtered time offset, samples 
   double time_correction = 0;      // The time correction to apply, seconds

   if (inbase)   // Not the first time through?
   {
      // Expected sample number of the second mark, based on our timebase
      double expected = sample_rate * our_srcal * 
                        timestamp_diff( timestamp_round(T_start), our_timebase);

      // Calculate timing error, units are number of samples
      time_error = incnt_this + secmark_index - expected - inbase;
      smoothed_time_error = time_smooth( time_error);

      // Time correction, units are seconds
      time_correction =  - smoothed_time_error * our_dt;
   }

   //
   //  Apply the timebase correction.
   //

   int timebase_valid = 0;
   static timestamp last_tstart = timestamp_ZERO;

   if (!inbase)  // First time through?
   {
      // First time through: initialise our timebase
      last_tstart = timestamp_round( T_start);
      our_timebase = timestamp_add( last_tstart, -secmark_index * our_dt);
      inbase = incnt_this;
   }
   else          // Subsequent times through
   {
      // Increment our timebase by about a second and add in the time correction
      int secs = round( timestamp_diff( T_start, last_tstart));
      int nadj = round( sample_rate * our_srcal * secs);
      inbase += nadj;

      our_timebase = timestamp_add( our_timebase, 
                                    nadj * our_dt + time_correction);

      // Add the timebase correction to every entry in the timing queue
      int i;
      for (i = 0; i < TQLEN; i++) tq[i].our_timebase = 
                   timestamp_add( tq[i].our_timebase, time_correction);

      last_tstart = timestamp_round( T_start);
      timebase_valid = 1;
   }

   //
   //  Append the current timebase to the timing queue.
   //

   tq[tqlp].our_timebase = our_timebase;
   tq[tqlp].inbase = inbase;
   tqlp = (tqlp + 1) % TQLEN;

   //
   //  Report a bunch of stuff to the log stream.
   //
   //  in_err is just for reporting the timing error of the input stream, it
   //  says how well ntpd is doing on the machine running vtcard.
   //

   double in_err = timestamp_diff( T_start, our_timebase)
                    + (inbase - incnt_this)  * our_dt;

   VT_report( 1, "st%d PPSpmr %.1f PPSmad %.3fuS in %.6fmS "
                    "out %7.3fuS "
                    "rate_err %6.3f inrate %.4f int %.4f tc %.3f", 
       state, pmr, 1e6 * ppsmad/sample_rate,
       1e3 * in_err,
       1e6 * smoothed_time_error/sample_rate,
       (double) rate_error,
       new_rate, (double) interval,
       time_correction/our_dt);

   //
   // If the output stream timestamp error reaches half a sample, then reset
   // everything and re-establish timing from scratch.  The rate error is
   // normally a small (< 0.05) fraction of a sample so 0.5 is pretty bad and
   // usually occurs when the soundcard sample rate or the PCs RTC is drifting
   // too fast to correct, or there has been a step change of timing.
   //

   if (fabs( time_correction/our_dt) > 0.5)
   {
      VT_report( 0, "input rate drifting too fast");
      if (state) reset_timebase();
      return;
   }

   //
   // Set up the timing baseline required for output resampling.   This occurs
   // only the first time through since the output timebase is constant.
   //

   if (!state && timebase_valid && fabs( rate_error) < 0.05)
   {
      our_timebase = timestamp_add( our_timebase, -time_error * our_dt);

      // Set the output to begin on the next second
      outbase = timestamp_add( timestamp_round(T_start), 1);
      VT_set_timebase( vtoutfile, outbase, 1.0);

      // Switch to run state.  The main loop will start to output samples
      // when the time reaches outbase.  outbase has been set to the next
      // second which means the main loop will discard the remainder of the
      // current second.

      state = 2;  
   }

   holdover_count = 0;
}

///////////////////////////////////////////////////////////////////////////////
//  PPS Waveform Filter                                                      //
///////////////////////////////////////////////////////////////////////////////

//
//  A low pass filter activated by lp= option in centroid mode.
//

static double pps_lp_freq = 0;    // Set by lp= option
static int pps_lp_poles = 0;      // Set by poles= option

static double a1, a2, b0, b1, b2; // Filter coefficients

#define MAX_SECT 6  // Max number of 2-pole filter sections

static struct SECT {
   double w1, w2;
} sects[MAX_SECT];   // State variables for each filter section

static double pps_lp_filter( double v)
{
   int NS = pps_lp_poles / 2;  // Number of 2-pole filter sections

   // Direct form 2 filter

   double filter( struct SECT *s, double x0)
   {
      double w0 = x0 - a1 * s->w1 - a2 * s->w2;
      double y0 = b0 * w0 + b1 * s->w1 + b2 * s->w2;
      s->w2 = s->w1; s->w1 = w0;
      return y0;
   }

   int i;
   for (i = 0; i < NS; i++) v = filter( sects + i, v);
   return v;
}

static void pps_filter_init( void)
{
   double T = 1/(double) sample_rate;

   double Wc = 2 * M_PI * pps_lp_freq;  // Requested cut-off frequency
   double Wa = 2/T * tan(Wc * T/2);     // Warped frequency
   double tau = 1/Wa;                   // Warped time constant

   double d = T*T + 4 * tau *T + 4 * tau * tau;
   b2 = T*T/d;
   b1 = 2*T*T/d;
   b0 = T*T/d;
   a1 = (2 * T*T - 8 * tau * tau)/d;
   a2 = (T*T - 4 * tau * T + 4 * tau * tau)/d;
}

//
//  A high pass filter.
//

static double hpRC = 0.01/(2 * M_PI);
static double pps_hp_filter( double v)
{
   static double x = 0;
   static double y = 0;

   double alpha = hpRC/(hpRC + 1.0/sample_rate);

   y = alpha * y + alpha * (v - x);
   x = v;
   return y;
}

static double pps_filter( double v)
{
   return pps_lp_filter( pps_hp_filter( v));
}

///////////////////////////////////////////////////////////////////////////////
//  Baseband PPS                                                             //
///////////////////////////////////////////////////////////////////////////////

//
//  Pulse buffer.
//

static double *cbuf = NULL;     // Capture buffer
static int cbuflen = 0;         // Capture buffer length, samples
static int pulse_sign = 1;      // +1 or -1, polarity expected of the PPS 

static double hw = 0.45;        // Half width of capture buffer
static double vmax = 0;         // Peak value of PPS channel
static int vmax_pos;            // Offset into capture buffer of peak location
static double pmr = 0;          // PPS peak to mean ratio

static int64_t incnt_this = 0;    // Frame count of first frame in capture buf
static timestamp T_start = timestamp_ZERO;   // Incoming timestamp of first
                                             // frame in capture buffer

static int load_buffer( timestamp T, double in)
{
   static int cbufp = 0;    // Capture buffer load index
   static int flag = FALSE; // TRUE while we're loading into the capture buffer
   static double vavg = 0;  // Average signal value of PPS channel

   if (!cbuflen)   // First time through?   Do some initialisation.
   {
      cbuflen = hw * 2 * sample_rate;
      cbuf = VT_malloc( sizeof( double) * cbuflen);
      our_srcal = in_srcal;
   }

   //
   //  As we reach within hw seconds of the second mark, start capturing 
   //

   double t = timestamp_frac( T);  // Fractional part of the input timestamp
 
   if (!flag && t > 1-hw)   // Time to start capturing?
   {
      //  Initialise capture buffer
      cbufp = 0;           // Capture buffer load index
      vmax = vavg = 0;     // Average and max accumulators
      incnt_this = incnt;  // Input sample number at start of capture buffer
      T_start = T;         // Input stream timestamp at start of capture buffer
      flag = TRUE;         // Capture buffer is loading
   }
  
   if (!flag) return FALSE;   // Not yet capturing?

   //
   //  Load this sample into the capture buffer, accumulate max and mean,
   //  and the location within the capture buffer of the peak.
   //

   double v = pulse_sign * in;

   if (v > vmax)
   {
      vmax = v;
      vmax_pos = cbufp;
   }

   vavg += fabs( v);
   cbuf[cbufp++] = v;

   if (cbufp < cbuflen) return FALSE;    // Capture buffer not yet full?

   flag = FALSE;      // End of capture, buffer is full
   vavg /= cbuflen;   // Average amplitude

   //
   //  Check pulse amplitude is sufficient, peak/mean ratio of buffer contents.
   //

   pmr = vmax / vavg;  // Peak to mean ratio
   if (!noppm && pmr < 30)
   {
      // Will occur if the PPS fails, or input timestamping is way out
      VT_report( 0, "insufficient PPS peak/mean ratio %.1f - skipped", pmr);
      holdover_count++;
      return FALSE;
   }

   return TRUE;
}

///////////////////////////////////////////////////////////////////////////////
//  Pulse Timing                                                             //
///////////////////////////////////////////////////////////////////////////////

//
//  Pulse timing: PPS pulse width less than or equal to two sound card
//  sample periods.  
//

static void timing_pulse( timestamp T, double in)
{
   int i;

   if (!load_buffer( T, in)) return;

   static int fftwid = 128;
   static fftw_plan ffp;
   static complex *X = NULL;
   static double *inbuf = NULL;
   static int b1, b2;   // First and last bins to use for phase slope
   static double *tsbuf = NULL;
   static double *window = NULL;

   if (!X)    // First time through initialisation
   {
      X = VT_malloc( sizeof( fftw_complex) * (fftwid/2 + 1));
      inbuf = VT_malloc( sizeof( double) * fftwid);
      ffp = fftw_plan_dft_r2c_1d( fftwid, inbuf, X, FFTW_ESTIMATE);
      b1 = fftwid/2 * 0.1;
      b2 = fftwid/2 * 0.6;
   
      tsbuf = VT_malloc( sizeof( double) * fftwid/2);
 
      window = VT_malloc( sizeof( double) * fftwid);
      for (i = 0; i < fftwid; i++)
         window[i] = pow( sin(M_PI * i/fftwid), 8);
   }

   if (pulse_width)
   {
      int i1 = vmax_pos - pulse_width * sample_rate;
      if (i1 < 0) i1 = 0;
      int i2 = vmax_pos + pulse_width * sample_rate;
      if (i2 > cbuflen) i2 = cbuflen;
      for (i = 0; i < i1; i++) cbuf[i] = 0;
      for (i = i2; i < cbuflen; i++) cbuf[i] = 0;
   }

   //
   //  Find the leading edge of the pulse in the capture buffer. Backtrack from
   //  the peak to find where the samples straddle 50% of the peak amplitude.
   //

   int vstart;
   for (vstart = vmax_pos;
        vstart >= 0 && cbuf[vstart] > 0.5 * vmax;
        vstart--) ;;;

   double dF = sample_rate/(double) fftwid;

   if (vstart < fftwid/2)
   {
      VT_report( 0, "pulse starts too early - skipped");
      holdover_count++;
      return;
   }

   //
   //  Extract the pulse for Fourier transform, centering the transform on
   //  the vstart sample.
   //

   for (i = 0; i < fftwid; i++) inbuf[i] = cbuf[vstart - fftwid/2 + i];
   diagnostic_td( inbuf, fftwid, fftwid/2);

   for (i = 0; i < fftwid; i++) inbuf[i] *= window[i];
   fftw_execute( ffp);

   for (i = 0; i < fftwid/2; i += 2) X[i] = -X[i];

   //
   //  Find the median phase difference between adjacent bins and take
   //  that as the phase slope per bin.
   //

   int n = 0;
   for (i = b1; i < b2; i++)
   {
      double dp = carg( X[i+1]/X[i]);
      tsbuf[n++] = dp;
   }

   qsort( tsbuf, n, sizeof( double), cmp_double);

   //
   //  Phase slope should be negative because the center of the FFT is placed
   //  at the sample before the pulse.
   //

   double slope = tsbuf[n/2]/(2 * M_PI);   // Cycles per bin

   if (slope > 0)
   {
      VT_report( 0, "invalid phase slope - skipped");
      holdover_count++;
      return;
   }

   double delay = slope/dF;           // Cycles per Hz

   VT_report( 2, "slope %.3e delay %.3f uS", slope, delay * 1e6);

   double secmark_index =
      vstart - (delay + calibration_offset) * sample_rate * our_srcal;
   revise_timebase( secmark_index, incnt_this, T_start, pmr);
}

///////////////////////////////////////////////////////////////////////////////
//  Centroid Timing                                                          //
///////////////////////////////////////////////////////////////////////////////

//
//  Pulse centroid timing:  Measures the centroid of a RC-shaped pulse.
//

static double set_w1 = 0, set_w2 = 0;

static void timing_centroid( timestamp T, double in)
{
   //
   //  Filter the PPS sample and load into the input buffer.
   //  Return if buffer not yet full.
   //

   in = pps_filter( in);
   if (!load_buffer( T, in)) return;

   //
   //  Determine centroid of the pulse whose peak is at index vmax_pos
   //  in the capture buffer.  aw1 and aw2 are buffer offsets which bound
   //  the centroid integration.
   //

   int aw1, aw2;
   if (set_w1 && set_w2)  // w1 and w2 given in the mode options?
   {
      aw1 = vmax_pos - set_w1 * sample_rate;
      aw2 = vmax_pos + set_w2 * sample_rate;
   }
   else
   if (!pulse_width_auto)    // Not using w=auto mode option?
   {
      // pulse_width (actually the half-width) is set by w= mode option
      aw1 = vmax_pos - pulse_width * sample_rate;
      aw2 = vmax_pos + pulse_width * sample_rate;
   }
   else  // Using w=auto mode option
   {
      aw1 = aw2 = vmax_pos;
      while (aw1 > 0 && cbuf[aw1] > vmax/1000) aw1--;
      while (aw2 < cbuflen-1 && cbuf[aw2] > vmax/1000) aw2++;
   }

   if (aw1 <= 0 || aw2 >= cbuflen-1)
   {
      // Input timing is so far out that the PPS doesn't sit completely in
      // the capture buffer - beginning or end of the pulse has been lost.
      VT_report( 0, "PPS at timing limit - skipped");
      holdover_count++;
      return;
   }

   diagnostic_td( cbuf + aw1, aw2 - aw1, vmax_pos - aw1);

   if (pulse_width_auto)
      VT_report( 2, "w1=%.3fmS w2=%.3fmS n=%d",
               (vmax_pos - aw1)/(double) sample_rate * 1000,
               (aw2 - vmax_pos)/(double) sample_rate * 1000, aw2 - aw1);

   //
   //  Centroid calculation.
   //

   double mprod = 0, msum = 0;
   int j;
   for (j = aw1; j < aw2; j++)
   {
      double v = cbuf[j];

      msum += v;
      mprod += v * j;
   }

   double centroid = mprod/msum;  // Index of centroid in capture buffer

   //
   //  Move back from the centroid to the timing mark of the pulse using the
   //  centroid offset from the command line.
   //

   double secmark_index =
      centroid - calibration_offset * sample_rate * our_srcal;

   revise_timebase( secmark_index, incnt_this, T_start, pmr);
}

///////////////////////////////////////////////////////////////////////////////
//  No Timing                                                                //
///////////////////////////////////////////////////////////////////////////////

//
//  With -m none, the program does interpolation to UT synchronous samples at
//  the exact sample rate, relying entirely on the incoming stream's timestamp.
//

static void timing_none( timestamp T, double in)
{
   static timestamp T_last = timestamp_ZERO;

   our_srcal = in_srcal;
   our_dt = 1/(sample_rate * our_srcal);

   if (!state)
   {
      our_timebase = T;
      inbase = 0;
      incnt = 0;
      outbase = timestamp_compose( timestamp_secs( our_timebase) + 2, 0);
      VT_set_timebase( vtoutfile, outbase, 1.0);
      state = 2;
   }
   else inbase++;

   if (timestamp_diff( T, T_last) >= 1)
   {
      tq[tqlp].our_timebase = T;
      tq[tqlp].inbase = inbase;
      tqlp = (tqlp + 1) % TQLEN;

      T_last = T;
   }
}

///////////////////////////////////////////////////////////////////////////////
//  Interpolation                                                            //
///////////////////////////////////////////////////////////////////////////////

//
//  Conversion to exact sample rate and UT synchronous samples is done by
//  sinc-weighted interpolation between the input samples.
//

static double **Iframe;                // Input signal buffer - array of frames
static timestamp *Istamp;                             // Input timestamp buffer

#define ILOG 8                         // Log base 2 of the input buffer length
#define Ilen  (1 << ILOG)                                // Input buffer length
#define Imask (Ilen - 1)                                  // Index cycling mask
static unsigned char Ip = 0;                              // Input buffer index 

// Macros to load incoming frames and timestamps
#define IFRAME_IN(A) Iframe[(Ip + (A) - (1<<(ILOG-1))) & Imask]
#define ISTAMP_IN(A) Istamp[(Ip + (A) - (1<<(ILOG-1))) & Imask]

// Macros to extract frames and timestamps.  The zero index is positioned
// half way along the buffer so that (A) can range -Ilen/2 to +Ilen/2
#define IFRAME(A) Iframe[(Ip + (A)) & Imask]
#define ISTAMP(A) Istamp[(Ip + (A)) & Imask]

static float *wsinc;                   // Interpolation kernel, a sinc function
static float **wsinck;                 // See setup_wsinc()

#define WNZ  36       // Half-length of sinc function, number of zero crossings
#define WNL  128            // Number of function points between zero crossings

static double *outframe;

//
//  Construct an output frame located in time between IFRAME(0) and IFRAME(1).
//  x (0 <= x < 1) is the fractional distance in time between the output
//  sample and IFRAME(0).   wsinc[] holds half of the sinc function, the other
//  half is just mirrored.   WNL is the number of sinc function entries per
//  output sampling period.
//

//        IFRAME(-1)     IFRAME(0)        IFRAME(1)       IFRAME(2)
// Input  ----X-------------X----------------X----------------X-------  ...
//                          |<  x >|
// Output ---------x---------------X---------------x----------------x--
//                 ^               ^               ^                ^
//             wsinc[-WNL]      wsinc[0]       wsinc[+WNL]      wsinc[+2*WNL]
// wsinc                    |< km >|
//            |             |                |                |
//      wsinc[-km-WNL]   wsinc[-km]       wsinc[-km+WNL]    wsinc[-km+2*WNL]

static inline void interpolate( int chans, double x)
{
   int i, j;

   int km = x * WNL;

   float *wm = wsinck[km];

   // This uses a lot more cpu if you swap the inner and outer loops
   for (i = 0; i < chans; i++)
   {
      double v = 0;

      for (j = -(WNZ-1); j < WNZ; j++) v += IFRAME(-j)[i] * wm[j];

      outframe[i] = v;
   }
}

//
//  Note from the above ASCII art, wsinc[] is used in steps of WNL all having
//  an offset of km.  Auxilliary array wsinck[km][i*WNL] reduces CPU cache
//  misses by listing wsinc[km + i*WNL] values these for all possible values
//  of km.
//
//  wsinck[km][0] is positioned in the center of the allocated space so that
//  signed values of i can be used.
//
 
static void setup_wsinc( void)
{
   int i, k, wsinc_len = WNZ * WNL + 1;

   wsinc = VT_malloc( sizeof(float) * wsinc_len);
   wsinc[0] = 1;
   for (i = 1; i < wsinc_len; i++)
   {
      double a = i/(double) WNL;
      wsinc[i]  = sin(M_PI * a)/(M_PI * a);
   }

   wsinck = VT_malloc( sizeof( float *) * (WNL+1));
   for (k = 0; k <= WNL; k++)
   {
      wsinck[k] = WNZ + (float *) VT_malloc( sizeof( float) * WNZ * 2);

      for (i = -WNZ; i < WNZ; i++)
      {
         int n = k + i * WNL;  if (n < 0) n = -n;
         wsinck[k][i] = wsinc[n];
      }
   }
}

///////////////////////////////////////////////////////////////////////////////
//  Parse Method Options                                                     //
///////////////////////////////////////////////////////////////////////////////

static void parse_method( char *s)
{
   // centroid+,w=width,c=offset
   // centroid-,w=width,c=offset
   // edge+,c=offset
   // edge-,c=offset

   int w_given = 0, c_given = 0;
  
   while (s && *s)
   {
      char *p = strchr( s, ','); if (p) *p++ = 0;

      if (!strcmp( s, "none"))
      {
         process_timing = timing_none;
      }
      else
      if (!strcmp( s, "ppsbase+") ||
          !strcmp( s, "centroid+"))
      {
         process_timing = timing_centroid;
         pulse_sign = 1;
      }
      else
      if (!strcmp( s, "ppsbase-") ||
          !strcmp( s, "centroid-"))
      {
         process_timing = timing_centroid;
         pulse_sign = -1;
      }
      else
      if (!strcmp( s, "edge+"))
      {
         process_timing = timing_pulse;
         pulse_sign = 1;
         noppm = TRUE;
      }
      else
      if (!strcmp( s, "edge-"))
      {
         process_timing = timing_pulse;
         pulse_sign = -1;
         noppm = TRUE;
      }
      else
      if (!strcmp( s, "pulse+"))
      {
         process_timing = timing_pulse;
         pulse_sign = 1;
      }
      else
      if (!strcmp( s, "pulse-"))
      {
         process_timing = timing_pulse;
         pulse_sign = -1;
      }
      else
      if (!strncmp( s, "c=", 2))
      {
         calibration_offset = atof( s+2);
         c_given = 1;
      }
      else
      if (!strncmp( s, "w=", 2))
      {
         if (!strncmp( s+2, "auto", 4)) pulse_width_auto = TRUE;
         else
         {
            pulse_width = atof( s+2);
            pulse_width_auto = FALSE;
            w_given = 1;
         }
      }
      else
      if (!strncmp( s, "lp=", 3))
      {
         pps_lp_freq = atof( s+3);
         if (!pps_lp_poles) pps_lp_poles = 6;
      }
      else
      if (!strncmp( s, "poles=", 6))
      {
         pps_lp_poles = atoi( s + 6);
      }
      else
      if (!strncmp( s, "w1=", 3)) set_w1 = atof( s+3);
      else
      if (!strncmp( s, "w2=", 3)) set_w2 = atof( s+3);
      else
      if (!strncmp( s, "noppm", 5)) noppm = TRUE;
      else
      if (!strncmp( s, "rs=", 3)) rs_TC = atof( s+3);
      else
      if (!strncmp( s, "ts=", 3)) ts_TC = atof( s+3);
      else
         VT_bailout( "unrecognised method option: %s", s);

      s = p;
   }

   // Apply some defaults with centroid timing
   if (process_timing == timing_centroid)
   {
      if (!w_given && c_given) pulse_width = 1.1 * calibration_offset;
      else
      if (!w_given)
      {
         pulse_width = 1e-3;
         pulse_width_auto = TRUE;
      }
   }

   if (pps_lp_poles % 2) VT_bailout( "poles= option must be even");

   if (pps_lp_poles && !pps_lp_freq)
      VT_bailout( "poles= given without lp=");

   if (pps_lp_freq < 0) VT_bailout( "invalid lp frequency");

   if (pps_lp_poles > MAX_SECT * 2)
   {
      VT_report( 0, "poles= value too high, reduced to 12");
      pps_lp_poles = 12;
   }
}

static void usage( void)
{
   fprintf( stderr,
        "usage: vttime [options] input output\n"
        "\n"
        "options:\n"
        "  -v       Increase verbosity\n"
        "  -B       Run in background\n"
        "  -L name  Specify logfile\n"
        "  -c chan  Channel containing timing signal\n"
        "           (default 1)\n"
        "\n"
        "  -m method,options  PPS timing method\n"
        "     -m centroid+,c=centroid,w=halfwidth\n"
        "     -m centroid-,c=centroid,w=halfwidth\n"
        "     -m edge+,c=offset\n"
        "     -m edge-,c=offset\n"
        "     -m pulse+,c=offset\n"
        "     -m pulse-,c=offset\n"
        "     -m none\n" 
     );

   exit( 1);
}

int main( int argc, char *argv[])
{
   VT_init( "vttime");

   int background = 0;

   while (1)
   {
      int c = getopt( argc, argv, "vdBc:m:L:h:?");

      if (c == 'v') VT_up_loglevel();
      else
      if (c == 'B') background = 1;
      else
      if (c == 'L') VT_set_logfile( "%s", optarg);
      else
      if (c == 'm') parse_method( strdup( optarg));
      else
      if (c == 'c') timing_chan = atoi( optarg);
      else
      if (c == 'd') diag = TRUE;
      else
      if (c == 'h') holdover_limit = atoi( optarg);
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

   if (!process_timing) VT_bailout( "no timing method specified: needs -m");

   VT_report( 1, "rate smoothing: %.1f seconds", rs_TC);
   VT_report( 1, "time smoothing: %.1f seconds", ts_TC);

   struct VT_CHANSPEC *chspec = VT_parse_chanspec( inname);
   vtinfile = VT_open_input( inname);
   if (!vtinfile)
      VT_bailout( "cannot open input %s: %s", inname, VT_error);

   sample_rate = VT_get_sample_rate( vtinfile);
   our_dt = 1/(double) sample_rate;   // Input sample interval
   long double DT = 1.0/sample_rate;  // Output sample interval - fixed

   VT_init_chanspec( chspec, vtinfile);

   VT_report( 1, "channels: %d, input rate: %d", chspec->n, sample_rate);
   VT_report( 1, "calibration offset: %.3e", calibration_offset);

   if (timing_chan < 1 || timing_chan > chspec->n)
      VT_bailout( "invalid channel %d given with -c", timing_chan);

   vtoutfile = VT_open_output( outname, chspec->n, 0, sample_rate);
   if (!vtoutfile) VT_bailout( "cannot open: %s", VT_error);

   if (pps_lp_poles) pps_filter_init();  // Set up PPS LP filter

   //
   //  Initialise interpolation buffer.
   //

   int i;
   Iframe = VT_malloc_zero( sizeof( double *) * Ilen);
   Istamp = VT_malloc_zero( sizeof( timestamp) * Ilen);

   for (i =0 ; i < Ilen; i++)
      Iframe[i] = VT_malloc_zero( sizeof( double) * chspec->n);
   setup_wsinc();    //  Prepare interpolation coefficients

   outframe = VT_malloc_zero( sizeof( double) * chspec->n);

   //
   //  Setup up the delay buffer.  This is a circular buffer which holds the
   //  input stream for a duration roughly equal to the response lag of the
   //  time/rate averaging.
   //

   int dbuflen = sample_rate * 10;
   double *dbuf = VT_malloc_zero( sizeof( double) * chspec->n * dbuflen);
   int dbufn = 0;                       // Number of frames in the delay buffer
   int dbufp = 0;                             // Delay buffer load/unload index

   //
   //  Main processing loop.
   //

   while (1)
   {
      //
      //  Get the input stream timestamp and check of end of input.
      //

      timestamp in_timebase = VT_get_timestamp( vtinfile);
      if (timestamp_is_NONE( in_timebase)) VT_exit( "end of input");
  
      in_srcal = VT_get_srcal( vtinfile); 

      if (vtinfile->rbreak)
         VT_report( 2, "timing break on input stream");

      //
      //  Read an incoming frame and call the timing processing
      //

      double *inframe = VT_get_frame( vtinfile);
      double pps_in = inframe[chspec->map[timing_chan - 1]];

      process_timing( in_timebase, pps_in);
      if (holdover_limit && holdover_count >= holdover_limit)
      {
         holdover_count = 0;
         VT_report( 0, "holdover limit reached, resetting");
         reset_timebase();
      }

      //
      //  Take out the oldest frame from the delay buffer and replace it with
      //  the new frame and move the buffer pointer along.   The old frame is
      //  moved to the head of the interpolation buffer.
      //

      double *dp = dbuf + dbufp * chspec->n;
      for (i = 0; i < chspec->n; i++)
      {
         IFRAME_IN(0)[i] = dp[i];
         dp[i] = inframe[chspec->map[i]];
      }
      dbufp = (dbufp + 1) % dbuflen;

      if (dbufn < dbuflen)   // Delay buffer not yet full?
      {
         dbufn++;
         incnt++;
         continue;
      }

      //
      //  Decide which entry in the timing queue to use.
      //

      int64_t ocnt = incnt - dbuflen;    // Output frame number
      int64_t offset = 0;
      while ((offset = ocnt - tq[tqup].inbase) > sample_rate/2 &&
             (i = (tqup+1) % TQLEN) != tqlp)  tqup = i;

      //
      //  Assign our timestamp to the head of the interpolation pipeline.
      //

      ISTAMP_IN(0) = timestamp_add( tq[tqup].our_timebase,
                                    our_dt * offset);

      //
      //  Generate output frames if we're in running state.
      //

      if (state == 2)               // In running state?
      {
         // Tout is the timestamp of the next output sample. 
         // ISTAMP(0) is the timestamp of the previous input sample,
         // ISTAMP(1) the next. Tout is somewhere in between the two.

         timestamp Tout = timestamp_add( outbase, vtoutfile->nft * DT);

         // Output as many frames as possible: 0, 1, or 2 frames are
         // interpolated between ISTAMP(0) and ISTAMP(1).

         double Tint = timestamp_diff(ISTAMP(1), ISTAMP(0));
         if (Tint > 0)
         {
            double x = timestamp_diff(Tout, ISTAMP(0))/Tint;
          
            while (x < 1)   // Until ISTAMP(1) is too old
            {
               interpolate( chspec->n, x < 0 ? 0 : x);
                         
               VT_insert_frame( vtoutfile, outframe);
   
               x += 1/(sample_rate * Tint);
            }
         }
      }

      Ip++;      // Shift the interpolation buffer
      incnt++;   // Input frame counter
   }
}


