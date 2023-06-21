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

//
// Data structure for each input channel
//

static struct CHAN {

   double *ibuf;                 // Input buffer
   double *obuf;                 // Output buffer
   fftw_plan fwd, rev;   // Forward and reverse transforms
   fftw_complex *X;              // Fourier transform data  
}
 *chans;

static int nchans = 0;   // Number of input channels

static VTFILE *vtinfile, *vtoutfile;
static char *inname = NULL;
static char *outname = NULL;

static int sample_rate;   // Nominal sample rate
static double srcal;       // Sample rate correction factor

static int fftwid;
static int overlap;

static int CFLAG = FALSE;   // Set by -c option to do convolution

static void usage( void)
{
   fprintf( stderr,
       "usage:  vtmatch [options] [input [output]]\n"
       "\n"
       "options:\n"
       "  -v        Increase verbosity\n"
       "  -B        Run in background\n"
       "  -L name   Specify logfile\n"
       "\n"
       "  -t filename  Template data file\n"
       "\n"
       "  -F freqspec  Select frequency range (default all)\n"
       "\n"
       "  -c        Convolve input stream with template\n"
       "            (default is matched filter)\n"
     );
   exit( 1);
}

static char *template_file = NULL;   // File which supplies the template
static complex double *template_FT = NULL;  // Transform of template

static double *template = NULL;
static int ntemplate = 0;
static int ntemplate_alloc = 0;

static double FS = 0, FE = 0;   // Set by -F option

static void load_template( void)
{
   if (!template_file) VT_bailout( "must specify template file with -t");

   FILE *f = fopen( template_file, "r");
   if (!f) VT_bailout( "cannot open %s, %s", template_file, strerror( errno));

   char *p, temp[100];
   while (fgets( temp, 99, f))
   {
      p = strchr( temp, '\n');  if (p) *p = 0;
      p = strchr( temp, '\r');  if (p) *p = 0;
      if (!temp[0]) continue;

      if (ntemplate == ntemplate_alloc)
      {
         ntemplate_alloc += 100;
         template = VT_realloc( template, ntemplate_alloc * sizeof( double));
      }

      template[ntemplate++] = atof( temp);
   }

   VT_report( 1, "template loaded: %d samples", ntemplate);
   fclose( f);

   // Make template size even
   if (ntemplate % 2) ntemplate++;
}

static void process_buffer( void)
{
   int i, j;

   for (i = 0; i < nchans; i++)
   {
      struct CHAN *cp = chans + i;
      fftw_execute( cp->fwd);

      for (j = 0; j < fftwid; j++) cp->X[j] *= template_FT[j];

      fftw_execute( cp->rev);

      // Retain overlap input samples from end of current transform, move to
      // start of next transform

      memmove( cp->ibuf, 
               cp->ibuf+fftwid-overlap,
               overlap * sizeof( double));
   }
}

int main( int argc, char *argv[])
{
   VT_init( "vtmatch");

   int background = 0;

   while (1)
   {
      int c = getopt( argc, argv, "vBcL:t:F:?");
      
      if (c == 'v') VT_up_loglevel();
      else
      if (c == 'B') background = 1;
      else
      if (c == 'L') VT_set_logfile( "%s", optarg);
      else
      if (c == 't') template_file = strdup( optarg);
      else
      if (c == 'c') CFLAG = TRUE;
      else
      if (c == 'F') VT_parse_freqspec( optarg, &FS, &FE);
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

   load_template();

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

   VT_init_chanspec( chspec, vtinfile);
   nchans = chspec->n;

   vtoutfile = VT_open_output( outname, nchans, 0, sample_rate);
   if (!vtoutfile)
      VT_bailout( "cannot open %s: %s", outname, VT_error);

   chans = (struct CHAN *) VT_malloc_zero( sizeof( struct CHAN) * nchans);

   // Choose FT width to give about 10 or so transforms per second and at
   // least four times the template length
   int i;
   for (i = 4; ; i *= 2)
   {
      fftwid = i * ntemplate;
      if (fftwid > sample_rate/10) break;
   }

   //
   // Extend template to fftwid with zeros
   //
   if (ntemplate_alloc < fftwid)
      template = VT_realloc( template, fftwid * sizeof( double));
   for (i = ntemplate; i < fftwid; i++) template[i] = 0;

   overlap = ntemplate;

   VT_report( 1, "template width %d fftwid %d overlap %d",
                 ntemplate, fftwid, overlap);
   if (!CFLAG)
   {
      // Scale the template
      double sumsq = 0;
      for ( i = 0; i < fftwid; i++) sumsq += template[i] * template[i];
      double s = sqrt( sumsq);
      for (i = 0; i < fftwid; i++) template[i] /= s;
   }

   //
   // Transform template
   //
   template_FT =
      (complex double *) VT_malloc_zero( sizeof( complex double) * 2 * fftwid);

   fftw_plan fp = fftw_plan_dft_r2c_1d(
                     fftwid, template, template_FT, FFTW_ESTIMATE);
   fftw_execute( fp);
   for (i = 0; i < fftwid; i++) template_FT[i] /= fftwid;

   if (!CFLAG)
      for (i = 0; i < fftwid; i++) template_FT[i] = conj( template_FT[i]);

   //
   // If -F given, zero the template outside of the requested range
   //
   if (FS >= 0 && FE == 0) FE = sample_rate/2.0;
   if (FE < FS) VT_bailout( "invalid frequency range given with -F");
   double dF = sample_rate / (double) fftwid;
   for (i = 0; i < FS/dF; i++) template_FT[i] = 0;
   for (i = FE/dF; i < fftwid/2; i++) template_FT[i] = 0;

   //
   //  Set up the buffers for each channel
   //

   for (i = 0; i < nchans; i++)
   {
      struct CHAN *cp = chans + i;
      cp->ibuf = (double *) VT_malloc( fftwid * sizeof( double));
      cp->obuf = (double *) VT_malloc_zero( fftwid * sizeof( double));
      cp->X = (complex double *) VT_malloc( 2*fftwid * sizeof( complex double));
      cp->fwd = fftw_plan_dft_r2c_1d( fftwid, cp->ibuf, cp->X,
                FFTW_ESTIMATE | FFTW_PRESERVE_INPUT);
      cp->rev = fftw_plan_dft_c2r_1d( fftwid, cp->X, cp->obuf, FFTW_ESTIMATE);
   }

   int nb = 0;         // Input buffer load offset
   int ob = 0;         // Output buffer unload offset
   int nfft = FALSE;   // TRUE after the first Fourier transform

   double *outframe = VT_malloc( nchans * sizeof( double));

   int64_t icnt = 0;   // Count of input samples
   int64_t ocnt = 0;   // Count of output samples

   int eof = FALSE;
   timestamp T = timestamp_ZERO;
   int64_t tcnt = 0;

   while (1)
   {
      int e;
      if (!eof && (e = VT_is_block( vtinfile)) < 0)
      {
         VT_report( 1, "end of input");
         eof = TRUE;
      }

      if (!eof)
      {     
         if (e)  // Starting a new input block?
         {
            srcal = VT_get_srcal( vtinfile);
            T = VT_get_timestamp( vtinfile);
            tcnt = icnt;
         }
   
         double *inframe = VT_get_frame( vtinfile);
         for (i = 0; i < nchans; i++)
            chans[i].ibuf[nb] = inframe[chspec->map[i]];
         icnt++;
      }
      else
      {
         //  Input has ended.  Push zeros into the FFT to flush
         for (i = 0; i < nchans; i++) chans[i].ibuf[nb] = 0;
      }

      // Defer output until the first FFT has happened
      if (nfft)
      {
         if (!vtoutfile->nfb)  // Starting a new output block?
         {
            double offset = (ocnt - tcnt)/(srcal * sample_rate);
            VT_set_timebase( vtoutfile, timestamp_add( T, offset), srcal);
         }

         for (i = 0; i < nchans; i++) outframe[i] = chans[i].obuf[ob];
         VT_insert_frame( vtoutfile, outframe);
         ob++;
         ocnt++;
      }

      if (ocnt == icnt) break;
      if (++nb == fftwid)
      {
         process_buffer();
         nb = overlap;
         ob = 0;
         nfft = TRUE;
      }
   }

   VT_release( vtoutfile);

   return 0;
}

