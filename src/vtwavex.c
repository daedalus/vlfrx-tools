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

static VTFILE *vtfile;
static char *inname = NULL;
static char *outname = NULL;

static int wav_sample_rate = 0;
static double inf1_sample_rate = 0;
static timestamp inf1_T = timestamp_NONE;
static FILE *fin = NULL;
static int chans = 0;
static uint16_t wav_format = 0;
static int wav_frame_size = 0;
static int bits_per_sample = 0;

static double override_rate = 0;
static double srcal = 1.0;               // Sample rate calibration coefficient
static unsigned int sample_rate = 0;                     // Nominal sample rate

static timestamp Tstart = timestamp_NONE;    // From -T option

static void (*cvt_frame)(uint8_t *, double *) = NULL; // Frame conversion function

//
//  Info on wav format:
//    http://www-mmsp.ece.mcgill.ca/Documents/AudioFormats/WAVE/WAVE.html
//

static void decode_fmt_chunk( char *buff, size_t n)
{
   VT_report( 1, "fmt chunk size %d", (int) n);

   wav_format = *(uint16_t *) (buff + 0);
   VT_report( 1, "format code %04X", wav_format);

   if (wav_format == 0x6) VT_bailout( "cannot read A-law wav files");
   if (wav_format == 0x7) VT_bailout( "cannot read u-law wav files");
   if (wav_format != 0x1 && wav_format != 0x3 && wav_format != 0xFFFE)
      VT_bailout( "unrecognised format code %04X", wav_format);

   chans = *(uint16_t *)(buff + 2);
   VT_report( 1, "channels %d", chans);

   wav_sample_rate = *(uint32_t *)(buff + 4);
   VT_report( 1, "wav sample rate %d", wav_sample_rate);

   uint32_t byte_rate = *(uint32_t *)(buff + 8);
   VT_report( 1, "byte rate %d", byte_rate);

   wav_frame_size = *(uint16_t *)(buff + 12);
   VT_report( 1, "frame size %d", wav_frame_size);

   bits_per_sample = *(uint16_t *)(buff + 14);
   VT_report( 1, "bits per sample %d", bits_per_sample);

   if (n > 16)
   {
      uint16_t ext_size = *(uint16_t *)(buff + 16);
      VT_report( 1, "ext size %d", ext_size);

      if (wav_format == 0xFFFE && ext_size == 22)
      {
         wav_format = *(uint16_t *)(buff + 24);
         VT_report( 1, "ext format %04X", wav_format);
      }
   }
}

static void open_wav( void)
{
   if (!strcmp( inname, "-")) fin = stdin;
   else
      fin = fopen( inname, "rb");

   if (!fin) VT_bailout( "cannot open %s, %s", inname, strerror( errno));

   char buff[500];

   if (fread( buff, 1, 12, fin) != 12)
      VT_bailout( "input file read error (1) %s", strerror(errno));

   if (strncmp( buff, "RIFF", 4))
      VT_bailout( "input file incorrect format (1)");

   if (strncmp( buff+8, "WAVE", 4))
      VT_bailout( "input file incorrect format (2)");

   // Parse following chunks looking for fmt, inf1 and data

   while (1)
   {
      // Read the chunk header
      if (fread( buff, 1, 8, fin) != 8)
         VT_bailout( "input file read error %s", strerror(errno));

      size_t n = *(uint32_t *)(buff + 4);

      if (!strncmp( buff, "fmt", 3))
      {
         if (fread( buff, 1, n, fin) != n)
            VT_bailout( "input file read error (2) %s", strerror(errno));
         decode_fmt_chunk( buff, n);
         continue;
      }

      //
      //  If the WAV file came from Spectrum Lab, it may have an 'inf1' chunk
      //
      if (!strncmp( buff, "inf1", 4))
      {
         if (fread( buff, 1, n, fin) != n)
            VT_bailout( "input file read error (3): %s", strerror( errno));
         buff[n] = 0;
         VT_report( 2, "inf1 [%s]", buff);

         char *p, *q;

         for (p = buff; p && *p; )
         {
            q = strchr( p, ' '); if (q) *q++ = 0;

            if (!strncmp( p, "ut=", 3)) inf1_T = VT_parse_timestamp( p + 3);
            else
            if (!strncmp( p, "sr=", 3)) inf1_sample_rate = atof( p + 3);
            p = q;
         }

         continue;
      }

      if (!strncmp( buff, "data", 4))
      {
         VT_report( 1, "found data chunk");
         break;
      }

      VT_report( 1, "ignoring chunk [%s] size %d", buff, (int) n);
      if (fread( buff, 1, n, fin) != n)
         VT_bailout( "input file read error (4): %s", strerror( errno));
   }
}

static void cvt_frame_uint8( uint8_t *buff, double *frame)
{
   int i;
   for (i = 0; i < chans; i++) frame[i] = (uint8_t) buff[i] - 127;
}

static void cvt_frame_int16( uint8_t *buff, double *frame)
{
   int i;
   for (i = 0; i < chans; i++)
      frame[i] = *(int16_t *) (buff + i*2)/(double)INT16_MAX;
}

static void cvt_frame_int32( uint8_t *buff, double *frame)
{
   int i;
   for (i = 0; i < chans; i++)
      frame[i] = *(int32_t *) (buff + i*4)/(double)INT32_MAX;
}

static void cvt_frame_float32( uint8_t *buff, double *frame)
{
   int i;
   for (i = 0; i < chans; i++) frame[i] = *(float *) (buff + i*4);
}

static void cvt_frame_float64( uint8_t *buff, double *frame)
{
   int i;
   for (i = 0; i < chans; i++) frame[i] = *(double *) (buff + i*8);
}

static void usage( void)
{
   fprintf( stderr,
       "usage:  vtwavex [options] [input [output]]\n"
       "\n"
       "options:\n"
       "  -v            Increase verbosity\n"
       "  -B            Run in background\n"
       "  -L name       Specify logfile\n"
       "\n"
       " -r sample/sec  Override WAV header sample rate\n"
       " -T timestamp   Specify start time (default WAV inf1 field)\n"
     );

   exit( 1);
}

int main( int argc, char *argv[])
{
   VT_init( "vtwavex");

   int background = 0;

   while (1)
   {
      int c = getopt( argc, argv, "vBL:r:T:?");

      if (c == 'v') VT_up_loglevel();
      else
      if (c == 'B') background = 1;
      else
      if (c == 'L') VT_set_logfile( "%s", optarg);
      else
      if (c == 'r') override_rate = atof( optarg);
      else
      if (c == 'T') Tstart = VT_parse_timestamp( optarg);
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

   open_wav();

   if (!timestamp_is_NONE( inf1_T))
   {
      char temp[50];
      VT_format_timestamp( temp, inf1_T);
      VT_report( 1, "inf1 timestamp %s", temp);
   }

   if (override_rate)
   {
      VT_report( 1, "override sample rate %.6f", override_rate);
      sample_rate = round( override_rate);
      if (sample_rate) srcal = override_rate/sample_rate;
      else
      {
         sample_rate = 1;
         srcal = override_rate;
      }
   }
   else
   if (inf1_sample_rate)
   {
      VT_report( 1, "inf1 sample rate %.6f", inf1_sample_rate);
      sample_rate = round( inf1_sample_rate);
      if (sample_rate) srcal = inf1_sample_rate/sample_rate;
      else
      {
         sample_rate = 1;
         srcal = inf1_sample_rate;
      }
   }
   else
   {
      sample_rate = wav_sample_rate;
      srcal = 1.0;
   }

   if (chans < 0) VT_bailout( "invalid number of channels: %d", chans);
   if (sample_rate <= 0) VT_bailout( "invalid or missing sample rate");

   VT_report( 1, "sample rate %d cal %.6f", sample_rate, srcal);

   vtfile = VT_open_output( outname, chans, 0, sample_rate);
   if (!vtfile) VT_bailout( "cannot open %s: %s", outname, VT_error);

   if (bits_per_sample == 32 && wav_format == 0x3)
      cvt_frame = cvt_frame_float32;
   else
   if (bits_per_sample == 64 && wav_format == 0x3)
      cvt_frame = cvt_frame_float64;
   else
   if (bits_per_sample == 8 && wav_format == 0x1)
      cvt_frame = cvt_frame_uint8;
   else
   if (bits_per_sample == 16 && wav_format == 0x1)
      cvt_frame = cvt_frame_int16;
   else
   if (bits_per_sample == 32 && wav_format == 0x1)
      cvt_frame = cvt_frame_int32;
   else
      VT_bailout( "cannot convert format %04X bits %d",
                   wav_format, bits_per_sample);

   if (background)
   {
      int flags = inname[0] == '-' ? KEEP_STDIN : 0;
      if (outname[0] == '-') flags |= KEEP_STDOUT;
      VT_daemonise( flags);
   }

   double *frame = VT_malloc( sizeof( double) * chans);
   unsigned char *buff = VT_malloc( wav_frame_size);

   timestamp T;
   if (!timestamp_is_NONE( Tstart)) T = Tstart;
   else
   if (!timestamp_is_NONE( inf1_T)) T = inf1_T;
   else
      T = VT_rtc_time();

   VT_set_timebase( vtfile, T, srcal);
   while (1)
   {
      int n = fread( buff, 1, wav_frame_size, fin);
      if (!n) break;   // End of file
      if (n != wav_frame_size) VT_bailout( "short read (%d)", n);

      cvt_frame( buff, frame);
      VT_insert_frame( vtfile, frame);
   }

   VT_release( vtfile);

   return 0;
}

