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

static int FASCII = 0;                              // -oa option: ASCII output
static int FWAV = 0;                                  // -ow option: WAV output
static int FBWF = FALSE;                            // -obwf option: BFW format
static int RFLAG = 0;                         // -r option: relative timestamps
static double gain = 1.0;                             // -g option: gain factor

static void usage( void)
{
   fprintf( stderr,
       "usage:  vtraw [options] [input]\n"
       "\n"
       "options:\n"
       "  -v        Increase verbosity\n"
       "  -B        Run in background\n"
       "  -L name   Specify logfile\n"
       "  -oa       ASCII output\n"
       "  -ob       Binary output, signed 16 bit\n"
       "  -ow       WAV output, signed 16 bit\n"
       "  -obwf     BWF format, signed 16 bit\n"
       "  -r        Relative time in ASCII output\n"
       "  -g factor Overall gain factor (default 1.0)\n"
       "\n"
       "  -m d='description'  Metadata description\n"
       "  -m a='author'       Metadata author\n"
       "  -m r='reference'    Metadata reference\n"
     );
   exit( 1);
}

static void parse_format_options( char *s)
{
   FASCII = FWAV = FBWF = FALSE;

   if (!strcmp( s, "a") ||
       !strcmp( s, "ascii")) { FASCII = TRUE; return; }
   if (!strcmp( s, "b")) { return; }
   if (!strcmp( s, "w") ||
       !strcmp( s, "wav")) { FWAV = TRUE; return; }
   if (!strcmp( s, "bwf")) { FWAV = TRUE; FBWF = TRUE; return; }

   VT_bailout( "unrecognised output format option: [%s]", s);
}

struct BWF {
   char desc[256];
   char orig[32];
   char orig_ref[32];
   char orig_date[10];
   char orig_time[8];
   uint32_t tref_low;
   uint32_t tref_high;
   uint16_t version;
   uint8_t umid[64];
   uint16_t loud_val;
   uint16_t loud_range;
   uint16_t mtpl;
   uint16_t mml;
   uint16_t mstl;
   uint8_t reserved[180];
} __attribute__((packed)) bwf;

static void output_bwf_chunk( VTFILE *vtfile)
{
   if (sizeof( bwf) != 602) VT_bailout( "incorrect BWF size: %d",
                                        (int) sizeof( bwf));

   uint32_t bb[2] = { 0x74786562, sizeof( struct BWF) };

   timestamp Tstart = VT_get_timestamp( vtfile); 

   char temp[50];
   time_t ts = timestamp_secs( Tstart);
   struct tm *tm = gmtime( &ts);
   sprintf( temp, "%04d-%02d-%02d_%02d:%02d:%02d",
           tm->tm_year + 1900, tm->tm_mon+1, tm->tm_mday,
           tm->tm_hour, tm->tm_min, tm->tm_sec);

   // 2018-06-01_00:00:00.000
   // 01234567890123456789
   memcpy( bwf.orig_date, temp, 10);
   memcpy( bwf.orig_time, temp+11, 8);

   bwf.loud_val = 0x7fff;
   bwf.loud_range = 0x7fff;
   bwf.mtpl = 0x7fff;
   bwf.mml = 0x7fff;
   bwf.mstl = 0x7fff;

   bwf.version = 1;

   uint32_t day = timestamp_secs( Tstart)/86400;
   double secs = timestamp_secs( Tstart) - day * 86400
                 + timestamp_frac( Tstart);
   uint64_t count = secs * VT_get_sample_rate( vtfile);;
   bwf.tref_low = count & 0xffffffffUL;
   bwf.tref_high = count >> 32;

   if (fwrite( bb, 4, 2, stdout) != 2 ||
       fwrite( &bwf, sizeof( bwf), 1, stdout) != 1)
      VT_bailout( "output failed: %s", strerror( errno));
}

static void output_wav_header( int chans, VTFILE *vtfile)
{
   uint32_t ua[3] = { 0x46464952, // RIFF
                      0xffffffff,
                      0x45564157};  // WAVE

   if (fwrite( ua, 4, 3, stdout) != 3)
      VT_bailout( "output failed: %s", strerror( errno));

   if (FBWF) output_bwf_chunk( vtfile);

   uint32_t ub[2] = { 0x20746d66, // 'fmt '
                      16};        // Header length

   uint16_t sa[2] = { 1, chans };

   int sample_rate = VT_get_sample_rate( vtfile);
   uint32_t uc[2] = { sample_rate, sample_rate * chans * 2 };

   uint16_t sb[2] = { 2 * chans, 16 }; 

   if (fwrite( ub, 4, 2, stdout) != 2 ||
       fwrite( sa, 2, 2, stdout) != 2 ||
       fwrite( uc, 4, 2, stdout) != 2 ||
       fwrite( sb, 2, 2, stdout) != 2)
      VT_bailout( "output failed: %s", strerror( errno));

   // 'data' and a field for the size        
   uint32_t ud[2] = { 0x61746164, 0xffffffff };
 
   if (fwrite( ud, 4, 2, stdout) != 2)
      VT_bailout( "output failed: %s", strerror( errno));
}

static void complete_wav_header( int chans, uint64_t frames)
{
   int offset = 40;
   if (FBWF) offset += 602 + 8;

   if (fseek( stdout, offset, SEEK_SET))
   {
      VT_report( 1, "cannot set wave header size, non-seekable");
      return;
   }

   uint32_t u1 = frames * 2 * chans;
   uint32_t u2 = u1 + 44 - 8;

   if (fwrite( &u1, 4, 1, stdout) != 1 ||
       fseek( stdout, 4, SEEK_SET) ||
       fwrite( &u2, 4, 1, stdout) != 1)
      VT_bailout( "output failed: %s", strerror( errno));
}

static void parse_meta( char *s)
{
   if (!strncmp( s, "d=", 2))
   {
      if (strlen( s+2) > 256) VT_bailout( "description field too long");
      strncpy( bwf.desc, s+2, 256);
   }
   else
   if (!strncmp( s, "a=", 2))
   {
      if (strlen( s+2) > 32) VT_bailout( "author field too long"); 
      strncpy( bwf.orig, s+2, 32);
   }
   else
   if (!strncmp( s, "r=", 2))
   {
      if (strlen( s+2) > 32) VT_bailout( "reference field too long"); 
      strncpy( bwf.orig_ref, s+2, 32);
   }
   else VT_bailout( "unrecognised -m argument [%s]", s);
}

int main( int argc, char *argv[])
{
   VT_init( "vtraw");

   int background = 0;
   
   while (1)
   {
      int c = getopt( argc, argv, "vBro:g:L:m:?");
      
      if (c == 'v') VT_up_loglevel();
      else
      if (c == 'B') background = 1;
      else
      if (c == 'L') VT_set_logfile( "%s", optarg);
      else
      if (c == 'o') parse_format_options( optarg);
      else
      if (c == 'g') gain = atof( optarg);
      else
      if (c == 'm') parse_meta( optarg);
      else
      if (c == 'r') RFLAG = 1;
      else
      if (c == -1) break;
      else
         usage(); 
   }  
  
   if (argc > optind + 1) usage();
   char *bname = strdup( optind < argc ? argv[optind] : "-");
 
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
   int chans = chspec->n;
   VT_report( 1, "channels: %d, sample_rate: %d",
                   chans, VT_get_sample_rate( vtfile));

   timestamp Tstart = VT_get_timestamp( vtfile); 
   double *frame;
   int ch;

   if (FASCII)
      while (1)
      {
         timestamp T = VT_get_timestamp( vtfile); 

         if ((frame = VT_get_frame( vtfile)) == NULL) break;

         if (RFLAG) printf( "%.7f", timestamp_diff( T, Tstart));
         else
         {
            char temp[30];   timestamp_string7( T, temp);
            printf( "%s", temp);
         }

         for (ch = 0; ch < chans; ch++)
         {
            double v = frame[chspec->map[ch]];
            if (printf( " %.5e", v) <= 0)
               VT_bailout( "output failed: %s", strerror( errno));
         }
         if (printf( "\n") <= 0)
            VT_bailout( "output failed: %s", strerror( errno));
      }
   else
   {
      if (FWAV) output_wav_header( chans, vtfile);
      uint64_t nout = 0;

      while (1)
      {
         if ((frame = VT_get_frame( vtfile)) == NULL) break;

         for (ch = 0; ch < chans; ch++)
         {
            double v = gain * frame[chspec->map[ch]] * 32767;
            short s;
            if (v > INT16_MAX) s = INT16_MAX;
            else
            if (v < INT16_MIN) s = INT16_MIN;
            else s = v;
   
            if (fwrite( &s, 2, 1, stdout) != 1)
               VT_bailout( "output failed: %s", strerror( errno));
         }

         nout++;
      }

      if (FWAV) complete_wav_header( chans, nout);
   }

   return 0;
}

