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

#include <sys/timex.h>  // XXX

#define UBX_CFG_TMODE2 0x3D06
#define UBX_CFG_RXM    0x1106
#define UBX_CFG_NMEA   0x1706
#define UBX_CFG_GNSS   0x3E06
#define UBX_CFG_ANT    0x1306
#define UBX_CFG_CFG    0x0906
#define UBX_CFG_SBAS   0x1606
#define UBX_CFG_TPS    0x3106

#define UBX_NAV_SVINFO 0x3001

#define UBX_MON_VER    0x040A
#define UBX_MON_HW     0x090A

static int XFLAG = FALSE;
static int background = FALSE;
static char *device = NULL;

struct UBX {

   int h;                    // File handle to the GPS
   pthread_t pid;            // Thread ID of the service task

   uint16_t wait_id;         // Class and ID that the application waits for
   volatile int wait_len;    // Length of the response payload
   uint8_t *wait_res;        // Place to store the response payload
   volatile int wait_done;   // Response complete flag
   volatile int wait_fail;   // Response failed flag

   // Position fix from NMEA GPGGA
   double lat, lon, alt, hdop;
   int nsat;
   char date[50], time[50];

   // Callback functions for NMEA sentences and fix strings
   void (*callback_nmea)(char *);
   void (*callback_fix)(char *);
};

static long double hires_rtc( void)
{
   struct timeval tv;
   gettimeofday( &tv, NULL);

   return tv.tv_sec + 1e-6 * (long double) tv.tv_usec;
}

static void ubx_read_byte( int h, uint8_t *u)
{
   long double t = hires_rtc();

   while (1)
   {
      int n = read( h, u, 1);

      if (n < 0)
      {
         if (errno == EAGAIN)
         {
            if (hires_rtc() - t > 1.0) VT_bailout( "gps comms timeout");
            usleep( 1000);
            continue;
         }

         VT_bailout( "gps read error: %s", strerror( errno));
      }

      break;
   }
}

//
//  Decode NMEA strings.
//

static double convert_gga_coordinate( char *coord, char *nsew)
{
   // [d][d]dmm.mmmmm
  
   double flop = atof( coord) / 100.0;  // d.mmmmmm
   int deg = floor( flop);              // d
   double mins = (flop - deg) * 100;    // mm.mmmm
   double angle = deg + mins/60;

   return *nsew == 'S' || *nsew == 'W' ? -angle : angle;
}

static void decode_nmea( struct UBX *ubx, char *s)
{
   //
   //  Discard bad framing and bad checksum. Crop to remove the checksum.
   //

   char *p = s + strlen(s) - 3;
   if (*p != '*') return;
   *p++ = 0;
   uint8_t u = 0;
   int i;
   for (i = 0; s[i]; i++) u ^= s[i];
   char csum[3];
   sprintf( csum, "%02X", u);
   if (strcmp( csum, p)) return;

   if (ubx->callback_nmea) ubx->callback_nmea( s);

   //
   //  Split sentence into fields.
   //

   char *fields[50];
   int nfields = 0;
   while (s && *s)
   {
      fields[nfields++] = s;
      s = strchr( s, ',');
      if (s) *s++ = 0;
   }

   //
   //  Pick out the info we want from valid sentences.
   //

   if ((!strcmp( fields[0], "GPRMC") || !strcmp( fields[0], "GNRMC")) &&
       nfields >= 10 &&
       strlen( fields[9]) == 6)
   {
      sprintf( ubx->date, "20%2.2s-%2.2s-%2.2s",
                          fields[9] + 4, fields[9] + 2, fields[9]);
   }

   if ((!strcmp( fields[0], "GPGGA") || !strcmp(fields[0], "GNGGA")) &&
       nfields >= 13 &&
       atoi( fields[6]) > 0 &&
       strlen( fields[1]) >= 6)
   {
      sprintf( ubx->time, "%2.2s:%2.2s:%2.2s", 
            fields[1], fields[1] + 2, fields[1] + 4);
      ubx->lat = convert_gga_coordinate( fields[2], fields[3]);
      ubx->lon = convert_gga_coordinate( fields[4], fields[5]);
      ubx->alt = atof( fields[9]);
      ubx->nsat = atoi( fields[7]);
      ubx->hdop = atof( fields[8]);

      if (ubx->callback_fix && ubx->date[0])
      {
         char temp[100];

         sprintf( temp, "%s %s", ubx->date, ubx->time);

         // 2017-11-04_07:34:58
         // 0123456789012345678
         temp[4] = temp[7] = temp[10] = temp[13] = temp[16] = 0;
         struct tm t; memset( &t, 0, sizeof( t));
         t.tm_sec = atoi( temp+17);
         t.tm_min = atoi( temp+14);
         t.tm_hour = atoi( temp+11);
         t.tm_mday = atoi( temp+8);
         t.tm_mon = atoi( temp+5) - 1;
         t.tm_year = atoi( temp) - 1900;
         time_t ut = mktime( &t);

         char *msg;
         if (asprintf( &msg, "%s_%s %lu %12.7f %11.7f %7.2f %2d %6.2f",
                              ubx->date, ubx->time, (unsigned long) ut,
                              ubx->lat, ubx->lon, ubx->alt,
                              ubx->nsat, ubx->hdop) > 0)
         {
            ubx->callback_fix( msg);
            free( msg);
         }
      }
   }
}

//
//  Decode a UBX reply and notify the application via the wait_done flag
//  if the app is waiting for it.
//

static void decode_ubx( struct UBX *ubx, uint8_t *buf, int n)
{
   int len = *(uint16_t *)(buf + 2);

   if (ubx->wait_done) return;

   if (n == 8 &&
       *(uint16_t *) buf == 0x0005 &&
       *(uint16_t *) (buf + 4) == ubx->wait_id)
   {
      VT_report( 2, "command failed: class 0x%02X id 0x%02X", buf[4], buf[5]);
      ubx->wait_fail = TRUE;
      return;
   }

   //
   //  Decode a reply with payload as response to a poll request.
   //

   if (ubx->wait_id &&
       *(uint16_t *) buf == ubx->wait_id)
   {
      int len = *(uint16_t *)(buf + 2);
      if (ubx->wait_res) memcpy( ubx->wait_res, buf + 4, len);
      ubx->wait_len = len;
      ubx->wait_done = TRUE;
   }

   //
   //  Receive an acknowledgement.  The class/id are the two bytes
   //  of the payload.
   //

   if (n == 8 &&
       ubx->wait_id &&
       *(uint16_t *) buf == 0x0105 &&
       *(uint16_t *) (buf + 4) == ubx->wait_id)
   {
      if (ubx->wait_res) memcpy( ubx->wait_res, buf + 4, len);
      ubx->wait_len = len;
      ubx->wait_done = TRUE;
   }
}

//
//  Calculate the checksum.
//

static uint16_t fletcher( uint8_t *s, int n)
{
   uint8_t ck_a = 0, ck_b = 0;
   int i;
   for (i = 0; i < n; i++)
   {
      ck_a += s[i];
      ck_b += ck_a;
   }

   return ck_b << 8 | ck_a;
}

//
//  Send buffer to device.
//

static int ubx_send( struct UBX *ubx, uint8_t *buf, int len)
{
   return write( ubx->h, buf, len) == len;
}

//
//  Background thread to receive messages from the GPS, asynchronously to the
//  application.
//

static void *ubx_service( void *arg)
{
   struct UBX *ubx = (struct UBX *) arg;
   int state = 0;
   int n = 0;
   int ubxlen = 0;
   uint8_t buff[500];

   while (1)
   {
      uint8_t u;
      ubx_read_byte( ubx->h, &u);

      switch (state)
      {
         case 0:  // Idle state
                  if (u == '$') state = 1, n = 0;
                  if (u == 0xb5) state = 2;
                  break;
         case 1:  // Receiving NMEA sentence
                  if (u == '\r' || u == '\n')
                  {
                     buff[n] = 0;
                     decode_nmea( ubx, (char *) buff);
                     state = 0;
                  }
                  buff[n++] = u;
                  break;

         case 2:  // Waiting for UBX sync char 2
                  if (u == 0x62) state = 3, n = 0;
                  else state = 0;
                  break;

         case 3:  // Receiving UBX packet
                  buff[n] = u;
                  if (n == 3) ubxlen = *(uint16_t *) (buff + 2);
                  n++;
                  if (n > 4 && n == ubxlen + 6)
                  {
                     if (fletcher( buff, n - 2) == 
                         *(uint16_t *)(buff + n - 2))
                        decode_ubx( ubx, buff, n);
                     state = 0;
                  }
                  break;
      }
   }

   return 0;
}

//
//  Make a UBX transaction.
//

static int ubx_transact( struct UBX *ubx, uint16_t id,
                         uint8_t *cmd, int cmdlen, uint8_t *res, int *reslen)
{
   if (reslen) *reslen = 0;
   uint8_t *buf = malloc( 8 + cmdlen);

   *(uint16_t *)(buf + 0) = 0x62b5;
   *(uint16_t *)(buf + 2) = id;
   *(uint16_t *)(buf + 4) = cmdlen;
   if (cmd) memcpy( buf + 6, cmd, cmdlen);
   *(uint16_t *)(buf + 6 + cmdlen) = fletcher( buf+2, cmdlen + 4);

   ubx->wait_id = id;
   ubx->wait_res = res;
   ubx->wait_len = 0;
   ubx->wait_done = FALSE;
   ubx->wait_fail = FALSE;

   int i;
   for (i = 0; i< cmdlen+8; i++) VT_report( 2, "send: %02x", buf[i]);

   if (!ubx_send( ubx, buf, cmdlen + 8))
   {
      free( buf);
      return FALSE;
   }

   free( buf);

   //
   //  Wait up to 2.0 seconds for the response.
   //

   long double t = hires_rtc() + 2.0;
   while (!ubx->wait_fail && hires_rtc() < t)
      if (ubx->wait_done)
      {
         if (reslen) *reslen = ubx->wait_len;
         return TRUE;
      }
      else usleep( 1000);

   return FALSE;
}

//
//  Open the device, set serial attributes and start service thread.
//

static struct UBX *ubx_open( char *dev)
{
   struct UBX *ubx = malloc( sizeof( struct UBX));

   if ((ubx->h = open( dev, O_RDWR | O_NOCTTY)) < 0)
      VT_bailout( "cannot open %s: %s", dev, strerror( errno));

   struct termios t;
   tcgetattr( ubx->h, &t);
   cfmakeraw( &t);
   t.c_cflag &= ~CSIZE;
   t.c_cflag |= CS8;
   t.c_cflag |= CREAD | CLOCAL | HUPCL;
   t.c_cflag &= ~CRTSCTS;
   t.c_iflag |= IGNBRK | IGNPAR;
   t.c_iflag &= ~(IXOFF | IMAXBEL);
   t.c_oflag &= ~(ONLCR | OCRNL | ONOCR | ONLRET);
   t.c_lflag &= ~ECHO;

   if (tcsetattr( ubx->h, TCSAFLUSH, &t) < 0)
      VT_bailout( "unable to set serial mode: %s", strerror( errno));

   if (pthread_create( &ubx->pid, NULL, ubx_service, (void *) ubx) < 0)
      VT_bailout( "unable to create service thread: %s", strerror( errno));

   return ubx;
}

static int ubx_series = 0;

static void ubx_get_series( struct UBX *ubx)
{
   uint8_t res[1024];
   int len;

   if (!ubx_transact( ubx, UBX_NAV_SVINFO, NULL, 0, res, &len))
      VT_bailout( "unable to determine device type");

   switch (res[5])  // globalFlags
   {
      case 0: ubx_series = 4; break;
      case 1: ubx_series = 5; break;
      case 2: ubx_series = 6; break;
      case 3: ubx_series = 7; break;
      case 4: ubx_series = 8; break;
   }
}

static void ubx_versions( struct UBX *ubx)
{
   ubx_get_series( ubx);
   printf( "u-blox series: %d\n", ubx_series);

   uint8_t res[1024];
   int len;

   if (!ubx_transact( ubx, UBX_MON_VER, NULL, 0, res, &len) || len < 40)
      VT_bailout( "unable to determine device versions");

   printf( "software version: %s\n", (char *) res);
   printf( "hardware version: %s\n", (char *) (res + 30));

   int i;
   for (i = 40; i < len; i += 30)
      printf( "ext: %s\n", res + i);
}

// Series 7:
//  1 to 32: GPS
//  120-158: SBAS
//  193-197: QZSS
//  65-96, 255: GLONASS

// Series 8:
//  1 to 32: GPS
//  120-158: SBAS
//  193-197: QZSS
//  211-246: Galileo
//  65-96,255: GLONASS
//  159-163,33-64: BeiDou
//  173-182: IMES

static void lookup_gnid( int svid, char *temp)
{
   char *prefix = "XX";
   int base = 0;

   if (svid >= 1 && svid <= 32) prefix = "GP", base = 1;
   else
   if (svid >= 120 && svid <= 158) prefix = "SB", base = 120;
   else
   if (svid >= 193 && svid <= 197) prefix = "QZ", base = 193;
   else
   if (svid >= 211 && svid <= 246) prefix = "GA", base = 211;
   else
   if (svid >= 65 && svid <= 96 || svid == 255) prefix = "GL", base = 65;
   if (svid == 255) prefix = "GL", base = 255 - 98;
   else
   if (svid >= 159 && svid <= 163) prefix = "BD", base = 32 + 159;
   else
   if (svid >= 33 && svid <= 64) prefix = "BD", base = 33;
   else
   if (svid >= 173 && svid <= 182) prefix = "IM", base = 173;

   sprintf( temp, "%s:%d", prefix, svid - base + 1);
}

static void ubx_list( struct UBX *ubx)
{
   uint8_t res[1024];
   int len;

   if (!ubx_transact( ubx, UBX_NAV_SVINFO, NULL, 0, res, &len))
      VT_bailout( "unable to determine device type");

   int chan;
   for (chan = 0; chan < res[4]; chan++)
   {
      uint8_t *p = res + 8 + 12 * chan;
      char temp[50];
      lookup_gnid( p[1], temp);
      printf( "chan %3u sv %-6s %3d dB elev %2d az %3d %s\n",
             p[0], temp, p[4], p[5], p[6], p[2] & 1 ? "used" : ""); 
   }
}

//
//  Set timepulse into frequency mode.
//  
//  tp:  0 or 1 to select one of the two timepulse outputs;
//  hz1: Frequency to use when GPS locked, Hz;
//  hz2: Frequency to use when GPS unlocked, Hz;
//  antdel; Antenna cable delay, nS;
//  rfdel;  RF group delay, nS;
//  usrdel: User configurable delay, nS;
//  duty:  Duty cycle, fraction, eg 0.5 to 50%;
//  pol: Polarity: 0 = falling edge on second, 1 = rising edge on second
//
//  Return: non-zero if successful, zero on error.
//

static void ubx_set_frequency( struct UBX *ubx, int tp,
                               uint32_t hz1, uint32_t hz2,
                               int16_t antdel, int16_t rfdel, int32_t usrdel,
                               double duty, int pol)
{
   uint8_t cmd[32];
   memset( cmd, 0, sizeof( cmd));
   cmd[0] = tp;
   *(int16_t *)(cmd+4) = antdel;
   *(int16_t *)(cmd+6) = rfdel;
   *(uint32_t *)(cmd + 8) = hz2;
   *(uint32_t *)(cmd + 12) = hz1;
   *(uint32_t *)(cmd + 16) = 
   *(uint32_t *)(cmd + 20) = (uint32_t) round( duty * ((uint64_t)1 << 32));
   *(int32_t *)(cmd + 24) = usrdel;

   uint32_t flags = 0;
   flags |= 0x1; // Active
   flags |= 0x2; // LockGpsFreq
   flags |= 0x4; // LockedOtherSet
   flags |= 0x8; // isFreq
   
   flags |= 0x20; // alignToTow
   if (pol == 1) flags |= 0x40; // polarity

   *(uint32_t *)(cmd + 28) = flags;

   if (!ubx_transact( ubx, UBX_CFG_TPS, cmd, 32, NULL, NULL))
      VT_bailout( "unable to set frequency");

   VT_report( 1, "timepulse %d: frequency set", tp+1);
}

//
//  Set timepulse into pulse mode.
//  
//  tp:  0 or 1 to select one of the two timepulse outputs;
//  period1: Period to use when GPS locked, uS;
//  period2: Period to use when GPS unlocked, uS;
//  antdel; Antenna cable delay, nS;
//  rfdel;  RF group delay, nS;
//  usrdel: User configurable delay, nS;
//  length:  Pulse length, uS;
//  pol: Polarity: 0 = falling edge on second, 1 = rising edge on second
//
//  Return: non-zero if successful, zero on error.
//

static void ubx_set_pulse( struct UBX *ubx, int tp,
                   uint32_t period1, uint32_t period2,
                   int16_t antdel, int16_t rfdel, int32_t usrdel,
                   uint32_t length, int pol)
{
   uint8_t cmd[32];
   memset( cmd, 0, sizeof( cmd));
   cmd[0] = tp;
   *(int16_t *)(cmd+4) = antdel;
   *(int16_t *)(cmd+6) = rfdel;
   *(uint32_t *)(cmd + 8) = period2;
   *(uint32_t *)(cmd + 12) = period1;
   *(uint32_t *)(cmd + 16) = 
   *(uint32_t *)(cmd + 20) = length;
   *(int32_t *)(cmd + 24) = usrdel;

   uint32_t flags = 0;
   flags |= 0x1; // Active
   flags |= 0x2; // LockGpsFreq
   flags |= 0x4; // LockedOtherSet
   flags |= 0x10; // isLength 
   flags |= 0x20; // alignToTow
   if (pol == 1) flags |= 0x40; // polarity
   flags |= 0x80;   // GPS grid

   *(uint32_t *)(cmd + 28) = flags;

   if (!ubx_transact( ubx, UBX_CFG_TPS, cmd, 32, NULL, NULL))
      VT_bailout( "cannot set pulse width");

   VT_report( 1, "timepulse %d: pulse set", tp+1);
}

static void ubx_get_sbas( struct UBX *ubx) 
{
   uint8_t res[1024];
   int len;


   if (!ubx_transact( ubx, UBX_CFG_SBAS, NULL, 0, res, &len))
      VT_bailout( "cannot retrieve SBAS settings");

   printf( "sbas enable: %s\n", res[0] & 1 ? "enabled" : "disabled");
   printf( "sbas channels: %u\n", res[2]);
}

//
//  Enable or disable SBAS.
//

static void ubx_set_sbas( struct UBX *ubx, int state) 
{
   uint8_t cmd[8];
   memset( cmd, 0, sizeof( cmd));
   cmd[0] = state ? 1 : 0;  // mode
   cmd[1] = 7;  // usage
   cmd[2] = 3;  // maxSBAS

   if (!ubx_transact( ubx, UBX_CFG_SBAS, cmd, 8, NULL, NULL))
      VT_bailout( "unable to set SBAS");

   VT_report( 1, "SBAS set %s", state ? "ON" : "OFF");
}

//
//  Save configuration to non-volatile memory.
//

static void ubx_save_config( struct UBX *ubx)
{
   uint8_t cmd[12];
   memset( cmd, 0, sizeof( cmd));
   *(uint32_t *)(cmd + 4) = 0x61f;  // saveMask

   if (!ubx_transact( ubx, UBX_CFG_CFG, cmd, 12, NULL, NULL))
      VT_bailout( "unable to save config");

   VT_report( 1, "configuration saved");
}

static void ubx_callback_nmea( struct UBX *ubx, void (*fcn)(char *))
{
   ubx->callback_nmea = fcn;
}

static void ubx_callback_fix( struct UBX *ubx, void (*fcn)(char *))
{
   ubx->callback_fix = fcn;
}

static void ubx_get_antenna( struct UBX *ubx)
{
   uint8_t res[1024];
   int len;

   if (!ubx_transact( ubx, UBX_CFG_ANT, NULL, 0, res, &len) ||
       len != 4) VT_bailout( "cannot get antenna settings");
 
   uint16_t flags = *(uint16_t *)(res + 0);

   printf( "svcs: %s\n", flags & 1 ? "enabled" : "disabled");
   printf( "scd: %s\n", flags & 2 ? "enabled" : "disabled");
   printf( "ocd: %s\n", flags & 4 ? "enabled" : "disabled");
   printf( "pd on scd: %s\n", flags & 8 ? "enabled" : "disabled");
   printf( "autorec scd: %s\n", flags & 16 ? "enabled" : "disabled");

   printf( "pins: %04X\n", *(uint16_t *)(res + 2));
}

static void ubx_get_status( struct UBX *ubx)
{
   uint8_t res[1024];
   int len;

   if (!ubx_transact( ubx, UBX_MON_HW, NULL, 0, res, &len) ||
       len < 60) VT_bailout( "cannot get hardware_status");
 
   printf( "noise: %u\n", *(uint16_t *)(res + 16));  // noisePerMS

   char *ant = "undefined";
   switch (res[20])
   {
      case 0: ant="init";     break;
      case 1: ant="unknown";  break;
      case 2: ant="ok";       break;
      case 3: ant="short";    break;
      case 4: ant="open";     break;
   }
   printf( "antenna: %s\n", ant);

   if (len == 68)
   {
      switch (res[21])
      {
         case 0: printf( "ant power: off\n");      break;
         case 1: printf( "ant power: on\n");       break;
         case 2: printf( "ant power: unknown\n");  break;
      }
   }
   uint8_t jam = (res[22] >> 2) & 0x3;
   switch (jam)
   {
      case 0: printf( "jamming: unknown\n");    break;
      case 1: printf( "jamming: ok\n");         break;
      case 2: printf( "jamming: present\n");    break;
      case 3: printf( "jamming: critical\n");    break;
   }

   if (jam)
   {
      if (len == 60) printf( "cwjam: %u\n", res[45]);   // jamInd
      if (len == 68) printf( "cwjam: %u\n", res[53]);   // jamInd
   }

   if (!ubx_transact( ubx, UBX_CFG_RXM, NULL, 0, res, &len) ||
       len < 2) VT_bailout( "cannot retrieve power save settings %d", len);

   printf( "power save mode: %d\n", res[1]);

}

static char *gnss_names[] = {
   "GPS", "SBAS", "Galileo", "BeiDou", "IMES", "QZSS", "GLONASS"
};

static void ubx_get_gnss( struct UBX *ubx)
{
   uint8_t res[1024];
   int len;

   if (!ubx_transact( ubx, UBX_CFG_GNSS, NULL, 0, res, &len) ||
       len < 4) VT_bailout( "cannot get GNSS settings");

   printf( "total channels: %u\n", res[1]);

   int i;
   for (i = 0; i < res[3]; i++)
   {
      uint8_t *p = res + 4 + i * 8;
      uint32_t flags = *(uint32_t *)(p+4);

      printf( "GNSS %s\n", gnss_names[p[0]]);
      printf( "  chans reserved: %u\n", p[1]);
      printf( "  chans used:     %u\n", p[2]);
      printf( "  enabled:        %s\n", flags & 1 ? "Yes" : "No");
   }
}

static void ubx_set_gnss( struct UBX *ubx, char *args)
{
   uint8_t res[1024];
   int len;

   if (!ubx_transact( ubx, UBX_CFG_GNSS, NULL, 0, res, &len) ||
       len < 4) VT_bailout( "failed to apply %s", args);

   char *arg = args;
   while (arg)
   {
      char *t = strchr( arg, ',');
      if (t) *t++ = 0;

      int enable = FALSE;
      switch (*arg)
      {
         case '+': enable = TRUE; break;
         case '-': enable = FALSE; break;
         default: VT_bailout( "invalid argument [%s]", arg);
      }
   
      arg++;
      int id;
      for (id = 0; id < 7; id++)
         if (!strncmp( arg, gnss_names[id], strlen( gnss_names[id]))) break;
      if (id == 7) VT_bailout( "unrecognised GNSS name [%s]", arg);
   
      VT_report( 1, "set gnss [%s] %d %s",
                     gnss_names[id], id, enable ? "Enable" : "Disable");
      char *p = strchr( arg, '=');
      int chans = p ? atoi( p+1) : 0;
  
      int i; 
      for (i = 0; i < res[3]; i++)
      {
         uint8_t *p = res + 4 + i * 8;
   
         if (p[0] != id) continue;
   
         if (enable) *(uint32_t *)(p+4) |= 1;
         else *(uint32_t *)(p+4) &= ~1;
   
         if (enable && chans) p[2] = chans;
      }

      arg = t;
   }

   if (!ubx_transact( ubx, UBX_CFG_GNSS, res, len, NULL, NULL)) 
      VT_bailout( "failed to apply %s", optarg);

}

static void ubx_get_nmea_version( struct UBX *ubx)
{
   uint8_t res[1024];
   int len;

   if (!ubx_transact( ubx, UBX_CFG_NMEA, NULL, 0, res, &len) ||
       len < 4) VT_bailout( "unable to get NMEA version");

   printf( "filter: %02X\n", res[0]);
   printf( "nmeaVersion: %02X\n", res[1]);
   printf( "numSV %d\n", res[2]);
   printf( "flags %02X\n", res[3]);
   printf( "gnssToFilter: %08X\n", *(uint32_t *)(res+4));
   printf( "svNumbering: %d\n", res[8]);
   printf( "gsvTalkerId: %d\n", res[10]);
}

static void ubx_set_nmea_version( struct UBX *ubx, uint8_t version)
{
   uint8_t res[1024];
   int len;

   if (!ubx_transact( ubx, UBX_CFG_NMEA, NULL, 0, res, &len) ||
       len < 4) VT_bailout( "failed to set NMEA version");

   res[1] = version;
   if (!ubx_transact( ubx, UBX_CFG_NMEA, res, len, NULL, NULL)) 
      VT_bailout( "failed to set NMEA version");

   VT_report( 1, "NMEA version %u set", version);
}

static void ubx_get_time_mode( struct UBX *ubx)
{
   uint8_t res[1024];
   int len;

   if (!ubx_transact( ubx, UBX_CFG_TMODE2, NULL, 0, res, &len) ||
       len < 28) VT_bailout( "unable to get time mode");

   switch (res[0])
   {
      case 0: printf( "time mode: off\n");        break;
      case 1: printf( "time mode: survey-in\n");  break;
      case 2: printf( "time mode: active\n");     break;
      default: printf( "time mode: unknown\n");   break;
   }

   if (res[2] & 1) // LLA data?
   {
      printf( "lat: %.7f deg\n", *(int32_t *) (res + 4) * 1e-7);
      printf( "lon: %.7f deg\n", *(int32_t *) (res + 8) * 1e-7);
      printf( "alt: %.7f m\n", *(int32_t *) (res + 12) * 1e-2);
   }
   else   // ECEF coordinates
   {
      printf( "ECEFx: %.7f m\n", *(int32_t *) (res + 4) * 1e-2);
      printf( "ECEFy: %.7f m\n", *(int32_t *) (res + 8) * 1e-2);
      printf( "ECEFz: %.7f m\n", *(int32_t *) (res + 12) * 1e-2);
   }
   
   printf( "accuracy: %.3f m\n", 1e-3 * *(uint32_t *) (res + 16));
   printf( "survey-in duration: %u secs\n", *(uint32_t *) (res + 20));
   printf( "survey-in accuracy: %.3f m\n", 1e-3 * *(uint32_t *) (res + 24));
}

static void ubx_start_survey_in( struct UBX *ubx, int secs, int h)
{
   uint8_t res[1024];
   int len;

   if (!ubx_transact( ubx, UBX_CFG_TMODE2, NULL, 0, res, &len) ||
       len < 28) VT_bailout( "unable to start survey-in");

   *(uint32_t *) (res + 20) = secs;
   *(uint32_t *) (res + 24) = h * 1000;
   res[0] = 1;  // Request survey-in

   if (!ubx_transact( ubx, UBX_CFG_TMODE2, res, len, NULL, NULL))
      VT_bailout( "unable to start survey-in");

   VT_report( 1, "survey-in started");
}

static void ubx_set_time_mode( struct UBX *ubx,
                              double lat, double lon, double alt)
{
   uint8_t res[1024];
   int len;

   if (!ubx_transact( ubx, UBX_CFG_TMODE2, NULL, 0, res, &len) ||
       len < 28) VT_bailout( "unable to set time mode");
   
   *(int32_t *) (res + 4) = lat * 1e7;
   *(int32_t *) (res + 8) = lon * 1e7;
   *(int32_t *) (res + 12) = alt * 1e2;

   res[0] = 2;   // Set time mode
   *(uint16_t *) (res + 2) = 1;  // Flags: use LLA

   if (!ubx_transact( ubx, UBX_CFG_TMODE2, res, len, NULL, NULL))
      VT_bailout( "unable to set time mode");
}

#ifdef __linux__
static double gps_delay = 0;
static int gpscnt = 0;
static long double diffsum = 0;

static void time_server( char *s)
{
   char temp[100];

   long double now;
   double lat, lon, alt, d;
   int ns;
   if (sscanf( s, "%s %Lg %lg %lg %lg %d %lg",
                   temp, &now, &lat, &lon, &alt, &ns, &d) != 7)
   {
      VT_report( 0, "bad gps string [%s]", s);
      return;
   }

   struct timeval tc;
   gettimeofday( &tc, NULL);

   now += gps_delay;

   long double diff =
      now - (int64_t) tc.tv_sec - 1e-6 * (int64_t) tc.tv_usec;

   diffsum += diff;
   gpscnt++;
   if (gpscnt < 10) return;

   diff = diffsum/gpscnt;
   gpscnt = 0;
   diffsum = 0;

   if (diff > -0.1 && diff < 0.1)
   {
      struct timex tx;

      tx.modes = ADJ_OFFSET_SINGLESHOT;
      tx.offset = (int) (0.05 * diff * 1000000);
      if (tx.offset < -500000) tx.offset = -500000;
      if (tx.offset > 500000) tx.offset = 500000;

      int r = adjtimex( &tx);
      if (r < 0) VT_bailout( "error return from adjtimex");

      VT_report( 0, "%.7f,%.7f,%.1f %d adjust %+.6Lf status %s",
                  lat, lon, alt, ns,
                  diff,
                  r < TIME_BAD ? "locked" : "unlocked");
   }
   else
   {
      struct timeval tv;
      tv.tv_sec = (time_t) now;
      tv.tv_usec = 1e6 * (now - tv.tv_sec);
      if (settimeofday( &tv, NULL) != 0) VT_bailout( "unable to set time");
      VT_report( 0, "step adjust %.6Lf secs", diff);
   }
}
#endif

static void usage( void)
{
   fprintf( stderr,
     "ubx [options] device\n"
     "\n"
     "options:\n"
     "\n"
     " -v             Increase verbosity\n"
     " -F [tp,]Hz     Set frequency output, Hertz\n"
     "                A negative frequency reverses the polarity\n"
     " -P [tp,]secs   Set pulse-per-second output, width in seconds\n"
     " -S             Save configuration to GPS non-volatile memory\n"
     " -b?            Query SBAS\n"
     " -b+            Enable SBAS\n"
     " -b-            Disable SBAS\n"
     " -rh            Report hardware status\n"
     " -ra            Report antenna settings\n"
     " -rv            Report GPS software and hardware versions\n"
     " -f             Output timestamp and position fix records\n"
     " -n             Output NMEA sentences\n"
     " -l             List satellites\n"
     " -N?            Query NMEA protocol version\n"
     " -N version     Set NMEA protocol version\n"
     " -g?            Query GNSS settings\n"
     " -g options     Setup GNSS\n"
     " -T?            Query time mode\n"
     " -T options     Set time mode\n"
     " -x             Time server\n"
   );

   exit( 1);
}

int main( int argc, char *argv[])
{
   VT_init( "vtubx");

   if (argc < 2 ||
       !strcmp( argv[1], "-?")) usage();

   device = argv[--argc];
   struct UBX *h = NULL; 

   while (1)
   {
      int c = getopt( argc, argv, "vSfnlBF:P:b:r:g:N:T:L:x:");

      if (c == 'v') VT_up_loglevel();
      else
      if (c == 'L') VT_set_logfile( "%s", optarg);
      else
      if (c == 'B') background = TRUE;
      else
      if (c == 'l')
      {
         if (!h) h = ubx_open( device);
         ubx_list( h);
      }
      else
      if (c == 'S')
      {
         if (!h) h = ubx_open( device);
         ubx_save_config( h);
      }
      else
      if (c == 'x')
      {
         XFLAG = TRUE;
         #ifndef __linux__
            VT_bailout( "-x option only supported under linux");
         #else
            gps_delay = atof( optarg);
         #endif
      }
      else
      if (c == 'f') 
      {
         void fixprint( char *s)
         {
            if (printf( "%s\n", s) <= 0 ||
                fflush( stdout) < 0) exit( 0);
         }
   
         if (!h) h = ubx_open( device);
         ubx_callback_fix( h, fixprint);
         while (1) usleep( 1000000);
      }
      else
      if (c == 'n') 
      {
         void nmeaprint( char *s)
         {
            if (printf( "%s\n", s) <= 0 ||
                fflush( stdout) < 0) exit( 0);
         }
   
         if (!h) h = ubx_open( device);
         ubx_callback_nmea( h, nmeaprint);
         while (1) usleep( 1000000);
      }
      else
      if (c == 'F')
      {
         char *s = strdup( optarg);
         char *p = strchr( s, ',');
         int tp;
         double freq;
         int polarity = 1;
         if (!p)
         {
            tp = 0;
            freq = atof( s);
         }
         else
         {
            *p++ = 0;
            tp = atoi( s) - 1;
            freq = atof( p);
         }

         if (freq < 0) freq = -freq, polarity = 0;
         if (!h) h = ubx_open( device);
         ubx_set_frequency( h, tp, freq, 0, 0, 0, 0, 0.5, polarity);
      }
      else
      if (c == 'P')
      {
         char *s = strdup( optarg);
         char *p = strchr( s, ',');
         int tp;
         double width;
         int polarity = 1;
         if (!p)
         {
            tp = 0;
            width = atof( s);
         }
         else
         {
            *p++ = 0;
            tp = atoi( s) - 1;
            width = atof( p);
         }

         width = round( width * 1e6);
         if (width < 0) width = -width, polarity = 0;
         if (!h) h = ubx_open( device);
         ubx_set_pulse( h, tp, 1000000, 1000000, 0, 0, 0, width, polarity);
      }
      else
      if (c == 'b')
      {
         if (!h) h = ubx_open( device);

         if (!strcmp( optarg, "?")) ubx_get_sbas( h);
         else
         if (!strcmp( optarg, "+")) ubx_set_sbas( h, TRUE);
         else
         if (!strcmp( optarg, "-")) ubx_set_sbas( h, FALSE);
         else
            VT_bailout( "invalid argument to -b");
      }
      else
      if (c == 'g')
      {
         if (!h) h = ubx_open( device);

         if (!strcmp( optarg, "?"))
         {
            ubx_get_gnss( h);
         }
         else
         {     
            ubx_set_gnss( h, strdup( optarg));
            VT_report( 1, "applied %s", optarg);
         }
      }
      else
      if (c == 'N')
      {
         if (!h) h = ubx_open( device);

         if (!strcmp( optarg, "?")) ubx_get_nmea_version( h);
         else
            ubx_set_nmea_version( h, atoi( optarg));
      }
      else
      if (c == 'T')
      {
         if (!h) h = ubx_open( device);

         if (!strcmp( optarg, "?"))
            ubx_get_time_mode( h);
         else
         {
            char *s = strdup( optarg);
            char *p = strchr( s, ',');
            if (!p) VT_bailout( "invalid argument for -T");
            *p++ = 0;
            char *q = strchr( p, ',');
            if (q)  // Three arguments?  set time mode
            {
               *q++ = 0;
               double lat = atof( s);
               double lon = atof( p);
               double alt = atof( q);
               ubx_set_time_mode( h, lat, lon, alt);
            }
            else    // Two arguments? start survey-in
            {
               int secs = atof( s) * 3600;
               int hm = atoi( p);
               ubx_start_survey_in( h, secs, hm);
            }
         }
      }
      else
      if (c == 'r')
      {
         if (!h) h = ubx_open( device);

         if (!strcmp( optarg, "a")) ubx_get_antenna( h);
         else
         if (!strcmp( optarg, "h")) ubx_get_status( h);
         else
         if (!strcmp( optarg, "v")) ubx_versions( h);
         else
            VT_bailout( "unrecognised -r option [%s]", optarg);
      }
      else
      if (c == -1) break;
      else
         usage();
   }

   if (background) VT_daemonise( 0);

   #ifdef __linux__
      if (XFLAG)
      {
         VT_report( 0, "starting time server");
         if (!h) h = ubx_open( device);
   
         ubx_callback_fix( h, time_server);
         while (1) usleep( 1000000);
      }
   #endif

   return 0;
}

