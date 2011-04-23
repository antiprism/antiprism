/*
   Copyright (c) 2003-2008, Adrian Rossiter

   Antiprism - http://www.antiprism.com
   
   Permission is hereby granted, free of charge, to any person obtaining a
   copy of this software and associated documentation files (the "Software"),
   to deal in the Software without restriction, including without limitation
   the rights to use, copy, modify, merge, publish, distribute, sublicense,
   and/or sell copies of the Software, and to permit persons to whom the
   Software is furnished to do so, subject to the following conditions:

      The above copyright notice and this permission notice shall be included
      in all copies or substantial portions of the Software.

  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
  IN THE SOFTWARE.
*/

/* \file timing.cc
   \brief Timing utilities
*/

#include "timing.h"

#ifdef HAVE_CONFIG_H
   #include "../config.h"
#endif


#if UTIMER == 1   // gettimeofday

#include <unistd.h>
#include <sys/time.h>
int timer::get_time(struct time_val *tp)
{
   return gettimeofday((timeval *)tp, 0);
}

#elif UTIMER == 2 // timeGetTimer

// Norman Vine ?
#include <windef.h>
#include <mmsystem.h>
int timer::get_time(struct time_val *tp)
{
   DWORD t = timeGetTime();
   tp->tv_sec = t/1000;
   tp->tv_usec = (t%1000)*1000;
   return 0; // success
}

#endif

#if USLEEP == 1 // usleep

void timer::u_sleep(unsigned long usecs)
{
   usleep(usecs);
}

#elif USLEEP == 2 // Sleep

#include <windows.h>
void timer::u_sleep(unsigned long usecs)
{   
   Sleep(usecs/1000);
}

#endif



time_val &tv_normalise(time_val &tv)
{
   long sec = tv.tv_usec/1000000;
   tv.tv_usec -= sec*1000000;
   if(tv.tv_usec<0) {
      tv.tv_usec += 1000000;
      tv.tv_sec -= 1;
   }
   tv.tv_sec += sec;
   return tv;
}

bool operator > (const time_val &t0, const time_val &t1)
   { return (t0.tv_sec>t1.tv_sec ||
            (t0.tv_sec==t1.tv_sec && t0.tv_usec>t1.tv_usec) ); }

time_val operator +(const time_val &t0, const time_val &t1)
{ 
   time_val ret;
   ret.tv_sec = t0.tv_sec + t1.tv_sec;
   ret.tv_usec = t0.tv_usec + t1.tv_usec;
   tv_normalise(ret);
   return ret;
}

time_val operator -(const time_val &t0, const time_val &t1)
{ 
   time_val ret;
   ret.tv_sec = t0.tv_sec - t1.tv_sec;
   ret.tv_usec = t0.tv_usec - t1.tv_usec;
   return tv_normalise(ret);
}

time_val to_time_val(double tm)
{ 
   time_val tv;
   tv.tv_sec=long(tm);
   tv.tv_usec=long((tm-tv.tv_sec)*1000000);
   return tv;
}

long to_long_usecs(time_val tv) { return tv.tv_sec*1000000 + tv.tv_usec; }

