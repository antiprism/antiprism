/*
   Copyright (c) 2003, Adrian Rossiter

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

/*
   Name: col_val.cc
   Description: representation of color values
   Project: Antiprism - http://www.antiprism.com
*/

#include <string.h>
#include "col_val.h"

using std::min;
using std::max;


col_val col_val::invisible(0,0,0,0);

bool col_val::operator ==(col_val c) const
{
   if(!is_set() && !c.is_set())
      return true;
   if(is_idx() && c.is_idx() && get_idx()==c.get_idx())
      return true;
   if(is_val() && c.is_val() && memcmp(rgba, c.rgba, 4)==0)
      return true;
   return false;
}

bool operator <(const col_val &c1, const col_val &c2)
{
   if(!c1.is_set())
      return c2.is_set();

   if(c1.is_idx())
      return c2.is_val() || c1.get_idx()<c2.get_idx();

   return c1.get_long()<c2.get_long();
}

// The following RGB / HSV functions are taken from
// http://www.cs.rit.edu/~ncs/color/t_convert.html

// r,g,b values are from 0 to 1
// h = [0,1], s = [0,1], v = [0,1]
//    if s == 0, then h = -1 (undefined)

void RGBtoHSV( double r, double g, double b, double *h, double *s, double *v )
{
   double *rgb[] = {&r, &g, &b};
   for(int i=0; i<3; i++) {
      if(*rgb[i]<0)
         *rgb[i]=0;
      else if(*rgb[i]>1)
         *rgb[i]=1;
   }

   double min, max, delta;

   if(r<=g && r<=b)
      min=r;
   else if(g<=r && g<=b)
      min=g;
   else
      min=b;
      
   if(r>=g && r>=b)
      max=r;
   else if(g>=r && g>=b)
      max=g;
   else
      max=b;
  
   *v = max;            // v

   delta = max - min;

   if( max > epsilon )     // avoid division by zero
      *s = delta / max;    // s
   else {
      // r = g = b = 0     // s = 0, v and h are not important
      *s = 0;
      *h = 0;
      return;
   }

   if( delta < epsilon)             // grey range, avoid division by zero
      *h = 0;
   else if( r == max )
      *h = ( g - b ) / delta;       // between yellow & magenta
   else if( g == max )
      *h = 2 + ( b - r ) / delta;   // between cyan & yellow
   else
      *h = 4 + ( r - g ) / delta;   // between magenta & cyan

   *h *= 60;            // degrees
   if( *h < 0 )
      *h += 360;
   *h /= 360;

}

void HSVtoRGB( double *r, double *g, double *b, double h, double s, double v )
{
   // bring h into range 0-360
   h = fmod(h, 1.0)*360;
   if(h<0)
      h += 360;

   if( s == 0 ) {
      // achromatic (grey)
      *r = *g = *b = v;
      return;
   }

   if(h>360-0.0001 || h<0.0001)
      h=0.01;
   h /= 60;       // sector 0 to 5
   int i = (int)floor( h );
   double f = h - i;        // factorial part of h
   double p = v * ( 1 - s );
   double q = v * ( 1 - s * f );
   double t = v * ( 1 - s * ( 1 - f ) );

   switch( i ) {
      case 0:
         *r = v;
         *g = t;
         *b = p;
         break;
      case 1:
         *r = q;
         *g = v;
         *b = p;
         break;
      case 2:
         *r = p;
         *g = v;
         *b = t;
         break;
      case 3:
         *r = p;
         *g = q;
         *b = v;
         break;
      case 4:
         *r = t;
         *g = p;
         *b = v;
         break;
      default:    // case 5:
         *r = v;
         *g = p;
         *b = q;
         break;
   }

}


void col_val::set_hsva(double hue, double sat, double val, double alpha)
{
   double *hsva[] = {&hue, &sat, &val, &alpha};
   for(int i=1; i<4; i++) {  // skip i=0 as hue can wrap
      if(*hsva[i] < 0)
         *hsva[i] = 0;
      else if(*hsva[i] > 1)
         *hsva[i] = 1;
   }

   double r, g, b;
   HSVtoRGB(&r, &g, &b, hue, sat, val);
   set_rgba(r, g, b, alpha);
}

void col_val::set_hsva(const vec4d &hsva)
{
   set_hsva(hsva[0], hsva[1], hsva[2], hsva[3]);
}


vec4d col_val::get_hsva() const
{
   vec4d hsva;
   vec4d rgba_d = get_vec4d();
   RGBtoHSV(rgba_d[0], rgba_d[1], rgba_d[2], &hsva[0], &hsva[1], &hsva[2]);
   hsva[3] = rgba_d[3];
   return hsva;
}

 
// HSL algorithms by Paul Bourke
// http://local.wasp.uwa.edu.au/~pbourke/texture_colour/convert/
/*
   Calculate HSL from RGB
   Hue is in degrees
   Lightness is between 0 and 1
   Saturation is between 0 and 1
*/
vec4d RGB2HSL(const col_val &c)
{
   vec4d c1(c[0],c[1],c[2],c[3]);
   
   double themin = 0.0;
   double themax = 0.0;
   double delta = 0.0;

   themin = min(c1[0],min(c1[1],c1[2]));
   themax = max(c1[0],max(c1[1],c1[2]));
   delta = themax - themin;
   double l = (themin + themax) / 2.0;
   l /= 255.0; // Antiprism
   double s = 0.0;
   if (l > 0.0 && l < 1.0)
      s = delta / (l < 0.5 ? (2.0*l) : (2.0-2.0*l));
   s /= 255.0; // Antiprism
   double h = 0.0;
   if (delta > 0.0) {
      if (themax == c1[0] && themax != c1[1])
         h += (c1[1] - c1[2]) / delta;
      if (themax == c1[1] && themax != c1[2])
         h += (2.0 + (c1[2] - c1[0]) / delta);
      if (themax == c1[2] && themax != c1[0])
         h += (4.0 + (c1[0] - c1[1]) / delta);
      h *= 60.0;
   }
   if (h<0.0) h+=360.0; // Antiprism
   return(vec4d(h/360.0,s,l,c[3]/255.0));
}

/*
   Calculate RGB from HSL, reverse of RGB2HSL()
   Hue is in degrees
   Lightness is between 0 and 1
   Saturation is between 0 and 1
*/
col_val HSL2RGB(vec4d c1)
{
   c1[0] = fmod(c1[0], 1.0)*360; // Antiprism
   if (c1[0] < 0)
      c1[0] += 360;
   
   vec4d c2,sat,ctmp;

   while (c1[0] < 0.0)
      c1[0] += 360.0;
   while (c1[0] > 360.0)
      c1[0] -= 360.0;

   if (c1[0] < 120.0) {
      sat[0] = (120.0 - c1[0]) / 60.0;
      sat[1] = c1[0] / 60.0;
      sat[2] = 0.0;
   } else if (c1[0] < 240.0) {
      sat[0] = 0.0;
      sat[1] = (240.0 - c1[0]) / 60.0;
      sat[2] = (c1[0] - 120.0) / 60.0;
   } else {
      sat[0] = (c1[0] - 240.0) / 60.0;
      sat[1] = 0.0;
      sat[2] = (360.0 - c1[0]) / 60.0;
   }
   sat[0] = min(sat[0],1.0);
   sat[1] = min(sat[1],1.0);
   sat[2] = min(sat[2],1.0);

   ctmp[0] = 2.0 * c1[1] * sat[0] + (1 - c1[1]);
   ctmp[1] = 2.0 * c1[1] * sat[1] + (1 - c1[1]);
   ctmp[2] = 2.0 * c1[1] * sat[2] + (1 - c1[1]);

   if (c1[2] < 0.5) {
      c2[0] = c1[2] * ctmp[0];
      c2[1] = c1[2] * ctmp[1];
      c2[2] = c1[2] * ctmp[2];
   } else {
      c2[0] = (1.0 - c1[2]) * ctmp[0] + 2.0 * c1[2] - 1.0;
      c2[1] = (1.0 - c1[2]) * ctmp[1] + 2.0 * c1[2] - 1.0;
      c2[2] = (1.0 - c1[2]) * ctmp[2] + 2.0 * c1[2] - 1.0;
   }
   
   return(col_val(c2[0],c2[1],c2[2],c1[3]));
}


void col_val::set_hsla(double hue, double sat, double lightness, double alpha)
{
   double *hsla[] = {&hue, &sat, &lightness, &alpha};
   for(int i=1; i<4; i++) {
      if(*hsla[i] < 0)
         *hsla[i] = 0;
      else if(*hsla[i] > 1)
         *hsla[i] = 1;
   }

   *this = HSL2RGB(vec4d(*hsla[0], *hsla[1], *hsla[2], *hsla[3]));
}

void col_val::set_hsla(const vec4d &hsla)
{
   set_hsla(hsla[0], hsla[1], hsla[2], hsla[3]);
}


vec4d col_val::get_hsla() const
{
   return RGB2HSL(*this);
}
