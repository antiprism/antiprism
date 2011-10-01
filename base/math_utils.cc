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

/* \file utils.cc
   \brief utility routines for maths operations.
*/

#ifdef HAVE_CONFIG_H
   #include "../config.h"
#endif

#include <algorithm>
#include <math.h>
#include "math_utils.h"

// Greatest Common Divisor
long gcd(long m, long n)
{
   if(m == 0)
      return n;
   if(n == 0)
      return m;
   return gcd(n,m%n);
}


int factorial(int num)
{
   int result=1;
   for (int i=1; i<=num; i++)
      result*=i;
   return result;
}

#include <stdio.h>
// http://ericlin2.tripod.com/mis/ration.html
double find_denom(double num, const double &eps, int steps)
{
   double dec = num - long(num);
   if(!steps)
      return 0; // less than eps to indic
   if(fabs(dec)<eps)
        return 1;
    num = 1/dec;
    return num*find_denom(num, eps, steps-1);
}

// Convert double to rational
void double2rational(double f, long &num, long &denom, double eps,
      int max_steps)
{
   int sign = (f>0) ? 1 : -1;
   f *= sign;
   denom = long(0.5+find_denom(f, eps, max_steps));
   num = sign*long(0.5 + denom*f);
}
   

static int signR(double Z)
{
   if (Z > 0.0)
      return 1;
   if (Z < 0.0)
      return -1;
   return 0;
}

static double CBRT(double Z)
{
   return fabs(pow(fabs(Z),1/3.0)) * signR(Z);
}


// Solution of a cubic equation
// Equations of lesser degree are solved by the appropriate formulas.
// The solutions are arranged in ascending order.
// 17-Oct-2004 / Raoul Rausch: Conversion from Fortran to C
int cubic(double coeffs[4], double sol[3])
{
   const double THIRD = 1./3.;
   double U[3],W, P, Q, DIS;
   int L;

   // ====determine the degree of the polynomial ====
   if (coeffs[3] != 0.0) {
      //cubic problem
      W = coeffs[2]/coeffs[3]*THIRD;
      P = pow((coeffs[1]/coeffs[3]*THIRD - pow(W,2)),3);
      Q = -.5*(2.0*pow(W,3)-(coeffs[1]*W-coeffs[0])/coeffs[3] );
      DIS = pow(Q,2)+P;
      if ( DIS < 0.0 ) {
         //three real solutions!
         double PHI = acos(safe_for_trig(Q/sqrt(-P)));
         P=2.0*pow((-P),(5.e-1*THIRD));
         for (int i=0; i<3; i++)
            U[i] = P*cos((PHI+2*((double)i)*M_PI)*THIRD)-W;
         sol[0] = std::min(U[0], std::min(U[1], U[2]));
         sol[1] = std::max(std::min(U[0], U[1]),
                           std::max( std::min(U[0], U[2]),
                                     std::min(U[1], U[2])  ));
         sol[2] = std::max(U[0], std::max(U[1], U[2]));
         L = 3;
      }
      else {
         // only one real solution!
         DIS = sqrt(DIS);
         sol[0] = CBRT(Q+DIS)+CBRT(Q-DIS)-W;
         L = 1;
      }
   }
   else if (coeffs[2] != 0.0) {  // quadratic problem
      P = 0.5*coeffs[1]/coeffs[2];
      DIS = pow(P,2)-coeffs[0]/coeffs[2];
      if (DIS > 0.0) {
         // 2 real solutions
         sol[0] = -P - sqrt(DIS);
         sol[1] = -P + sqrt(DIS);
         L = 2;
      }
      else // no real solution
         L = 0;
   }
   else if (coeffs[1] != 0.0) {  //linear equation
      sol[0] =coeffs[0]/coeffs[1];
      L = 1;
   }
   else //no equation
      L = 0;
   
   // perform one step of a newton iteration to minimize round-off errors
   for (int i=0;i<L;i++)
      sol[i] = sol[i] -
         (coeffs[0]+sol[i]*(coeffs[1]+sol[i]*(coeffs[2]+sol[i]*coeffs[3]))) /
         (coeffs[1]+sol[i]*(2.0*coeffs[2]+sol[i]*3.0*coeffs[3]));

   return L;
}


// Solution of a quartic equation
// ref.: J. E. Hacke, Amer. Math. Monthly, Vol. 48, 327-328, (1941)
// 17-Oct-2004 / Raoul Rausch: Conversion from Fortran to C
int quartic(double dd[5], double sol[4], double sol_i[4])
{
   double AA[4], z[3];
   double a, b, c, d, f, p, q, r, zsol, xK2, xL, xK, sqp, sqm;
   double tmpi[4];
   double *soli = sol_i ? sol_i : tmpi;
   int Nsol = 0;

   if (dd[4] == 0.0) // not a quartic
      return -1;

   a = dd[4];
   b = dd[3];
   c = dd[2];
   d = dd[1];
   f = dd[0];

   p = (-3.0*pow(b,2) + 8.0 *a*c)/(8.0*pow(a,2));
   q = (pow(b,3) - 4.0*a*b*c + 8.0 *d*pow(a,2)) / (8.0*pow(a,3));
   r = (-3.0*pow(b,4) + 16.0 *a*pow(b,2)*c - 64.0 *pow(a,2)*b*d + 256.0 *pow(a,3)*f)/(256.0*pow(a,4));
   
   // Solve cubic resolvent
   AA[3] = 8.0;
   AA[2] = -4.0*p;
   AA[1] = -8.0*r;
   AA[0] = 4.0*p*r - pow(q,2);

   int ncube = cubic(AA, z);
   
   zsol = -1.e99;
   for(int i=0;i<ncube;i++)
      zsol = std::max(zsol, z[i]); //Not sure C has max fct
   z[0] = zsol;
   xK2 = 2.0*z[0] -p;
   xK = sqrt(xK2);
   xL = q/(2.0*xK);
   sqp = xK2 - 4.0 * (z[0] + xL);
   sqm = xK2 - 4.0 * (z[0] - xL);

   for(int i=0;i<4;i++)
      soli[i] = 0.0;
   
   if( (sqp >= 0.0) && (sqm >= 0.0)) {
      sol[0] = 0.5 * (xK + sqrt(sqp));
      sol[1] = 0.5 * (xK - sqrt(sqp));
      sol[2] = 0.5 * (-xK + sqrt(sqm));
      sol[3] = 0.5 * (-xK - sqrt(sqm));
      Nsol = 4;
   }
   else if((sqp >= 0.0) && (sqm < 0.0)) {
      sol[0] = 0.5 * (xK + sqrt(sqp));
      sol[1] = 0.5 * (xK - sqrt(sqp));
      sol[2] = -0.5 * xK;
      sol[3] = -0.5 * xK;
      soli[2] =  sqrt(-.25 * sqm);
      soli[3] = -sqrt(-.25 * sqm);
      Nsol = 2;
   }
   else if((sqp < 0.0) && (sqm >= 0.0)) {
      sol[0] = 0.5 * (-xK + sqrt(sqm));
      sol[1] = 0.5 * (-xK - sqrt(sqm));
      sol[2] = 0.5 * xK;
      sol[3] = 0.5 * xK;
      soli[2] =  sqrt(-0.25 * sqp);
      soli[3] = -sqrt(-0.25 * sqp);
      Nsol = 2;
   }
   else if((sqp < 0.0) && (sqm < 0.0)) {
      sol[0] = -0.5 * xK;
      sol[1] = -0.5 * xK;
      soli[0] =  sqrt(-0.25 * sqm);
      soli[1] = -sqrt(-0.25 * sqm);
      sol[2] = 0.5 * xK;
      sol[3] = 0.5 * xK;
      soli[2] =  sqrt(-0.25 * sqp);
      soli[3] = -sqrt(-0.25 * sqp);
      Nsol = 0;
   }
   
   for(int i=0;i<4;i++)
      sol[i] -= b/(4.0*a);

   return Nsol;
}




