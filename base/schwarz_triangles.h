/*
   Copyright (c) 2012, Adrian Rossiter

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

#include "const.h"

//Schwarz triangles
static const int num_schwarz_tris = 44;

static int schwarz_triangles[num_schwarz_tris][6] =
{
   {2, 1, 3, 1, 3, 1}, //  0
   {2, 1, 3, 1, 3, 2}, //  1
   {2, 1, 3, 1, 4, 1}, //  2
   {2, 1, 3, 1, 4, 3}, //  3
   {2, 1, 3, 1, 5, 1}, //  4
   {2, 1, 3, 1, 5, 2}, //  5
   {2, 1, 3, 1, 5, 3}, //  6
   {2, 1, 3, 1, 5, 4}, //  7
   {2, 1, 3, 2, 3, 2}, //  8
   {2, 1, 3, 2, 4, 1}, //  9
   {2, 1, 3, 2, 4, 3}, // 10
   {2, 1, 3, 2, 5, 1}, // 11
   {2, 1, 3, 2, 5, 2}, // 12
   {2, 1, 3, 2, 5, 3}, // 13
   {2, 1, 3, 2, 5, 4}, // 14
   {2, 1, 5, 1, 5, 2}, // 15
   {2, 1, 5, 1, 5, 3}, // 16
   {2, 1, 5, 2, 5, 4}, // 17
   {2, 1, 5, 3, 5, 4}, // 18
   {3, 1, 3, 1, 3, 2}, // 19
   {3, 1, 3, 1, 5, 2}, // 20
   {3, 1, 3, 1, 5, 4}, // 21
   {3, 1, 3, 2, 5, 1}, // 22
   {3, 1, 3, 2, 5, 3}, // 23
   {3, 1, 4, 1, 4, 3}, // 24
   {3, 1, 5, 1, 5, 3}, // 25
   {3, 1, 5, 1, 5, 4}, // 26
   {3, 1, 5, 2, 5, 3}, // 27
   {3, 1, 5, 2, 5, 4}, // 28
   {3, 2, 3, 2, 3, 2}, // 29
   {3, 2, 3, 2, 5, 2}, // 30
   {3, 2, 3, 2, 5, 4}, // 31
   {3, 2, 4, 1, 4, 1}, // 32
   {3, 2, 4, 3, 4, 3}, // 33
   {3, 2, 5, 1, 5, 1}, // 34
   {3, 2, 5, 1, 5, 2}, // 35
   {3, 2, 5, 2, 5, 2}, // 36
   {3, 2, 5, 3, 5, 3}, // 37
   {3, 2, 5, 3, 5, 4}, // 38
   {3, 2, 5, 4, 5, 4}, // 39
   {5, 1, 5, 1, 5, 4}, // 40
   {5, 2, 5, 2, 5, 2}, // 41
   {5, 2, 5, 3, 5, 3}, // 42
   {5, 4, 5, 4, 5, 4}, // 43
};


static double schwarz_triangles_verts[num_schwarz_tris][9] = {
   { //  0
     0, 1, 0,
     1/sqrt_3, 1/sqrt_3, -1/sqrt_3,
     1/sqrt_3, 1/sqrt_3, 1/sqrt_3, },
   { //  1
     0, 1, 0,
     1/sqrt_3, -1/sqrt_3, 1/sqrt_3,
     1/sqrt_3, 1/sqrt_3, -1/sqrt_3, },
   { //  2
     1/sqrt_2, 1/sqrt_2, 0,
     1/sqrt_3, 1/sqrt_3, 1/sqrt_3,
     1, 0, 0, },
   { //  3
     1/sqrt_2, 1/sqrt_2, 0,
     -1/sqrt_3, -1/sqrt_3, -1/sqrt_3,
     1, 0, 0, },
   { //  4
     0.5/phi, phi/2, 0.5,
     1/sqrt_3, 1/sqrt_3, 1/sqrt_3,
     0, 1/sqrt_phi_plus_2, phi/sqrt_phi_plus_2, },
   { //  5
     0.5, -0.5/phi, phi/2,
     phi/sqrt_3, (phi-1)/sqrt_3, 0,
     0, 1/sqrt_phi_plus_2, phi/sqrt_phi_plus_2, },
   { //  6
     0, 1, 0,
     phi/sqrt_3, -(phi-1)/sqrt_3, 0,
     0, 1/sqrt_phi_plus_2, phi/sqrt_phi_plus_2, },
   { //  7
     0.5/phi, phi/2, 0.5,
     0, -phi/sqrt_3, -(phi-1)/sqrt_3,
     0, 1/sqrt_phi_plus_2, phi/sqrt_phi_plus_2, },
   { //  8
     0, -1, 0,
     1/sqrt_3, 1/sqrt_3, -1/sqrt_3,
     1/sqrt_3, 1/sqrt_3, 1/sqrt_3, },
   { //  9
     -1/sqrt_2, 1/sqrt_2, 0,
     -1/sqrt_3, 1/sqrt_3, -1/sqrt_3,
     1, 0, 0, },
   { // 10
     -1/sqrt_2, 0, -1/sqrt_2,
     1/sqrt_3, 1/sqrt_3, 1/sqrt_3,
     1, 0, 0, },
   { // 11
     0.5/phi, -phi/2, -0.5,
     1/sqrt_3, -1/sqrt_3, -1/sqrt_3,
     0, 1/sqrt_phi_plus_2, phi/sqrt_phi_plus_2, },
   { // 12
     0.5, 0.5/phi, -phi/2,
     phi/sqrt_3, -(phi-1)/sqrt_3, 0,
     0, 1/sqrt_phi_plus_2, phi/sqrt_phi_plus_2, },
   { // 13
     0, -1, 0,
     phi/sqrt_3, (phi-1)/sqrt_3, 0,
     0, 1/sqrt_phi_plus_2, phi/sqrt_phi_plus_2, },
   { // 14
     -0.5, -0.5/phi, -phi/2,
     1/sqrt_3, 1/sqrt_3, 1/sqrt_3,
     0, 1/sqrt_phi_plus_2, phi/sqrt_phi_plus_2, },
   { // 15
     0.5/phi, phi/2, 0.5,
     -1/sqrt_phi_plus_2, phi/sqrt_phi_plus_2, 0,
     0, 1/sqrt_phi_plus_2, phi/sqrt_phi_plus_2, },
   { // 16
     0.5/phi, phi/2, 0.5,
     1/sqrt_phi_plus_2, -phi/sqrt_phi_plus_2, 0,
     0, 1/sqrt_phi_plus_2, phi/sqrt_phi_plus_2, },
   { // 17
     0.5/phi, -phi/2, -0.5,
     0, 1/sqrt_phi_plus_2, phi/sqrt_phi_plus_2,
     -1/sqrt_phi_plus_2, -phi/sqrt_phi_plus_2, 0, },
   { // 18
     0.5/phi, -phi/2, -0.5,
     0, 1/sqrt_phi_plus_2, phi/sqrt_phi_plus_2,
     1/sqrt_phi_plus_2, phi/sqrt_phi_plus_2, 0, },
   { // 19
     1/sqrt_3, 1/sqrt_3, -1/sqrt_3,
     1/sqrt_3, -1/sqrt_3, 1/sqrt_3,
     1/sqrt_3, 1/sqrt_3, 1/sqrt_3, },
   { // 20
     1/sqrt_3, 1/sqrt_3, 1/sqrt_3,
     (phi-1)/sqrt_3, 0, phi/sqrt_3,
     0, 1/sqrt_phi_plus_2, phi/sqrt_phi_plus_2, },
   { // 21
     phi/sqrt_3, (phi-1)/sqrt_3, 0,
     -1/sqrt_3, -1/sqrt_3, 1/sqrt_3,
     0, 1/sqrt_phi_plus_2, phi/sqrt_phi_plus_2, },
   { // 22
     phi/sqrt_3, -(phi-1)/sqrt_3, 0,
     phi/sqrt_3, (phi-1)/sqrt_3, 0,
     0, 1/sqrt_phi_plus_2, phi/sqrt_phi_plus_2, },
   { // 23
     0, -phi/sqrt_3, -(phi-1)/sqrt_3,
     1/sqrt_3, 1/sqrt_3, 1/sqrt_3,
     0, 1/sqrt_phi_plus_2, phi/sqrt_phi_plus_2, },
   { // 24
     1/sqrt_3, 1/sqrt_3, -1/sqrt_3,
     0, 0, 1,
     0, 1, 0, },
   { // 25
     1/sqrt_3, 1/sqrt_3, 1/sqrt_3,
     0, -1/sqrt_phi_plus_2, phi/sqrt_phi_plus_2,
     0, 1/sqrt_phi_plus_2, phi/sqrt_phi_plus_2, },
   { // 26
     1/sqrt_3, 1/sqrt_3, 1/sqrt_3,
     -1/sqrt_phi_plus_2, -phi/sqrt_phi_plus_2, 0,
     0, 1/sqrt_phi_plus_2, phi/sqrt_phi_plus_2, },
   { // 27
     phi/sqrt_3, -(phi-1)/sqrt_3, 0,
     0, 1/sqrt_phi_plus_2, phi/sqrt_phi_plus_2,
     1/sqrt_phi_plus_2, phi/sqrt_phi_plus_2, 0, },
   { // 28
     1/sqrt_3, -1/sqrt_3, -1/sqrt_3,
     0, 1/sqrt_phi_plus_2, phi/sqrt_phi_plus_2,
     1/sqrt_phi_plus_2, phi/sqrt_phi_plus_2, 0, },
   { // 29
     1/sqrt_3, 1/sqrt_3, -1/sqrt_3,
     1/sqrt_3, -1/sqrt_3, 1/sqrt_3,
     -1/sqrt_3, 1/sqrt_3, 1/sqrt_3, },
   { // 30
     1/sqrt_3, -1/sqrt_3, -1/sqrt_3,
     (phi-1)/sqrt_3, 0, -phi/sqrt_3,
     0, 1/sqrt_phi_plus_2, phi/sqrt_phi_plus_2, },
   { // 31
     phi/sqrt_3, -(phi-1)/sqrt_3, 0,
     -1/sqrt_3, 1/sqrt_3, -1/sqrt_3,
     0, 1/sqrt_phi_plus_2, phi/sqrt_phi_plus_2, },
   { // 32
     1/sqrt_3, 1/sqrt_3, -1/sqrt_3,
     0, 1, 0,
     1, 0, 0, },
   { // 33
     1/sqrt_3, 1/sqrt_3, -1/sqrt_3,
     0, -1, 0,
     0, 0, 1, },
   { // 34
     1/sqrt_3, 1/sqrt_3, 1/sqrt_3,
     0, 1/sqrt_phi_plus_2, phi/sqrt_phi_plus_2,
     1/sqrt_phi_plus_2, phi/sqrt_phi_plus_2, 0, },
   { // 35
     1/sqrt_3, 1/sqrt_3, 1/sqrt_3,
     1/sqrt_phi_plus_2, -phi/sqrt_phi_plus_2, 0,
     0, 1/sqrt_phi_plus_2, phi/sqrt_phi_plus_2, },
   { // 36
     phi/sqrt_3, (phi-1)/sqrt_3, 0,
     0, 1/sqrt_phi_plus_2, phi/sqrt_phi_plus_2,
     1/sqrt_phi_plus_2, -phi/sqrt_phi_plus_2, 0, },
   { // 37
     1/sqrt_3, 1/sqrt_3, -1/sqrt_3,
     0, 1/sqrt_phi_plus_2, phi/sqrt_phi_plus_2,
     1/sqrt_phi_plus_2, -phi/sqrt_phi_plus_2, 0, },
   { // 38
     1/sqrt_3, -1/sqrt_3, -1/sqrt_3,
     0, 1/sqrt_phi_plus_2, phi/sqrt_phi_plus_2,
     0, 1/sqrt_phi_plus_2, -phi/sqrt_phi_plus_2, },
   { // 39
     1/sqrt_3, -1/sqrt_3, -1/sqrt_3,
     0, 1/sqrt_phi_plus_2, phi/sqrt_phi_plus_2,
     -1/sqrt_phi_plus_2, phi/sqrt_phi_plus_2, 0, },
   { // 40
     0, 1/sqrt_phi_plus_2, phi/sqrt_phi_plus_2,
     1/sqrt_phi_plus_2, -phi/sqrt_phi_plus_2, 0,
     0, -1/sqrt_phi_plus_2, phi/sqrt_phi_plus_2, },
   { // 41
     0, 1/sqrt_phi_plus_2, phi/sqrt_phi_plus_2,
     1/sqrt_phi_plus_2, phi/sqrt_phi_plus_2, 0,
     -1/sqrt_phi_plus_2, phi/sqrt_phi_plus_2, 0, },
   { // 42
     1/sqrt_phi_plus_2, -phi/sqrt_phi_plus_2, 0,
     0, 1/sqrt_phi_plus_2, phi/sqrt_phi_plus_2,
     1/sqrt_phi_plus_2, phi/sqrt_phi_plus_2, 0, },
   { // 43
     0, 1/sqrt_phi_plus_2, phi/sqrt_phi_plus_2,
     1/sqrt_phi_plus_2, -phi/sqrt_phi_plus_2, 0,
     0, 1/sqrt_phi_plus_2, -phi/sqrt_phi_plus_2, },
};


