/*
   Copyright (c) 1996-2000 Darko Kirovski, Miodrag Podkonjak and the Regents of
                           the University of California

   Contact author(s): darko@cs.ucla.edu, miodrag@cs.ucla.edu
   Affiliationss: UCLA, Computer Science
 
   Permission is hereby granted, free of charge, to any person obtaining 
   a copy of this software and associated documentation files (the "Software"),
   to deal in the Software without restriction, including without limitation 
   the rights to use, copy, modify, merge, publish, distribute, sublicense, 
   and/or sell copies of the Software, and to permit persons to whom the 
   Software is furnished to do so, subject to the following conditions:
 
   The above copyright notice and this permission notice shall be included
   in all copies or substantial portions of the Software.
 
   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, 
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
   OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
   IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, 
   DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR 
   OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
   THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

/*
   Name: color.h
   Description: minimal proper colouring
   Project: Antiprism - http://www.antiprism.com
   Changes: 16-04-04 Adrian Rossiter <adrian_r@terra.es>
      removed use of iostreams, graph initialised from a geometry object
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "geom.h"

#define max(x,y) (((x)>(y))?(x):(y))
#define min(x,y) (((x)<(y))?(x):(y))

enum Logic { _NO , _YES };
enum Status { _COLORED, _NEIGHBOR, _FREE, _SELECTED };

class Vertex;

struct Edge {
   Vertex *vert;
};

class Vertex {
   public:
      long COLOR;
      long EDGES;
      Edge *edge;
      Status STATUS;
      long NEIGHBORS;
      long SELECTED_NEIGHBORS;
      long COLORED_NEIGHBORS;
      float OBJECTIVE_VALUE;
      Vertex () {
         COLOR = -1; EDGES = 0; edge = 0; STATUS = _FREE; SELECTED_NEIGHBORS = 0;
         COLORED_NEIGHBORS = 0; OBJECTIVE_VALUE = 0.; NEIGHBORS = 0;
      }
      void SelectVertex ();
      void ColorVertex (long );
      Logic ConnectedToVertex (Vertex* );
};

class Graph {
   Vertex *vertex;
   long VERTICES;
   long RANDOM_1, RANDOM_2;
   float *rando;
   long LENGTH_RANDOM;
   long CURRENT_RANDOM;
   col_geom_v &geom;

   public:
   Graph (col_geom_v &dgeom);
   ~Graph () {
      for (long i = 0; i < VERTICES; i++)
         if (vertex[i].edge) delete vertex[i].edge;
      if (vertex) delete vertex;
      delete rando;
   }
   float Cost(Edge* , long );
   void SetValue(Vertex* , float* , float* );
   float LocalSearch(Edge* , long& , long* , Edge* , long );
   void GraphColoring (long* , long& );
   float Uniform();
   float Get() {
      CURRENT_RANDOM = (CURRENT_RANDOM + 1) % LENGTH_RANDOM;
      return(rando[CURRENT_RANDOM]);
   }
};

struct MaxIndSet {
   Edge *set;
   float COST;
   long SIZE;
   Logic FULL;
};

class Window {
   public:
      MaxIndSet *list;
      long ELEMENTS;
      float SMALLEST_COST;
      Window (long , long );
      ~Window ();
      void InsertMis (Edge* , long , float );
      void BestMis (Edge* , long& , float& );
      void ComputeCosts (Graph* );
};
