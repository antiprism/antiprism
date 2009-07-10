/*
   Copyright (c) 1996-2000 Darko Kirovski, Miodrag Podkonjak and the Regents of
                           the University of California

   Contact author(s) for original code: darko@cs.ucla.edu, miodrag@cs.ucla.edu
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

#include "prop_col.h"

Graph::Graph (col_geom_v &dgeom): geom(dgeom)
{
   VERTICES = geom.verts().size();
   vertex = new Vertex[VERTICES];
   
   // Structures
   struct ListEdge {
      Vertex *pointer;
      ListEdge *next;
   };
   struct ListVertex {
      ListEdge *first;
      long edges;
   } *arrayE = new ListVertex[VERTICES];
   long i;
   for (i = 0; i < VERTICES; i++) {
      arrayE[i].first = 0;
      arrayE[i].edges = 0;
   }  
   
   // Edge creation
   vector<vector<int> > edges;
   geom.get_impl_edges(edges);
   for (unsigned int k = 0; k < edges.size(); k++) {
      int i = edges[k][0];
      int j = edges[k][1];
      ListEdge *help = new ListEdge;
      help->pointer = &vertex[j];
      help->next = arrayE[i].first;
      arrayE[i].first = help;
      help = new ListEdge;
      help->pointer = &vertex[i];
      help->next = arrayE[j].first;
      arrayE[j].first = help;
      arrayE[i].edges++;
      arrayE[j].edges++;
   }
   
   // Creating dynamic arrays of edges
   for (i = 0; i < VERTICES; i++) {
      vertex[i].COLOR = -1;
      vertex[i].EDGES = arrayE[i].edges;
      if (arrayE[i].edges) vertex[i].edge = new Edge[arrayE[i].edges];
      long count = 0;
      while (arrayE[i].first) {
         vertex[i].edge[count++].vert = arrayE[i].first->pointer;
         ListEdge *help = arrayE[i].first;
         arrayE[i].first = arrayE[i].first->next;
         delete help;
      }
   }
   delete arrayE;
   RANDOM_1 = ((rand() >> 1) << 1) + 1;
   RANDOM_2 = ((rand() >> 1) << 1) + 1;
   LENGTH_RANDOM = 10 * VERTICES;
   rando = new float[LENGTH_RANDOM];
   for (i = 0; i < LENGTH_RANDOM; i++)
      rando[i] = Uniform();
   CURRENT_RANDOM = 0;
}


/************************* RANDOM *************************/
float Graph::Uniform() {
   long Z, k;
   k = RANDOM_1 / 53668;
   RANDOM_1 = 40014 * (RANDOM_1 - k * 53668) - k * 12211;
   if (RANDOM_1 < 0)
      RANDOM_1 = RANDOM_1 + 2147483563 ;
   k = RANDOM_2 / 52774;
   RANDOM_2 = 40692 * (RANDOM_2 - k * 52774) - k * 3791;
   if (RANDOM_2 < 0)
      RANDOM_2 = RANDOM_2 + 2147483399;
   Z = RANDOM_1 - RANDOM_2;
   if (Z < 1)
      Z = Z + 2147483562;
   return (float(Z) * (float) 4.656613E-10);
}
/************************* END RANDOM *************************/

/************************* OBJECTIVE *************************/
void Graph::SetValue(Vertex *vert, float *square, float *cube) {
   vert->OBJECTIVE_VALUE = 0.;
   for (long j = 0; j < vert->EDGES; j++)
      if (vert->edge[j].vert->STATUS == _FREE)
         vert->OBJECTIVE_VALUE += cube[vert->edge[j].vert->NEIGHBORS];
   vert->OBJECTIVE_VALUE *= square[vert->NEIGHBORS];
}

float Graph::Cost(Edge *MIS,
      long TOP_MIS) {
   float result = 0.;
   for (long i = 0; i <= TOP_MIS; i++)
      result += MIS[i].vert->OBJECTIVE_VALUE;
   return (result);
}
/************************* END OBJECTIVE *************************/

/*************************************/
/*************************************/
/***  G R A P H   C O L O R I N G  ***/
/*************************************/
/*************************************/
void Graph::GraphColoring (long *parameter, long& current_color) {
   current_color = 0;
   long number_of_colored = 0;
   //	float best_cost = 0.;
   Window FRAME(VERTICES, parameter[3]);
   Edge *GRAPH = new Edge[VERTICES];
   long TOP_GRAPH = -1;
   Edge *MIS = new Edge[VERTICES];
   long TOP_MIS = -1;
   Edge *FREES = new Edge[VERTICES];
   long TOP_FREES = -1;
   Edge *DEL = new Edge[VERTICES];
   long TOP_DEL = -1;
   float *SINGLES = new float[VERTICES + 1];
   float *SQUARES = new float[VERTICES + 1];
   float *CUBES = new float[VERTICES + 1];
   long i;
   for (i = 0; i <= VERTICES; i++) {
      SINGLES[i] = float(i);
      SQUARES[i] = SINGLES[i] * SINGLES[i];
      CUBES[i] = SQUARES[i] * SINGLES[i];
   }
   Logic FIRST_TRY = _YES;

   long EdgesLeft;

   /********** ADD NEW COLOR **********/

   do {

      /********** INITIALIZATION **********/

      EdgesLeft = 0;
      TOP_GRAPH = -1;
      long i;
      for (i = 0; i < VERTICES; i++)
         if (vertex[i].STATUS != _COLORED) {
            TOP_GRAPH++;
            GRAPH[TOP_GRAPH].vert = &vertex[i];
            vertex[i].STATUS = _FREE;
            vertex[i].SELECTED_NEIGHBORS = 0;
            vertex[i].NEIGHBORS = vertex[i].EDGES - vertex[i].COLORED_NEIGHBORS;
            EdgesLeft += vertex[i].NEIGHBORS;
         }
      EdgesLeft = EdgesLeft >> 1;
      for (i = 0; i <= TOP_GRAPH; i++)
         SetValue(GRAPH[i].vert, SQUARES, CUBES);
      FRAME.ComputeCosts(this);

      /**********************************/
      //if (current_color > 0) {
      //   cout << "\tC=" << current_color << " Clrd=" << number_of_colored << "                                        \n";
      //   cout.flush();
      //}
      /**********************************/

      /********** INITIAL INDEPENDENT SET **********/

      TOP_MIS = -1; float best_cost = 0.;
      for (i = 0; i <= TOP_GRAPH; i++)
         if (GRAPH[i].vert->STATUS == _FREE) {
            TOP_MIS++;
            MIS[TOP_MIS].vert = GRAPH[i].vert;
            GRAPH[i].vert->SelectVertex();
            best_cost += GRAPH[i].vert->OBJECTIVE_VALUE;
         }
      FRAME.InsertMis(MIS, TOP_MIS, best_cost);

      /**********************************/
      //cout << "\t\tInitial MIS=" << TOP_MIS + 1
      //   << " COST=" << best_cost << "\r";
      //cout.flush();
      /**********************************/

      if (VERTICES - number_of_colored - TOP_MIS - 1 > 0) {
         long limit = long(parameter[0] *
               (1. - SQUARES[number_of_colored] / SQUARES[VERTICES])
               + parameter[1]);
         if (FIRST_TRY == _YES) {
            limit *= parameter[3];
            FIRST_TRY = _NO;
         }

         /********** RANDOMIZED INDEPENDENT SETS **********/
         long counter = 0;
         long no_improvement = 0;
         float current_cost;
         if (limit > 0) do {

            /********** RANDOMIZED VERTEX EXCLUSION **********/
            TOP_DEL = TOP_FREES = -1;
            counter++;
            do {
               long index = long((TOP_MIS + 1) * Get());
               MIS[index].vert->STATUS = _FREE;
               MIS[index].vert->SELECTED_NEIGHBORS = 0;
               TOP_DEL++; DEL[TOP_DEL].vert = MIS[index].vert;
               long i;
               for (i = 0; i < MIS[index].vert->EDGES; i++)
                  if (MIS[index].vert->edge[i].vert->STATUS == _NEIGHBOR) {
                     MIS[index].vert->edge[i].vert->SELECTED_NEIGHBORS--;
                     if (MIS[index].vert->edge[i].vert->SELECTED_NEIGHBORS == 0) {
                        TOP_FREES++;
                        FREES[TOP_FREES].vert = MIS[index].vert->edge[i].vert;
                        MIS[index].vert->edge[i].vert->STATUS = _FREE;
                     }
                  }
               for (i = index; i < TOP_MIS; i++)
                  MIS[i].vert = MIS[i + 1].vert;
               TOP_MIS--;
            } while (TOP_FREES == -1);

            /********** RANDOMIZED FREES INCLUSION **********/
            Logic first = _YES;
            do {
               float maximum = -1.;
               long i, index = -1;
               for (i = 0; i <= TOP_FREES; i++)
                  if (FREES[i].vert->STATUS == _FREE) {
                     float temp = Get();
                     if (temp > maximum) { maximum = temp; index = i; }
                  }
                  else {
                     for (long j = i; j < TOP_FREES; j++)
                        FREES[j].vert = FREES[j + 1].vert;
                     TOP_FREES--;
                  }
               if (index > -1) {
                  FREES[index].vert->SelectVertex();
                  TOP_MIS++; MIS[TOP_MIS].vert = FREES[index].vert;
                  long i;
                  for (i = index; i < TOP_FREES; i++)
                     FREES[i].vert = FREES[i + 1].vert;
                  TOP_FREES--;
               }
               if (first == _YES) {
                  first = _NO;
                  long i;
                  for (i = 0; i <= TOP_DEL; i++)
                     if (DEL[i].vert->STATUS == _FREE) {
                        DEL[i].vert->STATUS = _FREE;
                        TOP_FREES++; FREES[TOP_FREES].vert = DEL[i].vert;
                     }
               }
            } while (TOP_FREES != -1);

            /********** MEMORIZE BEST INDEPENDENT SET **********/
            current_cost = Cost(MIS, TOP_MIS);
            if (current_cost > FRAME.SMALLEST_COST)
               FRAME.InsertMis(MIS, TOP_MIS, current_cost);
            if (current_cost > best_cost) {
               best_cost = current_cost;
               no_improvement = 0;

               /**********************************/
               //cout << "\t\tcounter=" << counter << " MIS=" << TOP_MIS + 1
               //   << " COST=" << current_cost << "\r";
               //cout.flush();
               /**********************************/

            }
            else no_improvement++;
         } while (no_improvement < limit);

         /********** LOCAL SEARCH **********/ 
         long i;
         for (i = 0; i < FRAME.ELEMENTS; i++)
            if (FRAME.list[i].FULL == _YES)
               FRAME.list[i].COST =
                  LocalSearch (&(FRAME.list[i].set[0]),
                        FRAME.list[i].SIZE,
                        parameter,
                        GRAPH,
                        TOP_GRAPH);
      }

      /********** COLOR THE BEST INDEPENDENT SET **********/
      FRAME.BestMis(MIS, TOP_MIS, best_cost);
      for(i = 0; i <= TOP_MIS; i++)
         MIS[i].vert->ColorVertex(current_color);
      number_of_colored += TOP_MIS + 1;
      current_color++;

   } while (number_of_colored < VERTICES);
   //cout << "\tC=" << current_color << " Clrd=" << number_of_colored << "                                        \n";
   //cout.flush();
   //cout << "ColorsUsed=" << current_color << "                                      \n";
   //cout.flush();
   for (i = 0; i < VERTICES; i++) {
      if (vertex[i].EDGES != vertex[i].COLORED_NEIGHBORS)
        fprintf(stderr, "\noff_color: warning: mistake in proper colouring!\n");
      //printf("%ld = %ld\n", i, vertex[i].COLOR);
      geom.set_v_col(i, vertex[i].COLOR);

/*      for (long j = 0; j < vertex[i].EDGES; j++)
         if (!(vertex[i].COLOR < current_color && vertex[i].COLOR >= 0) ||
               vertex[i].COLOR == vertex[i].edge[j].vert->COLOR) {
            fprintf(stderr, "\noff_color: error: while calculating proper coloring\n");
            exit (-1);
         }
*/
   }
   delete MIS; delete FREES; delete DEL; delete SQUARES;
   delete SINGLES; delete CUBES; delete GRAPH;
}

float Graph::LocalSearch(Edge *BEST_MIS,
      long& BEST_TOP_MIS,
      long *parameter,
      Edge *GRAPH,
      long TOP_GRAPH) {
   Edge *MIS = new Edge[VERTICES];
   long TOP_MIS = -1;
   Edge *FREES = new Edge[VERTICES];
   long TOP_FREES = -1;
   Edge *DEL = new Edge[VERTICES];
   long TOP_DEL = -1;
   long counter = 0, iterations = 0;
   float best_cost=0., current_cost;
   Logic done = _YES;

   /*
      cout << "Search starts at";
      for (long i = 0; i <= BEST_TOP_MIS; i++)
      cout << " " << BEST_MIS[i].vert - &vertex[0];
    */

   Logic improvement;
   do {
      improvement = _NO;
      counter = 0;
      iterations = 0;
      do {

         /************** INITIALIZATION **************/
         iterations++;
         counter++;
         current_cost = 0.;
         long i;
         for (i = 0; i <= TOP_GRAPH; i++) {
            GRAPH[i].vert->STATUS = _FREE;
            GRAPH[i].vert->SELECTED_NEIGHBORS = 0;
         }
         for (i = 0; i <= BEST_TOP_MIS; i++) {
            BEST_MIS[i].vert->SelectVertex();
            MIS[i].vert = BEST_MIS[i].vert;
            current_cost += BEST_MIS[i].vert->OBJECTIVE_VALUE;
         }
         TOP_MIS = BEST_TOP_MIS;
         if (done == _YES) { done = _NO; best_cost = current_cost; }

         /************** RANDOMIZED VERTEX EXCLUSION **************/
         TOP_DEL = TOP_FREES = -1;
         do {
            long index = long((TOP_MIS + 1) * Get());
            MIS[index].vert->STATUS = _FREE;
            MIS[index].vert->SELECTED_NEIGHBORS = 0;
            TOP_DEL++; DEL[TOP_DEL].vert = MIS[index].vert;
            long i;
            for (i = 0; i < MIS[index].vert->EDGES; i++)
               if (MIS[index].vert->edge[i].vert->STATUS == _NEIGHBOR) {
                  MIS[index].vert->edge[i].vert->SELECTED_NEIGHBORS--;
                  if (MIS[index].vert->edge[i].vert->SELECTED_NEIGHBORS == 0) {
                     TOP_FREES++;
                     FREES[TOP_FREES].vert = MIS[index].vert->edge[i].vert;
                     MIS[index].vert->edge[i].vert->STATUS = _FREE;
                  }
               }
            for (i = index; i < TOP_MIS; i++) MIS[i].vert = MIS[i + 1].vert;
            TOP_MIS--;
         } while (TOP_FREES == -1);

         /************** RANDOMIZED VERTEX INCLUSION **************/
         Logic first = _YES;
         do {
            float maximum = 0.; long index = -1;
            long i;
            for (i = 0; i <= TOP_FREES; i++)
               if (FREES[i].vert->STATUS == _FREE) {
                  float temp = Get();
                  if (temp > maximum) { maximum = temp; index = i; }
               }
               else {
                  for (long j = i; j < TOP_FREES; j++)
                     FREES[j].vert = FREES[j + 1].vert;
                  TOP_FREES--;
               }
            if (index > -1) {
               FREES[index].vert->SelectVertex();
               TOP_MIS++; MIS[TOP_MIS].vert = FREES[index].vert;
               long i;
               for (i = index; i < TOP_FREES; i++)
                  FREES[i].vert = FREES[i + 1].vert;
               TOP_FREES--;
            }
            if (first == _YES) {
               first = _NO;
               long i;
               for (i = 0; i <= TOP_DEL; i++)
                  if (DEL[i].vert->STATUS == _FREE) {
                     DEL[i].vert->STATUS = _FREE;
                     TOP_FREES++; FREES[TOP_FREES].vert = DEL[i].vert;
                  }
            }
         } while (TOP_FREES != -1);

         /************** MEMORIZE BEST INDEPENDENT SET **************/
         current_cost = Cost(MIS, TOP_MIS);
         if (current_cost > best_cost) {
            improvement = _YES;
            best_cost = current_cost;
            BEST_TOP_MIS = TOP_MIS; 
            counter = 0;
            long i;
            for (i = 0; i <= TOP_MIS; i++) 
               BEST_MIS[i].vert = MIS[i].vert;

            /**********************************/
            //cout << "\t\tlocal_counter=" << iterations << " MIS=" << TOP_MIS + 1 
            //   << " COST=" << current_cost << "\r";
            //cout.flush();
            /**********************************/

         }
      } while (counter < parameter[2]);
   } while (improvement == _YES);
   delete MIS;
   delete FREES;
   delete DEL;
   return (best_cost);
}

/****************/
/***  Vertex  ***/
/****************/
void Vertex::SelectVertex () {
   long i;
   for (i = 0; i < EDGES; i++) {
      if (edge[i].vert->STATUS == _FREE) edge[i].vert->STATUS = _NEIGHBOR;
      if (edge[i].vert->STATUS == _NEIGHBOR) edge[i].vert->SELECTED_NEIGHBORS++;
   }
   STATUS = _SELECTED;
}

void Vertex::ColorVertex (long color) {
   COLOR = color;
   STATUS = _COLORED;
   long i;
   for (i = 0; i < EDGES; i++) edge[i].vert->COLORED_NEIGHBORS++;
}

Logic Vertex::ConnectedToVertex(Vertex *guy) {
   long i;
   for (i = 0; i < EDGES; i++) 
      if (edge[i].vert == guy) return(_YES);
   return(_NO);
}

/****************/
/***  Window  ***/
/****************/
Window::Window(long size, 
      long elem) {
   ELEMENTS = elem;
   list = new MaxIndSet[ELEMENTS];
   SMALLEST_COST = -1.;
   long i;
   for (i = 0; i < ELEMENTS; i++) {
      list[i].COST = -1.;
      list[i].SIZE = -1;
      list[i].set = new Edge[size];
      list[i].FULL = _NO;
   }
}

Window::~Window() {
   long i;
   for (i = 0; i < ELEMENTS; i++) 
      delete list[i].set;
   delete list;
}

void Window::InsertMis (Edge *MIS,
      long TOP_MIS,
      float cost) {
   long i;
   for (i = 0; i < ELEMENTS; i++)
      if (list[i].FULL == _YES)
         if ((list[i].COST - cost) * (list[i].COST - cost) < 1.) return;
   long smallest_index=-1; // originally uninitialised - Adrian
   for (i = 0; i < ELEMENTS; i++)
      if (list[i].FULL == _NO) { smallest_index = i; break; }
      else if (list[i].COST == SMALLEST_COST) { smallest_index = i; break; }
      SMALLEST_COST = list[smallest_index].COST = cost;
      list[smallest_index].FULL = _YES;
      list[smallest_index].SIZE = TOP_MIS;
      for (i = 0; i <= TOP_MIS; i++)
         list[smallest_index].set[i].vert = MIS[i].vert;
      for (i = 0; i < ELEMENTS; i++)
         if (list[i].FULL == _YES && list[i].COST < SMALLEST_COST)
            SMALLEST_COST = list[i].COST;
         else if (list[i].FULL == _NO) { SMALLEST_COST = -1.; return; }
}

void Window::BestMis (Edge *MIS,
      long& TOP_MIS,
      float& cost) {
   float maximum = -1.;
   long index = -1;
   long i;
   for (i = 0; i < ELEMENTS; i++)
      if (list[i].FULL == _YES) {
         if (list[i].COST > maximum) { index = i; maximum = list[i].COST; }
         else if (maximum - list[i].COST < 0.005 * maximum &&
               list[i].SIZE > list[index].SIZE)
         { index = i; maximum = list[i].COST; }
      }
   for (i = 0; i <= list[index].SIZE; i++)
      MIS[i].vert = list[index].set[i].vert;
   TOP_MIS = list[index].SIZE;
   cost = list[index].COST;
   list[index].FULL = _NO;
   for (i = 0; i < ELEMENTS; i++)
      if (list[i].FULL == _YES) {
         Logic to_delete = _NO;
         for (long j = 0; j <= TOP_MIS; j++) {
            for (long k = 0; k <= list[i].SIZE; k++)
               if (list[i].set[k].vert == MIS[j].vert)
               { to_delete = _YES; break; }
            if (to_delete == _YES) break;
         }
         if (to_delete == _YES) list[i].FULL = _NO;
      }
   SMALLEST_COST = -1.;
}

void Window::ComputeCosts (Graph *pointer) {
   Logic check = _NO;
   long i;
   for (i = 0; i < ELEMENTS; i++)
      if (list[i].FULL == _YES)
         list[i].COST = pointer->Cost(&(list[i].set[0]), list[i].SIZE);
      else check = _YES;
      if (check == _NO) { fprintf(stderr, "\nError_1!\n"); exit (-1); }
      SMALLEST_COST = -1.;
}
