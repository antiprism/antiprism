/*
 * The author of this software is Michael Trick.  Copyright (c) 1994 by
 * Michael Trick.
 *
 * Permission to use, copy, modify, and distribute this software for any
 * purpose without fee is hereby granted, provided that this entire notice
 * is included in all copies of any software which is or includes a copy
 * or modification of this software and in all copies of the supporting
 * documentation for such software.
 * THIS SOFTWARE IS BEING PROVIDED "AS IS", WITHOUT ANY EXPRESS OR IMPLIED
 * WARRANTY.  IN PARTICULAR, NEITHER THE AUTHOR DOES NOT MAKE ANY
 * REPRESENTATION OR WARRANTY OF ANY KIND CONCERNING THE MERCHANTABILITY
 * OF THIS SOFTWARE OR ITS FITNESS FOR ANY PARTICULAR PURPOSE.
 */

/*
   COLOR.C: Easy code for graph Coloring
   Author: Michael A. Trick, Carnegie Mellon University, trick+@cmu.edu
   Last Modified: November 2, 1994
   http://mat.gsia.cmu.edu/COLOR/solvers/trick.c

Graph is input in a file.  First line contains the number of nodes and
edges.  All following contain the node numbers (from 1 to n) incident to
each edge.  Sample:

4 4
1 2
2 3
3 4
1 4

represents a four node cycle graph.

Code is probably insufficiently debugged, but may be useful to some people.

For more information on this code, see Anuj Mehrotra and Michael A. Trick,
"A column generation approach to graph Coloring", GSIA Technical report series.

*/

/*
   Name: prop_col.cc
   Description: minimal proper colouring
   Project: Antiprism - http://www.antiprism.com
   Changes: 112/09/11 Adrian Rossiter <adrian@antiprism.com>
      convert to class and use STL conatiners rather than arrays
*/

#include <set>
#include <utility>
#include <vector>

#include "private_prop_col.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

using std::pair;
using std::set;
using std::vector;

int ProperColor::find_colors()
{
  prob_count = 0;
  visit_cnt.resize(num_node + 1, 0);
  ColorAdj.resize(num_node, vector<int>(num_node + 1, 0));
  set<pair<int, int>>::const_iterator si;
  for (si = adj.begin(); si != adj.end(); ++si) {

    ColorAdj[si->first][0]++;
    ColorAdj[si->second][0]++;
  }

  ColorCount.resize(num_node, 0);
  ColorClass.resize(num_node, 0);
  Handled.resize(num_node, false);
  Order.resize(num_node, 0);
  BestColoring = num_node + 1;

  int valid[num_node], clique[num_node];
  for (int i = 0; i < num_node; i++)
    valid[i] = true;

  best_clique = 0;
  num_prob = 0;
  max_prob = 10000;

  lb = max_w_clique(valid, clique, 0, num_node);

  int place = 0;

  for (int i = 0; i < num_node; i++) {
    if (clique[i]) {
      Order[place] = i;
      Handled[i] = true;
      place++;
      assign_color(i, place);
      for (int j = 0; j < num_node; j++)
        if ((i != j) && clique[j] && (!is_adj(i, j)))
          fprintf(stderr,
                  "warning: proper colouring, result is not a clique\n");
    }
  }

  if (color(place, place) == 0) {
    ColorClass.clear();
    Order.clear();
    Handled.clear();
    ColorAdj.clear();
    ColorCount.clear();
    visit_cnt.clear();

    long parameter[] = {1000, 10, 50, 5};
    long colours;
    Graph g(*this);
    adj.clear();
    g.GraphColoring(parameter, colours);
  }

  if ((int)BestColorClass.size() < num_node) {
    BestColorClass.resize(num_node);
    for (int i = 0; i < num_node; i++)
      BestColorClass[i] = i + 1;
  }

  return 0;
}

int ProperColor::greedy_clique(int *valid, int *clique)
{
  for (int i = 0; i < num_node; i++)
    clique[i] = 0;

  vector<int> order(num_node + 1, 0);
  int place = 0;
  for (int i = 0; i < num_node; i++) {
    if (valid[i]) {
      order[place] = i;
      place++;
    }
  }

  vector<int> weight(num_node, 0);
  for (int i = 0; i < num_node; i++) {
    if (!valid[i])
      continue;
    for (int j = 0; j < num_node; j++) {
      if (!valid[j])
        continue;
      if (is_adj(i, j))
        weight[i]++;
    }
  }

  bool done = false;
  while (!done) {
    done = true;
    for (int i = 0; i < place - 1; i++) {
      int j = order[i];
      int k = order[i + 1];
      if (weight[j] < weight[k]) {
        order[i] = k;
        order[i + 1] = j;
        done = false;
      }
    }
  }

  clique[order[0]] = true;
  for (int i = 1; i < place; i++) {
    int j = order[i];
    int k;
    for (k = 0; k < i; k++)
      if (clique[order[k]] && !is_adj(j, order[k]))
        break;
    if (k == i)
      clique[j] = true;
    else
      clique[j] = false;
  }
  int max = 0;
  for (int i = 0; i < place; i++)
    if (clique[order[i]])
      max++;

  return max;
}

/* Target is a goal value:  once a clique is found with value target
   it is possible to return

   Lower is a bound representing an already found clique:  once it is
   determined that no clique exists with value better than lower, it
   is permitted to return with a suboptimal clique.

   Note, to find a clique of value 1, it is not permitted to just set
   the lower to 1:  the recursion will not work.  Lower represents a
   value that is the goal for the recursion.
   */

int ProperColor::max_w_clique(int *valid, int *clique, int lower, int target)
{
  /*  printf("entered with lower %d target %d\n",lower,target);*/
  num_prob++;
  if (num_prob > max_prob)
    return -1;

  for (int j = 0; j < num_node; j++)
    clique[j] = 0;

  int total_left = 0;
  for (int i = 0; i < num_node; i++)
    if (valid[i])
      total_left++;
  if (total_left < lower)
    return 0;

  int incumb = greedy_clique(valid, clique);
  if (incumb >= target)
    return incumb;
  if (incumb > best_clique) {
    best_clique = incumb;
    /*    printf("Clique of size %5d found.\n",best_clique);*/
  }
  /*  printf("Greedy gave %f\n",incumb);*/

  vector<int> order(num_node, 0);
  int place = 0;
  for (int i = 0; i < num_node; i++) {
    if (clique[i]) {
      order[place] = i;
      total_left--;
      place++;
    }
  }
  int start = place;
  for (int i = 0; i < num_node; i++) {
    if (!clique[i] && valid[i]) {
      order[place] = i;
      place++;
    }
  }

  vector<int> value(num_node, 0);
  int finish = place;
  for (place = start; place < finish; place++) {
    int i = order[place];
    value[i] = 0;
    for (int j = 0; j < num_node; j++) {
      if (valid[j] && is_adj(i, j))
        value[i]++;
    }
  }

  bool done = false;
  while (!done) {
    done = true;
    for (place = start; place < finish - 1; place++) {
      int i = order[place];
      int j = order[place + 1];
      if (value[i] < value[j]) {
        order[place] = j;
        order[place + 1] = i;
        done = false;
      }
    }
  }

  for (place = start; place < finish; place++) {
    if (incumb + total_left < lower)
      return 0;
    int j = order[place];
    total_left--;

    if (clique[j])
      continue;

    int valid1[num_node];
    for (int place1 = 0; place1 < num_node; place1++)
      valid1[place1] = false;
    for (int place1 = 0; place1 < place; place1++) {
      int k = order[place1];
      if (valid[k] && is_adj(j, k))
        valid1[k] = true;
      else
        valid1[k] = false;
    }
    int clique1[num_node];
    int new_weight = max_w_clique(valid1, clique1, incumb - 1, target - 1);
    if (new_weight + 1 > incumb) {
      /*      printf("Taking new\n");*/
      incumb = new_weight + 1;
      for (int k = 0; k < num_node; k++)
        clique[k] = clique1[k];
      clique[j] = true;
      if (incumb > best_clique) {
        best_clique = incumb;
        /*	printf("Clique of size %5d found.\n",best_clique);*/
      }
    }

    // else fprintf(stderr, "Taking incumb\n");
    if (incumb >= target)
      break;
  }
  return (incumb);
}

void ProperColor::assign_color(int node, int color)
{
  // fprintf(stderr, "  %d color +%d\n",node,color);
  ColorClass[node] = color;
  for (int node1 = 0; node1 < num_node; node1++) {
    if (node == node1)
      continue;
    if (is_adj(node, node1)) {
      if (ColorAdj[node1][color] == 0)
        ColorCount[node1]++;
      ColorAdj[node1][color]++;
      ColorAdj[node1][0]--;
      if (ColorAdj[node1][0] < 0)
        fprintf(stderr, "warning: proper colouring, error setting colour\n");
    }
  }
}

void ProperColor::remove_color(int node, int color)
{
  // fprintf(stderr, "  %d color -%d\n",node,color);
  ColorClass[node] = 0;
  for (int node1 = 0; node1 < num_node; node1++) {
    if (node == node1)
      continue;
    if (is_adj(node, node1)) {
      ColorAdj[node1][color]--;
      if (ColorAdj[node1][color] == 0)
        ColorCount[node1]--;
      if (ColorAdj[node1][color] < 0)
        fprintf(stderr, "warning: proper colouring, error setting colour\n");
      ColorAdj[node1][0]++;
    }
  }
}

int ProperColor::color(int i, int current_color)
{
  visit_cnt[i]++;
  // fprintf(stderr, "entering BestColoring = %d, visit_cnt[%d]=%d\n",
  // BestColoring, i, visit_cnt[i]);
  if (visit_cnt[i] > max_num_visits) {
    //   fprintf(stderr, "too many visits\n");
    return 0;
  }
  prob_count++;
  if (current_color >= BestColoring)
    return (current_color);
  if (BestColoring <= lb)
    return (BestColoring);
  if (i >= num_node)
    return (current_color);
  /*  printf("Node %d, num_color %d\n",i,current_color);*/

  /* Find node with maximum color_adj */
  int max = -1;
  int place = -1;
  for (int k = 0; k < num_node; k++) {
    if (Handled[k])
      continue;
    if ((ColorCount[k] > max) ||
        ((ColorCount[k] == max) && (ColorAdj[k][0] > ColorAdj[place][0]))) {
      // fprintf(stderr, "Best now at %d\n",k);
      max = ColorCount[k];
      place = k;
    }
  }
  // Adrian: Disconnected graphs haven't triggered this. Original code
  // exited. I added 'return BestColoring' as a guess of what to do if
  // it does get triggered. Maybe some cleanup is also required to
  // avoid a segfault
  if (place == -1) {
    fprintf(stderr, "Graph is disconnected.  This code needs to be updated for "
                    "that case.\n");
    return BestColoring;
  }

  Order[i] = place;
  Handled[place] = true;

  // fprintf(stderr, "Using node %d at level %d\n",place,i);
  for (int j = 1; j <= current_color; j++) {
    if (!ColorAdj[place][j]) {
      ColorClass[place] = j;
      assign_color(place, j);
      int new_val = color(i + 1, current_color);
      if (new_val < BestColoring) {
        BestColoring = new_val;
        BestColorClass = ColorClass;
      }
      remove_color(place, j);
      if (BestColoring < current_color) {
        Handled[place] = false;
        return (BestColoring);
      }
    }
  }

  if (current_color + 1 < BestColoring) {
    ColorClass[place] = current_color + 1;
    assign_color(place, current_color + 1);
    int new_val = color(i + 1, current_color + 1);
    if (new_val < BestColoring) {
      BestColoring = new_val;
      BestColorClass = ColorClass;
    }

    remove_color(place, current_color + 1);
  }
  Handled[place] = false;
  // fprintf(stderr, "BestColoring = %d\n", BestColoring);
  return (BestColoring);
}

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
   Name: prop_color.cc
   Description: minimal proper colouring
   Project: Antiprism - http://www.antiprism.com
   Changes: 16-04-04 Adrian Rossiter <adrian_r@terra.es>
      removed use of iostreams, graph initialised from a geometry object
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <algorithm>
using std::max;
using std::min;
//#define max(x,y) (((x)>(y))?(x):(y))
//#define min(x,y) (((x)<(y))?(x):(y))

Graph::Graph(ProperColor &pr_col) : prop(pr_col)
{
  VERTICES = prop.num_node;
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
    arrayE[i].first = nullptr;
    arrayE[i].edges = 0;
  }

  // Edge creation
  set<pair<int, int>>::const_iterator si;
  for (si = prop.adj.begin(); si != prop.adj.end(); ++si) {

    int i = si->first;
    int j = si->second;
    if (i < j) {
      auto *help = new ListEdge;
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
  }

  // Creating dynamic arrays of edges
  for (i = 0; i < VERTICES; i++) {
    vertex[i].COLOR = -1;
    vertex[i].EDGES = arrayE[i].edges;
    if (arrayE[i].edges)
      vertex[i].edge = new Edge[arrayE[i].edges];
    long count = 0;
    while (arrayE[i].first) {
      vertex[i].edge[count++].vert = arrayE[i].first->pointer;
      ListEdge *help = arrayE[i].first;
      arrayE[i].first = arrayE[i].first->next;
      delete help;
    }
  }
  delete[] arrayE;
  RANDOM_1 = ((rand() >> 1) << 1) + 1;
  RANDOM_2 = ((rand() >> 1) << 1) + 1;
  LENGTH_RANDOM = 10 * VERTICES;
  rando = new float[LENGTH_RANDOM];
  for (i = 0; i < LENGTH_RANDOM; i++)
    rando[i] = Uniform();
  CURRENT_RANDOM = 0;
}

/************************* RANDOM *************************/
float Graph::Uniform()
{
  long Z, k;
  k = RANDOM_1 / 53668;
  RANDOM_1 = 40014 * (RANDOM_1 - k * 53668) - k * 12211;
  if (RANDOM_1 < 0)
    RANDOM_1 = RANDOM_1 + 2147483563;
  k = RANDOM_2 / 52774;
  RANDOM_2 = 40692 * (RANDOM_2 - k * 52774) - k * 3791;
  if (RANDOM_2 < 0)
    RANDOM_2 = RANDOM_2 + 2147483399;
  Z = RANDOM_1 - RANDOM_2;
  if (Z < 1)
    Z = Z + 2147483562;
  return (float(Z) * (float)4.656613E-10);
}
/************************* END RANDOM *************************/

/************************* OBJECTIVE *************************/
void Graph::SetValue(Vertex *vert, float *square, float *cube)
{
  vert->OBJECTIVE_VALUE = 0.;
  for (long j = 0; j < vert->EDGES; j++)
    if (vert->edge[j].vert->STATUS == _FREE)
      vert->OBJECTIVE_VALUE += cube[vert->edge[j].vert->NEIGHBORS];
  vert->OBJECTIVE_VALUE *= square[vert->NEIGHBORS];
}

float Graph::Cost(Edge *MIS, long TOP_MIS)
{
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
void Graph::GraphColoring(long *parameter, long &current_color)
{
  current_color = 0;
  long number_of_colored = 0;
  //	float best_cost = 0.;
  Window FRAME(VERTICES, parameter[3]);
  auto *GRAPH = new Edge[VERTICES];
  long TOP_GRAPH = -1;
  auto *MIS = new Edge[VERTICES];
  long TOP_MIS = -1;
  auto *FREES = new Edge[VERTICES];
  long TOP_FREES = -1;
  auto *DEL = new Edge[VERTICES];
  long TOP_DEL = -1;
  auto *SINGLES = new float[VERTICES + 1];
  auto *SQUARES = new float[VERTICES + 1];
  auto *CUBES = new float[VERTICES + 1];
  long i;
  for (i = 0; i <= VERTICES; i++) {
    SINGLES[i] = float(i);
    SQUARES[i] = SINGLES[i] * SINGLES[i];
    CUBES[i] = SQUARES[i] * SINGLES[i];
  }
  Logic FIRST_TRY = _YES;

  // unneeded AR: long EdgesLeft;

  /********** ADD NEW COLOR **********/

  do {

    /********** INITIALIZATION **********/

    // unneeded AR: EdgesLeft = 0;
    TOP_GRAPH = -1;
    long i;
    for (i = 0; i < VERTICES; i++)
      if (vertex[i].STATUS != _COLORED) {
        TOP_GRAPH++;
        GRAPH[TOP_GRAPH].vert = &vertex[i];
        vertex[i].STATUS = _FREE;
        vertex[i].SELECTED_NEIGHBORS = 0;
        vertex[i].NEIGHBORS = vertex[i].EDGES - vertex[i].COLORED_NEIGHBORS;
        // unneeded AR: EdgesLeft += vertex[i].NEIGHBORS;
      }
    // unneeded AR: EdgesLeft = EdgesLeft >> 1;
    for (i = 0; i <= TOP_GRAPH; i++)
      SetValue(GRAPH[i].vert, SQUARES, CUBES);
    FRAME.ComputeCosts(this);

    /**********************************/
    // if (current_color > 0) {
    //   cout << "\tC=" << current_color << " Clrd=" << number_of_colored << "
    //   \n";
    //   cout.flush();
    //}
    /**********************************/

    /********** INITIAL INDEPENDENT SET **********/

    TOP_MIS = -1;
    float best_cost = 0.;
    for (i = 0; i <= TOP_GRAPH; i++)
      if (GRAPH[i].vert->STATUS == _FREE) {
        TOP_MIS++;
        MIS[TOP_MIS].vert = GRAPH[i].vert;
        GRAPH[i].vert->SelectVertex();
        best_cost += GRAPH[i].vert->OBJECTIVE_VALUE;
      }
    FRAME.InsertMis(MIS, TOP_MIS, best_cost);

    /**********************************/
    // cout << "\t\tInitial MIS=" << TOP_MIS + 1
    //   << " COST=" << best_cost << "\r";
    // cout.flush();
    /**********************************/

    if (VERTICES - number_of_colored - TOP_MIS - 1 > 0) {
      long limit = long(
          parameter[0] * (1. - SQUARES[number_of_colored] / SQUARES[VERTICES]) +
          parameter[1]);
      if (FIRST_TRY == _YES) {
        limit *= parameter[3];
        FIRST_TRY = _NO;
      }

      /********** RANDOMIZED INDEPENDENT SETS **********/
      long counter = 0;
      long no_improvement = 0;
      float current_cost;
      if (limit > 0)
        do {

          /********** RANDOMIZED VERTEX EXCLUSION **********/
          TOP_DEL = TOP_FREES = -1;
          counter++;
          do {
            long index = long((TOP_MIS + 1) * Get());
            MIS[index].vert->STATUS = _FREE;
            MIS[index].vert->SELECTED_NEIGHBORS = 0;
            TOP_DEL++;
            DEL[TOP_DEL].vert = MIS[index].vert;
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
                if (temp > maximum) {
                  maximum = temp;
                  index = i;
                }
              }
              else {
                for (long j = i; j < TOP_FREES; j++)
                  FREES[j].vert = FREES[j + 1].vert;
                TOP_FREES--;
              }
            if (index > -1) {
              FREES[index].vert->SelectVertex();
              TOP_MIS++;
              MIS[TOP_MIS].vert = FREES[index].vert;
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
                  TOP_FREES++;
                  FREES[TOP_FREES].vert = DEL[i].vert;
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
            // cout << "\t\tcounter=" << counter << " MIS=" << TOP_MIS + 1
            //   << " COST=" << current_cost << "\r";
            // cout.flush();
            /**********************************/
          }
          else
            no_improvement++;
        } while (no_improvement < limit);

      /********** LOCAL SEARCH **********/
      long i;
      for (i = 0; i < FRAME.ELEMENTS; i++)
        if (FRAME.list[i].FULL == _YES)
          FRAME.list[i].COST =
              LocalSearch(&(FRAME.list[i].set[0]), FRAME.list[i].SIZE,
                          parameter, GRAPH, TOP_GRAPH);
    }

    /********** COLOR THE BEST INDEPENDENT SET **********/
    FRAME.BestMis(MIS, TOP_MIS, best_cost);
    for (i = 0; i <= TOP_MIS; i++)
      MIS[i].vert->ColorVertex(current_color);
    number_of_colored += TOP_MIS + 1;
    current_color++;

  } while (number_of_colored < VERTICES);
  // cout << "\tC=" << current_color << " Clrd=" << number_of_colored << " \n";
  // cout.flush();
  // cout << "ColorsUsed=" << current_color << " \n";
  // cout.flush();
  for (i = 0; i < VERTICES; i++) {
    if (vertex[i].EDGES != vertex[i].COLORED_NEIGHBORS)
      fprintf(stderr, "\noff_color: warning: mistake in proper colouring!\n");
    // printf("%ld = %ld\n", i, vertex[i].COLOR);
    prop.set_color(i, vertex[i].COLOR);

    /*      for (long j = 0; j < vertex[i].EDGES; j++)
             if (!(vertex[i].COLOR < current_color && vertex[i].COLOR >= 0) ||
                   vertex[i].COLOR == vertex[i].edge[j].vert->COLOR) {
                fprintf(stderr, "\noff_color: error: while calculating proper
       Coloring\n");
                exit (-1);
             }
    */
  }
  delete[] MIS;
  delete[] FREES;
  delete[] DEL;
  delete[] SQUARES;
  delete[] SINGLES;
  delete[] CUBES;
  delete[] GRAPH;
}

float Graph::LocalSearch(Edge *BEST_MIS, long &BEST_TOP_MIS, long *parameter,
                         Edge *GRAPH, long TOP_GRAPH)
{
  auto *MIS = new Edge[VERTICES];
  long TOP_MIS = -1;
  auto *FREES = new Edge[VERTICES];
  long TOP_FREES = -1;
  auto *DEL = new Edge[VERTICES];
  long TOP_DEL = -1;
  long counter = 0, iterations = 0;
  float best_cost = 0., current_cost;
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
      if (done == _YES) {
        done = _NO;
        best_cost = current_cost;
      }

      /************** RANDOMIZED VERTEX EXCLUSION **************/
      TOP_DEL = TOP_FREES = -1;
      do {
        long index = long((TOP_MIS + 1) * Get());
        MIS[index].vert->STATUS = _FREE;
        MIS[index].vert->SELECTED_NEIGHBORS = 0;
        TOP_DEL++;
        DEL[TOP_DEL].vert = MIS[index].vert;
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

      /************** RANDOMIZED VERTEX INCLUSION **************/
      Logic first = _YES;
      do {
        float maximum = 0.;
        long index = -1;
        long i;
        for (i = 0; i <= TOP_FREES; i++)
          if (FREES[i].vert->STATUS == _FREE) {
            float temp = Get();
            if (temp > maximum) {
              maximum = temp;
              index = i;
            }
          }
          else {
            for (long j = i; j < TOP_FREES; j++)
              FREES[j].vert = FREES[j + 1].vert;
            TOP_FREES--;
          }
        if (index > -1) {
          FREES[index].vert->SelectVertex();
          TOP_MIS++;
          MIS[TOP_MIS].vert = FREES[index].vert;
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
              TOP_FREES++;
              FREES[TOP_FREES].vert = DEL[i].vert;
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
        // cout << "\t\tlocal_counter=" << iterations << " MIS=" << TOP_MIS + 1
        //   << " COST=" << current_cost << "\r";
        // cout.flush();
        /**********************************/
      }
    } while (counter < parameter[2]);
  } while (improvement == _YES);
  delete[] MIS;
  delete[] FREES;
  delete[] DEL;
  return (best_cost);
}

/****************/
/***  Vertex  ***/
/****************/
void Vertex::SelectVertex()
{
  long i;
  for (i = 0; i < EDGES; i++) {
    if (edge[i].vert->STATUS == _FREE)
      edge[i].vert->STATUS = _NEIGHBOR;
    if (edge[i].vert->STATUS == _NEIGHBOR)
      edge[i].vert->SELECTED_NEIGHBORS++;
  }
  STATUS = _SELECTED;
}

void Vertex::ColorVertex(long color)
{
  COLOR = color;
  STATUS = _COLORED;
  long i;
  for (i = 0; i < EDGES; i++)
    edge[i].vert->COLORED_NEIGHBORS++;
}

Logic Vertex::ConnectedToVertex(Vertex *guy)
{
  long i;
  for (i = 0; i < EDGES; i++)
    if (edge[i].vert == guy)
      return (_YES);
  return (_NO);
}

/****************/
/***  Window  ***/
/****************/
Window::Window(long size, long elem)
{
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

Window::~Window()
{
  long i;
  for (i = 0; i < ELEMENTS; i++)
    delete[] list[i].set;
  delete[] list;
}

void Window::InsertMis(Edge *MIS, long TOP_MIS, float cost)
{
  long i;
  for (i = 0; i < ELEMENTS; i++)
    if (list[i].FULL == _YES)
      if ((list[i].COST - cost) * (list[i].COST - cost) < 1.)
        return;
  long smallest_index = -1; // originally uninitialised - Adrian
  for (i = 0; i < ELEMENTS; i++)
    if (list[i].FULL == _NO) {
      smallest_index = i;
      break;
    }
    else if (list[i].COST == SMALLEST_COST) {
      smallest_index = i;
      break;
    }
  SMALLEST_COST = list[smallest_index].COST = cost;
  list[smallest_index].FULL = _YES;
  list[smallest_index].SIZE = TOP_MIS;
  for (i = 0; i <= TOP_MIS; i++)
    list[smallest_index].set[i].vert = MIS[i].vert;
  for (i = 0; i < ELEMENTS; i++)
    if (list[i].FULL == _YES && list[i].COST < SMALLEST_COST)
      SMALLEST_COST = list[i].COST;
    else if (list[i].FULL == _NO) {
      SMALLEST_COST = -1.;
      return;
    }
}

void Window::BestMis(Edge *MIS, long &TOP_MIS, float &cost)
{
  float maximum = -1.;
  long index = -1;
  long i;
  for (i = 0; i < ELEMENTS; i++)
    if (list[i].FULL == _YES) {
      if (list[i].COST > maximum) {
        index = i;
        maximum = list[i].COST;
      }
      else if (maximum - list[i].COST < 0.005 * maximum &&
               list[i].SIZE > list[index].SIZE) {
        index = i;
        maximum = list[i].COST;
      }
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
          if (list[i].set[k].vert == MIS[j].vert) {
            to_delete = _YES;
            break;
          }
        if (to_delete == _YES)
          break;
      }
      if (to_delete == _YES)
        list[i].FULL = _NO;
    }
  SMALLEST_COST = -1.;
}

void Window::ComputeCosts(Graph *pointer)
{
  Logic check = _NO;
  long i;
  for (i = 0; i < ELEMENTS; i++)
    if (list[i].FULL == _YES)
      list[i].COST = pointer->Cost(&(list[i].set[0]), list[i].SIZE);
    else
      check = _YES;
  if (check == _NO) {
    fprintf(stderr, "\nError_1!\n");
    exit(-1);
  }
  SMALLEST_COST = -1.;
}
