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

#ifndef PRIVATE_PROP_COL_H
#define PRIVATE_PROP_COL_H

/*
   Name: prop_col.h
   Description: minimal proper colouring
   Project: Antiprism - http://www.antiprism.com
   Changes: 12/09/11 Adrian Rossiter <adrian@antiprism.com>
      convert to class and use STL containers rather than arrays
*/

// Adrian Rossiter: converted code in class

#include <map>
#include <set>
#include <utility>
#include <vector>

class ProperColor {
private:
  std::set<std::pair<int, int>> adj;
  int BestColoring;
  std::vector<int> ColorClass;
  std::vector<int> BestColorClass;
  int prob_count;
  std::vector<int> Order;
  std::vector<bool> Handled;
  std::vector<std::vector<int>> ColorAdj;
  std::vector<int> ColorCount;
  std::vector<int> visit_cnt;
  const int max_num_visits; // chosen as indicator of non-completion

  int lb;
  int num_prob, max_prob;

  int best_clique;

  int num_node;

  int greedy_clique(int *valid, int *clique);
  int max_w_clique(int *valid, int *clique, int lower, int target);
  void assign_color(int node, int color);
  void remove_color(int node, int color);
  int color(int i, int current_color);
  bool is_adj(int i, int j)
  {
    if (i > j)
      std::swap(i, j);
    return adj.find(std::make_pair(i, j)) != adj.end();
  }
  void set_color(int i, int col) { BestColorClass[i] = col + 1; }

public:
  ProperColor(int nodes, int max_visits = 100)
      : max_num_visits(max_visits), num_node(nodes){};
  void set_adj(int i, int j)
  {
    if (i > j)
      std::swap(i, j);
    adj.insert(std::make_pair(i, j));
  }
  int find_colors();
  int get_color(int i) { return BestColorClass[i] - 1; }

  friend class Graph;
};

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

enum Logic { _NO, _YES };
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
  Vertex()
  {
    COLOR = -1;
    EDGES = 0;
    edge = nullptr;
    STATUS = _FREE;
    SELECTED_NEIGHBORS = 0;
    COLORED_NEIGHBORS = 0;
    OBJECTIVE_VALUE = 0.;
    NEIGHBORS = 0;
  }
  void SelectVertex();
  void ColorVertex(long);
  Logic ConnectedToVertex(Vertex *);
};

class Graph {
  Vertex *vertex;
  long VERTICES;
  long RANDOM_1, RANDOM_2;
  float *rando;
  long LENGTH_RANDOM;
  long CURRENT_RANDOM;
  ProperColor &prop;

public:
  Graph(ProperColor &pr_col);
  ~Graph()
  {
    for (long i = 0; i < VERTICES; i++)
      if (vertex[i].edge)
        delete vertex[i].edge;
    if (vertex)
      delete vertex;
    delete rando;
  }
  float Cost(Edge *, long);
  void SetValue(Vertex *, float *, float *);
  float LocalSearch(Edge *, long &, long *, Edge *, long);
  void GraphColoring(long *, long &);
  float Uniform();
  float Get()
  {
    CURRENT_RANDOM = (CURRENT_RANDOM + 1) % LENGTH_RANDOM;
    return (rando[CURRENT_RANDOM]);
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
  Window(long, long);
  ~Window();
  void InsertMis(Edge *, long, float);
  void BestMis(Edge *, long &, float &);
  void ComputeCosts(Graph *);
};

#endif // PRIVATE_PROP_COL_H
