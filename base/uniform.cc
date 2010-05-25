/*
   Copyright (c) 2007-2009, Adrian Rossiter

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
   Name: uniform.cc
   Description: Uniforn polyhedra in OFF format 
   Project: Antiprism - http://www.antiprism.com
*/



#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <unistd.h>

#include <string>
#include <vector>
#include <set>
#include <algorithm>

#include "std_polys.h"


using std::string;
using std::vector;
using std::set;
using std::swap;


int uni_poly::get_poly(geom_if &geom, int sym_no)
{
   char sym[6];
   snprintf(sym, 6, "#%d", sym_no+1);
   poly = kaleido(sym, uniform, last_uniform);
   if(!poly)
      return 0;
   
   vector<vec3d> &verts = *geom.get_verts();
   vector<vector<int> > &faces = *geom.get_faces();
   
   faces.resize(poly->F);
   vector<int> edges;
   for (int i=0; i<poly->V; i++) {
      verts.push_back(vec3d(poly->v[i].x, poly->v[i].y, poly->v[i].z));

      for (int j=0; j<poly->M; j++) {
         faces[poly->incid[j][i]].push_back(i);
         if(i<poly->adj[j][i]) {
            edges.push_back(i);
            edges.push_back(poly->adj[j][i]);
         }
      }
   }

   for(unsigned int i=0; i<faces.size(); i++) {
      set<int> vs;
      for(unsigned int j=0; j<faces[i].size(); j++)
         vs.insert(faces[i][j]);
      
      vector<int> lns;
      for(unsigned int j=0; j<edges.size(); j+=2)
         if(vs.find(edges[j])!=vs.end() && vs.find(edges[j+1])!=vs.end() ) {
            lns.push_back(edges[j]);
            lns.push_back(edges[j+1]);
         }
                  
      int num_edges = faces[i].size();
      faces[i].clear();
      faces[i].push_back(lns[0]);
      int pt = lns[1];
      for(int j=1; j<num_edges; j++) {
         faces[i].push_back(pt);
         for(int k=j; k<num_edges; k++) {
            if(lns[k*2]==pt) {
               pt=lns[k*2+1];
               swap(lns[j*2], lns[k*2+1]);
               swap(lns[j*2+1], lns[k*2]);
               break;
            }
            else if(lns[k*2+1]==pt) {
               pt=lns[k*2];
               swap(lns[j*2], lns[k*2]);
               swap(lns[j*2+1], lns[k*2+1]);
               break;
            }
         }
      }
   }
   geom.orient();
   geom.sym_align();
   return 1;
}



int uni_poly::lookup_sym_no(string sym, int is_dual)
{
   // remove double spaces and spaces at beginning and end
   string sym_norm;
   bool ignore_if_space = true;
   for(unsigned int i=0; i<sym.length(); i++) {
      if(sym[i]==' ') {
         if(ignore_if_space)
            continue;
         else
            ignore_if_space = true;
      }
      else
         ignore_if_space = false;
      sym_norm+=sym[i];
   }

   if(sym_norm[sym_norm.size()-1]==' ')
      sym_norm.resize(sym_norm.size()-1);
   

   // remove spaces either side of a punctuation mark (for Wythoff)
   string sym_norm2;
   for(unsigned int i=0; i<sym_norm.length(); i++) {
      if(sym_norm[i]==' ' &&
          ( (i>0 && ispunct(sym_norm[i-1])) ||
            (i<sym_norm.length() && ispunct(sym_norm[i+1])) ) )
         continue;
      sym_norm2+=sym_norm[i];
   }

   // sym_norm2 is now normalised
  
   // is it blank
   if(sym_norm2=="")
      return -1;

   
   // is it a number
   int offset = (sym_norm[0]=='u'||sym_norm[0]=='U');
   char *endptr;
   int idx = strtol(sym_norm2.c_str()+offset, &endptr, 10);
   
   if(!*endptr)     // all of string is an integer
      return idx-1;

   if(!is_dual) {
      // is it a Wythoff symbol
      for(int i=0; i<last_uniform; i++) {
         if(sym_norm2==uniform[i].Wythoff)
            return i;
      }
   }

   idx= -1;
 
   // is it a poly name
   for(unsigned int i=0; i<sym_norm2.size(); i++)
      if(isalpha(sym_norm2[i]))
         sym_norm2[i] = tolower(sym_norm2[i]);
   // remove any space after an -akis
   if(is_dual) {
      size_t kis_pos = sym_norm2.find("akis ");
      if(kis_pos != string::npos)
         sym_norm2.erase(kis_pos+4, 1);
   }
   //fprintf(stderr, "sym_name = '%s'\n", sym_norm2.c_str());
   for(int i=0; i<last_uniform; i++) {
      const char *name = (is_dual)?uniform[i].dual:uniform[i].name;
      if(sym_norm2==name)
         return i;
      
      if(idx<0 && strncmp(sym_norm2.c_str(), name, sym_norm2.size())==0)
         idx = i;
   }

   return idx;
}




