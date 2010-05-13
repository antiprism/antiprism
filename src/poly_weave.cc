/*
   Copyright (c) 2003-2009, Adrian Rossiter

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
   Name: weave.cc
   Description: make a weave based on edge mid-points
   Project: Antiprism - http://www.antiprism.com
*/

#include <string.h>
#include <math.h>
#include <unistd.h>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <functional>
#include "../base/antiprism.h"

using std::string;
using std::vector;
using std::swap;
using std::not_equal_to;

class wv_opts: public prog_opts
{
   private:
   
   public:
      //vec3d centre;

      string ifile;
      string ofile;

      wv_opts() : prog_opts("poly_weave") //, centre(vec3d(0, 0, 0))
                  {}

      void process_command_line(int argc, char **argv);
      void usage();
};

void wv_opts::usage()
{
   fprintf(stdout,
"\n"
"Usage: %s [options] [input_file]\n"
"\n"
"Read a file in OFF format and make a weave based on the mid-points of the\n"
"edges. If input_file is not given the program reads from standard input.\n"
"\n"
"Options\n"
"%s"
"  -o <file> write output to file (default: write to standard output)\n"
"\n"
"\n", prog_name(), help_ver_text);
}

     
void wv_opts::process_command_line(int argc, char **argv)
{
   //char errmsg[MSG_SZ];
   opterr = 0;
   char c;
   
   handle_long_opts(argc, argv);

   while( (c = getopt(argc, argv, ":ho:")) != -1 ) {
      if(common_opts(c, optopt))
         continue;

      switch(c) {
         /*case 'c':
            if(!centre.read(optarg, errmsg))
               error(errmsg, c);
            break;
           */ 
         case 'o':
            ofile = optarg;
            break;
         
         default:
            error("unknown command line error");
      }
   }

   if(argc-optind > 1) {
      error("too many arguments");
      exit(1);
   }
   
   if(argc-optind == 1)
      ifile=argv[optind];

}

void weave(geom_if &geom, geom_if &wv_geom)
{
   vector<vec3d> &verts = *geom.get_verts();
   vector<vector<int> > &faces = *geom.get_faces();
   vector<vec3d> &wv_verts = *wv_geom.get_verts();

   map<vector<int>, int> edge_vno;
   map<vector<int>, vector<int> > edges;
   geom.get_edge_face_pairs(edges, true);
   map<vector<int>, vector<int> >::iterator mi;
   for(mi=edges.begin(); mi!=edges.end(); mi++) {
      edge_vno[mi->first] = wv_verts.size();
      wv_verts.push_back((verts[mi->first[0]] + verts[mi->first[1]])/2.0);
   }

   vector<int> w0(2);
   vector<int> w1(2);
   vector<int>::iterator vi;
   vector<int> v_prs;
   for(mi=edges.begin(); mi!=edges.end(); mi++) {
      if(mi->first.size() <2 || mi->first[0] ==-1 || mi->first[1] ==-1 ||
         mi->second.size()<2 || mi->second[0]==-1 || mi->second[1]==-1)
         continue;
      vector<int> *face = &faces[mi->second[0]];
      vi = find(face->begin(), face->end(), mi->first[1]);
      //if(vi==face->end())  // may happen if not oriented?
      //   continue;
      w1[1] = (vi+1 != face->end()) ? *(vi+1) : face->front();
      w1[0] = mi->first[1];
      if(w1[0] > w1[1])
         swap(w1[0], w1[1]);
      int v1 = edge_vno.find(w1)->second;
            
      face = &faces[mi->second[1]];
      vi = find(face->begin(), face->end(), mi->first[0]);
      //if(vi==face->end())  // may happen if not oriented?
      //   continue;
      w0[0] = (vi+1 != face->end()) ? *(vi+1) : face->front();
      w0[1] = mi->first[0];
      if(w0[0] > w0[1])
         swap(w0[0], w0[1]);
      int v0 = edge_vno.find(w0)->second;

      v_prs.push_back(v1);
      v_prs.push_back(v0);
   }

   
   while((vi = find_if(v_prs.begin(), v_prs.end(), bind2nd(not_equal_to<int>(), -1))) != v_prs.end()) {
      vector<int> face;
      face.push_back(*vi);
      face.push_back(*(vi+1));
      *vi = *(vi+1) = -1;
      while((vi = find(v_prs.begin(), v_prs.end(), face.back())) != v_prs.end()) {
         if(is_even(vi-v_prs.begin())) {
            face.push_back(*(vi+1));
            *vi = *(vi+1) = -1;
         }
         else {
            face.push_back(*(vi-1));
            *vi = *(vi-1) = -1;
         }
      }
      face.resize(face.size()-1);
      wv_geom.get_faces()->push_back(face);
   }
}





int main(int argc, char *argv[])
{
   wv_opts opts;
   opts.process_command_line(argc, argv);

   char errmsg[MSG_SZ];
   col_geom_v geom;
   if(!geom.read(opts.ifile, errmsg))
      opts.error(errmsg);
   if(*errmsg)
      opts.warning(errmsg);

   geom_info info(geom);
   if(!info.is_oriented())
      opts.warning("base polyhedron is not oriented");

   geom_v wv_geom;
   weave(geom, wv_geom);

   if(!wv_geom.write(opts.ofile, errmsg))
      opts.error(errmsg);

   return 0;
}
   

