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
   Name: minmax.cc
   Description: minimise the maximum edgelength on a sphere
   Project: Antiprism - http://www.antiprism.com
*/

#include <string.h>
#include <math.h>
#include <ctype.h>
#include <stdlib.h>
#include <string>
#include <vector>

#include "../base/antiprism.h"

using std::string;
using std::vector;



class mm_opts: public prog_opts
{
   public:
      int num_iters;
      char algm;
      char placement;
      double shorten_by;
      double lengthen_by;
      vec4d ellipsoid;
      
      string ifile;
      string ofile;

      mm_opts(): prog_opts("minmax"), num_iters(1000), algm('v'),
                 placement('n'), shorten_by(1.0), lengthen_by(0.0)
                 {}

      void process_command_line(int argc, char **argv);
      void usage();
};


void mm_opts::usage()
{
   fprintf(stdout,
"\n"
"Usage: %s [options] [input_file]\n"
"\n"
"Read a file in OFF format containing a graph of a polyhedron, with or\n"
"without vertex coordinates, and try to create a spherical or ellipsoidal\n"
"tesselation where the maximum edge is a minimum length. If input_file is\n"
"not given the program reads from standard input.\n"
"\n"
"Options\n"
"%s"
"  -n <itrs> number of iterations (default 1000)\n" 
"  -s <perc> percentage to shorten the maximum edge by (default 1)\n" 
"  -l <perc> percentage to lengthen the minimum edge by (default 1)\n" 
"  -a <alg>  length changing algorithm\n"
"              v - shortest and longest edges attached to a vertex (default)\n"
"              a - shortest and longest of all edges\n"
"  -p <mthd> method of placement onto a unit sphere:\n"
"              n - project onto the sphere (default)\n"
"              r - random placement\n"
"              u - unscramble: place a small polygon on one side and the\n"
"                  rest of the vertices at a point on the other side\n"
"  -E <prms> use ellipsoid, three numbers separated by commas are the\n"
"            axis lengths (for a superellipsoid an optional fourth number\n"
"            gives the power)\n"
"  -o <file> write output to file (default: write to standard output)\n"
"\n"
"\n", prog_name(), help_ver_text);
}


     
void mm_opts::process_command_line(int argc, char **argv)
{
   char errmsg[MSG_SZ];
   opterr = 0;
   char c;
   vector<double> nums;
   
   handle_long_opts(argc, argv);

   while( (c = getopt(argc, argv, ":hn:s:l:a:p:E:o:")) != -1 ) {
      if(common_opts(c, optopt))
         continue;

      switch(c) {
         case 'o':
            ofile = optarg;
            break;

         case 'n':
            if(!read_int(optarg, &num_iters, errmsg))
               error(errmsg, c);
            if(num_iters < 0)
               error("number of iterations must be greater than 0", c);
            break;

         case 's':
            if(!read_double(optarg, &shorten_by, errmsg))
               error(errmsg, c);
            if(shorten_by <= 0 || shorten_by >=100)
               warning("not inside range 0 to 100", c);
            break;

         case 'l':
            if(!read_double(optarg, &lengthen_by, errmsg))
               error(errmsg, c);
            if(shorten_by <= 0 || lengthen_by >=100) {
               warning("not inside range 0 to 100", c);
            }
            break;

         case 'a':
            if(strlen(optarg) > 1 || !strchr("avu", *optarg))
               error("method is '"+string(optarg)+"' must be a or v");
            algm = *optarg;
            break;

         case 'p':
            if(strlen(optarg) > 1 || !strchr("nur", *optarg))
               error("method is '"+string(optarg)+"' must be n, u, or r");
            placement = *optarg;
            break;

         case 'E': 
            if(!read_double_list(optarg, nums, errmsg))
               error(errmsg, c);
            if(nums.size()<3 || nums.size()>4)
               error("must give three or four numbers only", c);
            ellipsoid = vec4d(nums[0], nums[1], nums[2],
                                              (nums.size()==4) ? nums[3]: 2);
            if(ellipsoid[3]<=0)
               error("superellipsoid power must be greater than zero", c);
            break;

         default:
            error("unknown command line error");
      }
   }

   if(argc-optind > 1)
      error("too many arguments");
   
   if(argc-optind == 1)
      ifile=argv[optind];

}


void to_ellipsoid(vec3d &v, vec4d ellipsoid)
{
   if(!ellipsoid.is_set())
      v.to_unit();
   else if (v.mag2()>epsilon) {
      double scale=0;
      for(unsigned int i=0; i<3; i++)
         scale += pow(fabs(v[i]/ellipsoid[i]), ellipsoid[3]);
      v *= pow(scale, -1/ellipsoid[3]);
   }
}

void initial_placement(geom_if &geom, char placement, vec4d ellipsoid)
{
   vector<int> face = (*geom.get_faces())[0];
   switch(placement) {
      case 'n':
         for(unsigned int i=0; i<geom.get_verts()->size(); i++)
            to_ellipsoid((*geom.get_verts())[i], ellipsoid);
         break;
      
      case 'u':
         for(unsigned int i=0; i<geom.get_verts()->size(); i++)
            (*geom.get_verts())[i] = vec3d(-1,0,0);
        
         for(unsigned int i=0; i<face.size(); i++) {
            (*geom.get_verts())[face[i]] = vec3d(1,
                    0.01*cos(i*2*M_PI/face.size()),
                    0.01*sin(i*2*M_PI/face.size()));
            (*geom.get_verts())[face[i]].to_unit();
         }
         break;
         
      case 'r':
      {
         rand_gen rnd;
         rnd.time_seed();
         for(unsigned int i=0; i<geom.get_verts()->size(); i++) {
            (*geom.get_verts())[i] = vec3d::random(rnd);
            to_ellipsoid((*geom.get_verts())[i], ellipsoid);
         }
         break;
      }
   }
}

void minmax_a(geom_v &geom, double shorten_factor, double lengthen_factor,
      int n=1000, vec4d ellipsoid=vec4d())
{
   vector<vec3d> &verts = *geom.get_verts();
   vector<vector<int> > &edges = *geom.get_edges();
   int max_edge=0, min_edge=0, p0, p1;
   double dist, max_dist=0, min_dist=1e100;
   for(int cnt=0; cnt<n; cnt++) {
      max_dist = 0;
      min_dist = 1e100;
      for(unsigned int i=0; i<edges.size(); i++) {
         dist = (verts[edges[i][1]] - verts[edges[i][0]]).mag2();
         if(dist > max_dist) {
            max_dist = dist;
            max_edge = i;
         }
         if(dist < min_dist) {
            min_dist = dist;
            min_edge = i;
         }
      }
      p0 = edges[max_edge][0];
      p1 = edges[max_edge][1];
      vec3d diff = verts[p1] - verts[p0];
      verts[p0] += diff*shorten_factor;
      to_ellipsoid(verts[p0], ellipsoid);
      verts[p1] -= diff*shorten_factor;
      to_ellipsoid(verts[p1], ellipsoid);

      p0 = edges[min_edge][0];
      p1 = edges[min_edge][1];
      diff = verts[p1] - verts[p0];
      verts[p0] -= diff*lengthen_factor;
      to_ellipsoid(verts[p0], ellipsoid);
      verts[p1] += diff*lengthen_factor;
      to_ellipsoid(verts[p1], ellipsoid);
      
      if((cnt+1)%100 == 0)
         fprintf(stderr, ".");
      if((cnt+1)%1000 == 0)
         fprintf(stderr, "\n%-15d %12.10g  (%12.10g) ", cnt+1, max_dist, min_dist);
   }
   if(n%1000 != 0)
      fprintf(stderr, "\n%-15d %12.10g  (%12.10g) ", n, max_dist, min_dist);
   fprintf(stderr, "\n");
}

/*
bool get_dist_sum(vector<vec3d> &verts, vec3d vert, vector<int> &eds, double &dist)
{
   bool valid = true;
   dist = 0;
   for(unsigned int i=0; i<eds.size(); i++) {
      double inc = (vert - verts[eds[i]]).mag();
      if(inc<1)
         valid = false;
      dist += inc;
   }

   return valid;
}
      
void minmax_v(geom_v &geom, vector<vector<int> > &eds, double shorten_factor, double, int n=1000)
{
   vector<vec3d> &verts = *geom.get_verts();
   double perim=24;
   rand_gen rnd;
   rnd.time_seed();
   for(int cnt=0; cnt<n; cnt++) {
      for(unsigned int v=0; v<verts.size(); v++) {
         double cur_sum, test_sum;
         get_dist_sum(verts, verts[v], eds[v], cur_sum);
         for(int j=0; j<20;j++) {
            vec3d test_vert =(verts[v] + vec3d::random(rnd)*shorten_factor).unit();
            if( get_dist_sum(verts, test_vert, eds[v], test_sum) &&
                test_sum<cur_sum ) {
               cur_sum = test_sum;
               verts[v] = test_vert;
               break;
            }
         }
      }
     
      perim = 0.0;
      for(unsigned int v=0; v<verts.size(); v++) {
         double sum;
         get_dist_sum(verts, verts[v], eds[v], sum);
         perim += sum;
      }
      if((cnt+1)%100 == 0)
         fprintf(stderr, ".");
      if((cnt+1)%1000 == 0)
         fprintf(stderr, "\n%-15d  perim=%12.10g", cnt+1, perim);
   }
   if(n%1000 != 0)
      fprintf(stderr, "\n%-15d  perim=%12.10g", n, perim);
   fprintf(stderr, "\n");
}
*/
 
void minmax_v(geom_v &geom, vector<vector<int> > &eds, double shorten_factor,
      double lengthen_factor, int n=1000, vec4d ellipsoid=vec4d())
{
   vector<vec3d> &verts = *geom.get_verts();
   int max_edge=0, min_edge=0, p0, p1;
   double dist, max_dist=0, min_dist=1e100, g_max_dist=0, g_min_dist=1e100;
   for(int cnt=0; cnt<n; cnt++) {
      g_max_dist = 0;
      g_min_dist = 1e100;
      for(unsigned int v=0; v<verts.size(); v++) {
         if(eds[v].size()==0)
            continue;
         max_dist = 0;
         min_dist = 1e100;
         for(unsigned int i=0; i<eds[v].size(); i++) {
            dist = (verts[v] - verts[eds[v][i]]).mag2();
            //fprintf(stderr, "dist=%g\n", dist);
            if(dist > max_dist) {
               max_dist = dist;
               max_edge = i;
            }
            if(dist > g_max_dist)
               g_max_dist = dist;
            if(dist < min_dist) {
               min_dist = dist;
               min_edge = i;
            }
            if(dist < g_min_dist)
               g_min_dist = dist;
         }
         
         p0 = v;
         p1 = eds[v][max_edge];
         vec3d diff = verts[p1] - verts[p0];
         verts[p0] += diff*shorten_factor;
         to_ellipsoid(verts[p0], ellipsoid);

         p0 = v;
         p1 = eds[v][min_edge];
         diff = verts[p1] - verts[p0];
         verts[p0] -= diff*lengthen_factor;
         to_ellipsoid(verts[p0], ellipsoid);
      }
      
      if((cnt+1)%100 == 0)
         fprintf(stderr, ".");
      if((cnt+1)%1000 == 0)
         fprintf(stderr, "\n%-15d %12.10g  (%12.10g) ", cnt+1, g_max_dist, g_min_dist);
   }
   if(n%1000 != 0)
      fprintf(stderr, "\n%-15d %12.10g  (%12.10g) ", n, g_max_dist, g_min_dist);
   fprintf(stderr, "\n");
}



void minmax_unit(geom_v &geom, vector<vector<int> > &eds, double shorten_factor, int n=1000)
{
   vector<vec3d> &verts = *geom.get_verts();
   vector<vector<int> > diags (verts.size(), vector<int>());
   vector<vector<int> > &faces = *geom.get_faces();
   for(unsigned int i=0; i<faces.size(); i++) {
      int fsz = faces[i].size();
      for(int j=0; j<fsz/2; j++)  {
         if(fsz==3)
            continue;
         //fprintf(stderr, "adding %d - %d\n", faces[i][j], faces[i][j+fsz/2]);
         diags[faces[i][j]].push_back(faces[i][j+fsz/2]);
         diags[faces[i][j+fsz/2]].push_back(faces[i][j]);
      }
   }
         
   double dist, max_dist=0, min_dist=1e100, g_max_dist=0, g_min_dist=1e100;
   for(int cnt=0; cnt<n; cnt++) {
      g_max_dist = 0;
      g_min_dist = 1e100;
      for(unsigned int v=0; v<verts.size(); v++) {
         max_dist = 0;
         min_dist = 1e100;
         vec3d off_sum(0,0,0);
         for(unsigned int i=0; i<eds[v].size(); i++) {
            vec3d off = verts[v] - verts[eds[v][i]];
            dist = off.mag();
            off_sum += off.unit()*(1-dist)*shorten_factor;
            
            //fprintf(stderr, "dist=%g\n", dist);
            if(dist > max_dist) {
               max_dist = dist;
            }
            if(dist > g_max_dist)
               g_max_dist = dist;
            if(dist < min_dist) {
               min_dist = dist;
            }
            if(dist < g_min_dist)
               g_min_dist = dist;
         }
         
         verts[v] += off_sum;
         off_sum = vec3d(0,0,0);
     
         vector<vec3d> vs(verts.size(), vec3d(0,0,0));
         for(unsigned int f=0; f<faces.size(); f++) {
            if(faces[f].size()==3)
               continue;
            vec3d norm = face_norm(verts, faces[f]).unit();
            vec3d f_cent = centroid(verts, faces[f]);
            if(vdot(norm, f_cent)<0)
               norm *= -1.0;
            for(unsigned int v=0; v<faces[f].size(); v++)
               vs[faces[f][v]] +=
                  vdot(shorten_factor*norm, f_cent - verts[faces[f][v]]) *norm;
         }

         // adjust vertices post-loop
         for(unsigned int i=0; i<vs.size(); i++)
            verts[i] += vs[i];

         for(unsigned int i=0; i<diags[v].size(); i++) {
            vec3d off = verts[v] - verts[diags[v][i]];
            dist = off.mag();
            off_sum += off.unit()*(sqrt(2)-dist)*shorten_factor;

            if(dist > max_dist) {
               max_dist = dist;
            }
            if(dist > g_max_dist)
               g_max_dist = dist;
            if(dist < min_dist) {
               min_dist = dist;
            }
            if(dist < g_min_dist)
               g_min_dist = dist;
         }
      verts[v] += off_sum;
      }

      
      if((cnt+1)%100 == 0)
         fprintf(stderr, ".");
      if((cnt+1)%1000 == 0)
         fprintf(stderr, "\n%-15d %12.10g  (%12.10g) ", cnt+1, g_max_dist, g_min_dist);
   }
   if(n%1000 != 0)
      fprintf(stderr, "\n%-15d %12.10g  (%12.10g) ", n, g_max_dist, g_min_dist);
   fprintf(stderr, "\n");
}

 

int main(int argc, char *argv[])
{
   mm_opts opts;
   opts.process_command_line(argc, argv);

   char errmsg[MSG_SZ] = "";
   col_geom_v geom;
   if(!geom.read(opts.ifile, errmsg))
   if(!geom)
      opts.error(errmsg);
   if(*errmsg)
      opts.warning(errmsg);
      
   if(!geom.get_edges()->size())
      geom.add_missing_impl_edges();
   
   if(geom.get_edges()->size()) {
      if(opts.algm!='u')
         initial_placement(geom, opts.placement, opts.ellipsoid);
      if(opts.algm=='a')
         minmax_a(geom, opts.shorten_by/200, opts.lengthen_by/200,
               opts.num_iters, opts.ellipsoid);
      else{
         vector<vector<int> > &edges = *geom.get_edges();
         vector<vector<int> > eds(geom.get_verts()->size());
         for(unsigned int i=0; i<edges.size(); i++) {
            eds[edges[i][0]].push_back(edges[i][1]);
            eds[edges[i][1]].push_back(edges[i][0]);
         }
         if(opts.algm=='v')
            minmax_v(geom, eds, opts.shorten_by/200, opts.lengthen_by/200,
               opts.num_iters, opts.ellipsoid);
         else if(opts.algm=='u')
            minmax_unit(geom, eds, opts.shorten_by/200, opts.num_iters);
      }
   }
   else
      opts.warning("input file contains no edges");

   if(!geom.write(opts.ofile, errmsg))
      opts.error(errmsg);

   return 0;
}
   

