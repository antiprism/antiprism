/*
   Copyright (c) 2014, Roger Kaufman

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
   Name: symmetro.cc
   Description: Make symmetrohedra
   Project: Antiprism - http://www.antiprism.com
*/

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <float.h>

#include <string>
#include <vector>
#include <algorithm>

#include "../base/antiprism.h"

using std::string;
using std::vector;
using std::fill;
using std::swap;


class symmetro_opts: public prog_opts {
   public:
      int sym;
      vector<int> multipliers;
      double rotation;
      vector<int> scale_direction;
      double scale;
      int convex_hull;
      bool unitize;
      col_val axis_0_color;
      col_val axis_1_color;
      col_val axis_2_color;
      col_val convex_hull_color;
      bool verbose;
      string ofile;

      symmetro_opts(): prog_opts("symmetro"),
                       sym(0),
                       rotation(0.0),
                       scale(0.0),
                       convex_hull(4),
                       unitize(false),
                       axis_0_color(col_val(255,0,0,255)), // red
                       axis_1_color(col_val(0,0,255,255)), // blue
                       axis_2_color(col_val(0,100,0,255)), // darkgreen
                       convex_hull_color(col_val(255,255,0,255)), // yellow
                       verbose(false) {}
      
      void process_command_line(int argc, char **argv);
      void usage();
};

void symmetro_opts::usage()
{
   fprintf(stdout,
"\n"
"Usage: %s [options]\n"
"\n"
"This program is inspired by a study done by Craig S. Kaplan and George W. Hart\n"
"(http://www.georgehart.com). Project information can be found at\n"
"http://www.cgl.uwaterloo.ca/~csk/projects/symmetrohedra\n"
"\n"
"Symmetrohedra are created by placing equilateral polygons centered on the\n"
"symmetry axes of Icosahedral, Octahedral, or Tetrahedral symmetry. The number of\n"
"sides of the polygons will be a multiple number of the axis reflection number.\n"
"Axes are in order as 0, 1 and 2 corresponding to icosahedral {5,3,2},\n"
"octahedral {4,3,2}, or tetrahedral {3,3,2}. Only two multipliers may be given.\n"
"\n"
"It is possible to generate models such that the polygons intersect. If a\n"
"collision is detected, convex hull will be suppressed\n" 
"\n"
"Options\n"
"%s"
"  -s <type> symmetry type of Symmetrohedra. sets {p,q}\n"
"               icosahedral {5,3}, octahedral {4,3}, or tetrahedral {3,3}\n"
"  -m <vals> multipliers for axis polygons. Given as three integers\n"
"               separated by commas. One multiplier must be 0. e.g. 2,3,0\n"
"  -r <val>  A value which is an amount of rotation given to polygons. A value\n"
"               of zero turns polygons on vertex. A value of 1 turns them on edge\n"
"               or give a face rotation type: vertex=0, edge=1  (default: 0)\n"
"  -S <s,n,m> scale s, from axis n polygon applied to axis m polygon\n"
"               if n and m are not specified, implies first 2 polygons in order\n"
"               e.g. 0.5,1,0  (default: calculated for unit edge length)\n"
"  -C <mode> convex hull. polygons=1, suppress=2, force=3, normal=4\n"
"               (default: 4)\n"
"  -u        make the average edge length 1 unit (performed before convex hull)\n"
"  -V        verbose output\n"
"  -o <file> write output to file (default: write to standard output)\n"
"\nColoring Options (run 'off_util -H color' for help on color formats)\n"
"              and for options below, key word: none - sets no color\n"
"  -x <col>  color for axis 0 polygons (default: red)\n"
"  -y <col>  color for axis 1 polygons (default: blue)\n"
"  -z <col>  color for axis 2 polygons (default: darkgreen)\n"
"  -w <col>  color for polygons resulting from convex hull (default: yellow)\n"
"\n"
"\n", prog_name(), help_ver_text);
}

void symmetro_opts::process_command_line(int argc, char **argv)
{
   opterr = 0;
   char c;
   char errmsg[MSG_SZ];
   
   string id;
   
   handle_long_opts(argc, argv);

   while( (c = getopt(argc, argv, ":hs:m:r:S:C:ux:y:z:w:Vo:")) != -1 ) {
      if(common_opts(c, optopt))
         continue;

      switch(c) {
         case 's': // symmetry
            id = get_arg_id(optarg, "icosahedral=1|octahedral=2|tetrahedral=3", argmatch_add_id_maps, errmsg);
            if(id=="")
               error(errmsg);
            sym = atoi(id.c_str());
            break;
            
         case 'm': // multiplier
            if(!read_int_list(optarg, multipliers, errmsg, true, 3))
               error(errmsg, c);
            if( (int)multipliers.size() != 3 )
               error("multipliers must be specified as three integers", c);
            if( multipliers[0]+multipliers[1]+multipliers[2] == 0 )
               error("at least one multiplier must be positive", c);
            if( multipliers[2] == 1 )
               //error("multiplier for axis 2 cannot be 1", c);
               warning("model will contain digons");
            break;
            
         case 'r': // rotate from side to point
            id = get_arg_id(optarg, "vertex=0|edge=1", argmatch_add_id_maps, errmsg);
            if ( id == "0" )
               rotation = 0.0;
            else
            if ( id == "1" )
               rotation = 1.0;
            else
            {
               //error(errmsg);
               if(!read_double(optarg, &rotation, errmsg))
                  error(errmsg, c);
            }
            break;
           
         case 'S': // scale direction and scale
         {
            vector<double> scale_direction_tmp;
            if(!read_double_list(optarg, scale_direction_tmp, errmsg, 3))
               error(errmsg, c);
               
            // pull out ratio
            scale = scale_direction_tmp[0];
            // if zero, make a minimum scale
            // a little lower than built in epsilon
            if ( scale == 0.0 )
               scale = epsilon/10.0;
            
            if ( scale_direction_tmp.size() == 2 ) {
               error("scale takes 1 or 3 arguments",c);
            }
            else
            if ( scale_direction_tmp.size() == 3 ) {           
               // place integer portions of possible doubles
               for( int i=0; i<2; i++ ) {
                  int a = (int)floor(scale_direction_tmp[i+1]);
                  //if ( scale_direction_tmp[i] - a > 0.0 )
                  //   error(msg_str("axis numbers must be specified by an integer: '%g'", scale_direction_tmp[i+1]), c);
                  scale_direction.push_back(a);
               }
                  
               scale_direction_tmp.clear();
                  
               for( int i=0; i<(int)scale_direction.size(); i++ ) {
                  if( scale_direction[i] > 2 )
                     error("ratio direction should be 0, 1 or 2", c);
                  for( int j=i+1; j<(int)scale_direction.size(); j++ )
                     if( scale_direction[i] == scale_direction[j] )
                         error(msg_str("an axis number is specified more than once: '%d'", scale_direction[j]), c);
               }
            }      
            break;
         }
            
         case 'C':
            id = get_arg_id(optarg, "polygons=1|suppress=2|force=3|normal=4", argmatch_add_id_maps, errmsg);
            if(id=="")
               error(errmsg,c);
            convex_hull = atoi(id.c_str());
            break;
            
         case 'u':
            unitize = true;
            break;
            
         case 'x':
            if(!axis_0_color.read(optarg, errmsg))
               error(errmsg, c);
            break;

         case 'y':
            if(!axis_1_color.read(optarg, errmsg))
               error(errmsg, c);
            break;

         case 'z':
            if(!axis_2_color.read(optarg, errmsg))
               error(errmsg, c);
            break;
           
         case 'w':
            if(!convex_hull_color.read(optarg, errmsg))
               error(errmsg, c);
            break;
                                    
         case 'V':
            verbose = true;
            break;

         case 'o':
            ofile = optarg;
            break;

         default:
            error("unknown command line error");
      }
   }
   
   if(argc-optind > 0)
      error("too many arguments");
}

class symmetro
{
public:
	symmetro()
	{
		fill( mult, mult + 3, 0 );
		
		for( int i=0; i<3; i++ )
		   sym_vec.push_back(vec3d(0.0,0.0,0.0));
	}
	
	void debug();
	void setSym( const int &psym, const int &qsym );
	void setMult( const int &a, const int &m );
	int getOrder( const int &a );
	int getN( const int &a );
	
	double getAngleBetweenAxes( const int &axis1, const int &axis2 );
	void fill_sym_vec();
	
   double angle( const int &n, const int &d );
   double circumradius( const int &n, const int &d );

   vector<col_geom_v> makePolygons( const symmetro_opts &opts );
	
	~symmetro()
	{}
	
private:
	int 		p;
	int 		q;

	int 		mult[3];

   vector<vec3d> sym_vec;
};

void symmetro::debug()
{
   fprintf(stderr,"\nsymmetry = {%d,%d,2}\n\n", p, q);
   
   for( int i=0; i<3; i++ )
      fprintf(stderr,"axis %d: mult = %d\n", i, mult[i]);
   fprintf(stderr,"\n");
   
   for( int i=0; i<3; i++ )
      if ( mult[i] )
         fprintf(stderr,"axis %d polygon: %d-gon\n", i, getN(i));
   fprintf(stderr,"\n");
}

void symmetro::setSym( const int &psym, const int &qsym )
{
   p = psym;
   q = qsym;
}

void symmetro::setMult( const int &a, const int &m )
{
	mult[a] = m;
}

int symmetro::getOrder( const int &a )
{
   switch( a ) {
   case 0: return p;
   case 1: return q;
   case 2: return 2;
   default: return 0;
   }
}

int symmetro::getN( const int &a )
{
	return ( getOrder( a ) * mult[ a ] );
}

double symmetro::getAngleBetweenAxes( const int &axis1, const int &axis2 ) {
   return ( acos(vdot(sym_vec[axis1].unit(), sym_vec[axis2].unit())) );
}

// RK - fill symvec here for angle_between_axes
void symmetro::fill_sym_vec() {
   // unit vectors on symmetry axes on centers of polygons
   if ( ( p == 5 && q == 3 ) || ( p == 3 && q == 5 ) ) {
      //double phi = (1 + sqrt(5))/2;
      sym_vec[0] = vec3d( 0.0, 1.0, phi );
      sym_vec[1] = vec3d( 1.0, 1.0, 1.0 );
      sym_vec[2] = vec3d( 1.0, phi+1, phi );
   }
   else 
   if ( ( p == 4 && q == 3 ) || ( p == 3 && q == 4 ) ) {
      sym_vec[0] = vec3d( 0.0, 0.0, 1.0 );
      sym_vec[1] = vec3d( 1.0, 1.0, 1.0 );
      sym_vec[2] = vec3d( 0.0, 1.0, 1.0 );
   }
   else
   if ( p == 3 && q == 3 ) {
      sym_vec[0] = vec3d(  1.0, 1.0, 1.0 );
      sym_vec[1] = vec3d( -1.0, 1.0, 1.0);
      sym_vec[2] = vec3d(  0.0, 0.0, 1.0);
   }
   
   //for( int i=0; i<3; i++ )
   //   fprintf(stderr,"%.17lf %.17lf %.17lf\n", sym_vec[i][0], sym_vec[i][1], sym_vec[i][2]);
   //fprintf(stderr,"\n");
}

double symmetro::angle( const int &n, const int &d )
{
   return ( (2.0*M_PI*double(d)/double(n)) );
}

double symmetro::circumradius( const int &n, const int &d )
{
//fprintf(stderr,"angle = %.17lf\n",angle(n,d));
   return ( 0.5/sin(angle(n,d)/2.0) );
}

vector<col_geom_v> symmetro::makePolygons( const symmetro_opts &opts )
{
   vector<double> scales(2);
   scales[0] = 1.0;
   scales[1] = 1.0;
   if ( opts.scale )
      scales[( opts.scale_direction[0] < opts.scale_direction[1] ) ? 0 : 1] = opts.scale;

   vector<int> axis; 
   for( int i=0; i<(int)opts.multipliers.size(); i++ ) {
      if ( opts.multipliers[i] ) {
         axis.push_back(i);
      }            
   }
   
   double r0 = scales[0] * circumradius( getN(axis[0]), 1.0 );
   double r1 = scales[1] * circumradius( getN(axis[1]), 1.0 );

   double angle_between_axes = getAngleBetweenAxes( axis[0], axis[1] );
   if ( opts.verbose )
      fprintf(stderr,"angle_between_axes = %.17lf\n",rad2deg(angle_between_axes));
   mat3d rot = mat3d::rot(vec3d(0, 1, 0), angle_between_axes);
   mat3d rot_inv = mat3d::rot(vec3d(0, 1, 0), -angle_between_axes);

   double ang = opts.rotation;
   ang = ang * (2.0*M_PI*1.0/getN(axis[0]))/2.0;
   
   vec3d V = mat3d::rot(vec3d(0, 0, 1), ang) * vec3d(r0, 0, 0);
   vec3d q = rot * V;
   vec3d u = rot * vec3d(0, 0, 1);
   
   double a = u[0]*u[0] + u[1]*u[1];
   double b = 2*(q[0]*u[0] + q[1]*u[1]);
   double c = q[0]*q[0] + q[1]*q[1] - r1*r1;

   double disc = b*b - 4*a*c;
   if (disc < -epsilon)
     opts.error("model is not geometrically constructible");
   else
   if (disc < 0)
     disc = 0;

   double t = (-b - sqrt(disc))/(2*a);

   vec3d P = V + vec3d(0, 0, t);
   vec3d Q = rot * P;
   
   // there can only be 2, one polygon will not be filled
   vector<col_geom_v> pgeom(3);
   
   for( int i=0; i<2; i++ ) {
      int j = axis[i];
	   int n = getN(j);
	   if( (n > 0) && (scales[i] > epsilon) ) {
         for( int idx = 0; idx < n; ++idx ) {
            if ( i == 0 ) {
	            pgeom[j].add_vert( mat3d::rot(vec3d(0, 0, 1), idx * angle(n,1.0)) * P );
	         }
	         else
	         if ( i == 1 ) {
	            pgeom[j].add_vert( rot_inv * mat3d::rot(vec3d(0, 0, 1), idx * angle(n,1.0)) * Q );
	         }
         }

	      pgeom[j].transform( mat3d::refl(vec3d(0, 0, 1)) );
	      pgeom[j].transform( mat3d::alignment(vec3d(0, 0, 1), vec3d(1, 0, 0), sym_vec[axis[0]], sym_vec[axis[1]]) ); 
         
         vector<int> face;
         for( int idx = 0; idx < n; ++idx )
            face.push_back( idx );
         pgeom[j].add_face( face );
      }
      
      // epsilon size faces are because scale was set at 0
      if ( fabs( scales[i] ) <= epsilon ) {
         pgeom[j].clear_all();
      }
   }
   
   return pgeom;   
}

bool is_point_on_polygon_edges(const geom_if &polygon, const vec3d &P, const double &eps)
{
   const vector<int> &face = polygon.faces()[0];
   const vector<vec3d> &verts = polygon.verts();

   bool answer = false;

   int fsz = face.size();
   for(int i=0; i<fsz; i++) {
      vec3d v1 = verts[face[i]];
      vec3d v2 = verts[face[(i+1)%fsz]];
      if ((point_in_segment(P, v1, v2, eps)).is_set()) {
         answer = true;
         break;
      }
   }

   return answer;
}

bool detect_collision( const col_geom_v &geom )
{
   const vector<vector<int> > &faces = geom.faces();
   const vector<vec3d> &verts = geom.verts();
   
   // don't test digons
   for( int i=0; i<(int)faces.size(); i++) {
      vector<int> face0 = faces[i];
      // digons won't work in plane intersection
      if ( face0.size() < 3 )
         continue;
      for( int j=i+1; j<(int)faces.size(); j++) {
         vector<int> face1 = faces[j];
         if ( face1.size() < 3 )
            continue;

         vec3d P, dir;
         if ( two_plane_intersect(  centroid(verts, face0), face_norm(verts, face0),
                                    centroid(verts, face1), face_norm(verts, face1),
                                    P, dir, epsilon ) ) {
            if ( !P.is_set() )
               continue;               
            // if two polygons intersect, see if intersection point is inside polygon
            vector<int> face_idxs;
            face_idxs.push_back(i);
            col_geom_v polygon = faces_to_geom(geom, face_idxs);
            int winding_number = 0;
            // get winding number, if not zero, point is on a polygon
            wn_PnPoly( polygon, P, 2, winding_number, epsilon );
            // if point in on an edge set winding number back to zero
            if ( winding_number ) {
               if ( is_point_on_polygon_edges( polygon, P, epsilon ) )
                  winding_number = 0;
            }
            if ( winding_number ) {
               return true;
            }
         }
      }
   }
   
   return false;
}

col_geom_v build_geom(vector<col_geom_v> &pgeom, const symmetro_opts &opts)
{
   col_geom_v geom;
   
   for( int i=0; i<3; i++ ) {
      // skip empty geoms
      if ( !pgeom[i].verts().size() )
         continue;
         
      // if not polygon, repeat for symmetry type
      if ( opts.convex_hull > 1 ) {
         sch_sym sym; 
         if (opts.sym == 1)
            sym.init(sch_sym::I);
         else
         if (opts.sym == 2)
            sym.init(sch_sym::O);
         else
         if (opts.sym == 3)
            sym.init(sch_sym::T);
         sym_repeat(pgeom[i], pgeom[i], sym, ELEM_FACES);
      }
      
      coloring clrng(&pgeom[i]);
      col_val c;
      if (i == 0)
         c = opts.axis_0_color;
      else if (i == 1)
         c = opts.axis_1_color;
      else if (i == 2)
         c = opts.axis_2_color; 
      // if color is unset we need to keep track of it because of convex hull coloring
      // use index of INT_MAX as a place holder
      if (!c.is_set())
         c = INT_MAX;
      clrng.f_one_col(c);
      
      geom.append(pgeom[i]);
   }
   
   if ( opts.convex_hull > 1 )
      sort_merge_elems(geom, "vf", epsilon);
   
   // check for collisions
   bool collision = false;
   if ( opts.convex_hull > 2 )
      collision = detect_collision( geom );
      if ( collision ) {
         fprintf(stderr,"collision detected. convex hull is suppressed\n");
   }
   
   // if making unit edges, to it before convex hull
   if (opts.unitize) {
      geom_info info(geom);
      if (info.num_iedges() > 0) {
         double val = info.iedge_lengths().sum/info.num_iedges();
         geom.transform(mat3d::scale(1/val));
      }
   }      
   
   if ( ( !collision && opts.convex_hull == 4 ) || ( opts.convex_hull == 3 ) ) {
      char errmsg[MSG_SZ];
      int ret = geom.add_hull("",errmsg);
      if(!ret)
         if (opts.verbose)
            fprintf(stderr,"%s\n",errmsg);

      // merged faces will retain RGB color
      sort_merge_elems(geom, "f", epsilon);

      // after sort merge, only new faces from convex hull will be uncolored
      for( int i=0; i<(int)geom.faces().size(); i++ ) {
         col_val c = geom.get_f_col(i);
         if (!c.is_set())
            geom.set_f_col(i, opts.convex_hull_color);
      }
   }
   
   // unset the colors for real
   for( int i=0; i<(int)geom.faces().size(); i++ ) {
      col_val c = geom.get_f_col(i);
      if ( c.is_idx() && c == INT_MAX )
         geom.set_f_col(i, col_val());
   }
   
   geom.orient();
   
   return geom;
}

int main(int argc, char *argv[])
{
   symmetro_opts opts;
   opts.process_command_line(argc, argv);
 
   // some pre-checking  
   int num_multipliers = 0;
   for( int i=0; i<(int)opts.multipliers.size(); i++ ) {
      if ( opts.multipliers[i]>0 )
         num_multipliers++;
   }
   
   if ( num_multipliers == 3 )
      opts.error("one multiplier must be 0",'m');
      
   // 1 polygon is done by scale to 0
   if ( num_multipliers == 1 )
   {
      int second_direction = 0;
      bool done = false;
    
      opts.scale_direction.clear();  
      for( int i=0; i<3; i++ ) {
         if ( opts.multipliers[i] != 0 )
            second_direction = i;
         else {
            if (!done) {
               opts.multipliers[i] = 1;
               opts.scale = epsilon/10.0;
               opts.scale_direction.push_back(i);
               done = true;
            }
         }
      }
      
      opts.scale_direction.push_back(second_direction);
   }
      
   symmetro s;
   
   if ( opts.sym == 1 ) {
      s.setSym( 5, 3 );
   }
   else
   if ( opts.sym == 2 )
      s.setSym( 4, 3 );
   else
   if ( opts.sym == 3 )
      s.setSym( 3, 3 );
   else
      opts.error("symmetry type not set",'s');
   
   for( int i=0; i<(int)opts.multipliers.size(); i++ ) {
      s.setMult( i, opts.multipliers[i] );
   }
   
   if ( !opts.scale_direction.size() ) {
      if ( opts.scale && opts.multipliers.size() == 1 )
         opts.error("cannot scale when only one polygon is generated",'S');

      int j = 0;
      for( int i=0; i<(int)opts.multipliers.size(); i++ ) {
         if ( j == 2 )
            continue;
         if ( opts.multipliers[i] ) {
            opts.scale_direction.push_back(i);
            j++;
         }            
      }
   }   

   for( int i=0; i<(int)opts.scale_direction.size(); i++ ) {
      if ( opts.multipliers[opts.scale_direction[i]] == 0 )
         opts.error(msg_str("polygon '%d' is not generated so cannot be used for scaling", opts.scale_direction[i]), 'S');
   }

   // calculate axes here      
   s.fill_sym_vec();
   
   vector<col_geom_v> pgeom = s.makePolygons( opts );   
   
   if ( opts.verbose ) {
      s.debug();
      
      double edge_length[3];
      for( int i=0; i<(int)pgeom.size(); i++ ) {
         geom_info info(pgeom[i]);
         if (info.num_iedges() > 0) {
            edge_length[i] = info.iedge_lengths().sum/info.num_iedges();
            fprintf(stderr,"Edge length of polygon %d = %.17lf %s\n", i, edge_length[i], (opts.unitize ? "(before unit edges)" : ""));
         }
      }
      
      fprintf(stderr,"\n");
      for( int i=0; i<3; i++ ) {
         for( int j=0; j<3; j++ ) {
            if (i==j)
               continue;
            if ( edge_length[i] > epsilon && edge_length[j] > epsilon )
               fprintf(stderr,"edge length ratio of polygon %d to %d = %.17lf\n", i, j, edge_length[i] / edge_length[j] );
         }
      }
      
      fprintf(stderr,"\n");
   }
   
   col_geom_v geom;
   geom = build_geom(pgeom, opts);
   
   geom_write_or_error(geom, opts.ofile, opts);
   
   return 0;
}

