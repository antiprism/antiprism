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
      bool reverse;
      vector<int> multipliers;
      vector<bool> theta;
      char rotation_method;
      vector<double> extra_rotation;
      vector<int> scale_direction;
      double scale;
      int convex_hull;
      bool unitize;
      int patch;
      col_val axis_0_color;
      col_val axis_1_color;
      col_val axis_2_color;
      col_val convex_hull_color;
      bool verbose;
      string ofile;

      symmetro_opts(): prog_opts("symmetro"),
                       sym(0),
                       reverse(false),
                       rotation_method('v'),
                       scale(0.0),
                       convex_hull(4),
                       unitize(false),
                       patch(0),
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
"This program is derived from a program written by Craig S. Kaplan\n"
"and George W. Hart (http://www.georgehart.com). Project information\n"
"can be found at http://www.cgl.uwaterloo.ca/~csk/projects/symmetrohedra\n"
"\n"
"Symmetrohedra are created by placing equilateral polygons centered on\n"
"the symmetry axes of Icosahedral, Octahedral, or Tetrahedral symmetry.\n"
"The number of sides of the polygons will be a multiple number of the\n"
"axis reflection number. Axes are in order as 0, 1 and 2 corresponding\n"
"to icosahedral {5,3,2}, octahedral {4,3,2}, or tetrahedral {3,3,2}\n"
"The end result vertices are all 1 unit from the polyhedron center. Note\n"
"that when all three multipliers are used the solution will generally\n"
"not yield polygons with equal edge length. The one exception is the\n"
"truncated octahedron\n"
"\n"
"These types of polyhedra will either have all polygons touching on edge\n"
"or all on vertices. In a case where a vertex meets an edge, a warning\n"
"will be displayed and the model will be generated. In this case, convex\n"
"hull will be suppressed. Try using -r to rotate the polygons to possible\n"
"valid combinations. It is also possible to size polygons such that they\n"
"intersect. If a collision is detected, convex hull will be suppressed\n" 
"\n"
"Options\n"
"%s"
"  -s <type> symmetry type of Symmetrohedra. sets {p,q}\n"
"               icosahedral {5,3}, octahedral {4,3}, or tetrahedral {3,3}\n"
"  -R        reverse p and q. e.g. {5,3,2} becomes {3,5,2}\n"
"  -m <vals> multipliers for axis polygons. Given as three integers\n"
"               separated by commas. e.g. 2,3,0\n"
"  -r <vals> which axis polygons are rotated to edge or on point. Up to\n"
"               three values from 0, 1 and 2 separated by commas. e.g. 0,1\n"
"               or use v - connect on vertex  e - connect on edge  (default: v)\n"
"  -q <vals> angles in degrees to add rotation. Given as three floating\n"
"               point numbers separated by commas. e.g. 45,60,45\n"
"  -S <a,b,s> scale s, applied from axis a polygon applied to axis b polygon\n"
"               e.g. 0,1,0.5  (default: calculated for unit edge length)\n"
"  -C <mode> convex hull. polygons=1, suppress=2, force=3, normal=4\n"
"               (default: 4)\n"
"  -u        make the average edge length 1 unit (performed before convex hull)\n"
//"  -g <val>  edge=1, vertex=2  force meeting at edge or vertex to true. patch\n"
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
   vector<double> scale_direction_tmp;
   
   handle_long_opts(argc, argv);

   while( (c = getopt(argc, argv, ":hs:Rm:r:q:S:C:ug:x:y:z:w:Vo:")) != -1 ) {
      if(common_opts(c, optopt))
         continue;

      switch(c) {
         case 's': // symmetry
            id = get_arg_id(optarg, "icosahedral=1|octahedral=2|tetrahedral=3", argmatch_add_id_maps, errmsg);
            if(id=="")
               error(errmsg);
            sym = atoi(id.c_str());
            break;
            
         case 'R':
            reverse = true;
            break;
            
         case 'm': // multiplier
            if(!read_int_list(optarg, multipliers, errmsg, true, 3))
               error(errmsg, c);
            if((int)multipliers.size() != 3)
               error("multipliers must be specified as three integers", c);
            if(multipliers[0]+multipliers[1]+multipliers[2] == 0)
               error("at least one multiplier must be positive", c);
            if(multipliers[2] == 1)
               //error("multiplier for axix 2 cannot be 1", c);
               warning("model will contain digons");
            break;
            
         case 'r': // theta (rotate from side to point)
            for( int i=0; i<3; i++ )
               theta.push_back(false);
               
            if(strchr("e", *optarg))
               rotation_method = *optarg;
            else
            if(strchr("v", *optarg))
               rotation_method = *optarg;
            else {
               vector<int> theta_input;
               if(!read_int_list(optarg, theta_input, errmsg, true, 3))
                  error(errmsg, c);
               for( int i=0; i<(int)theta_input.size(); i++ ) {
                  if( theta_input[i] > 2 )
                     error("axes specified for rotation should be 0, 1 or 2", c);
                  for( int j=i+1; j<(int)theta_input.size(); j++ )
                     if( theta_input[i] == theta_input[j] )
                         error(msg_str("an axis number is specified more than once: '%d'", theta_input[j]), c);
                  theta[ theta_input[i] ] = true;
               }
            }  
            break;
            
         case 'q': // extra rotation to add to theta
            if(!read_double_list(optarg, extra_rotation, errmsg, 3))
               error(errmsg, c);
            if((int)extra_rotation.size() != 3)
               error("extra rotation values must be specified as three floating point numbers", c);
            break;
           
         case 'S': // scale direction and scale
            if(!read_double_list(optarg, scale_direction_tmp, errmsg, 3))
               error(errmsg, c);
            
            // place integer portions of possible doubles
            for( int i=0; i<2; i++ ) {
               double a = (int)floor(scale_direction_tmp[i]);
               if ( scale_direction_tmp[i] - a > 0.0 )
                  error(msg_str("axis numbers must be specified by an integer: '%g'", scale_direction_tmp[i]), c);
               scale_direction.push_back((int)a);
            }
               
            // pull out ratio
            scale = scale_direction_tmp[2];
            if ( scale == 0.0 )
               scale = DBL_MIN;
               
            scale_direction_tmp.clear();
               
            for( int i=0; i<(int)scale_direction.size(); i++ ) {
               if( scale_direction[i] > 2 )
                  error("ratio direction should be 0, 1 or 2", c);
               for( int j=i+1; j<(int)scale_direction.size(); j++ )
                  if( scale_direction[i] == scale_direction[j] )
                      error(msg_str("an axis number is specified more than once: '%d'", scale_direction[j]), c);
            }      
            break;
            
         case 'C':
            id = get_arg_id(optarg, "polygons=1|suppress=2|force=3|normal=4", argmatch_add_id_maps, errmsg);
            if(id=="")
               error(errmsg,c);
            convex_hull = atoi(id.c_str());
            break;
            
         case 'u':
            unitize = true;
            break;
            
         case 'g':
            id = get_arg_id(optarg, "edge=1|vertex=2", argmatch_add_id_maps, errmsg);
            if(id=="")
               error(errmsg,c);
            patch = atoi(id.c_str()); 
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
		fill( scale, scale + 3, 0.0 );
		fill( theta, theta + 3, 0.0 );
		
		// RK - added items
		for( int i=0; i<3; i++ )
		   sym_vec.push_back(vec3d(0.0,0.0,0.0));
		
		force_edge = false;
		force_vertex = false;
		
		old_tie = false;
	}
	
	// RK - added methods
	void setSym( const int &psym, const int &qsym )
	{
	   p = psym;
	   q = qsym;
	   
		// fill fund[] here
		ssym_repl( p, q );
	}
	
	double getScale( const int &a )
	{
	   return ( scale[a] );
	}
	
	double getTheta( const int &a )
	{
	   return ( theta[a] );
	}
	
	void setPatch( const int &a )
	{
	   if ( a == 1 )
	      force_edge = true;
	   else
	   if ( a == 2 )
	      force_vertex = true;
	}
	
   // RK - return polygon number of sides
	int getN( const int &a ) {
   	return ( getOrder( a ) * mult[ a ] );
	}
	
	// RK - return angle_between_axes
	double getAngleBetweenAxes( const int &axis1, const int &axis2 ) {
	   return ( acos(vdot(sym_vec[axis1].unit(), sym_vec[axis2].unit())) );
	}
	
	// RK - fill symvec here for angle_between_axes (used to be in makePolygons)
	void fill_sym_vec( const bool &reverse ) {
      // unit vectors on symmetry axes on centers of polygons
      if ( ( p == 5 && q == 3 ) || ( p == 3 && q == 5 ) ) {
         double a = sqrt((10.0-sqrt(20.0))/20.0);
         double b = sqrt((10.0+sqrt(20.0))/20.0);
         double c = sqrt((3.0-sqrt(5.0))/6.0);
         double d = sqrt((3.0+sqrt(5.0))/6.0);
         sym_vec[0] = vec3d( 0.0, a, b );
         sym_vec[1] = vec3d( c, 0.0, d );
         sym_vec[2] = vec3d( 0.0, 0.0, 1.0 );
      }
      else 
      if ( ( p == 4 && q == 3 ) || ( p == 3 && q == 4 ) ) {
         double a = sqrt(1.0/3.0);
         double b = sqrt(0.5);
         sym_vec[0] = vec3d( 0.0, 0.0, 1.0 );
         sym_vec[1] = vec3d( a, a, a );
         sym_vec[2] = vec3d( 0.0, b, b );
      }
      else
      if ( p == 3 && q == 3 ) {
         double a = sqrt(1.0/3.0);
         sym_vec[0] = vec3d( a, a, a );
         sym_vec[1] = vec3d( -a, a, a);
         sym_vec[2] = vec3d( 0.0, 0.0, 1.0);
      }
      if ( reverse )
         swap( sym_vec[0], sym_vec[1] );
      
      //for( int i=0; i<3; i++ )
      //   fprintf(stderr,"%.17lf %.17lf %.17lf\n", sym_vec[i][0], sym_vec[i][1], sym_vec[i][2]);
      //fprintf(stderr,"\n");
	}
	
	// RK - wrapper for only one polygon
	// unfortunately needs some hacks for this
	void tie1()
	{
	   // tie will perform as though indexes 1 and 2 were never swapped (patch)
		old_tie = true;
	
	   // hack for platonics as multiple of 1 does not give the correct answer
	   if ( mult[ 0 ] + mult[ 1 ] + mult[ 2 ] == 1 ) {
	      int idx = mult[0] ? 0 : 1;
	      if ( idx == 0 )
            setTheta( 0, false );
         else
         if ( idx == 1 )
            setTheta( 1, true );
            
	      double face_radius = 0.0;   
         if ( ( p == 5 && q == 3 ) || ( p == 3 && q == 5 ) )
            face_radius = acos(sqrt((5.0+sqrt(20.0))/15.0));
         else 
         if ( ( p == 4 && q == 3 ) || ( p == 3 && q == 4 ) )
            face_radius = acos(sqrt(1.0/3.0));
         else 
         if ( p == 3 && q == 3 )
            face_radius = acos(1.0/3.0);
	      
	      scale[ idx ] = face_radius;
	   }
	   else
	   if ( mult[ 0 ] || mult[ 1 ] ) {
	      // axis 0 alone returns the wrong value from tie()
	      // move axis 0 to axis 1 for correct calculation from tie()
	      // also move theta value along with it
	      if ( mult[ 0 ] ) {
	         mult[ 1 ] = mult[ 0 ];
	         if ( theta[ 0 ] )
	            setTheta( 1, true );
	         mult[ 0 ] = 0;
	         setTheta( 0, false );
	      }
	      
	      // temporarily set 2-fold symmetry to 1, creating digon
	      mult[ 2 ] = 1;
	      
	      // do tie for (a,b,c)
	      tie();
	      
	      // set scale so digon will not be generated
	      scale[ 2 ] = 0.0;
	      // reset 2-fold multiplier
	      mult[ 2 ] = 0;
	   }  
	   else
	   if ( mult[ 2 ] ) {
	      if ( theta[ 2 ] ) {
            // create a triangle we will not use
            mult[ 1 ] = 1;
            setTheta( 1, true );
            
            if ( !is_even( mult[ 2 ] ) )
               setTheta( 2, false );
            
            if ( !is_even( mult[ 2 ] ) )
               // force edge meeting to true
               force_edge = true;
             
            tie(1, 2);
            
            // discard triangle
            scale[ 1 ] = 0.0;
            // reset multiplier
            mult[ 1 ] = 0;
         }
         else {
            double face_radius = 0.0;
            if ( ( p == 5 && q == 3 ) || ( p == 3 && q == 5 ) ) {
	            face_radius = acos(sqrt((3.0+sqrt(5.0))/6.0));
	            if ( !is_even( mult[ 2 ] ) )
                  setTheta( 2, true );
	         }
            else 
            if ( ( p == 4 && q == 3 ) || ( p == 3 && q == 4 ) ) {
               // for even mult[ 2 ] edge-wise pairing is possible
               if ( !is_even( mult[ 2 ] ) )
                  face_radius = acos(sqrt(0.5));
               else {
                  // alternative form
                  // use multiplier 1 for correct distance
                  mult[ 1 ] = mult[ 2 ] / 2;
                  
                  setTheta( 1, true );
	               
	               // temporarily set 2-fold symmetry to 1, creating digon
	               mult[ 2 ] = 1;
	               
	               // do tie for (a,b,c)
	               tie();

	               face_radius = scale[ 1 ];
	               
	               // restore axis 2
	               mult[ 2 ] = mult[ 1 ] * 2;
	               setTheta( 2, true );
	               
	               // discard axis 1
                  scale[ 1 ] = 0.0;
	               // reset multiplier
	               mult[ 1 ] = 0;
	            }
            }
            else 
            if ( p == 3 && q == 3 )
               face_radius = acos(sqrt(1.0/3.0));
                     
            scale[ 2 ] = face_radius;
         }
	   }
	}
	
	// RK - error checking wrapper for tie()
	int tie3( const int &a, const int &b, const int &c, char *errmsg)
	{
	   int ret = 0;
	   if( isEdgeOn( a, b ) && isEdgeOn( b, a ) && isEdgeOn( a, c ) && isEdgeOn( c, a ) && isEdgeOn( b, c ) && isEdgeOn( c, b ) )
	      ret = 1;
	   else
	   if( isVertexOn( a, b ) && isVertexOn( b, a ) && isVertexOn( a, c ) && isVertexOn( c, a ) && isVertexOn( b, c ) && isVertexOn( c, b ) )
	      ret = 1;
	      
	   tie();
	   
	   if (!ret)
	      strcpy(errmsg, "to tie three polygons correctly, they must all meet at edges or vertices");
	   
	   return ret;
	}
	
	//
	// scale, theta calculations need this code
	//
	
	// from class ssystem
	bool isEdgeOn( const int &a, const int &b );
	bool isVertexOn( const int &a, const int &b );
	int getOrder( const int &a );
	double getAlpha( const int &a, const int &b );
	double getEdgeLength( const int &a, const int &b );
   double getLengthOn( const int &a, const int &b );
   void setLengthOn( const int &a, const int &b, const double &s );

   void tieTo( const int &a, const int &b );
	int tie( const int &a, const int &b ); // RK - added success/fail return flag
	void tie();
	
	void debug();
	
   void pin( const int &a );
	void setMult( const int &a, const int &m );
   void setScale( const int &a, const double &s );
   void setTheta( const int &a, const double &t );
   void setTheta( const int &a, const bool &b );
   
   //
   // code needed to calculate fund array
   //
   
   // from class stransform (antiprism: mat3d)
   mat3d toPole( const vec3d &p );
   mat3d fromPole( const vec3d &p );
   mat3d rotateAboutPole( const vec3d &p, double theta );
   
   // from class spoint (antiprism: vec3d)
   vec3d moveZ( const double &dist );
   vec3d moveZ( const double &dist, const double &theta );
   
   // from class ssystem
	void ssym_repl( const int &p, const int &q );
	
	//
	// code needed to make polygons
	//
   
   // from class stransform
   mat3d match( const int &idx, const vec3d &p, const bool &reverse );
   
   // from class rawreceiver
   vector<col_geom_v> makePolygons( const bool &reverse );
	
	~symmetro()
	{}
	
private:
	int 		p;
	int 		q;

	int 		mult[3];
	double   scale[3];
	double	theta[3];
	
	vec3d 	fund[3];
	
	// RK - added items
   vector<vec3d> sym_vec;
	
	bool force_edge;
	bool force_vertex;
	bool old_tie;
};

// patch with force_edge when the function fails
bool symmetro::isEdgeOn( const int &a, const int &b )
{
   if ( force_edge )
      return true;
   else
   if ( force_vertex )
      return false;

	int pm = -1;

	if( b == ((a+1)%3) ) {
		pm = 1;
	}

	if( (pm == 1) || ((pm == -1) && ((mult[a]%2)==0)) ) {
		if( fabs( theta[a] ) > epsilon ) { // was 1e-7
			return true;
		} else {
			return false;
		}
	} else {
		if( fabs( theta[a] ) > epsilon ) { // was 1e-7
			return false;
		} else {
			return true;
		}
	}
}

bool symmetro::isVertexOn( const int &a, const int &b )
{
   if ( force_vertex )
      return true;
   else
   if ( force_edge )
      return false;
      
   return !isEdgeOn( a, b );
}

// RK 1 and 2 were reversed so that we can have p,q,2
int symmetro::getOrder( const int &a )
{
	switch( a ) {
	case 0: return p;
	case 1: return q;
	case 2: return 2;
	default: return 0;
	}
}

double symmetro::getAlpha( const int &a, const int &b )
{
	if( isVertexOn( a, b ) ) {
		return 1.0;
	} else {
		double n = double( getN( a ) );
		return cos( M_PI / n );
	}
}
double symmetro::getEdgeLength( const int &a, const int &b )
{
//fprintf(stderr,"edge length = %.17lf\n", vdot( fund[a], fund[b] ));
	return acos( vdot( fund[a], fund[b] ) );
}

double symmetro::getLengthOn( const int &a, const int &b )
{
	if( isVertexOn( a, b ) ) {
		return scale[a];
	} else {
		double n = double( getN( a ) );
		/*
		double al = sin( M_PI / n );
		return atan( al * tan( scale[a] ) );
		*/
		double al = cos( M_PI / n );
		return atan( al * tan( scale[a] ) );
	}
}

void symmetro::setLengthOn( const int &a, const int &b, const double &s )
{
	if( mult[a] == 0 ) {
		return;
	}

	if( isVertexOn( a, b ) ) {
		scale[a] = s;
	} else {
		double n = double( getN( a ) );
		double al = cos( M_PI / n );
		scale[a] = atan( tan( s ) / al );
	}
}

void symmetro::tieTo( const int &a, const int &b )
{
	double s = getLengthOn( b, a );
	if( s > 0.0 ) {
		setLengthOn( a, b, getEdgeLength( a, b ) - s );
	}
}

// RK - tie can fail when vertex meets edge. set return flag if it does
int symmetro::tie( const int &a, const int &b )
{
	double delta = getEdgeLength( a, b );
	double na = double( getN( a ) );
	double nb = double( getN( b ) );

	if( isEdgeOn( a, b ) && isEdgeOn( b, a ) ) {
		double tpa = 1.0 / tan( M_PI / na );
		double tpb = 1.0 / tan( M_PI / nb );

		double tanx = sin( delta ) / 
			sqrt( tpa*tpa + tpb*tpb + 2.0*cos(delta) * tpa*tpb );
		double el = atan( tanx );

		scale[a] = asin( sin( el ) / sin( M_PI / na ) );
		scale[b] = asin( sin( el ) / sin( M_PI / nb ) );
		
		return 1;
	} else if( isVertexOn( a, b ) && isVertexOn( b, a ) ) {
		double spa = 1.0 / sin( M_PI / na );
		double spb = 1.0 / sin( M_PI / nb );

		double sinx = sin( delta ) / 
			sqrt( spa*spa + spb*spb + 2.0*cos(delta) * spa*spb );
		double el = asin( sinx );

		scale[a] = asin( sin( el ) / sin( M_PI / na ) );
		scale[b] = asin( sin( el ) / sin( M_PI / nb ) );
		
		return 1;
	}
	
	return 0;
}

void symmetro::tie()
{
   // RK - because I swapped indexes 1 and 2, bb and cc will be swapped unless patch
   int aa = 0;
   int bb = 2;
   int cc = 1;
   
   // for single polygons it works in the original form
   if ( old_tie )
      swap( bb, cc );
   
	double D1 = tan( getEdgeLength( aa, bb ) ); // was 0 1
	double D2 = tan( getEdgeLength( bb, cc ) ); // was 1 2
	double D3 = tan( getEdgeLength( cc, aa ) ); // was 2 0

	double a12 = getAlpha( aa, bb ); // was 0 1
	double a21 = getAlpha( bb, cc ); // was 1 0
	double a13 = getAlpha( aa, cc ); // was 0 2
	double a31 = getAlpha( cc, aa ); // was 2 0
	double a23 = getAlpha( bb, cc ); // was 1 2
	double a32 = getAlpha( cc, bb ); // was 2 1
	
	/*
   fprintf(stderr,"D1 = %.17lf\n",D1);
   fprintf(stderr,"D2 = %.17lf\n",D2);
   fprintf(stderr,"D3 = %.17lf\n",D3);
   
   fprintf(stderr,"a12 = %.17lf\n",a12);
   fprintf(stderr,"a21 = %.17lf\n",a21);
   fprintf(stderr,"a13 = %.17lf\n",a13);
   fprintf(stderr,"a31 = %.17lf\n",a31);
   fprintf(stderr,"a23 = %.17lf\n",a23);
   fprintf(stderr,"a32 = %.17lf\n",a32);
   */

/*
	cout << "a12 = " << a12 << ", a21 = " << a21 << ", a13 = " << a13 << endl;
	cout << "a31 = " << a31 << ", a23 = " << a23 << ", a32 = " << a32 << endl;
	cout << "D1 = " << D1 << ", D2 = " << D2 << ", D3 = " << D3 << endl;
*/

	double a = 
		-D2*a21*a23*a32*a13 +
		-a13*D1*a12*a21*a23 +
		-D1*D2*D3*a12*a21*a23*a32 +
		D3*a13*a21*a23*a31;
		// D3 - D1 - D2 - D1*D2*D3;
	double b = 
		-a13*a21*a32 + 
		D1*D2*a13*a23*a32 +
		-a31*a12*a23 +
		D1*D2*a31*a12*a21 +
		-D2*D3*a12*a23*a32 +
		-D1*D3*a12*a21*a32 +
		-D1*D3*a12*a31*a23 +
		-D2*D3*a13*a31*a21;
		// 2.0*(D1*D2 - 1.0 - D2*D3 - D1*D3);
	double c =
		D1*a13*a32 +
		a31*a12*D2 +
		-D3*a12*a32 +
		a13*a31*D1*D2*D3;
		// D1 + D2 - D3 + D1*D2*D3;
	
	double disc = sqrt( b*b - 4.0*a*c );
	double B = (-b + disc) / (2.0*a);
	if( B < 0.0 ) {
		B = (-b - disc) / (2.0*a);
	}
	double A = (D1-a21*B)/(a12+D1*a12*a21*B);
	double C = (D2-a23*B)/(a32+D2*a23*a32*B);

/*
	cout << "A = " << A << ", B = " << B << ", C = " << C << endl;
	cout << "MAGIC: " << 
		(B*B*(-D2-D1-D1*D2*D3+D3)+B*(-1.0+D1*D2-1.0+D1*D2-D2*D3-D1*D3-D1*D3-D2*D3)+(D1+D2-D3+D1*D2*D3)) << endl;

	cout << "a12*A+a21*B = " << (a12*A + a21*B) << endl;
	cout << "D1(1-a12*a21*A*B) = " << (D1*(1.0-a12*a21*A*B)) << endl << endl;
	cout << "tan(l12+l21) = " << ((a12*A+a21*B)/(1.0-a12*a21*A*B)) <<endl<<endl;
	cout << "l12+l21 = " << atan((a12*A+a21*B)/(1.0-a12*a21*A*B)) <<endl<<endl;
	cout << "atan(D1) = " << atan(D1) << endl;

	double l12 = atan( a12 * A );
	double l21 = atan( a21 * B );

	cout << "l12 = " << l12 << endl;
	cout << "l21 = " << l21 << endl;

	cout << "t+t = " << ((tan(l12)+tan(l21))/(1.0-tan(l12)*tan(l21))) << endl;
*/

   // RK - because I swapped indexes 1 and 2
	scale[ aa ] = atan( A );
	scale[ bb ] = atan( B );
	scale[ cc ] = atan( C );

/*
	if( scale[0] < 0.0 ) {
		scale[0] += M_PI;
	}
	if( scale[1] < 0.0 ) {
		scale[1] += M_PI;
	}
	if( scale[2] < 0.0 ) {
		scale[2] += M_PI;
	}
*/

/*
	double d1 = getEdgeLength( 0, 1 );
	double d2 = getEdgeLength( 1, 2 );
	double d3 = getEdgeLength( 2, 0 );
	scale[0] = 0.5*( d1+d3-d2 );
	scale[1] = d1 - scale[0];
	scale[2] = d3 - scale[0];
*/
}

void symmetro::debug()
{
   fprintf(stderr,"\nsymmetry = {%d,%d,2}\n\n", p, q);
   
   for( int i=0; i<3; i++ )
      fprintf(stderr,"axis %d: mult = %d, scale = %.17lf, theta = %.17lf (%.17lf degrees)\n", i, mult[i], scale[i], theta[i], rad2deg(theta[i]));
   fprintf(stderr,"\n");
   
   for( int i=0; i<3; i++ )
      if ( mult[i] )
         fprintf(stderr,"axis %d: face center radius = cosine(scale) = %.17lf\n", i, cos( scale[i] ) );
   fprintf(stderr,"\n");
   
   for( int i=0; i<3; i++ )
      if ( mult[i] )
         fprintf(stderr,"axis %d polygon: %d-gon\n", i, getN(i));
   fprintf(stderr,"\n");
}

void symmetro::pin( const int &a )
{
	double d1 = getEdgeLength( a, (a+1)%3 );
	double d2 = getEdgeLength( a, (a+2)%3 );

	if( d1 < d2 ) {
		setLengthOn( a, (a+1)%3, d1 );
	} else {
		setLengthOn( a, (a+2)%3, d2 );
	}
}

void symmetro::setMult( const int &a, const int &m )
{
	mult[a] = m;
}

void symmetro::setScale( const int &a, const double &s )
{
	scale[a] = s;
}

// not used in symshell
void symmetro::setTheta( const int &a, const double &t )
{
	theta[a] = t;
}

void symmetro::setTheta( const int &a, const bool &b )
{
	if( !b ) {
		theta[a] = 0.0;
	} else {
		int order = getN( a );

		if( order > 0 ) {
			theta[a] = M_PI / double(order);
		}
	}
}

// Get the rotation that maps the positive Z axis to the given pole.
mat3d symmetro::toPole( const vec3d &p )
{
    double x = p[0];
    double y = p[1];
    double z = p[2];

    double d = sqrt( y*y + z*z );

    mat3d Ryi(
        vec3d(d, 0.0, x),
        vec3d(0.0, 1.0, 0.0),
        vec3d(-x, 0.0, d ));

    if( d > 0.0 ) {
        mat3d Rxi(
            vec3d(1.0, 0.0, 0.0),
            vec3d(0.0, z/d, y/d),
            vec3d(0.0, -y/d, z/d ));

        return Rxi * Ryi;
    } else {
        return Ryi;
    }
}

// Get the rotation that maps the given pole to the positive Z axis.
mat3d symmetro::fromPole( const vec3d &p )
{
    double x = p[0];
    double y = p[1];
    double z = p[2];

    double d = sqrt( y*y + z*z );

    mat3d Ry(
        vec3d(d, 0.0, -x),
        vec3d(0.0, 1.0, 0.0),
        vec3d(x, 0.0, d ));

    if( d > 0.0 ) {
        mat3d Rx(
            vec3d(1.0, 0.0, 0.0),
            vec3d(0.0, z/d, -y/d),
            vec3d(0.0, y/d, z/d ));

        return Ry * Rx;
    } else {
        return Ry;
    }
}

mat3d symmetro::rotateAboutPole( const vec3d &p, double theta )
{
	double ct = cos( theta );
	double st = sin( theta );

	mat3d rotatez( 
		vec3d(ct, -st, 0.0),
		vec3d(st, ct, 0.0),
		vec3d(0.0, 0.0, 1.0 ));

	return toPole( p ) * rotatez * fromPole( p );
}

vec3d symmetro::moveZ( const double &dist )
{
	return vec3d( sin( dist ), 0.0, cos( dist ) );
}

vec3d symmetro::moveZ( const double &dist, const double &theta )
{
	double s = sin( dist );
	return vec3d( s * cos( theta ), s * sin( theta ), cos( dist ) );
}

void symmetro::ssym_repl( const int &p, const int &q )
{
	double pp = M_PI / double( p );
	double pq = M_PI / double( q );
	
	// Fund[] is required for the algorithm
	// RK - fund 1 and 2 were reverse so we can have p,q,2
	fund[0] = vec3d( 0.0, 0.0, 1.0 );
	fund[1] = rotateAboutPole( vec3d(0.0,0.0,1.0), pp ) * moveZ( acos( 1.0 / (tan( pp ) * tan( pq )) ) );
	fund[2] = moveZ( acos( cos( pq ) / sin( pp ) ) );
	
//for( int i=0; i<3; i++ )
//   fprintf(stderr,"%.17lf %.17lf %.17lf\n", fund[i][0], fund[i][1], fund[i][2]);
//fprintf(stderr,"\n");
}

// RK - match() distills to this because we only do the first three polygons
mat3d symmetro::match( const int &idx, const vec3d &pseg, const bool &reverse )
{
   double ang = 0.0;
   if ( ( p == 5 && q == 3 ) || ( p == 3 && q == 5 ) ) {
      if (idx == 0 )
         ang = ( !reverse ) ? 18.0 : 120.0;
      else
      if (idx == 1 )
         ang = ( !reverse ) ? 60.0 : 54.0;
      else
      if (idx == 2 )
         ang = ( !reverse ) ? 90.0 : 0.0;
   }
   else 
   if ( ( p == 4 && q == 3 ) || ( p == 3 && q == 4 ) ) {
      if (idx == 0 )
         ang = ( !reverse ) ? 45.0 : 120.0;
      else
      if (idx == 1 )
         ang = ( !reverse ) ? 60.0 : 90.0;
      else
      if (idx == 2 )
         ang = ( !reverse ) ? 90.0 : 0.0;
   } 
   else 
   if ( p == 3 && q == 3 ) {
      if (idx == 0 )
         ang = ( !reverse ) ? 180.0 : 0.0;
      else
      if (idx == 1 )
         ang = ( !reverse ) ? 60.0 : 0.0;
      else
      if (idx == 2 )
         ang = ( !reverse ) ? 45.0 : -45.0;
   }

	double c = cos( deg2rad(ang) );
	double s = sin( deg2rad(ang) );

	return toPole( pseg ) * mat3d( vec3d(c, -s, 0.0), 
	                               vec3d(s, c, 0.0),
	                               vec3d(0.0, 0.0, 1.0) );
}

// distilled from constructFaces()
// RK - added non-zero scale qualifier for making polygon
// RK - added reverse so we can also reverse {3,3}
// find the first 3 polygons. Antiprism takes it from there...
vector<col_geom_v> symmetro::makePolygons( const bool &reverse )
{
   vector<col_geom_v> pgeom(3);
   
   for( int i=0; i<3; i++ ) {
	   int n = getN( i );
	   // RK - added scale check. don't make polygon if scale is zero or nan
	   if( n > 0 && ( scale[i] != 0.0 && !isnan(scale[i]) ) ) {
	      // RK - tried to transform built in polygon but couldn't figure out the transforms
		   // build polygon here works ok and is accurate
		   mat3d T = match( i, sym_vec[i], reverse );
         for( int idx = 0; idx < n; ++idx ) {
	         double ang = theta[i] - 2.0*M_PI*(double)idx/(double)n;
	         pgeom[i].add_vert( T * moveZ( scale[i], ang ) );
         }
         vector<int> face;
         for( int idx = 0; idx < n; ++idx )
            face.push_back( idx );
         pgeom[i].add_face( face );
      }
      
      // epsilon size faces are because scale was set at 0
      if ( fabs( scale[i] ) < epsilon ) {
         pgeom[i].clear_all();
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

bool detect_collision( col_geom_v &geom )
{
   const vector<vector<int> > &faces = geom.faces();
   const vector<vec3d> &verts = geom.verts();
   
   for( int i=0; i<(int)faces.size(); i++) {
      vector<int> face0 = faces[i];
      for( int j=i; j<(int)faces.size(); j++) {
         vector<int> face1 = faces[j];
         vec3d P, dir;
         if ( two_plane_intersect(  centroid(verts, face0), face_norm(verts, face0),
                                    centroid(verts, face1), face_norm(verts, face1),
                                    P, dir, epsilon ) ) {
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

   char errmsg[MSG_SZ];
 
   // some pre-checking  
   bool scale_zero = false;
   int num_multipliers = 0;
   int total_multiples = 0;
   vector<int> axes;
   for( int i=0; i<(int)opts.multipliers.size(); i++ ) {         
      if ( opts.multipliers[i]>0 )
         axes.push_back(i);
      total_multiples += opts.multipliers[i];
   }
   num_multipliers = (int)axes.size();
   
   if ( num_multipliers == 0 )
      opts.error("no multipliers specified",'m');
      
   // theta must exist
   if ( (int)opts.theta.size() == 0 )
      for( int i=0; i<3; i++ )
         opts.theta.push_back(false);

   // control reverse when generating only one polygon on axis 0 and 1
   if ( ( num_multipliers == 1 ) && ( total_multiples > 1 ) ) {
      if ( opts.multipliers[ 2 ] && opts.sym == 2 && ( !(opts.theta[ 2 ] || opts.rotation_method == 'v') ) ) {
         opts.reverse = is_even( opts.multipliers[ 2 ] ) ? true : false;
      }
      else 
      if ( !opts.multipliers[ 2 ] ) {
         opts.reverse = false;
         if ( opts.multipliers[0] )
            opts.reverse = true;
      }
   }

   if ( fabs( opts.scale ) > 0.0 && fabs( opts.scale ) < epsilon )
      scale_zero = true;

   // legacy code as object
   symmetro s;
   
   if (opts.sym == 1) {
      if (!opts.reverse)
         s.setSym( 5, 3 );
      else
         s.setSym( 3, 5 );
   }
   else
   if (opts.sym == 2)
      if (!opts.reverse)
         s.setSym( 4, 3 );
      else
         s.setSym( 3, 4 );
   else
   if (opts.sym == 3)
      s.setSym( 3, 3 );
   else
      opts.error("symmetry type not set",'s');
   
   bool zero_found = false;   
   for( int i=0; i<(int)opts.multipliers.size(); i++ ) {
      if ( opts.multipliers[i] == 0 )
         zero_found = true;
      s.setMult( i, opts.multipliers[i] );
   }
   if ( !zero_found )
      opts.warning("when all three multipliers are used, polygons will not be of equal edge length",'m');
   
   // for edge on model or vertex on model
   if ( opts.rotation_method ) {
      if ( num_multipliers == 1 ) {
         int place = opts.multipliers[ 0 ] ? 0 : ( opts.multipliers[ 1 ] ? 1 : 2 );
         if (opts.rotation_method == 'e' )
            opts.theta[ place ] = true;
         // patch for octahedral, is reversed when even
         if ( opts.sym == 2 && is_even( opts.multipliers[ 2 ] ) )
            opts.theta[ 2 ] = !opts.theta[ 2 ];
      }
      else
      if ( num_multipliers > 1 ) {
         bool edge_on0 = s.isEdgeOn( axes[0], axes[1] );
         bool edge_on1 = s.isEdgeOn( axes[1], axes[0] );
         if ( ( opts.rotation_method == 'e' && !edge_on0 ) || ( opts.rotation_method == 'v' && edge_on0 ) )
            opts.theta[ axes[0] ] = true;
         if ( ( opts.rotation_method == 'e' && !edge_on1 ) || ( opts.rotation_method == 'v' && edge_on1 ) )
            opts.theta[ axes[1] ] = true;
         // only when 2-fold axis is used
         if ( num_multipliers == 3 ) {
            bool edge_on2 = s.isEdgeOn( axes[0], axes[2] );
            if ( ( opts.rotation_method == 'e' && !edge_on2 ) || ( opts.rotation_method == 'v' && edge_on2 ) )
               opts.theta[ axes[2] ] = true;
         }
      }
   }   
   for( int i=0; i<(int)opts.theta.size(); i++ )
      s.setTheta( i, opts.theta[i] );
   // RK - I think if used as a third multiplier, the 2-fold axis has greater than a square and edge on, this is alway true
   if ( !opts.patch && ( num_multipliers == 3 ) && ( opts.multipliers[2] * 2 > 4 ) ) {
      bool edge_on2 = s.isEdgeOn( axes[0], axes[2] ) && s.isEdgeOn( axes[2], axes[0] );
      if ( edge_on2 && opts.convex_hull > 3 ) {
         opts.warning("when axis 2 multiplier is greater than 2, model will not close");
         opts.warning("convex hull is suppressed");
         opts.convex_hull = 2;
      }
   }
   
   // can't set patch until here
   if ( opts.patch )
      s.setPatch( opts.patch );
      
   // do tie completely automatic based on multipliers
   vector<int> auto_tie;
   for( int i=0; i<(int)opts.multipliers.size(); i++ )
      if ( opts.multipliers[i]>0 )
         auto_tie.push_back(i);
         
   if( (int)auto_tie.size() == 3 ) {
      if ( !s.tie3( 0, 1, 2, errmsg ) ) {
         opts.warning(errmsg);
         if ( opts.convex_hull > 3 ) {
            opts.warning("convex hull is suppressed");
            opts.convex_hull = 2;
         }
      }
   }
   else
   if( (int)auto_tie.size() == 2 ) {
      if ( !s.tie( auto_tie[0], auto_tie[1] ) ) {
         // RK - catch if tie(a,b) failed. Try to tie all three
         s.tie3( 0, 1, 2, errmsg ); // will be in error. just for unit vector solution
         if ( !scale_zero ) { // silence warnings
            opts.warning("tie(a,b) FAILED because edge meets a vertex. Trying tie(a,b,c)");
            if ( opts.convex_hull > 3 ) {
               opts.warning("convex hull is suppressed");
               opts.convex_hull = 2;
            }
         }
      }
   }
   else
   if( (int)auto_tie.size() == 1 ) {
      s.tie1();
   }
   else
      opts.error("automatic tie did not work"); // probably can't get here
   
   // does a propeller like operation
   for( int i=0; i<(int)opts.extra_rotation.size(); i++ )
      if ( s.getScale( i ) )
         s.setTheta( i, s.getTheta( i ) + deg2rad(opts.extra_rotation[i]) );
   
   // calculate axes here      
   s.fill_sym_vec( opts.reverse );
   
   // change ratio
   if ( opts.scale ) {
      // edge scale math furnished by Adrian Rossiter
      double angle_between_axes = s.getAngleBetweenAxes( opts.scale_direction[0], opts.scale_direction[1] );
      if ( opts.verbose )
         fprintf(stderr,"angle_between_axes = %.17lf\n",rad2deg(angle_between_axes));
      
      double p_ang0 = (2.0*M_PI/s.getN(opts.scale_direction[0]))/2.0;
      double r0 = opts.scale*0.5/sin(p_ang0);
      double p_ang1 = (2.0*M_PI/s.getN(opts.scale_direction[1]))/2.0;
      double r1 = 0.5/sin(p_ang1);
      double d = sqrt(r0*r0 + r1*r1 + 2.0*r0*r1*cos(angle_between_axes));
      double a0 = acos((r1*r1 + d*d - r0*r0)/(2.0*r1*d));
      double a1 = acos((r0*r0 + d*d - r1*r1)/(2.0*r0*d));
   
      if ( (int)opts.scale_direction.size() == 0 )
         opts.error("ratio direction not set",'d');
      for( int i=0; i<(int)opts.scale_direction.size(); i++ ) {
         if ( s.getScale( opts.scale_direction[i] ) == 0.0 )
            opts.error(msg_str("scale of axis polygon '%d' is zero and cannot be used for scaling", opts.scale_direction[i]), 'd');
      }
      if ( s.isEdgeOn( opts.scale_direction[0], opts.scale_direction[1] ) )
         opts.error(msg_str("polygon '%d' and '%d' are not vertex connected", opts.scale_direction[0], opts.scale_direction[1]), 'd');
         
      s.setScale( opts.scale_direction[0], a0 );
      s.setScale( opts.scale_direction[1], a1 );
      
      // old method only varies face center radii
      //s.setScale( opts.scale_direction[0], s.getScale( opts.scale_direction[0] ) * opts.scale  );
      //s.tieTo( opts.scale_direction[1], opts.scale_direction[0] );

      // scale third axis if present
      if ( num_multipliers == 3 ) {
         int third_axis = -1;
         int sum = opts.scale_direction[0] + opts.scale_direction[1];
         if ( sum == 1 ) // 0,1
            third_axis = 2;
         else
         if ( sum == 2 ) // 0,2
            third_axis = 1;
         else
         if ( sum == 3 ) // 1,2
            third_axis = 0;
         
         if ( s.isEdgeOn( third_axis, opts.scale_direction[1] ) )
            opts.error(msg_str("polygon '%d' and '%d' are not vertex connected", third_axis, opts.scale_direction[1]), 'd');
         
         s.tieTo( third_axis, opts.scale_direction[1] );
      }
   }
   
   vector<col_geom_v> pgeom = s.makePolygons( opts.reverse );
   
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

