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
   Name: off_trans.cc
   Description: transformations for OFF files
   Project: Antiprism - http://www.antiprism.com
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <string>
#include <vector>

#include "../base/antiprism.h"

using std::string;
using std::vector;
using std::map;


class trans_opts : public prog_opts {
   public:
      mat3d trans_m;
      col_geom_v geom;
      
      string ifile;
      string ofile;

      trans_opts(): prog_opts("off_trans") {}
      void process_command_line(int argc, char **argv);
      void usage();
};


void trans_opts::usage()
{
   fprintf(stdout,
"\n"
"Usage: %s [options] [input_file]\n"
"\n"
"Read a file in OFF format and apply transformations to it. If\n"
"input_file is not given the program reads from standard input.\n"
"\n"
"Options\n"
"%s"
"  -T <tran> translate, three numbers separated by commas which are\n"
"            used as the x, y and z displacements\n"
"  -R <rot>  rotate about an axis, three, four or six numbers separated by\n"
"            commas. If three numbers these are angles (degrees) to rotate\n"
"            about the x, y and z axes. If four numbers, the first three\n"
"            are a direction vector for the axis, the last number is the\n"
"            angle (degrees) to rotate. If six numbers, these are two\n"
"            vectors (from,to) and rotate to carry the first to the second.\n"
"            If twelve numbers these are four vectors (from1,from2,to1,to2)\n"
"            and rotate to carry the first onto the third then rotate around\n"
"            the third to carry the second onto the fourth\n"
"  -M <norm> reflect in a plane, three numbers separated by commas which\n"
"            give a vector normal to the plane of reflection.\n"
"  -S <scal> scale, one, three or four numbers separated by commas. If one\n"
"            number then scale by this factor in all directions. If three\n"
"            numbers these are the factors to scale along the x, y and\n"
"            z axes. If four numbers, the first three are a direction\n"
"            vector for the scaling, the last number is the factor to scale\n"
"  -I        inversion\n"
"  -A <crds> transformation that will align two sets of three points\n"
"            (18 numbers coordinates of from1,from2,from3,to1,to2,to3)\n"
"  -a <angs> transformation that makes particular angles between the\n"
"            mapped axes, angles in degrees in form yz_ang,zx_ang,xy_ang\n"
"            (corresponding to the angles opposite the x-, y- and z-axis)\n"
"  -X <mtrx> transformation matrix of 9 or 12 values, given left-to-right\n"
"            top-to-bottom, used to premultipy each coordinate\n"
"  -C        translation that carries the centroid to the origin\n"
"  -y        align geometry with the standard alignment for a symmetry type,\n"
"            up to three comma separated parts: symmetry subgroup (Schoenflies\n"
"            notation) or 'full', conjugation type (integer), realignment\n"
"            (colon separated list of an integer then decimal numbers)\n"
"  -Y        align standard alignment of a symmetry type with its position\n"
"            as a subgroup of another symmetry type, up to four comma\n"
"            separated parts: containing symmetry (Schoenflies notation) or\n"
"            file name, subgroup *Schoenflies notation), conjugation type\n"
"            (integer), realignment (colon separated list of an integer\n"
"            then decimal numbers)\n"
"  -s <type> relative scaling, scale so a measure has a value of 1.\n"
"            VAa need an oriented polyhedron, V needs a closed polyhedron.\n"
"            A - area                       a - average face area\n"
"            E - perimeter (sum of edges)   e - average edge length\n"
"            V - volume                     r - radius, from centroid\n"
"                                               to furthest vertex\n"
"  -i        replace the current combined transformation by its inverse\n"
"  -o <file> write output to file (default: write to standard output)\n"
"\n"
"\n", prog_name(), help_ver_text);
}


bool rel_scale_val(geom_if &geom, char rel_scale, double *scale, char *errmsg)
{
   geom_info info(geom);
   *errmsg = '\0';
   if(strchr("VAa", rel_scale)) {
      if(!geom.get_faces()->size()) {
         strcpy(errmsg, "scaling type is V, A or a but there is no face data");
         return false;
      }
      if(!info.is_oriented())
         strcpy(errmsg,"scaling type is V, A or a and polyhedron is not"
               "oriented");
   }
   if(strchr("V", rel_scale) && !info.is_closed())
      strcpy(errmsg, "scaling type is V and polyhedron is not closed");
   if(strchr("Ee", rel_scale) && (!info.num_edges() && !info.num_iedges()) ) {
      strcpy(errmsg, "scaling type is E or e but there is no edge data");
      return false;
   }

   switch(rel_scale) {
      case 'V':
         *scale = pow(fabs(info.volume()), 1.0/3.0);
         break;
      case 'A':
         *scale = sqrt(info.face_areas().sum);
         break;
      case 'a':
         *scale = sqrt(fabs(info.face_areas().sum)/info.num_faces());
         break;
      case 'E':
         if(info.num_edges())  // use explicit edges, if present
            *scale = info.edge_lengths().sum;
         else
            *scale = info.iedge_lengths().sum;
         break;
      case 'e':
         if(info.num_edges())  // use explicit edges, if present
            *scale = info.edge_lengths().sum/info.num_edges();
         else
            *scale = info.iedge_lengths().sum/info.num_iedges();
         break;
      case 'r':
         info.set_center(centroid(*geom.get_verts()));
         *scale = info.vert_dists().max;
         break;
   }

   return true;
}


void trans_opts::process_command_line(int argc, char **argv)
{
   char errmsg[MSG_SZ];
   opterr = 0;
   char c;
   vector<pair<char, char *> > args;
   
   handle_long_opts(argc, argv);

   while( (c = getopt(argc, argv, ":hT:R:M:S:IX:A:a:CY:y:s:io:")) != -1 ) {
      if(common_opts(c, optopt))
         continue;

      switch(c) {  // Keep switch for consistency/maintainability
         default:
            args.push_back(pair<char, char *>(c, optarg));
      }
   }

   if(argc-optind > 1)
      error("too many arguments");
   
   if(argc-optind == 1)
      ifile=argv[optind];

   if(!geom.read(ifile, errmsg))
      error(errmsg);
   if(*errmsg)
      warning(errmsg);

   
   vector<double> nums;
   mat3d trans_m2;
   
   for(unsigned int i=0; i<args.size(); ++i) {
      c = args[i].first;
      char *optarg = args[i].second;
      switch(c) {
         case 'R':
            if(!read_double_list(optarg, nums, errmsg))
               error(errmsg, c);
            if(nums.size()==3)
               trans_m2 = mat3d::rot(deg2rad(nums[0]), deg2rad(nums[1]),
                     deg2rad(nums[2]));
            else if(nums.size()==4)
               trans_m2 = mat3d::rot(vec3d(nums[0], nums[1], nums[2]),
                     deg2rad(nums[3]));
            else if(nums.size()==6)
               trans_m2 = mat3d::rot(vec3d(nums[0], nums[1], nums[2]),
                     vec3d(nums[3], nums[4], nums[5]));
            else if(nums.size()==12)
               trans_m2 = mat3d::alignment(vec3d(nums[0], nums[1], nums[2]),
                                           vec3d(nums[3], nums[4], nums[5]),
                                           vec3d(nums[6], nums[7], nums[8]),
                                           vec3d(nums[9], nums[10], nums[11]));
            else
               error(msg_str( "must give 3, 4, 6 of 12 numbers (%lu were given)", (unsigned long)nums.size()), c);
            trans_m = trans_m2 * trans_m;
            break;

         case 'S':
            if(!read_double_list(optarg, nums, errmsg))
               error(errmsg, c);
            if(nums.size()==1)
               trans_m2=mat3d::scale(nums[0]);
            else if(nums.size()==3)
               trans_m2=mat3d::scale(nums[0], nums[1], nums[2]);
            else if(nums.size()==4)
               trans_m2=mat3d::scale(vec3d(nums[0], nums[1], nums[2]),nums[3]);
            else
               error(msg_str("must give 1, 3 or 4 numbers (%lu were given)",
                        (unsigned long)nums.size()), c);
            
            trans_m = trans_m2 * trans_m;
            break;

         case 'T':
            if(!read_double_list(optarg, nums, errmsg))
               error(errmsg, c);
            if(nums.size()!=3)
               error(msg_str("must give exactly three numbers (%lu were "
                        "given)" , (unsigned long)nums.size()), c);
            trans_m2 =mat3d::transl(vec3d(nums[0], nums[1], nums[2]));
            trans_m = trans_m2 * trans_m;
            break;

         case 'M':
            if(!read_double_list(optarg, nums, errmsg))
               error(errmsg, c);
            if(nums.size()!=3)
               error(msg_str("must give exactly three numbers (%lu were "
                        "given)", (unsigned long)nums.size()), c);
            
            trans_m2 = mat3d::refl(vec3d(nums[0], nums[1], nums[2]));
            trans_m = trans_m2 * trans_m;
            break;

         case 'I':
            trans_m = mat3d::inversion() * trans_m;
            break;

         case 'A':
            if(!read_double_list(optarg, nums, errmsg))
               error(errmsg, c);
            if(nums.size()==18)
               trans_m2 = mat3d::alignment(vec3d(nums[0], nums[1], nums[2]),
                                           vec3d(nums[3], nums[4], nums[5]),
                                           vec3d(nums[6], nums[7], nums[8]),
                                           vec3d(nums[9], nums[10], nums[11]),
                                           vec3d(nums[12], nums[13], nums[14]),
                                           vec3d(nums[15], nums[16], nums[17]));
            else
               error(msg_str("must give 18 numbers (%lu were given)",
                        (unsigned long)nums.size()), c);
            trans_m = trans_m2 * trans_m;
            break;

         case 'a':
            if(!read_double_list(optarg, nums, errmsg))
               error(errmsg, c);
            if(nums.size()!=3)
               error(msg_str("must give exactly three numbers (%lu were "
                        "given)", (unsigned long)nums.size()), c);
            bool valid;
            trans_m2 = mat3d::trans_by_angles(deg2rad(nums[0]),
                  deg2rad(nums[1]), deg2rad(nums[2]), &valid);
            if(!valid)
               error("the sum of any two angles must be greater than the third",                     c);

            trans_m = trans_m2 * trans_m;
            break;

         case 'X':
            if(!read_double_list(optarg, nums, errmsg))
               error(errmsg, c);
            trans_m2 = mat3d::unit();
            if(nums.size()==9) {
               trans_m2[0]=nums[0]; trans_m2[1]=nums[1]; trans_m2[2]=nums[2];
               trans_m2[4]=nums[3]; trans_m2[5]=nums[4]; trans_m2[6]=nums[5];
               trans_m2[8]=nums[6]; trans_m2[9]=nums[7]; trans_m2[10]=nums[8];
            }
            else if(nums.size()==12) {
               for(int i=0; i<12; i++)
                  trans_m2[i]=nums[i];
            }
            else
               error(msg_str("must give 9 or 12 numbers (%lu were given)",
                        (unsigned long)nums.size()), c);
            
            trans_m = trans_m2 * trans_m;
            break;

         case 'C':
            trans_m = mat3d::transl(-(trans_m*geom.centroid())) * trans_m;
            break;
         
         case 'Y':
         {
            vector<char *> parts;
            split_line(optarg, parts, ",");
            if(parts.size()<1 || parts.size()>4)
               error("argument should have 1-4 comma separated parts", c);
           
            sch_sym sym;
            col_geom_v sgeom;
            if(sgeom.read(parts[0], errmsg))
               sym.init(sgeom);
            else if(!sym.init(parts[0], mat3d(), errmsg))
               error(msg_str("invalid filename or symmetry type name: %s",
                        errmsg), c);
            
            sch_sym sub_sym = sym;
            if(parts.size()>1 && !sub_sym.init(parts[1], mat3d(), errmsg))
               error(msg_str("sub-symmetry type: %s", errmsg), c);
               
            int sub_sym_conj = 0;
            if(parts.size()>2) {
               if(!read_int(parts[2], &sub_sym_conj, errmsg))
                  error(msg_str("sub-symmetry conjugation number: %s",
                           errmsg), c);
            }
 
            sch_sym final_sym = sym.get_sub_sym(sub_sym, sub_sym_conj, errmsg);
            if(final_sym.get_sym_type() == sch_sym::unknown)
               error(msg_str("sub-symmetry: %s", errmsg), c);

            if(parts.size()>3) {
               if(!final_sym.get_autos().set_realignment(parts[3], errmsg))
                  error(msg_str("sub-symmetry realignment: %s",
                           errmsg), c);
            }

            trans_m = mat3d::inverse(final_sym.get_autos().get_realignment() *
                  final_sym.get_to_std()) * trans_m;
            break;
         }

         case 'y':
         {
            geom_v geom_cur = geom;
            geom_cur.transform(trans_m);
            sch_sym full_sym(geom_cur);

            vector<char *> parts;
            split_line(optarg, parts, ",");
            if(parts.size()==0 || parts.size()>3)
               error("argument should have 1-3 comma separated parts", c);
           
            sch_sym sub_sym;
            if(strncmp(parts[0], "full", strlen(optarg))==0) {
               sub_sym = full_sym;
            }
            else if(!sub_sym.init(parts[0], mat3d(), errmsg))
               error(msg_str("sub-symmetry type: %s", errmsg), c);
               
            int sub_sym_conj = 0;
            if(parts.size()>1) {
               if(!read_int(parts[1], &sub_sym_conj, errmsg))
                  error(msg_str("sub-symmetry conjugation number: %s",
                           errmsg), c);
            }

            sch_sym sym = full_sym.get_sub_sym(sub_sym, sub_sym_conj, errmsg);
            if(sym.get_sym_type() == sch_sym::unknown)
               error(msg_str("sub-symmetry: %s", errmsg), c);

            if(parts.size()>2) {
               if(!sym.get_autos().set_realignment(parts[2], errmsg))
                  error(msg_str("sub-symmetry realignment: %s",
                           errmsg), c);
            }

            trans_m = sym.get_autos().get_realignment() *
                      sym.get_to_std() * trans_m;
            break;
         }

         case 's': {
            char rel_scale = '\0';
            if(strlen(optarg)==1 && strchr("VAaEer", *optarg))
               rel_scale=*optarg;
            else
               error("relative scale must be V,A,a,E,e or r", c);
            geom_v geom_cur = geom;
            geom_cur.transform(trans_m);
            double scale;
            if(!rel_scale_val(geom_cur, rel_scale, &scale, errmsg))
               error(errmsg, 's');
            if(*errmsg)
               warning(errmsg, 's');
            trans_m = trans_m * mat3d::scale(1/scale);
            break;
         }

         case 'i':
            trans_m.set_inverse();
            break;

         case 'o':
            ofile = optarg;
            break;

         default:
            error("unknown command line error");
      }
   }

}


int main(int argc, char *argv[])
{
   char errmsg[MSG_SZ];
   trans_opts opts;
   opts.process_command_line(argc, argv);
    
   opts.geom.transform(opts.trans_m);

   if(!opts.geom.write(opts.ofile, errmsg))
      opts.error(errmsg);

    return 0;
}


