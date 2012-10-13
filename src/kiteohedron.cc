/*
   Copyright (c) 2012, Adrian Rossiter

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
   Name: kitohedron.cc
   Description: make a kite faced polyhedron
   Project: Antiprism - http://www.antiprism.com
*/


#include "../base/antiprism.h"

using std::string;
using std::vector;

int read_fraction(char *frac_str, int &num, int &denom, char *errmsg)
{
   if(errmsg)
      *errmsg = '\0';

   char errmsg2[MSG_SZ];
   denom=1;
   char *p = strchr(frac_str, '/');
   if(p!=0) {
      *p++='\0';
      if(!read_int(p, &denom, errmsg2)) {
         sprintf(errmsg, "denominator, %s", errmsg2);
         return false;
      }
   }
   if(!read_int(frac_str, &num, errmsg2)) {
      sprintf(errmsg, "numerator, %s", errmsg2);
      return false;
   }
   return true;
}

int read_triangle_fraction(char *frac_str, int &num, int &denom, char *errmsg)
{
   if(!read_fraction(frac_str, num, denom, errmsg))
      return false;

   if(num<2) {
      sprintf(errmsg, "numerator must be 2 or greater");
      return false;
   }
   if(denom==0) {
      sprintf(errmsg, "denominator cannot be 0");
      return false;
   }
   if(denom>=num) {
      sprintf(errmsg, "denominator must be less than the numerator");
      return false;
   }
   return true;
}


class kt_opts: public prog_opts
{
   private:

   public:
      vector<int> triangle;
      char height_type;
      double height;
      bool kite_only;
      int color_by_vertex;

      string ofile;

      kt_opts():  prog_opts("kiteohedron"),
                  height_type('x'),
                  kite_only(false),
                  color_by_vertex(-1)
                  {}

      void process_command_line(int argc, char **argv);
      void usage();
};

void kt_opts::usage()
{
   fprintf(stdout,
"\n"
"Usage: %s [options] frac_A frac_B [frac_C]\n"
"\n"
"The three fractions A, B, C (default 2) specifying the Schwarz triangle\n"
"vertices, where n/d corresponds to an angle PI*d/n. The kite has vertices\n"
"along A, C, B, and C reflected in AB. The kite is repeated with the\n"
"symmetry corresponding to the base Schwarz triangle.\n"
"\n"
"Options\n"
"  -h        this help message\n"
"  -l ht     height of kite vertex on OA (cannot use with -p))\n"
"  -p ht     height of pivot, midpoint of C and C reflected in AB (cannot\n"
"            use with -l))\n"
"  -k        output a single kite (colours not applied)\n"
"  -c type   colour the faces around each vertex of a type, from AaBbCc\n"
"            (colouring by value/index for upper/lower case) using a\n"
"            different colour for each set\n"
"  -o <file> write output to file (default: write to standard output)\n"
"\n"
"\n", prog_name());
}


void kt_opts::process_command_line(int argc, char **argv)
{
   char errmsg[MSG_SZ];
   opterr = 0;
   char c;

   handle_long_opts(argc, argv);

   while( (c = getopt(argc, argv, ":hl:p:c:ko:")) != -1 ) {
      if(common_opts(c, optopt))
         continue;

      switch(c) {
         case 'l':
            if(height_type!='x' && height_type!='l')
               error("option -l: cannot be used with option -p");
            if(!read_double(optarg, &height, errmsg))
               error(errmsg, c);
            height_type = 'l';
            break;

         case 'p':
            if(height_type!='x' && height_type!='p')
               error("option -p: cannot be used with option -l");
            if(!read_double(optarg, &height, errmsg))
               error(errmsg, c);
            height_type = 'p';
            break;

         case 'k':
            kite_only = true;;
            break;

         case 'c':
            if(strlen(optarg)>1)
               error(msg_str("colour letter has extra letters"));
            else {
               const char *v_letters = "AaBbCc";
               const char *p = strchr(v_letters, *optarg);
               if(p)
                  color_by_vertex = p - v_letters;
               else
                  error(msg_str("colour letter is '%c', must be from AaBbCc",
                        *optarg));
            }
            break;

         case 'o':
            ofile = optarg;
            break;

         default:
            error("unknown command line error");
      }
   }

   triangle.resize(6);
   int num_args = argc-optind;
   if(argc-optind > 3)
      error(msg_str("should give 2 or 3 fractions (%d given)", num_args));

   if(!read_triangle_fraction(argv[optind], triangle[0], triangle[1], errmsg))
      error(msg_str("fraction A is '%s': %s", argv[optind], errmsg));
   optind++;

   if(!read_triangle_fraction(argv[optind], triangle[2], triangle[3], errmsg))
      error(msg_str("fraction B is '%s': %s", argv[optind], errmsg));
   optind++;

   if(num_args==2) {
      triangle[4] = 2;
      triangle[5] = 1;
   }
   else if(!read_triangle_fraction(argv[optind], triangle[4], triangle[5],
            errmsg))
      error(msg_str("fraction C is '%s': %s", argv[optind], errmsg));


}


bool triangle_to_kite(const vector<vec3d> &tri, vector<vec3d> &kite,
      char ht_type, double ht)
{
   kite.resize(4);
   vec3d A  = tri[0].unit();
   vec3d B  = tri[1].unit();
   vec3d n  = vcross(A, B);    // angle AOB is not 0 or PI in Schwarz triangle
   vec3d P0 = tri[2];
   if(ht_type=='p')
      P0 *= ht;
   vec3d P1 = mat3d::refl(n) * P0;

   if(ht_type=='l')
      A *= ht;

   vec3d pivot = (P0 + P1)/2;

   B = lines_intersection(A, pivot, vec3d(0,0,0), B);
   for(int i=0; i<3; i++)
      if(isnan(B[i]))
         return false;

   kite[0] = A;
   kite[1] = P0;
   kite[2] = B;
   kite[3] = P1;

   return true;
}

void repeat_kite_with_color(col_geom_v &out_geom, col_geom_v kite_geom,
      sch_sym sym, const vector<int> &fracs, int col_vert_type)
{
   int vert_no = col_vert_type/2;
   int v_map[] = {0, 2, 1};  // B is at 3 and C is at 2
   vec3d axis = kite_geom.verts(v_map[vert_no]);
   sch_sym c_sym;
   c_sym.init(sch_sym::C, fracs[2*vert_no], mat3d::rot(axis, vec3d::Z));

   col_geom_v vert_kites_geom;
   sym_repeat(vert_kites_geom, kite_geom, c_sym);
   coloring clrngs[3];
   if(col_vert_type%2==0)     // value letter is followed by index letter
      clrngs[2].add_cmap(init_color_map("spread"));

   t_set min_ts;
   min_ts.min_set(sym.get_trans(), c_sym.get_trans());

   col_geom_v comp_geom;
   sym_repeat(out_geom, vert_kites_geom, min_ts, ELEM_FACES, clrngs);
}




int main(int argc, char *argv[])
{
   kt_opts opts;
   opts.process_command_line(argc, argv);

   char errmsg[MSG_SZ];
   vector<vec3d> schwarz_verts;
   sch_sym sym;
   if(!get_schwarz_tri_verts(opts.triangle, schwarz_verts, &sym))
      opts.error("fractions do not correspond to a Schwarz triangle");

   col_geom_v kite_geom;
   if(!triangle_to_kite(schwarz_verts, kite_geom.raw_verts(),
         opts.height_type, opts.height))
      opts.error("could not calculate model with this height");
   kite_geom.add_face(0, 1, 2, 3, -1);

   col_geom_v out_geom;
   if(!opts.kite_only) {
      if(opts.color_by_vertex<0)
         sym_repeat(out_geom, kite_geom, sym);
      else
         repeat_kite_with_color(out_geom, kite_geom, sym,
               opts.triangle, opts.color_by_vertex);
      sort_merge_elems(out_geom, "v", epsilon);
   }
   else
      out_geom = kite_geom;

   if(!out_geom.write(opts.ofile, errmsg))
      opts.error(errmsg);

   return 0;
}


