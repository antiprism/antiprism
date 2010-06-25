/*
   Copyright (c) 2008-2009, Roger Kaufman

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
   Name: iso_delta.cc
   Description: Isohedral Deltahedra
                Based on a paper by G. C. Shephard
                Periodica Mathematica Hungarica, Volume 39, Numbers 1-3, 2000 , pp. 83-106(24).
                with enhancements by Jim McNeill (http://www.orchidpalms.com/polyhedra)
                and Adrian Rossiter (http://www.antiprism.com)
   Project: Antiprism - http://www.antiprism.com
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <ctype.h>
#include <unistd.h>

#include <string>
#include <vector>

#include "../base/antiprism.h"
#include "../base/transforms.h"

using std::string;
using std::vector;


struct IsoDeltaItem {
   const char *name;
   const char *sym_type;
   const char *symbol;
	const char *comment;
};

struct IsoDeltaVector {
   const char *name;
   vec3d A;
   vec3d B;
   vec3d C;
   double alpha;
   double beta;
   double gamma;
};

IsoDeltaItem iso_delta_item_list[] = {
   {"T1(1)",   "Oh", "[1/3,2/3,1/3]", "Augmented Cube (=O5(1))"},
   {"T1(2)",   "Oh", "[2/3,2/3,2/3]", "Excavated Cube (=O5(2))"},
   {"T1(3)",   "Td", "[1/3,2/3,1/3]", "Additonal Isomer of T1(1) (Mobius)"},
   {"T2(1)",   "Td", "[1/3,1/3,2/3]", "Augmented Tetrahedron"},
   {"T2(2)",   "Td", "[2/3,1/3,1/3]", "Tetrahedron"},
   {"O1",      "Oh", "[1/2,1/2,2/3]", "COMPOUND of 3 Octahedra (UD08)"},
   {"O2(1)",   "Oh", "[1/2,1/3,1/4]", "Relaxed Excavated Rhombic Dodecahedron (Mobius)"},
   {"O2(2)",   "Oh", "[1/2,1/3,1/4]", "Relaxed Augmented Rhombic Dodecahedron (Mobius)"},
   {"O3",      "Oh", "[1/2,1/4,2/4]", "COMPOUND of 3 8/3 Star Bipyramids"},
   {"O4(1)",   "Oh", "[1/3,1/3,2/3]", "COMPOUND of 2 Augmented Tetrahedra T2(1) (UC54_d)"},
   {"O4(2)",   "Oh", "[2/3,1/3,1/3]", "COMPOUND of 2 Tetrahedra (Stella Octangula UC04)"},
   {"O5(1)",   "Oh", "[1/3,2/3,1/3]", "Augmented Cube (=T1(1))"},
   {"O5(2)",   "Oh", "[2/3,2/3,2/3]", "Excavated Cube (=T1(2))"},
   {"O5(3)",   "Oh", "[1/3,2/3,1/3]", "COMPOUND of 2 T1(3)"},
   {"O6(1)",   "Oh", "[2/3,1/4,1/4]", "Augmented Octahedron"},
   {"O6(2)",   "Oh", "[2/3,3/4,3/4]", "Excavated Octahedron"},
   {"O7",      "Oh", "[2/4,2/4,2/4]", "Octahedron"},
   {"I1(1)",   "Ih", "[1/2,1/3,1/5]", "Relaxed Excavated Rhombic Triacontahedron (Mobius)"},
   {"I1(2)",   "Ih", "[1/2,1/3,1/5]", "Relaxed Augmented Rhombic Triacontahedron (Mobius)"},
   {"I2(1)",   "Ih", "[1/2,2/5,4/5]", "Relaxed Augmented Did (U36) *"},
   {"I2(2)",   "Ih", "[1/2,3/5,4/5]", "Relaxed Excavated Did (U36) *"},
   {"I3(1)",   "Ih", "[1/2,1/3,2/5]", "Relaxed Augmented Gid (U54) *"},
   {"I3(2)",   "Ih", "[1/2,1/3,2/5]", "Relaxed Excavated Gid (U54) *"},
   {"I3(3)",   "Ih", "[1/2,2/3,3/5]", "Relaxed Dual of Gaquatid (UD73) *"},
   {"I3(4)",   "Ih", "[1/2,2/3,2/5]", "Additional Isomer of Augmented RT I1(2) **"},
   {"I4",      "Ih", "[1/2,1/2,1/2]", "COMPOUND of 5 Octahedra (UC17)"},
   {"I5(1)",   "Ih", "[2/3,3/5,4/5]", "Relaxed Augmented Dual of Ditdid (UD41)"},
   {"I5(2)",   "Ih", "[1/3,2/5,4/5]", "Relaxed Excavated Dual of Ditdid (UD41)"},
   {"I6(1)",   "Ih", "[2/3,1/5,1/5]", "Augmented Icosahedron"},
   {"I6(2)",   "Ih", "[2/3,4/5,4/5]", "Excavated Icosahedron"},
   {"I7(1)",   "Ih", "[1/3,2/5,3/5]", "Excavated Gike (U53)"},
   {"I7(2)",   "Ih", "[2/3,2/5,3/5]", "Augmented Gike (U53)"},
   {"I7(3)",   "Ih", "[1/3,2/5,3/5]", "Additional Isomer of Augmented Icosahedron I6(1)"},
   {"I8(1)",   "Ih", "[1/3,1/3,4/5]", "Augmented Gissid (U52)"},
   {"I8(2)",   "Ih", "[1/3,1/3,4/5]", "Excavated Gissid (U52)"},
   {"I9(1)",   "Ih", "[1/3,1/3,2/5]", "Augmented Dodecahedron"},
   {"I9(2)",   "Ih", "[1/3,1/3,2/5]", "Excavated Dodecahedron"},
   {"I9(3)",   "Ih", "[1/3,1/3,2/5]", "Relaxed Excavated Dual of Sidtid (UD30)"},
   // Redundant 
   //{"I9(4)",   "Ih", "[1/3,1/3,2/5]", "Relaxed Excavated Dual of Sidtid (UD30) (Repeat)"},
   {"I10(1)",  "Ih", "[4/5,4/5,4/5]", "Great Icosahedron or Gike (U53)"},
   {"I10(2)",  "Ih", "[1/5,1/5,4/5]", "Augmented Sissid (U34)"},
   // Redundant 
   //{"I10(3)",  "Ih", "[1/5,1/5,4/5]", "Augmented Sissid (U34) (Repeat)"},
   {"I11(1)",  "Ih", "[2/5,2/5,2/5]", "Icosahedron"},
   {"I11(2)",  "Ih", "[2/5,3/5,3/5]", "Excavated Gad (U35)"},
};

#define phi ((1+sqrt(5))/2)
#define iphi (1/phi)
#define phi2 (phi*phi)
#define term1 (sqrt(5.0/2+sqrt(5)/2))
// = 1.90211303259030714... also sqrt((10+sqrt(20))/4) or sqrt(phi)*pow(5, 0.25)
#define term2 sqrt(5.0/8+sqrt(5)/8)
// = 0.95105651629515357... also sqrt((10+sqrt(20))/16) or 0.5*sqrt(0.5*(5+sqrt(5)))
#define term3 sqrt(5.0/8-sqrt(5)/8)
// = 0.58778525229247312... also sqrt((10-sqrt(20))/16) or 0.5*sqrt(0.5*(5-sqrt(5)))

// notes added as to what changed from what was listed in the paper
IsoDeltaVector iso_delta_vector_list[] = {
   {"T1(1)",   vec3d(1,-1,1)/sqrt(3), vec3d(1,0,0), vec3d(1,-1,-1)/sqrt(3),
               sqrt(3)/2, 1.207106781186550, sqrt(3)/2},
   {"T1(2)",   vec3d(1,-1,1)/sqrt(3), vec3d(-1,0,0), vec3d(1,-1,-1)/sqrt(3),
               sqrt(3)/2, 0.207106781186547, sqrt(3)/2},
   {"T1(3)",   vec3d(1,-1,1)/sqrt(3), vec3d(1,0,0), vec3d(1,-1,-1)/sqrt(3),
               0.237248933904040, 1.118033988749890,1.053745514831770},
   // B changed
   {"T2(1)",   vec3d(1,-1,1)/sqrt(3), vec3d(1,-1,-1)/sqrt(3), vec3d(-1,-1,-1)/sqrt(3),
               sqrt(6)/4, 1.020620726159660, sqrt(6)/4},
   {"T2(2)",   vec3d(1,-1,1)/sqrt(3), vec3d(-1,1,1)/sqrt(3), vec3d(-1,-1,-1)/sqrt(3),
               sqrt(6)/4, sqrt(6)/4, sqrt(6)/4},
   {"O1",      vec3d(1,0,1)/sqrt(2), vec3d(1,0,-1)/sqrt(2), vec3d(0,1,0),
               sqrt(2)/2, sqrt(2)/2, sqrt(2)/2},
   // alpha changed
   {"O2(1)",   vec3d(1,0,0), vec3d(1,0,1)/sqrt(2), vec3d(1,1,1)/sqrt(3),
               1.080488239084420,0.118828664898468, 1.094667047615130},
   {"O2(2)",   vec3d(1,0,0), vec3d(1,0,1)/sqrt(2), vec3d(1,1,1)/sqrt(3),
               1.204738778767500, 1.375617671088820,0.515547971500924},
   {"O3",      vec3d(1,0,1)/sqrt(2), vec3d(-1,0,0), vec3d(0,1,0),
               0.541196100146197, 0.541196100146197, 0.840896415253715},
   {"O4(1)",   vec3d(1,-1,1)/sqrt(3), vec3d(1,-1,-1)/sqrt(3), vec3d(-1,-1,-1)/sqrt(3),
               sqrt(6)/4, 1.020620726159660, sqrt(6)/4},
   // B changed
   {"O4(2)",   vec3d(-1,1,-1)/sqrt(3), vec3d(1,-1,-1)/sqrt(3), vec3d(1,1,1)/sqrt(3),
               sqrt(6)/4, sqrt(6)/4, sqrt(6)/4},
   {"O5(1)",   vec3d(1,-1,1)/sqrt(3), vec3d(1,0,0), vec3d(1,-1,-1)/sqrt(3),
               sqrt(3)/2, 1.207106781186550, sqrt(3)/2},
   {"O5(2)",   vec3d(1,-1,1)/sqrt(3), vec3d(-1,0,0), vec3d(1,-1,-1)/sqrt(3),
               sqrt(3)/2, 0.207106781186547, sqrt(3)/2},
   {"O5(3)",   vec3d(1,-1,1)/sqrt(3), vec3d(1,0,0), vec3d(1,-1,-1)/sqrt(3),
               0.237248933904040, 1.118033988749890, 1.053745514831770},
   {"O6(1)",   vec3d(1,-1,1)/sqrt(3), vec3d(1,0,0), vec3d(0,-1,0),
               sqrt(6)/2, sqrt(2)/2, sqrt(2)/2},
   {"O6(2)",   vec3d(-1,1,-1)/sqrt(3), vec3d(1,0,0), vec3d(0,-1,0),
               sqrt(6)/6, sqrt(2)/2, sqrt(2)/2},
   {"O7",      vec3d(1,0,0), vec3d(0,1,0), vec3d(0,0,1),
               sqrt(2)/2, sqrt(2)/2, sqrt(2)/2},
   // alpha changed
   {"I1(1)",   vec3d(phi,phi2,1)/(2*phi), vec3d(1,1,1)/sqrt(3), vec3d(0,phi,1)/term1,
               0.674181480291942, 1.600435268943570, 1.508572470741690},
   {"I1(2)",   vec3d(phi,phi2,1)/(2*phi), vec3d(1,1,1)/sqrt(3), vec3d(0,phi,1)/term1,
               1.901892201462340, 1.042221422390510, 1.602608615469960},
   // the following six added to the table (not in the paper)
   // I2(1) thru I3(3) added by Jim McNeill
   {"I2(1)",   vec3d(-phi,phi2,-1)/(2*phi), vec3d(-1,0,-phi)/term1, vec3d(0,phi,-1)/term1,
               1.035310785747180, 1.017991958789870,0.041794129091506},
   {"I2(2)",   vec3d(phi,phi2,1)/(2*phi), vec3d(-1,0,-phi)/term1, vec3d(0,-phi,-1)/term1,
               0.069590240312534, 0.961660565811817, 0.940133523007424},
   {"I3(1)",   vec3d(-phi,phi2,1)/(2*phi), vec3d(1,1,1)/sqrt(3), vec3d(0,phi,-1)/term1,
               1.067540004877190,0.454812301059250,0.979984198587005},
   {"I3(2)",   vec3d(phi,-phi2,-1)/(2*phi), vec3d(1,1,1)/sqrt(3), vec3d(0,phi,-1)/term1,
               0.365481361740863, 0.809498075533857, 0.758298681854245},
   {"I3(3)",   vec3d(-phi,phi2,1)/(2*phi), vec3d(1,1,1)/sqrt(3), vec3d(0,phi,-1)/term1,
               1.037389492284540,0.123497623343534, 1.015782486189710},
   // I3(4) result added from programmatic search by Adrian Rossiter
   {"I3(4)",   vec3d(-phi,phi2,1)/(2*phi), vec3d(1,1,1)/sqrt(3), vec3d(0,-phi,1)/term1,
               0.840436103166278, 0.919239756651211, 0.257365254781315},
   // A, B changed
   {"I4",      vec3d(phi,-phi2,-1)/(2*phi), vec3d(-1,-phi,phi2)/(2*phi), vec3d(phi2,1,phi)/(2*phi),
               1/sqrt(2), 1/sqrt(2), 1/sqrt(2)},
   // C changed
   {"I5(1)",   vec3d(1,1,1)/sqrt(3), vec3d(-1,0,-phi)/term1, vec3d(phi,-1,0)/term1,
               0.740118744863774, 0.305243296610152, 0.825499999909858},
   // alpha changed
   {"I5(2)",   vec3d(1,1,1)/sqrt(3), vec3d(-1,0,-phi)/term1, vec3d(-phi,1,0)/term1,
               0.095586833624572, 0.922356501413017, 0.977651218274709},
   {"I6(1)",   vec3d(1,1,1)/sqrt(3), vec3d(1,0,phi)/term1, vec3d(phi,1,0)/term1,
               1.572257895003900,term2, term2},
   {"I6(2)",   vec3d(1,1,1)/sqrt(3), vec3d(-1,0,-phi)/term1, vec3d(-phi,-1,0)/term1,
               0.060735266851555, term2, term2},
   // C changed
   {"I7(1)",   vec3d(1,1,1)/sqrt(3), vec3d(0,-phi,1)/term1, vec3d(-phi,1,0)/term1,
               0.706232491219458, term3, term3},
   {"I7(2)",   vec3d(1,1,1)/sqrt(3), vec3d(0,phi,-1)/term1, vec3d(phi,-1,0)/term1,
               0.926760670635994, term3, term3},
   {"I7(3)",   vec3d(1,1,1)/sqrt(3), vec3d(0,-phi,1)/term1, vec3d(phi,-1,0)/term1,
               0.330792269124804, 0.883687468829836, 1.007795749176520},
   // C, beta changed
   {"I8(1)",   vec3d(1,1,1)/sqrt(3), vec3d(0,-iphi,-phi)/sqrt(3), vec3d(phi,-1,0)/term1,
               0.535233134659635, 0.535233134659635, term2},
   // B changed
   {"I8(2)",   vec3d(1,1,1)/sqrt(3), vec3d(0,-iphi,-phi)/sqrt(3), vec3d(-phi,1,0)/term1,
               0.535233134659635, 0.535233134659635, 0.750245100408926},
   {"I9(1)",   vec3d(1,1,1)/sqrt(3), vec3d(1,0,phi)/term1, vec3d(phi,0,iphi)/sqrt(3),
               1.401258538444070, 1.639247476530740, 1.401258538444070},
   {"I9(2)",   vec3d(1,1,1)/sqrt(3), vec3d(1,0,phi)/term1, vec3d(phi,0,iphi)/sqrt(3),
               1.401258538444070,term3, 1.401258538444070},
   {"I9(3)",   vec3d(1,1,1)/sqrt(3), vec3d(1,0,phi)/term1, vec3d(phi,0,iphi)/sqrt(3),
               1.056261606816530, 1.606723212103540, 1.497317965649610},
   // Redundant
   //{"I9(4)",   vec3d(1,1,1)/sqrt(3), vec3d(1,0,phi)/term1, vec3d(phi,0,iphi)/sqrt(3),
   //            1.497317965649610, 1.606723212103540, 1.056261606816530},
   {"I10(1)",  vec3d(1,0,phi)/term1, vec3d(0,-phi,-1)/term1, vec3d(-phi,1,0)/term1,
               term3, term3, term3},
   // C changed
   {"I10(2)",  vec3d(1,0,phi)/term1, vec3d(0,-phi,-1)/term1, vec3d(phi,-1,0)/term1,
               term3, term3, 1.113516364411610},
   // Redundant
   //{"I10(3)",  vec3d(1,0,phi)/term1, vec3d(0,phi,1)/term1, vec3d(phi,-1,0)/term1,
   //            1.113516364411610,term3, term3},
   {"I11(1)",  vec3d(1,0,phi)/term1, vec3d(0,phi,1)/term1, vec3d(phi,1,0)/term1,
               term2, term2, term2},
   {"I11(2)",  vec3d(1,0,phi)/term1, vec3d(0,-phi,-1)/term1, vec3d(phi,1,0)/term1,
               term2, 0.100405707943114, term2},
};

#undef phi
#undef iphi
#undef phi2
#undef term1
#undef term2
#undef term3

class id_poly
{
   private:
      int last_iso_delta;
      IsoDeltaItem* iso_delta_items;
      IsoDeltaVector* iso_delta_vectors;

   public:
      id_poly();
      void list_poly(int idx, FILE *fp=stderr);
      void list_polys(FILE *fp=stderr);
      int lookup_sym_no(string sym);
      int get_last_iso_delta() { return last_iso_delta; }

      const char *get_sym_type(int i);
      vec3d A(int i);
      vec3d B(int i);
      vec3d C(int i);
      double alpha(int i);
      double beta(int i);
      double gamma(int i);
};

id_poly::id_poly()
{ 
   iso_delta_items = iso_delta_item_list;
   iso_delta_vectors = iso_delta_vector_list;
   last_iso_delta = sizeof (iso_delta_item_list) / sizeof (iso_delta_item_list[0]);
}

void id_poly::list_poly(int idx, FILE *fp)
{
   fprintf(fp, "%2d) %-7s %-2s %13s %-45s\n",
      idx+1, iso_delta_items[idx].name, iso_delta_items[idx].sym_type, iso_delta_items[idx].symbol,
         iso_delta_items[idx].comment);
}

void id_poly::list_polys(FILE *fp)
{
   for(int i=0; i<last_iso_delta; i++)
      list_poly(i, fp);
   fprintf(fp,"\n");
   fprintf(fp,"36 Isohedral Deltahedra + 6 Compounds = 42 Total\n\n");
   fprintf(fp,"*  added to the table by Jim McNeill (http://www.orchidpalms.com/polyhedra)\n");
   fprintf(fp,"** added to the table by Adrian Rossiter (http://www.antiprism.com)\n");
}

int id_poly::lookup_sym_no(string sym)
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
         
   // remove spaces either side of a punctuation mark
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
   
   // is it the list order number
   char *endptr;
   int idx = strtol(sym_norm2.c_str(), &endptr, 10);
   if(!*endptr)     // all of string is an integer
      return idx-1;

   idx= -1;
   
   // is it a poly name
   for(unsigned int i=0; i<sym_norm2.size(); i++)
      if(isalpha(sym_norm2[i]))
         sym_norm2[i] = toupper(sym_norm2[i]);
   for(int i=0; i<last_iso_delta; i++) {
      if(sym_norm2==iso_delta_item_list[i].name)
         return i;
      
      if(idx<0 && strncmp(sym_norm2.c_str(), iso_delta_item_list[i].name, sym_norm2.size())==0)
         idx = i;
   }

   return idx;
}

const char *id_poly::get_sym_type(int i)
{
   return iso_delta_item_list[i].sym_type;
}

vec3d id_poly::A(int i)
{
   return iso_delta_vectors[i].A;
}

vec3d id_poly::B(int i)
{
   return iso_delta_vectors[i].B;
}

vec3d id_poly::C(int i)
{
   return iso_delta_vectors[i].C;
}

double id_poly::alpha(int i)
{
   return iso_delta_vectors[i].alpha;
}

double id_poly::beta(int i)
{
   return iso_delta_vectors[i].beta;
}

double id_poly::gamma(int i)
{
   return iso_delta_vectors[i].gamma;
}

class id_opts: public prog_opts {
   public:
      string ifile;
      string ofile;

      bool list_polys;
      string poly;

      bool make_dipyramid;
      bool triangle_only;
      bool verbose;
      bool allow_angles;
      double angle;
      int n;
      int d;
      int k;
      char coloring_method;
      int face_opacity;

      color_map_multi map;
      string case_type;

      id_opts(): prog_opts("iso_delta"),
                 list_polys(false),
                 make_dipyramid(false),
                 triangle_only(false),
                 verbose(false),
                 allow_angles(false),
                 angle(INFINITY),
                 n(0),
                 d(1),
                 k(0),
                 coloring_method('c'),
                 face_opacity(255)
             {}

      void process_command_line(int argc, char **argv);
      void usage();
};

void id_opts::usage()
{
   fprintf(stdout,
"\n"
"Usage: %s [options] polyhedron\n"
"\n"
"Make Isohedral Deltahedra in OFF format. The polyhedron may be specified\n"
"by its list number, or the (start of the) name of the polyhedron.\n"
"Additional examples can be generated using -d or -c paramters.\n"
"Based on a paper by G. C. Shephard\n"
"Periodica Mathematica Hungarica, Volume 39, Numbers 1-3, 2000 , pp. 83-106(24).\n"
"with enhancements by Jim McNeill (http://www.orchidpalms.com/polyhedra)\n"
"and Adrian Rossiter (http://www.antiprism.com)\n"
"\n"
"Options\n"
"%s"
"  -l        display the list of Isohedral Deltahedra 1 thru 42\n"
"  -t        generate triangle only (Isohedral Deltahedra 1 to 42 and option -d)\n"
"  -v        verbose output (Isohedral Deltahedra 1 thru 42 and option -d)\n"
"  -o <file> write output to file (default: write to standard output)\n"
// undocumented switch
//"  -w        allow angle on b,e,f,m,n,o,p,q creates bi-hedral forms\n"
"\nIsohedral Deltahedra Options\n"
"  -a <ang>  angle\n"
"  -n <n/d>  n/d (d is optional)\n"
"              note: for option -d and compound cases c, g, h:\n"
"              n and d must be such that (2 < n/d < 6)\n"
"  -k <k>    k (for various use)\n"
"\nIsohedral Deltahedra Special Cases\n"
"  -d        dipyramid of n/d using -n n/d (infinite set)\n"
"  -c <type> compound cases (a thru f from Shephard's paper)\n"
"              a - tetrahedron repeated k times, evenly spaced\n"
"                     k=1 tetrahedron, k=2 Stella Octangula\n"
"                     Uniform Compound Set UC23 when n/d is 2/1\n"
"              b - 5 or 10 tetrahedra\n"
"                     k=1 icosahedral, k=2 with horizontal reflection\n"
"                     Uniform Compounds UC05 and UC06\n"
"              c - 2 dipyramids of n/d using -a angle (default: calculated)\n"
"                     relaxed dual of Uniform Compound Set UC20 and k=1\n"
"              d - 6 octahedra using -a angle (default: 22.5)\n"
"                     At 45.0 degrees is 3 Octahedra\n"
"                     dual of Uniform Compound Set UC07\n"
"              e - 4 or 8 triangular dipyramids\n"
"                     k=1 octahedral, k=2 with horizontal reflection\n"
"                     relaxed duals of Uniform Compounds UC30 & UC31\n"
"              f - 6 or 12 5/1 pentagonal or 5/2 star dipyramids\n"
"                     5/1: k=1 icosahedral, k=2 with horizontal reflection\n"
"                          relaxed duals of Uniform Compounds UC34 & UC35\n"
"                     5/2: k=3 icosahedral, k=4 with horizontal reflection\n"
"                          relaxed duals of Uniform Compounds UC36 & UC37\n"
"\n         additional cases:\n"
"              g - 2 tetrahedra using -a angle (default: 45.0)\n"
"                     At 45.0 degrees is Uniform Compound UC04\n"
"                     Uniform Compound Set UC23 when n/d is 2/1 and k=1\n"
"              h - 2 tetrahedra repeated k times, evenly spaced\n"
"                     using -a angle (default: 1.0)\n"
"                     Uniform Compound Set UC23 when n/d is 2/1 for any k\n"
"              i - 6 tetrahedra using -a angle (default: 45.0)\n"
"                     Uniform Compound UC01. At 45.0 degrees is UC03\n"
"              j - 12 tetrahedra using -a angle (default: 30.0)\n"
"                     Uniform Compound UC02. At 45.0 degress is UC03\n"
"              k - 2 dipyramids of n/d repeated k times, evenly spaced\n"
"                     using -a angle (default: 1.0)\n"
"                     relaxed dual of Uniform Compound Set UC20\n"
"              l - k dipyramids of n/d using -n n/d, evenly spaced\n"
"                     relaxed dual of Uniform Compound Set UC21\n"
"              m - 10 or 20 triangular dipyramids\n"
"                     k=1 icosahedral, k=2 with horizontal reflection\n"
"                     relaxed duals of Uniform Compounds UC32 & UC33\n"
"              n - 6 10/3 star dipyramids\n"
"                     relaxed dual of Uniform Compounds UC41\n"
"              o - 5 or 10 Augmented Tetrahedra T2(1)\n"
"                     k=1 icosahedral, k=2 with horizontal reflection\n"
"                     relaxed duals of Uniform Compounds UC55 & UC56\n"
"              p - 5 Augmented Octahedra O6(1)\n"
"                     relaxed dual of Uniform Compounds UC57\n"
"              q - 5 Excavated Octahedra O6(2)\n"
"                     relaxed dual of Uniform Compounds UC58\n"
"\nColoring Options\n"
"  -f <opt> compound coloring\n"
"              key word: none - sets no color (default: c)\n"
"              c - unique coloring for each constituent\n"
"              s - symmetric colouring (should always be one color)\n"
"  -T <tran> face transparency for color values. valid range from 0 to 255\n"
"               0 - invisible  255 - opaque (default 255)\n"
"  -m <maps> color maps for all elements to be tried in turn (default: compound)\n"
"\n"
"\n", prog_name(), help_ver_text);
}

void id_opts::process_command_line(int argc, char **argv)
{
   opterr = 0;
   char c;
   char errmsg[MSG_SZ];

   string map_file;
   
   handle_long_opts(argc, argv);

   while( (c = getopt(argc, argv, ":hldtvwc:n:a:k:f:T:m:o:")) != -1 ) {
      if(common_opts(c, optopt))
         continue;

      switch(c) {
         case 'l':
            list_polys = true;
            break;

         case 'd':
            make_dipyramid = true;
            break;

         case 't':
            triangle_only = true;
            break;

         case 'v':
            verbose = true;
            break;

         // undocumented switch
         case 'w':
            allow_angles = true;
            warning("using -w to allow angles creates bihedral forms");
            break;

         case 'c':
            if(strspn(optarg, "abcdefghijklmnopq") != strlen(optarg) ||
                                                         strlen(optarg)>1)
               error(msg_str("case type is '%s' must be only one of a thru q",
                        optarg), c);
            case_type=optarg;
            break;

         case 'a':
            if(!read_double(optarg, &angle, errmsg))
               error(errmsg, c);
            angle = deg2rad(angle);
            break;
            
         case 'n': {
            char *p;
            p = strchr(optarg, '/');
            if(p!=0) {
               *p++='\0';
               if(!read_int(p, &d, errmsg))
                  error(errmsg, "n/d (d part)");
            }

            if(!read_int(optarg, &n, errmsg))
               error(errmsg, "n/d (n part)");
            if(n<1)
               error("must be an integer 1 or greater", "n/d (n part)");
            if(d < 1)
               error("d must be 1 or greater", "n/d (d part)");
            break;
         }

         case 'k':
            if(!read_int(optarg, &k, errmsg))
               error(errmsg, c);
            if(k<1)
               error("must be an integer 1 or greater", c);
            break;
            
         case 'f':
            if(!strcasecmp(optarg,"none"))
               coloring_method = '\0';
            else
            if(strspn(optarg, "cs") != strlen(optarg) || strlen(optarg)>1)
               error(msg_str("invalid coloring method '%s'", optarg), c);
            else
               coloring_method = *optarg;
            break;

         case 'T':
            if(!read_int(optarg, &face_opacity, errmsg))
               error(errmsg, c);
            if(face_opacity < 0 || face_opacity > 255) {
               error("face transparency must be between 0 and 255", c);
            }
            break;

         case 'm':
            map_file = optarg;
            break;

         case 'o':
            ofile = optarg;
            break;

         default:
            error("unknown command line error");
      }
   }


   if (case_type == "c" || case_type == "k" || case_type == "l" || make_dipyramid) {
      if (n <= 0)
         error("n must be greater than 0");

      if(d >= n)
         error("d must be less than number of sides", "n/d (d part)");
      if(((double)n/d <= 2.0 ) || ((double)n/d >= 6.0))
         error("2 < n/d < 6 is enforced");
   }
   else
   if (n != 0)
      warning("n/d is ignored");


   if ((case_type == "a" || make_dipyramid) || 
       (!allow_angles && 
        (case_type == "b" || case_type == "e" || case_type == "f" || case_type == "l" ||
         case_type == "m" || case_type == "n" || case_type == "o" || case_type == "p" || case_type == "q"))) {
      if (angle != INFINITY) {
         warning("-a angle is ignored");
         angle = INFINITY;
      }
   }

   // check k
   if (case_type == "a" || case_type == "h" || case_type == "k" || case_type == "l") {
      if (k <= 0)
         error("k must be greater then 0");
   }
   else
   if (case_type == "b" || case_type == "e" || case_type == "m" || case_type == "o") {
      if (k != 1 && k !=2)
         error("k must 1 or 2");
   }
   else
   if (case_type == "f") {
      if (k < 1 || k > 4)
         error("k must 1, 2, 3 or 4");
   }
   else
   if (k != 0)
      warning("k is ignored");
    
   // if one from the table or -d then force color by constituent    
   if (!case_type.length() && (coloring_method == 'T' || coloring_method == 't')) {
      if (coloring_method == 'T')
         coloring_method = 'C';
      else
      if (coloring_method == 't')
         coloring_method = 'c';
   }

   if(argc==optind && !list_polys && !case_type.length() && !make_dipyramid)
      error("no polyhedron specified", "polyhedron");
   
   while(optind < argc) {
      poly += argv[optind];
      poly += " ";
      optind++;
   }

   if(poly.length() != 0 && (case_type.length() || make_dipyramid))
      warning("polyhedron specifier ignored if using -c or -d");

   if (!map_file.size())
      map_file = "compound";
   if(!map.init(map_file.c_str(), errmsg))
      error(errmsg, c);
}

void verbose_output(vec3d A, vec3d B, vec3d C, double alpha, double beta, double gamma)
{
   A.dump("\nVector A",stderr);
   B.dump("Vector B",stderr);
   C.dump("Vector C",stderr);

   fprintf(stderr,"\nalpha = %1.15f\n", alpha);
   fprintf(stderr,"beta = %1.15f\n", beta);
   fprintf(stderr,"gamma = %1.15f\n\n", gamma);
}

void make_triangle(col_geom_v &geom, vec3d A, vec3d B, vec3d C, double alpha, double beta, double gamma)
{
   geom.add_vert(A*alpha);
   geom.add_vert(B*beta);
   geom.add_vert(C*gamma);

   vector<int> face;
   face.push_back(0);
   face.push_back(1);
   face.push_back(2);
   geom.add_face(face);
}

void make_poly(col_geom_v &geom, string sym_type, bool triangle_only, bool verbose,
               vec3d A, vec3d B, vec3d C, double alpha, double beta, double gamma)
{
   if (verbose)
      verbose_output(A, B, C, alpha, beta, gamma);

   // make poly called twice from dipyramid code, don't remake triangle
   if (geom.verts().size() == 0)
      make_triangle(geom, A, B, C, alpha, beta, gamma);

   if(sym_type[0]=='D' || sym_type[0]=='I' )
      geom.transform(mat3d::rot(0, M_PI/2,0));

   if (!triangle_only) {
      sym_repeat(geom, geom, sym_type);
      sort_merge_elems(geom, "vef", 1e-8);
   }
}

void make_delta_dipyramid(col_geom_v &geom, int n, int d, bool triangle_only = false, bool verbose = false)
{
   char buf1[MSG_SZ];
   char buf2[MSG_SZ];

   sprintf(buf1,"D%dh",n);
   sprintf(buf2,"[%d/%d,1/2,1/2]",d,n);
   fprintf(stderr, "Dihedral Group:  %s  %s  %d/%d %sdipyramid\n", buf1, buf2, n, d, ((d == 1) ? "" : "star "));

   string sym_type = buf1;

   vec3d A = vec3d(1,0,0);
   vec3d B = vec3d(cos((2.0*M_PI*d)/n),sin((2.0*M_PI*d)/n),0);
   vec3d C = vec3d(0,0,1);
   double alpha = 1/(2.0*sin((M_PI*d)/n));
   double beta = alpha;
   double gamma = sqrt(1.0 - alpha*alpha);

   // force triangle_only, verbose to false
   make_poly(geom, sym_type, true, false, A, B, C, alpha, beta, gamma);

   if (!triangle_only)
      make_poly(geom, sym_type, false, verbose, A, B, C, alpha, beta, gamma);
   else if (verbose)
      verbose_output(A, B, C, alpha, beta, gamma);
}

/*
void unit_tetrahedron(col_geom_v &geom)
{
   tetrahedron(geom);
   geom.transform(mat3d::scale((1.0/4)*sqrt(2)));
}

void unit_octahedron(col_geom_v &geom)
{
   octahedron(geom);
   geom.transform(mat3d::scale(1.0/sqrt(2)));
}
*/

void tet_to_dihedral(col_geom_v &geom, string sym_from, int k, mat3d pos=mat3d())
{
   char sym_to[MSG_SZ];
   sprintf(sym_to,"D%d%s",k,((k%4 == 0) ? "h" : "v"));
   transform_and_repeat(geom, sym_to, sym_from, pos);
}


void case_a_star_tetrahedron(col_geom_v &geom, int k)
{
   geom.read_resource("u1");
   tet_to_dihedral(geom, "Td", 2*k);
}

void case_b_5_or_10_tetrahedra(col_geom_v &geom, double angle, int k)
{
   if (angle == INFINITY)
      angle = 0;

   geom.read_resource("u1");
   
   // to construct in one statement for Ih
   // transform_and_repeat(geom, (k == 1 ? "I" : "Ih"), "Td", mat3d::rot(0,angle,0));

   if (k == 1)
      transform_and_repeat(geom, "I", "Td", mat3d::rot(0,angle,0));
   else
   if (k == 2) {
      transform_and_repeat(geom, "Oh", "Td", mat3d::rot(0,angle,0));
      transform_and_repeat(geom, "I", "Oh");
   }
}

void case_c_2_dipyramids(col_geom_v &geom, double angle, int n, int d)
{
   if (angle == INFINITY) {
      angle = (M_PI/2)/n; // 90/n degrees
      fprintf(stderr,"angle calculated is %g\n",rad2deg(angle));
   }

   fprintf(stderr,"Using: ");
   make_delta_dipyramid(geom, n, d);

   char sym_from[MSG_SZ];
   sprintf(sym_from,"D%dh",n);
   char sym_to[MSG_SZ];
   sprintf(sym_to,"D%dh",2*n);

   // advance angle so that angle = 0 is coincident constituents
   geom.transform(mat3d::rot(0,0,angle+M_PI/(2*n)));
   transform_and_repeat(geom, sym_to, sym_from);
}

void case_d_6_octahedra(col_geom_v &geom, double angle)
{
   if (angle == INFINITY)
      angle = (M_PI/8); // 22.5 degrees

   geom.read_resource("u5");

   // at 0 degrees, produced 3 coincident octahedra
   transform_and_repeat(geom, "D2h", "Oh", mat3d::rot(0,0,angle));
   transform_and_repeat(geom, "T", "D2h");
}

void case_e_4_or_8_triangular_dipyramids(col_geom_v &geom, double angle, int k)
{
   if (angle == INFINITY)
      angle = 0;

   fprintf(stderr,"Using: ");
   make_delta_dipyramid(geom, 3, 1);

   // to construct in one statement for Oh
   // transform_and_repeat(geom, (k == 1 ? "O" : "Oh"), "D3h",
   //    mat3d::rot(vec3d(0,0,1),vec3d(1,1,1)) * mat3d::rot(0,0,angle+M_PI/12));
   
   if (k == 1)
      transform_and_repeat(geom, "O", "D3h",
         mat3d::rot(vec3d(0,0,1),vec3d(1,1,1)) * mat3d::rot(0,0,angle+M_PI/12));
   else
   if (k == 2) {
      transform_and_repeat(geom, "D6h", "D3h");
      transform_and_repeat(geom, "O", "D6h",
         mat3d::rot(vec3d(0,0,1),vec3d(1,1,1)) * mat3d::rot(0,0,angle+M_PI/12));
   }
}

void case_f_6_or_12_pentagonal_dipyramids(col_geom_v &geom, double angle, int k)
{
   if (angle == INFINITY)
      angle = 0;

   fprintf(stderr,"Using: ");
   if (k == 1 || k == 2)
      make_delta_dipyramid(geom, 5, 1);
   else
   if (k == 3 || k == 4)
      make_delta_dipyramid(geom, 5, 2);
      
   double phi = (1 + sqrt(5))/2;

   // to construct in one statement for Ih
   // transform_and_repeat(geom, ((k == 1 || k == 3) ? "I" : "Ih"), "D5h",
   //    mat3d::rot(vec3d(0,0,1),vec3d(0,1,phi)) * mat3d::rot(0,0,angle+M_PI/5));

   if (k == 1 || k == 3)
      transform_and_repeat(geom, "I", "D5h",
         mat3d::rot(vec3d(0,0,1),vec3d(0,1,phi)) * mat3d::rot(0,0,angle+M_PI/5));
   else
   if (k == 2 || k == 4) {
      transform_and_repeat(geom, "D10h", "D5h");
      transform_and_repeat(geom, "I", "D10h",
         mat3d::rot(vec3d(0,0,1),vec3d(0,1,phi)) * mat3d::rot(0,0,angle+M_PI/5));
   }
}

void case_g_2_tetrahedra(col_geom_v &geom, double angle)
{
   if (angle == INFINITY)
      angle = (M_PI/4); // 45 degrees

   geom.read_resource("u1");

   // advance angle so that angle = 0 is coincident constituents
   geom.transform(mat3d::rot(0,0,angle));
   tet_to_dihedral(geom, "S4", 2); // 2*k=4
}

void case_h_2k_tetrahedra(col_geom_v &geom, double angle, int k)
{
   if (angle == INFINITY)
      angle = deg2rad(1.0);

   case_g_2_tetrahedra(geom, angle);
   tet_to_dihedral(geom, "D2v", 2*k);
}

void case_i_6_tetrahedra(col_geom_v &geom, double angle)
{
   geom.read_resource("u1");
   
   transform_and_repeat(geom, "T", "T", mat3d::rot(0,0,angle));
}

void case_j_12_tetrahedra(col_geom_v &geom, double angle)
{
   if (angle == INFINITY)
      angle = (M_PI/6); // 30 degrees;
   
   geom.read_resource("u1");
   
   transform_and_repeat(geom, "S2", "T");
   transform_and_repeat(geom, "T", "Oh", mat3d::rot(0,0,angle));
}

void case_k_2k_dipyramids(col_geom_v &geom, double angle, int k, int n, int d)
{
   if (angle == INFINITY)
      angle = deg2rad(1.0);

   fprintf(stderr,"Using: ");
   make_delta_dipyramid(geom, n, d);

   char sym_from[MSG_SZ];
   sprintf(sym_from,"D%dh",n);   
   char sym_to[MSG_SZ];
   sprintf(sym_to,"D%dh",k*n);

   transform_and_repeat(geom, sym_from, sym_from, mat3d::rot(0,0,angle+M_PI));
   transform_and_repeat(geom, sym_to, sym_from);
}

void case_l_k_dipyramids(col_geom_v &geom, int k, int n, int d)
{
   fprintf(stderr,"Using: ");
   make_delta_dipyramid(geom, n, d);

   char sym_from[MSG_SZ];
   sprintf(sym_from,"D%dh",n);
   char sym_to[MSG_SZ];
   sprintf(sym_to,"D%dh",k*n);

   transform_and_repeat(geom, sym_to, sym_from);
}

void case_m_10_or_20_triangular_dipyramids(col_geom_v &geom, double angle, int k)
{
   if (angle == INFINITY)
      angle = 0;

   fprintf(stderr,"Using: ");
   make_delta_dipyramid(geom, 3, 1);

   double phi = (1 + sqrt(5))/2;
   double iphi = 1/phi;

   // to construct in one statement for Ih
   // transform_and_repeat(geom, (k == 1 ? "I" : "Ih"), "D3h",
   //   mat3d::rot(vec3d(0,0,1), vec3d(iphi,0,phi)) * mat3d::rot(0,0,angle+M_PI/6));

   if (k == 1)
      transform_and_repeat(geom, "I", "D3h",
         mat3d::rot(vec3d(0,0,1), vec3d(iphi,0,phi)) * mat3d::rot(0,0,angle+M_PI/6));
   else
   if (k == 2) {
      transform_and_repeat(geom, "D6h", "D3h");
      transform_and_repeat(geom, "I", "D6h",
         mat3d::rot(vec3d(0,0,1), vec3d(iphi,0,phi)) * mat3d::rot(0,0,angle+M_PI/6));
   }
}

void case_n_6_10_3_star_dipyramids(col_geom_v &geom, double angle)
{
   if (angle == INFINITY)
      angle = 0;

   fprintf(stderr,"Using: ");
   make_delta_dipyramid(geom, 10, 3);

   double phi = (1 + sqrt(5))/2;

   transform_and_repeat(geom, "I", "D10h",
      mat3d::rot(vec3d(0,0,1), vec3d(0,1,phi)) * mat3d::rot(0,0,angle+M_PI/10));
}

void case_o_5_or_10_augmented_tetrahedra(col_geom_v &geom, double angle, int k)
{
   if (angle == INFINITY)
      angle = 0;

   // to construct in one statement for Ih
   // transform_and_repeat(geom, (k == 1 ? "I" : "Ih"), "Td", mat3d::rot(0,angle,0));

   if (k == 1)
      transform_and_repeat(geom, "I", "Td", mat3d::rot(0,angle,0));
   else
   if (k == 2) {
      transform_and_repeat(geom, "Oh", "Td", mat3d::rot(0,angle,0));
      transform_and_repeat(geom, "I", "Oh");
   }
}

void case_p_5_augmented_octahedra(col_geom_v &geom, double angle)
{
   if (angle == INFINITY)
      angle = 0;

   transform_and_repeat(geom, "I", "Oh", mat3d::rot(0,angle,0));
}

void case_q_5_excavated_octahedra(col_geom_v &geom, double angle)
{
   if (angle == INFINITY)
      angle = 0;

   transform_and_repeat(geom, "I", "Oh", mat3d::rot(0,angle,0));
}

void compound_coloring(col_geom_v &geom, char coloring_method, color_map_multi &map, int face_opacity)
{
   // color by sub-symmetry  as map indexes happened by default in sym_repeat()
   if (!coloring_method) {
      // no color, strip colors
      geom.clear_f_cols();
      geom.clear_edges();
      geom.clear_v_cols();
   }
   else {
      coloring clrng;
      clrng.add_cmap(map.clone());
      clrng.set_geom(&geom);

      if (coloring_method == 'c') {
         // color by constituents
         clrng.f_parts(true);
      }
      else
      if (coloring_method == 's') {
         sch_sym sym;
         vector<vector<set<int> > > sym_equivs;
         sym.init(geom, &sym_equivs);
         clrng.f_sets(sym_equivs[2], true);
      }

      // blend edges
      geom.add_missing_impl_edges();
      clrng.e_face_color();
      clrng.v_face_color();

      // transparency
      if (face_opacity != 255) {
         for (unsigned int i=0;i<geom.faces().size();i++) {
            col_val col = geom.get_f_col(i);
            if (col.is_val())
               col = col_val(col[0],col[1],col[2],face_opacity);
            geom.set_f_col(i,col);
         }
      }
   }
}

int main(int argc, char *argv[])
{
   id_opts opts;
   opts.process_command_line(argc, argv);
   id_poly id_polys;

   if(opts.list_polys) {
      id_polys.list_polys();
      exit(0);
   }

   col_geom_v geom;
   if (opts.make_dipyramid)
      make_delta_dipyramid(geom, opts.n, opts.d, opts.triangle_only, opts.verbose);
   else
   if (opts.case_type.length()) {
      if (opts.case_type == "a")
         case_a_star_tetrahedron(geom, opts.k);
      else
      if (opts.case_type == "b")
         case_b_5_or_10_tetrahedra(geom, opts.angle, opts.k);
      else
      if (opts.case_type == "c")
         case_c_2_dipyramids(geom, opts.angle, opts.n, opts.d);
      else
      if (opts.case_type == "d")
         case_d_6_octahedra(geom, opts.angle);
      else
      if (opts.case_type == "e")
         case_e_4_or_8_triangular_dipyramids(geom, opts.angle, opts.k);
      else
      if (opts.case_type == "f")
         case_f_6_or_12_pentagonal_dipyramids(geom, opts.angle, opts.k);
      else
      if (opts.case_type == "g")
         case_g_2_tetrahedra(geom, opts.angle);
      else
      if (opts.case_type == "h")
         case_h_2k_tetrahedra(geom, opts.angle, opts.k);
      else
      if (opts.case_type == "i")
         case_i_6_tetrahedra(geom, opts.angle);
      else
      if (opts.case_type == "j")
         case_j_12_tetrahedra(geom, opts.angle);
      else
      if (opts.case_type == "k")
         case_k_2k_dipyramids(geom, opts.angle, opts.k, opts.n, opts.d);
      else
      if (opts.case_type == "l")
         case_l_k_dipyramids(geom, opts.k, opts.n, opts.d);
      else
      if (opts.case_type == "m")
         case_m_10_or_20_triangular_dipyramids(geom, opts.angle, opts.k);
      else
      if (opts.case_type == "n")
         case_n_6_10_3_star_dipyramids(geom, opts.angle);
      else
      if (opts.case_type == "o") {
         int sym_no = 3;
         string sym_type = id_polys.get_sym_type(sym_no);
         make_poly(geom, sym_type, false, false,
                   id_polys.A(sym_no), id_polys.B(sym_no), id_polys.C(sym_no),
                   id_polys.alpha(sym_no), id_polys.beta(sym_no), id_polys.gamma(sym_no));
         case_o_5_or_10_augmented_tetrahedra(geom, opts.angle, opts.k);
      }
      else
      if (opts.case_type == "p") {
         int sym_no = 14;
         string sym_type = id_polys.get_sym_type(sym_no);
         make_poly(geom, sym_type, false, false,
                   id_polys.A(sym_no), id_polys.B(sym_no), id_polys.C(sym_no),
                   id_polys.alpha(sym_no), id_polys.beta(sym_no), id_polys.gamma(sym_no));
         case_p_5_augmented_octahedra(geom, opts.angle);
      }
      else
      if (opts.case_type == "q") {
         int sym_no = 15;
         string sym_type = id_polys.get_sym_type(sym_no);
         make_poly(geom, sym_type, false, false,
                   id_polys.A(sym_no), id_polys.B(sym_no), id_polys.C(sym_no),
                   id_polys.alpha(sym_no), id_polys.beta(sym_no), id_polys.gamma(sym_no));
         case_q_5_excavated_octahedra(geom, opts.angle);
      }
   }
   else {
      int sym_no = id_polys.lookup_sym_no(opts.poly);
      if(sym_no >= id_polys.get_last_iso_delta())
         opts.error("polyhedron number '"+opts.poly+"' out of range");
      if(sym_no <0)
         opts.error("unknown polyhedron '"+opts.poly+"'");

      id_polys.list_poly(sym_no);

      string sym_type = id_polys.get_sym_type(sym_no);
      make_poly(geom, sym_type, opts.triangle_only, opts.verbose,
                id_polys.A(sym_no), id_polys.B(sym_no), id_polys.C(sym_no),
                id_polys.alpha(sym_no), id_polys.beta(sym_no), id_polys.gamma(sym_no));
   }

   geom.orient();

//fprintf(stderr,"coloring_method = >%c<\n",opts.coloring_method);

   compound_coloring(geom, opts.coloring_method, opts.map, opts.face_opacity);
   
   char errmsg[MSG_SZ]="";
   if(!geom.write(opts.ofile, errmsg))
      opts.error(errmsg);

   return 0;
}
