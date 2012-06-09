/*
   Copyright (c) 2003, Adrian Rossiter

   Antiprism - http://www.antiprism.com

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
   Name: std_polys.cc
   Description: Uniform Polyhedra, Uniform Compounds and Johnson Solids
   Project: Antiprism - http://www.antiprism.com
*/

#include <stdio.h>
#include <stdlib.h>
#include <algorithm>

#include "polygons.h"
#include "std_polys.h"
#include "math_utils.h"
#include "utils.h"
#include "info.h"
#include "coloring.h"

using std::swap;


void normalised_face_list(geom_if &geom)
{
   geom.orient();
   if(vdot(geom.face_norm(0), geom.face_cent(0))<0)
      geom.orient_reverse();
   sort(geom.raw_faces().begin(), geom.raw_faces().end());
}


static void rh_dodecahedron(geom_if &geom)
{
   geom.clear_all();
   col_geom_v geom2;
   geom2.read_resource("std_cube");
   geom.add_verts(geom2.verts());
   geom2.read_resource("std_oct");
   geom2.transform(mat3d::scale(2));
   geom.add_verts(geom2.verts());
   geom.add_hull();
   normalised_face_list(geom);
}


static void rh_triacontahedron(geom_if &geom)
{
   geom.clear_all();
   col_geom_v geom2;
   geom2.read_resource("std_ico");
   geom.add_verts(geom2.verts());
   geom2.read_resource("std_dod");
   geom.add_verts(geom2.verts());
   geom.add_hull();
   normalised_face_list(geom);
}


static void rh_enneacontahedron(geom_if &geom)
{
   geom.clear_all();
   col_geom_v geom2;
   geom2.read_resource("std_dod");
   geom2.transform(mat3d::scale(1/(2*phi*phi)));
   geom.set_zono(geom2.verts());
   normalised_face_list(geom);
}


int make_resource_std_poly(geom_if &geom, string name, char *errmsg=0)
{
   if(name.size()<5 || name.substr(0,4)!="std_" || name.find('.')!=string::npos)
      return -1; // not std poly name (the "." indicates a likely local file)
                 // so the name is not handled

   map<string, model_func> models;
   models["rhombic_dodecahedron"]        = rh_dodecahedron;
   models["rh_dodecahedron"]             = rh_dodecahedron;
   models["rh_dod"]                      = rh_dodecahedron;
   models["rd"]                          = rh_dodecahedron;
   models["rhombic_triacontahedron"]     = rh_triacontahedron;
   models["rh_triacontahedron"]          = rh_triacontahedron;
   models["rh_tri"]                      = rh_triacontahedron;
   models["rt"]                          = rh_triacontahedron;
   models["rhombic_enneacontahedron"]    = rh_enneacontahedron;
   models["rh_enneacontahedron"]         = rh_enneacontahedron;
   models["rh_ennea"]                    = rh_enneacontahedron;
   models["re"]                          = rh_enneacontahedron;

   map<string, model_func>::iterator mi = models.find(name.substr(4));
   if(mi != models.end())
      mi->second(geom);
   else {
      if(errmsg)
         snprintf(errmsg, MSG_SZ, "invalid standard polyhedron name");
      return 1; // fail
   }
   return 0; // name found
}


string expand_abbrevs(const string &name, const char *abbrevs[][2], size_t last)
{
   string expanded;
   char name_cpy[MSG_SZ];
   strncpy(name_cpy, name.c_str(), MSG_SZ);
   vector<char *> parts;
   split_line(name_cpy, parts, RES_SEPARATOR);
   for(unsigned int i=1; i<parts.size(); i++) {
      size_t j;
      for(j=0; j<last; j++) {
         if(strcmp(parts[i], abbrevs[j][0])==0)
            break;
      }
      if(i>1)
         expanded += ' ';
      expanded += (j<last) ? abbrevs[j][1] : parts[i];
   }
   return expanded;
}

void set_resource_polygon_color(geom_if &geom)
{
   if(col_geom_v *cg = dynamic_cast<col_geom_v *>(&geom)) {
      coloring clrng(cg);
      color_map *cmap = init_color_map("uniform");
      clrng.add_cmap(cmap);
      clrng.f_avg_angle(true);
      clrng.v_avg_angle(true);
   }
}

static void make_resource_dual(geom_if &geom, bool is_std = false)
{
      vec3d cent = geom.centroid();
      geom_info info(geom);
      info.set_center(cent);
      double rad = info.impl_edge_dists().sum/info.num_iedges();
      col_geom_v dual;
      const double inf = 1200;
      get_dual(geom, dual, rad, cent, inf);
      add_extra_ideal_elems(dual, cent, 0.95*inf); // limit closer than inf
      geom.clear_all();
      geom.append(dual);
      if (!is_std)
         set_resource_polygon_color(geom);
}


const char *alt_names[][2] = {
   {"tet", "u1" }, 
   {"tetrahedron", "u1" }, 
   {"truncated_tetrahedron", "u2"}, 
   {"tr_tetrahedron", "u2"}, 
   {"tr_tet", "u2"}, 
   {"octahedron", "u5" }, 
   {"oct", "u5" }, 
   {"cube", "u6" }, 
   {"cuboctahedron", "u7" }, 
   {"cubo", "u7" }, 
   {"truncated_octahedron", "u8" }, 
   {"tr_octahedron", "u8" }, 
   {"tr_oct", "u8" }, 
   {"truncated_cube", "u9" }, 
   {"tr_cube", "u9" }, 
   {"rhombicuboctahedron", "u10" }, 
   {"rhombicubo", "u10" }, 
   {"rh_cubo", "u10" }, 
   {"truncated_cuboctahedron", "u11" }, 
   {"tr_cuboctahedron", "u11" }, 
   {"tr_cubo", "u11" }, 
   {"snub_cube", "u12" }, 
   {"sn_cube", "u12" }, 
   {"icosahedron", "u22" }, 
   {"ico", "u22" }, 
   {"icosa", "u22" }, 
   {"dodecahedron", "u23" }, 
   {"dod", "u23" }, 
   {"icosidodecahedron", "u24" }, 
   {"icosid", "u24" }, 
   {"truncated_icosahedron", "u25" }, 
   {"tr_icosahedron", "u25" }, 
   {"tr_icosa", "u25" }, 
   {"tr_ico", "u25" }, 
   {"truncated_dodecahedron", "u26" }, 
   {"tr_dodecahedron", "u26" }, 
   {"tr_dod", "u26" }, 
   {"rhombicosidodecahedron", "u27" }, 
   {"rhombicosid", "u27" }, 
   {"rh_icosid", "u27" }, 
   {"truncated_icosidodecahedron", "u28" }, 
   {"tr_icosidodecahedron", "u28" }, 
   {"tr_icosid", "u28" }, 
   {"snub_dodecahedron", "u29" },
   {"snub_dod", "u29" },
   {"sn_dod", "u29" },
   {"small_stellated_dodecahedron", "u34"},
   {"sm_st_dodecahedron", "u34"},
   {"sm_st_dod", "u34"},
   {"great_dodecahedron", "u35"},
   {"gr_dodecahedron", "u35"},
   {"gr_dod", "u35"},
   {"great_stellated_dodecahedron", "u52"},
   {"gr_st_dodecahedron", "u52"},
   {"gr_st_dod", "u52"},
   {"great_icosahedron", "u53"},
   {"gr_icosahedron", "u53"},
   {"gr_ico", "u53"},
   {"triakistetrahedron", "u2_d"},
   {"triakis_tetrahedron", "u2_d"},
   {"tri_tet", "u2_d"},
   {"rhombic_dodecahedron", "u7_d"},
   {"rh_dodecahedron", "u7_d"},
   {"rh_dod", "u7_d"},
   {"rd", "u7_d"},
   {"triakisoctahedron", "u9_d"},
   {"triakis_octahedron", "u9_d"},
   {"tri_oct", "u9_d"},
   {"tetrakishexahedron", "u8_d"},
   {"tetrakis_hexahedron", "u8_d"},
   {"tetr_hex", "u8_d"},
   {"disdyakishexahedron", "u8_d"},
   {"disdyakis_hexahedron", "u8_d"},
   {"disd_hex", "u8_d"},
   {"deltoidal_icositetrahedron", "u10_d"},
   {"delt_icosit", "u10_d"},
   {"trapezoidal_icositetrahedron", "u10_d"},
   {"trap_icosit", "u10_d"},
   {"hexakisoctahedron", "u11_d"},
   {"hexakis_octahedron", "u11_d"},
   {"hex_oct", "u11_d"},
   {"disdyakisdodecahedron", "u11_d"},
   {"disdyakis_dodecahedron", "u11_d"},
   {"disd_dod", "u11_d"},
   {"pentagonal_icositetrahedron", "u12_d"},
   {"pen_icosit", "u12_d"},
   {"rhombic_triacontahedron", "u24_d"},
   {"rh_triacontahedron", "u24_d"},
   {"rh_tri", "u24_d"},
   {"rt", "u24_d"},
   {"triakisicosahedron", "u26_d"},
   {"triakis_icosahedron", "u26_d"},
   {"tri_ico", "u26_d"},
   {"pentakisdodecahedron", "u25_d"},
   {"pentakis_dodecahedron", "u25_d"},
   {"pent_dod", "u25_d"},
   {"deltoidal_hexacontahedron", "u27_d"},
   {"delt_hexac", "u27_d"},
   {"trapezoidal_hexacontahedron", "u27_d"},
   {"trap_hexac", "u27_d"},
   {"hexakisicosahedron", "u28_d"},
   {"hexakis_icosahedron", "u28_d"},
   {"hex_ico", "u28_d"},
   {"disdyakistriacontahedron", "u28_d"},
   {"disdyakis_triacontahedron", "u28_d"},
   {"disd_tri", "u28_d"},
   {"pentagonal_hexacontahedron", "u29_d"},
   {"pen_hexac", "u29_d"},

   // Wenninger Numbers which are not Uniforms
   {"w19", "uc4"},
   {"w23", "uc17"},
   {"w24", "uc5"},
   {"w25", "uc6"},
   {"w26", "u30_d"},
   {"w34", "u47_d"},

   // Bowers short names for Uniforms
   {"tet", "u1"},
   {"tut", "u2"},
   {"oho", "u3"},
   {"thah", "u4"},
   {"oct", "u5"},
   {"cube", "u6"},
   {"co", "u7"},
   {"toe", "u8"},
   {"tic", "u9"},
   {"sirco", "u10"},
   {"girco", "u11"},
   {"snic", "u12"},
   {"socco", "u13"},
   {"gocco", "u14"},
   {"cho", "u15"},
   {"cotco", "u16"},
   {"querco", "u17"},
   {"sroh", "u18"},
   {"quith", "u19"},
   {"quitco", "u20"},
   {"groh", "u21"},
   {"ike", "u22"},
   {"doe", "u23"},
   {"id", "u24"},
   {"ti", "u25"},
   {"tid", "u26"},
   {"srid", "u27"},
   {"grid", "u28"},
   {"snid", "u29"},
   {"sidtid", "u30"},
   {"siid", "u31"},
   {"seside", "u32"},
   {"saddid", "u33"},
   {"sissid", "u34"},
   {"gad", "u35"},
   {"did", "u36"},
   {"tigid", "u37"},
   {"raded", "u38"},
   {"sird", "u39"},
   {"siddid", "u40"},
   {"ditdid", "u41"},
   {"gidditdid", "u42"},
   {"sidditdid", "u43"},
   {"ided", "u44"},
   {"idtid", "u45"},
   {"sided", "u46"},
   {"gidtid", "u47"},
   {"giid", "u48"},
   {"seihid", "u49"},
   {"siddy", "u50"},
   {"sidhid", "u51"},
   {"gissid", "u52"},
   {"gike", "u53"},
   {"gid", "u54"},
   {"tiggy", "u55"},
   {"ri", "u56"},
   {"gosid", "u57"},
   {"quitsissid", "u58"},
   {"quitdid", "u59"},
   {"isdid", "u60"},
   {"gaddid", "u61"},
   {"sidhei", "u62"},
   {"giddy", "u63"},
   {"gisdid", "u64"},
   {"gidhei", "u65"},
   {"quitgissid", "u66"},
   {"qrid", "u67"},
   {"gaquatid", "u68"},
   {"gisid", "u69"},
   {"gidhid", "u70"},
   {"geihid", "u71"},
   {"sirsid", "u72"},
   {"gird", "u73"},
   {"girsid", "u74"},
   {"gidrid", "u75"},
   {"pip", "u76"},
   {"pap", "u77"},
   {"stip", "u78"},
   {"stap", "u79"},
   {"starp", "u80"},

   // Bowers short names for Prisms and Antiprisms
   {"trip", "pri3"},
   {"hip", "pri6"},
   {"hep", "pri7"},
   {"ship", "pri7/2"},
   {"giship", "pri7/3"},
   {"op", "pri8"},
   {"stop", "pri8/3"},
   {"ep", "pri9"},
   {"step", "pri9/2"},
   {"gistep", "pri9/4"},
   {"dip", "pri10"},
   {"stiddip", "pri10/3"},

   {"squap", "ant4"},
   {"hap", "ant6"},
   {"heap", "ant7"},
   {"sthap", "ant7/2"},
   {"gisthap", "ant7/3"},
   {"gisthirp", "ant7/4"},
   {"oap", "ant8"},
   {"stoap", "ant8/3"},
   {"storp", "ant8/5"},
   {"eap", "ant9"},
   {"steap", "ant9/2"},
   {"gisteap", "ant9/4"},
   {"gisterp", "ant9/5"},
   {"dap", "ant10"},
   {"stiddap", "ant10/3"},

   // Bowers short names for Uniform Compounds
   {"sis", "uc1"},
   {"dis", "uc2"},
   {"snu", "uc3"},
   {"so", "uc4"},
   {"ki", "uc5"},
   {"e", "uc6"},
   {"risdoh", "uc7"},
   {"rah", "uc8"},
   {"rhom", "uc9"},
   {"dissit", "uc10"},
   {"doso", "uc11"},
   {"sno", "uc12"},
   {"addasi", "uc13"},
   {"dasi", "uc14"},
   {"gissi", "uc15"},
   {"si", "uc16"},
   {"se", "uc17"},
   {"hirki", "uc18"},
   {"sapisseri", "uc19"},
   {"gadsid", "uc26"},
   {"gassid", "uc27"},
   {"gidasid", "uc28"},
   {"gissed", "uc29"},
   {"ro", "uc30"},
   {"dro", "uc31"},
   {"kri", "uc32"},
   {"dri", "uc33"},
   {"kred", "uc34"},
   {"dird", "uc35"},
   {"gikrid", "uc36"},
   {"giddird", "uc37"},
   {"griso", "uc38"},
   {"rosi", "uc39"},
   {"rassid", "uc40"},
   {"grassid", "uc41"},
   {"gassic", "uc42"},
   {"gidsac", "uc43"},
   {"sassid", "uc44"},
   {"sadsid", "uc45"},
   {"siddo", "uc46"},
   {"sne", "uc47"},
   {"presipsido", "uc48"},
   {"presipsi", "uc49"},
   {"passipsido", "uc50"},
   {"passipsi", "uc51"},
   {"sirsido", "uc52"},
   {"sirsei", "uc53"},
   {"tisso", "uc54"},
   {"taki", "uc55"},
   {"te", "uc56"},
   {"harie", "uc57"},
   {"quahri", "uc58"},
   {"arie", "uc59"},
   {"gari", "uc60"},
   {"iddei", "uc61"},
   {"rasseri", "uc62"},
   {"rasher", "uc63"},
   {"rahrie", "uc64"},
   {"raquahri", "uc65"},
   {"rasquahr", "uc66"},
   {"rasquahpri", "uc67"},
   {"disco", "uc68"},
   {"dissid", "uc69"},
   {"giddasid", "uc70"},
   {"gidsid", "uc71"},
   {"gidrissid", "uc72"},
   {"disdid", "uc73"},
   {"idisdid", "uc74"},
   {"desided", "uc75"},

   // Bowers short names for Johnson Solids
   {"squippy", "j1"},
   {"peppy", "j2"},
   {"tricu", "j3"},
   {"squicu", "j4"},
   {"pecu", "j5"},
   {"pero", "j6"},
   {"etripy", "j7"},
   {"esquippy", "j8"},
   {"epeppy", "j9"},
   {"gyesp", "j10"},
   {"gyepip", "j11"},
   {"tridpy", "j12"},
   {"pedpy", "j13"},
   {"etidpy", "j14"},
   {"esquidpy", "j15"},
   {"epedpy", "j16"},
   {"gyesqidpy", "j17"},
   {"etcu", "j18"},
   {"escu", "j19"},
   {"epcu", "j20"},
   {"epro", "j21"},
   {"gyetcu", "j22"},
   {"gyescu", "j23"},
   {"gyepcu", "j24"},
   {"gyepro", "j25"},
   {"gybef", "j26"},
   {"tobcu", "j27"},
   {"squobcu", "j28"},
   {"squigybcu", "j29"},
   {"pobcu", "j30"},
   {"pegybcu", "j31"},
   {"pocuro", "j32"},
   {"pegycuro", "j33"},
   {"pobro", "j34"},
   {"etobcu", "j35"},
   {"etigybcu", "j36"},
   {"esquigybcu", "j37"},
   {"epobcu", "j38"},
   {"epigybcu", "j39"},
   {"epocuro", "j40"},
   {"epgycuro", "j41"},
   {"epobro", "j42"},
   {"epgybro", "j43"},
   {"gyetibcu", "j44"},
   {"gyesquibcu", "j45"},
   {"gyepibcu", "j46"},
   {"gyepcuro", "j47"},
   {"gyepabro", "j48"},
   {"autip", "j49"},
   {"bautip", "j50"},
   {"tautip", "j51"},
   {"aupip", "j52"},
   {"baupip", "j53"},
   {"auhip", "j54"},
   {"pabauhip", "j55"},
   {"mabauhip", "j56"},
   {"tauhip", "j57"},
   {"aud", "j58"},
   {"pabaud", "j59"},
   {"mabaud", "j60"},
   {"taud", "j61"},
   {"mibdi", "j62"},
   {"teddi", "j63"},
   {"auteddi", "j64"},
   {"autut", "j65"},
   {"autic", "j66"},
   {"bautic", "j67"},
   {"autid", "j68"},
   {"pabautid", "j69"},
   {"mabautid", "j70"},
   {"tautid", "j71"},
   {"gyrid", "j72"},
   {"pabgyrid", "j73"},
   {"mabgyrid", "j74"},
   {"tagyrid", "j75"},
   {"dirid", "j76"},
   {"pagydrid", "j77"},
   {"magydrid", "j78"},
   {"bagydrid", "j79"},
   {"pabidrid", "j80"},
   {"mabidrid", "j81"},
   {"gybadrid", "j82"},
   {"tedrid", "j83"},
   {"snadow", "j84"},
   {"snisquap", "j85"},
   {"waco", "j86"},
   {"auwaco", "j87"},
   {"wamco", "j88"},
   {"hawmco", "j89"},
   {"dawci", "j90"},
   {"bilbiro", "j91"},
   {"thawro", "j92"},
};

const char *u_abbrevs[][2] = {
   {"tr",     "truncated"},
   {"sm",     "small"},
   {"gr",     "great"},
   {"st",     "stellated"},
   {"sn",     "snub"},
   {"tet",    "tetrahedron"},
   {"ico",    "icosahedron"},
   {"icosa",  "icosahedron"},
   {"dod",    "dodecahedron"},
   {"oct",    "octahedron"},
   {"cubo",   "cuboctahedron"},
   {"icosid", "icosidodecahedron"}
};

const char *ud_abbrevs[][2] = {
   {"sm",      "small"},
   {"gr",      "great"},
   {"st",      "stellated"},
   {"inv",     "inverted"},
   {"delt",    "deltoidal"},
   {"pen",     "pentagonal"},
   {"med",     "medial"},
   {"triam",   "triambic"},
   {"ditrig",  "ditrigonal"},
   {"tri",     "triakis"},
   {"tetr",    "tetrakis"},
   {"pent",    "pentakis"},
   {"hex",     "hexakis"},
   {"disd",    "disdyakis"},
   {"tet",     "tetrahedron"},
   {"ico",     "icosahedron"},
   {"icosa",   "icosahedron"},
   {"dod",     "dodecahedron"},
   {"oct",     "octahedron"},
   {"cubo",    "cuboctahedron"},
   {"hexa",    "hexahedron"},
   {"hexec",   "hexecontahedron"},
   {"icositet","icositetrahedron"}
};

const char *uc_abbrevs[][2] = {
   {"tr",     "truncated"},
   {"sm",     "small"},
   {"gr",     "great"},
   {"st",     "stellated"},
   {"sn",     "snub"},
   {"tet",    "tetrahedra"},
   {"ico",    "icosahedra"},
   {"icosa",  "icosahedra"},
   {"dod",    "dodecahedra"},
   {"oct",    "octahedra"},
   {"cubo",   "cuboctahedra"},
   {"icosid", "icosidodecahedra"},
   {"pri",    "prisms"},
   {"ant",    "antiprisms"},
   {"rot",    "rotational"}
};

const char *j_abbrevs[][2] = {
   {"tri", "triangular"},
   {"sq",  "square"},
   {"squ", "square"},
   {"pe",  "pentagonal"},
   {"pen", "pentagonal"},
   {"el",  "elongated"},
   {"ge",  "gyroelongated"},
   {"tr",  "truncated"},
   {"au",  "augmented"},
   {"ba",  "biaugmenbted"},
   {"ta",  "triaugmented"},
};


int make_resource_uniform(geom_if &geom, string name, bool is_std, char *errmsg=0)
{
   if(name.size()<2 || !strchr("uUkKcCwW", name[0]) || name.find('.')!=string::npos)
      return -1; // not uniform name (the "." indicates a likely local file)
                 // so the name is not handled

   bool is_dual = strchr("dD", name[1]);
   uni_poly uni;
   int sym_no;
   if(read_int(name.c_str()+1+is_dual, &sym_no)) {
      // Uniform, Kaleido, Coxeter, or Wenninger number
      if (is_dual)
         name.erase(1,1);
      sym_no = uni.lookup_sym_no(name, is_dual);
      if(sym_no == -1) {
         if(errmsg)
            snprintf(errmsg, MSG_SZ, "invalid uniform polyhedron %sname",
                  (is_dual)?"dual ":"");
         return 1; // fail
      }
   }
   // if name starts with "u_" prefix
   else if(strchr(RES_SEPARATOR, name[1+is_dual])) {
      string expanded;
      if(is_dual)
         expanded = expand_abbrevs(name, ud_abbrevs,
               sizeof(ud_abbrevs)/sizeof(ud_abbrevs[0]));
      else
         expanded = expand_abbrevs(name, u_abbrevs,
               sizeof(u_abbrevs)/sizeof(u_abbrevs[0]));
      sym_no = uni.lookup_sym_no(expanded.c_str(), is_dual);
      if(sym_no == -1) {
         if(errmsg)
            snprintf(errmsg, MSG_SZ, "invalid uniform polyhedron %sname",
                  (is_dual)?"dual ":"");
         return 1; // fail
      }
   }
   else
      return -1; // not uniform name

   
   uni.get_poly(geom, sym_no);

   if (!is_std) {
      double e_len = geom.edge_vec(geom.faces(0,0), geom.faces(0,1)).mag();
      geom.transform(mat3d::scale(1/e_len));
      set_resource_polygon_color(geom);
   }

   // this is here if 'd' was the second character
   if(is_dual)
      make_resource_dual(geom, is_std);

   return 0; // name found
}


int make_resource_uniform_compound(geom_if &geom, string name, bool is_std, char *errmsg=0)
{
   // lower case the name
   transform(name.begin(), name.end(), name.begin(), ::tolower);

   if(name.size()<3 || name.substr(0,2)!="uc")
      return -1; // not uniform compound name

   double angle = INFINITY;
   int n = -1;
   int d = -1;
   int k = -1;

   uc_poly uniform_compounds;
      
   // if complex parms, parse them
   if(strpbrk(name.c_str(),RES_SEPARATOR)) {
      char errmsg2[MSG_SZ];
      uniform_compounds.parse_uc_args(name,angle,n,d,k,errmsg2);
   }

   int sym_no;
   if(read_int(name.c_str()+2, &sym_no)) {
      sym_no--;
      if(sym_no<0 || sym_no >= uniform_compounds.get_last_uc()) {
         if(errmsg)
            snprintf(errmsg, MSG_SZ, "uniform compound number out of range");
         return 1; // fail
      }
   }
   else if(strchr(RES_SEPARATOR, name[2])) {
      string expanded = expand_abbrevs(name, uc_abbrevs,
            sizeof(uc_abbrevs)/sizeof(uc_abbrevs[0]));
      sym_no = uniform_compounds.lookup_sym_no(expanded.c_str());
      if(sym_no == -1)
         return 1; // fail
   }
   else
      return -1; // not uniform compound name
     
   char errmsg2[MSG_SZ];
   int ret = uniform_compounds.set_uc_args(sym_no+1,angle,n,d,k,errmsg2);
   if (ret>0)
      return 1; // fail

   if (angle != INFINITY || n != -1 || d != -1 || k != -1) {
      char angle_str[MSG_SZ] = "";
      if (angle != INFINITY)
         snprintf(angle_str, MSG_SZ, "a%g",rad2deg(angle));
      char n_str[MSG_SZ] = "";
      if (n != -1)
         snprintf(n_str, MSG_SZ, "n%d",n);
      char d_str[MSG_SZ] = "";
      if (d > 1)
         snprintf(d_str, MSG_SZ, "/%d",d);
      char k_str[MSG_SZ] = "";
      if (k != -1)
         snprintf(k_str, MSG_SZ, "k%d",k);
      //fprintf(stderr,"UC%d_%s%s%s%s\n",sym_no+1,angle_str,n_str,d_str,k_str);
   }

   uniform_compounds.get_poly(geom,sym_no,angle,n,d,k,is_std);
   
   col_geom_v *cg = dynamic_cast<col_geom_v *>(&geom); 
   if(cg && !is_std) {      
      coloring clrng(cg);
      color_map *cmap = init_color_map("compound");
      clrng.add_cmap(cmap);

      clrng.v_apply_cmap();
      clrng.e_apply_cmap();
      clrng.f_apply_cmap();
   }

   return 0; // name found
}

int make_resource_johnson(geom_if &geom, string name, bool is_std, char *errmsg=0)
{
   if(name.size()<2 || !strchr("jJ", name[0]) || name.find('.')!=string::npos)
      return -1; // not johnson name (the "." indicates a likely local file)
                 // so the name is not handled

   j_poly json;
   int sym_no;
   if(read_int(name.c_str()+1, &sym_no)) {
      sym_no--;
      if(sym_no<0 || sym_no >= json.get_last_J()) {
         if(errmsg)
            snprintf(errmsg, MSG_SZ, "Johnson polyhedron number out of range");
         return 1; // fail
      }
   }
   else if(strchr(RES_SEPARATOR, name[1])) {
      string expanded = expand_abbrevs(name, j_abbrevs,
            sizeof(j_abbrevs)/sizeof(j_abbrevs[0]));
      sym_no = json.lookup_sym_no(expanded.c_str());
      if(sym_no == -1) {
         if(errmsg)
            snprintf(errmsg, MSG_SZ, "invalid Johnson polyhedron name");
         return 1; // fail
      }
   }
   else
      return -1; // not johnson name
   
   json.get_poly(geom, sym_no);

   if (!is_std)
      set_resource_polygon_color(geom);

   return 0; // name found
}

class res_coloring : public coloring
{
   public:
      res_coloring(col_geom_v *geo=0) : coloring(geo) {}
      void f_all_angles(bool as_values);
      void f_max_angle(bool as_values);
};


void res_coloring::f_all_angles(bool as_values)
{
   geom_info info(*get_geom());
   int faces_sz = get_geom()->faces().size();
   for(int i=0; i<faces_sz; i++) {
      vector<double> f_angs;
      info.face_angles_lengths(i, f_angs);
      sort(f_angs.rbegin(), f_angs.rend());
      double ang_sum = 0;
      for(unsigned int j=0; j<f_angs.size(); j++)
         ang_sum += j*f_angs[j];
      int idx = (int)(rad2deg(50*ang_sum/(f_angs.size()*(f_angs.size()+1)/2))+0.5);
      if(as_values)
         get_geom()->set_f_col(i, get_col(idx));
      else
         get_geom()->set_f_col(i, idx);
   }
}


int make_resource_geodesic(geom_if &geom, string name, bool is_std, char *errmsg)
{
   if(name.size()<5 || name.substr(0,4)!="geo_" || name.find('.')!=string::npos)
      return -1; // not geodesic name (the "." indicates a likely local file)
                 // so the name is not handled

   col_geom_v base;
   int offset = 4;
   if(name[offset]=='t') {
      base.read_resource("std_tet");
      offset++;
   }
   else if(name[offset]=='o') {
      base.read_resource("std_oct");
      offset++;
   }
   else if(name[offset]=='i') {
      base.read_resource("std_ico");
      offset++;
   }
   else if(name[offset]=='T') {
      double X = 0.25;
      double Y = sqrt(3.0)/12;
      double Z = 0.8;
      base.add_vert(vec3d(-X, -Y,  Z)); // 0
      base.add_vert(vec3d( X, -Y,  Z)); // 1
      base.add_vert(vec3d( 0,2*Y,  Z)); // 3

      vector<int> face(3);
      for(int i=0; i<3; i++)
            face[i] = i;
      base.add_face(face);
      offset++;
   }
   else if(isdigit(name[offset])) {
      base.read_resource("std_ico");
   }
   else {
      if(errmsg)
         snprintf(errmsg, MSG_SZ, "geodesic base polyhedron not i, o, t or T");
      return 1; // fail
   }

   char str[MSG_SZ];
   strncpy(str, name.c_str()+offset, MSG_SZ-1);
   str[MSG_SZ-1] = '\0';
   vector<char *> num_strs;
   split_line(str, num_strs, "_");
   int nums[2] = {1, 0};
   if(num_strs.size()>2) {
      if(errmsg)
         snprintf(errmsg, MSG_SZ, "geodesic division: must be given by 1 or 2 integers");
      return 1; // fail
   }

   if(num_strs.size()>0 && !read_int(num_strs[0], &nums[0])) {
      if(errmsg)
         snprintf(errmsg, MSG_SZ, "geodesic division:%s division number not an integer", (num_strs.size()>1) ? " first" : "");
      return 1; // fail
   }

   if(num_strs.size()>1 && !read_int(num_strs[1], &nums[1])) {
      if(errmsg)
         snprintf(errmsg, MSG_SZ, "geodesic division: second division number not an integer");
      return 1; // fail
   }

   if(!geom.set_geodesic_sphere(base, nums[0], nums[1])) {
      if(errmsg)
         snprintf(errmsg, MSG_SZ, "geodesic division: invalid division");
      return 1; // fail
   }

   col_geom_v *cg = dynamic_cast<col_geom_v *>(&geom); 
   if(cg && !is_std) {
      res_coloring clrng(cg);
      color_map *cmap = init_color_map("rng256_R1:0.7:1:0.5:1G1:0.4:1:0.3:1B1:0.2:1:0.7:1%");
      clrng.add_cmap(cmap);
      clrng.f_all_angles(true);
   }
   
   return 0; // name found
}

// if name ends in "_d" then make dual of resource 
void process_for_dual(string &name, bool &dual_flag)
{
   int name_sz = name.size();
   if(name_sz>2 && name[name_sz-2]=='_' && name[name_sz-1]=='d') {
      name.resize(name_sz-2);
      dual_flag = !dual_flag;
   }
}


static void get_arrow(col_geom_v &geom, const sch_sym &sym)
{
   geom.clear_all();
   geom.add_vert(vec3d(0.0, 0.0, 0.0));
   geom.add_vert(vec3d(3.0, 0.0 ,0.0));
   geom.add_vert(vec3d(2.0, 1.0, 0.0));
   geom.add_vert(vec3d(2.0, 0.5, 0.0));
   geom.add_vert(vec3d(0.0, 0.5, 0.0));
   geom.add_vert(vec3d(0.0, 0.0, 0.5));
   geom.add_vert(vec3d(3.0, 0.0, 0.5));
   geom.add_vert(vec3d(2.0, 0.5, 0.25));
   geom.add_vert(vec3d(0.0, 0.5, 0.25));
   geom.set_f_col(geom.add_face(1, 2, 3, 4, 0, -1), col_val(0.9, 0.9, 0.7));
   geom.set_f_col(geom.add_face(6, 2, 7, 8, 5, -1), col_val(0.0, 0.2, 0.6));
   geom.set_f_col(geom.add_face(0, 1, 6, 5, -1), col_val(0.9, 0.9, 0.7));
   geom.set_f_col(geom.add_face(1, 2, 6, -1), col_val(0.9, 0.9, 0.7));
   geom.set_f_col(geom.add_face(2, 3, 7, -1), col_val(0.9, 0.9, 0.7));
   geom.set_f_col(geom.add_face(3, 4, 8, 7, -1), col_val(0.9, 0.9, 0.7));
   geom.set_f_col(geom.add_face(4, 0, 5, 8, -1), col_val(0.9, 0.9, 0.7));

   mat3d trans = mat3d::transl(vec3d(0.1, 0.2, 0.4));
   int fold = sym.get_nfold();

   switch(sym.get_sym_type()) {
      case sch_sym::C1:
      case sch_sym::Cs:
      case sch_sym::Ci:
         break;

      case sch_sym::C:
      case sch_sym::Ch:
      case sch_sym::S:
         trans = mat3d::rot(vec3d::Z, -M_PI/2) *
                 mat3d::transl(vec3d(-1.6,-0.2,0.0)) * trans;
         if(fold>1)
            trans = mat3d::transl(vec3d(1.6/sin(M_PI/fold), 0, 0)) * trans;
         break;

      case sch_sym::Cv:
      case sch_sym::D:
      case sch_sym::Dh:
         trans = mat3d::rot(vec3d::Z, -M_PI/2) * trans;
         if(fold>1)
            trans = mat3d::transl(vec3d(3.2/tan(M_PI/fold), 0, 0)) * trans;
         break;

      case sch_sym::Dv:
         if(fold>1)
            trans = mat3d::rot(vec3d::Z, 0.5*M_PI/fold) *
                    mat3d::transl(vec3d(3.2/tan(M_PI/fold), 0, 0)) *
                    mat3d::rot(vec3d::Z, -M_PI/2) * trans;
         break;

      case sch_sym::T:
      case sch_sym::Td:
      case sch_sym::Th:
         trans = mat3d::alignment(vec3d::Z, vec3d::X,
                                  vec3d(1,1,1), vec3d(0,-1,1)) *
                 mat3d::transl(vec3d(0, 0.5, 3)) *
                 mat3d::rot(vec3d::Z, M_PI/2) * trans;
         break;

      case sch_sym::O:
      case sch_sym::Oh:
         //trans = mat3d::transl(vec3d(0, 3, 4)) * trans;
         trans = mat3d::transl(vec3d(1, 0, 4)) * trans;
         break;

      case sch_sym::I:
      case sch_sym::Ih:
         trans = mat3d::rot(vec3d::Z, vec3d(0,1,1.618)) *
                 mat3d::transl(vec3d(0, 1, 5.5)) *
                 mat3d::rot(vec3d::Z, M_PI/2) * trans;
         break;

   }

   geom.transform(trans);
}


int make_resource_sym(geom_if &geom, string name, char *errmsg=0) {
   if(name.size()<5 || name.substr(0,4)!="sym_" || name.find('.')!=string::npos)
      return -1; // not sym_ name (the "." indicates a likely local file)
                 // so the name is not handled
      
   sch_sym sym;
   if(!sym.init(name.substr(4), mat3d(), errmsg))
      return 1;

   col_geom_v arrow;
   get_arrow(arrow, sym);
   sym_repeat(geom, arrow, sym.get_trans());
   return 0;
}


int make_resource_pgon(geom_if &geom, string name, bool is_std, char *errmsg)
{
   if(name.find('.')!=string::npos)
      return -1; // not polygon res name (the "." indicates a likely local file)
                 // so the name is not handled

   char pnam[MSG_SZ];
   strncpy(pnam, name.c_str(), MSG_SZ);
   int num_sides;
   int step=1;
   char *pnum = pnam+3;
   char *p = strchr(pnum, '/');
   if(p!=0) {
      *p++='\0';
      if(!read_int(p, &step))
         return -1;
   }
   if(!read_int(pnum, &num_sides))
      return -1;
   if(num_sides<2)
      return -1;
   if(step<1)
      return -1;
   if(step%num_sides==0)
      return -1;
   
   *pnum = '\0';

   polygon pgon(num_sides, step);
   polygon *poly;
   if(strcasecmp("pri", pnam)==0)
      poly = new prism(pgon);
   else if(strcasecmp("ant", pnam)==0)
      poly = new antiprism(pgon);
   else if(strcasecmp("pyr", pnam)==0)
      poly = new pyramid(pgon);
   else if(strcasecmp("dip", pnam)==0)
      poly = new dipyramid(pgon);
   else if(strcasecmp("cup", pnam)==0)
      poly = new cupola(pgon);
   else if(strcasecmp("ort", pnam)==0)
      poly = new orthobicupola(pgon);
   else if(strcasecmp("gyr", pnam)==0)
      poly = new gyrobicupola(pgon);
   else if(strcasecmp("snu", pnam)==0)
      poly = new snub_antiprism(pgon);
   else
      return -1;

   // check can have unit edge
   int ret;
   poly->set_edge(1);
   if(!poly->set_edge2(1)) {
      if(errmsg)
         strcpy(errmsg, "polyhedron cannot have unit edges");
      ret = 1;
   }
   else {
      poly->make_poly(geom);
      ret = 0;
   }

   delete poly;

   if (!is_std)
      set_resource_polygon_color(geom);

   return ret;
}

int make_resource_schwarz(geom_if &geom, string name, bool is_std,
      char *errmsg=0)
{
   if(name.size()<(7+5) || name.substr(0,7)!="schwarz" ||
         name.find('.')!=string::npos)
      return -1; // not schwarz name (the "." indicates a likely local file)
                 // so the name is not handled

   if(name.find(':')!=string::npos || name.find('|')!=string::npos) {
      if(errmsg)
         snprintf(errmsg, MSG_SZ, "schwarz name: cannot contain bar character");
      return 1; // fail
   }

   string symbol_str = name.substr(7);
   for(unsigned int i=0; i<symbol_str.size(); i++) {
      if(symbol_str[i] == '_')      // allow _ in place of space
         symbol_str[i] = ' ';
   }

   bool as_poly = false;
   if(symbol_str.size() && symbol_str[symbol_str.size()-1]=='p') {
      as_poly = true;
      symbol_str.resize(symbol_str.size()-1);
   }

   char errmsg2[MSG_SZ];
   symbol_str += "|"; // add to make sure symbol will be valid for wythoff
   wythoff_poly wyt(symbol_str.c_str(), errmsg2);
   if(!wyt.is_set()) {
      if(errmsg)
         snprintf(errmsg, MSG_SZ, "schwarz name: invalid symbol: %s",
               errmsg2);
      return 1; // fail
   }

   if(as_poly)
      wyt.make_tri_poly(geom);
   else
      wyt.make_tri(geom);

   if(!is_std) {
      col_geom_v *cg = dynamic_cast<col_geom_v *>(&geom); 
      if(cg) {
         res_coloring clrng(cg);
         color_map *cmap = init_color_map("map_grey20:ivory");
         clrng.add_cmap(cmap);
         clrng.f_apply_cmap();
      }
   }

   return 0; // name found
}


int make_resource_wythoff(geom_if &geom, string name, bool is_std,
      char *errmsg=0)
{
   if(name.size()<(7+5) || name.substr(0,7)!="wythoff" ||
         name.find('.')!=string::npos)
      return -1; // not wythoff name (the "." indicates a likely local file)
                 // so the name is not handled

   string symbol_str = name.substr(7);
   for(unsigned int i=0; i<symbol_str.size(); i++) {
      if(symbol_str[i] == '_')     // allow _ in place of space
         symbol_str[i] = ' ';
      else if(symbol_str[i] == ':')     // allow : in place of |
         symbol_str[i] = '|';
   }

   char errmsg2[MSG_SZ];
   wythoff_poly wyt(symbol_str.c_str(), errmsg2);
   if(!wyt.is_set()) {
      if(errmsg)
         snprintf(errmsg, MSG_SZ, "wythoff name: invalid symbol: %s",
               errmsg2);
      return 1; // fail
   }

   wyt.make_poly(geom, errmsg2);
   if(*errmsg2)
      fprintf(stderr, "warning: wythoff name: %s\n", errmsg2);

   if(!is_std) {
      double e_len = geom.edge_vec(geom.faces(0,0), geom.faces(0,1)).mag();
      geom.transform(mat3d::scale(1/e_len));
      col_geom_v *cg = dynamic_cast<col_geom_v *>(&geom); 
      if(cg) {
         res_coloring clrng(cg);
         color_map *cmap = init_color_map("map_red:blue:yellow:ivory");
         clrng.add_cmap(cmap);
         clrng.f_apply_cmap();
      }
   }

   return 0; // name found
}


bool make_resource_geom(geom_if &geom, string name, char *errmsg)
{
   if(errmsg)
      *errmsg = '\0';

   geom.clear_all();
   
   if(!name.size())
      return false;

   bool make_dual = false;
   process_for_dual(name, make_dual);

   char errmsg2[MSG_SZ];
   bool geom_ok = false;

   // RK: now check for standard polys before any work is done on the name
   // if found, then end here. if not fall through
   if(!geom_ok) {
      int ret = make_resource_std_poly(geom, name, errmsg2);
      if(ret==0) {
         if(make_dual)
            make_resource_dual(geom, true);
         return true;
      }
   }

   // RK: standard keyword "std_" now used for all built-in polyhedra
   // truncate "std_" but keep track of it
   bool is_std = (name.size()>3 && name.substr(0,4)=="std_");
   if (is_std)
      name = name.substr(4);
   
   // Look for an internal alternative name
   char alt_name[MSG_SZ];
   to_resource_name(alt_name, name.c_str());

   for(unsigned int i=0; i<sizeof(alt_names)/sizeof(alt_names[0]); i++) {
      if(strcmp(alt_name, alt_names[i][0])==0) {
         name = alt_names[i][1];  // set name to the usual name for the model
         process_for_dual(name, make_dual);
         break;
      }
   }
  
   if(!geom_ok) {
      int ret = make_resource_pgon(geom, name, is_std, errmsg2);
      if(ret==0)
         geom_ok = true;
      else if(ret > 0) {
         if(errmsg)
            strcpy(errmsg, errmsg2);
         return false;
      }
   }

   if(!geom_ok) {
      int ret = make_resource_uniform(geom, name, is_std, errmsg2);
      if(ret==0)
         geom_ok = true;
      else if(ret > 0) {
         if(errmsg)
            strcpy(errmsg, errmsg2);
         return false;
      }
   }

   if(!geom_ok) {
      int ret = make_resource_uniform_compound(geom, name, is_std, errmsg2);
      if(ret==0)
         geom_ok = true;
      else if(ret > 0) {
         if(errmsg)
            strcpy(errmsg, errmsg2);
         return false;
      }
   }

   if(!geom_ok) {
      int ret = make_resource_johnson(geom, name, is_std, errmsg2);
      if(ret==0)
         geom_ok = true;
      else if(ret > 0) {
         if(errmsg)
            strcpy(errmsg, errmsg2);
         return false;
      }
   }

   if(!geom_ok) {
      int ret = make_resource_geodesic(geom, name, is_std, errmsg2);
      if(ret==0)
         geom_ok = true;
      else if(ret > 0) {
         if(errmsg)
            strcpy(errmsg, errmsg2);
         return false;
      }
   }

   if(!geom_ok) {
      int ret = make_resource_wythoff(geom, name, is_std, errmsg2);
      if(ret==0)
         geom_ok = true;
      else if(ret > 0) {
         if(errmsg)
            strcpy(errmsg, errmsg2);
         return false;
      }
   }

   if(!geom_ok) {
      int ret = make_resource_schwarz(geom, name, is_std, errmsg2);
      if(ret==0)
         geom_ok = true;
      else if(ret > 0) {
         if(errmsg)
            strcpy(errmsg, errmsg2);
         return false;
      }
   }

   if(!geom_ok) {
      int ret = make_resource_sym(geom, name, errmsg2);
      if(ret==0)
         geom_ok = true;
      else if(ret > 0) {
         if(errmsg)
            strcpy(errmsg, errmsg2);
         return false;
      }
   }

   if(make_dual)
      make_resource_dual(geom, is_std);

   return geom_ok;
}

