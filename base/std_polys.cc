/*
   Copyright (c) 2003-2021, Adrian Rossiter

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

#include "coloring.h"
#include "geometryinfo.h"
#include "geometryutils.h"
#include "mathutils.h"
#include "polygon.h"
#include "private_std_polys.h"
#include "utils.h"

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <regex>
#include <set>
#include <string>
#include <unordered_map>
#include <vector>

using std::map;
using std::set;
using std::string;
using std::vector;

using namespace ::anti;

/// A basic function type that makes a model that has a normal form
// (e.g. with unit edges) and a standard form (e.g with commonly used
// coordinates)
typedef void (*std_model_func)(Geometry &, bool is_std);

void normalized_face_list(Geometry &geom)
{
  geom.orient(1);
  sort(geom.raw_faces().begin(), geom.raw_faces().end());
}

namespace {

/// Convert to a normalised resource name
/** Remove leading and trailing space, convert any whitespace to a
 *  single space
 * \param from the string to convert
 * \return the converted string. */
string to_resource_name(string from)
{
  string res_name(from);
  clear_extra_whitespace(res_name);
  return res_name;
}

string expand_abbrevs(const string &name, const char *abbrevs[][2], size_t last)
{
  string expanded;
  Split parts(name, RES_SEPARATOR);
  for (unsigned int i = 1; i < parts.size(); i++) {
    size_t j;
    for (j = 0; j < last; j++) {
      if (strcmp(parts[i], abbrevs[j][0]) == 0)
        break;
    }
    if (i > 1)
      expanded += ' ';
    expanded += (j < last) ? abbrevs[j][1] : parts[i];
  }
  return expanded;
}

static void set_resource_polygon_color(Geometry &geom)
{
  Coloring clrng(&geom);
  ColorMap *cmap = colormap_from_name("uniform");
  clrng.add_cmap(cmap);
  clrng.f_avg_angle(true);
  clrng.v_avg_angle(true);
}

void set_resource_unique_color(Geometry &geom)
{
  Coloring clrng(&geom);
  ColorMap *cmap = colormap_from_name("spread");
  clrng.add_cmap(cmap);
  clrng.f_unique(true);
  clrng.v_unique(true);
}

static void make_resource_dual(Geometry &geom, bool is_std = false)
{
  Vec3d cent = geom.centroid();
  GeometryInfo info(geom);
  info.set_center(cent);
  double rad = info.iedge_dist_lims().sum / info.num_iedges();
  Geometry dual;
  const double inf = 1200;
  get_dual(dual, geom, rad, cent, inf);
  add_extra_ideal_elems(dual, cent, 0.95 * inf); // limit closer than inf
  geom.clear_all();
  geom.append(dual);
  if (!is_std)
    set_resource_polygon_color(geom);
}

// clang-format off
std::unordered_map<string, string> alt_names = {
    {"tet", "u1"},
    {"tetrahedron", "u1"},
    {"truncated_tetrahedron", "u2"},
    {"tr_tetrahedron", "u2"},
    {"tr_tet", "u2"},
    {"octahedron", "u5"},
    {"oct", "u5"},
    {"cube", "u6"},
    {"cuboctahedron", "u7"},
    {"cubo", "u7"},
    {"truncated_octahedron", "u8"},
    {"tr_octahedron", "u8"},
    {"tr_oct", "u8"},
    {"truncated_cube", "u9"},
    {"tr_cube", "u9"},
    {"rhombicuboctahedron", "u10"},
    {"rhombicubo", "u10"},
    {"rh_cubo", "u10"},
    {"truncated_cuboctahedron", "u11"},
    {"tr_cuboctahedron", "u11"},
    {"tr_cubo", "u11"},
    {"snub_cube", "u12"},
    {"sn_cube", "u12"},
    {"icosahedron", "u22"},
    {"ico", "u22"},
    {"icosa", "u22"},
    {"dodecahedron", "u23"},
    {"dod", "u23"},
    {"icosidodecahedron", "u24"},
    {"icosid", "u24"},
    {"truncated_icosahedron", "u25"},
    {"tr_icosahedron", "u25"},
    {"tr_icosa", "u25"},
    {"tr_ico", "u25"},
    {"truncated_dodecahedron", "u26"},
    {"tr_dodecahedron", "u26"},
    {"tr_dod", "u26"},
    {"rhombicosidodecahedron", "u27"},
    {"rhombicosid", "u27"},
    {"rh_icosid", "u27"},
    {"truncated_icosidodecahedron", "u28"},
    {"tr_icosidodecahedron", "u28"},
    {"tr_icosid", "u28"},
    {"snub_dodecahedron", "u29"},
    {"snub_dod", "u29"},
    {"sn_dod", "u29"},
    {"sn_dodecahedron", "u29"},
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

    // Uniform duals
    {"triakistetrahedron", "u2_d"},
    {"triakis_tetrahedron", "u2_d"},
    {"tri_tet", "u2_d"},
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
    {"triakisicosahedron", "u26_d"},
    {"triakis_icosahedron", "u26_d"},
    {"tri_ico", "u26_d"},
    {"pentakisdodecahedron", "u25_d"},
    {"pentakis_dodecahedron", "u25_d"},
    {"pent_dod", "u25_d"},
    {"deltoidal_hexecontahedron", "u27_d"},
    {"delt_hexec", "u27_d"},
    {"trapezoidal_hexecontahedron", "u27_d"},
    {"trap_hexec", "u27_d"},
    {"hexakisicosahedron", "u28_d"},
    {"hexakis_icosahedron", "u28_d"},
    {"hex_ico", "u28_d"},
    {"disdyakistriacontahedron", "u28_d"},
    {"disdyakis_triacontahedron", "u28_d"},
    {"disd_tri", "u28_d"},
    {"pentagonal_hexecontahedron", "u29_d"},
    {"pen_hexec", "u29_d"},

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

// alternative names for uniform polyhedra (for names following u_)
// alternative name is lowercase, no abbreviations, words separated by spaces
std::unordered_map<string, string> u_alt_names = {
    // Jim McNeill http://www.orchidpalms.com/polyhedra/uniform/uniform.html
    {"snub cuboctahedron","u12"},
    {"small rhombicube","u18"},
    {"stellatruncated cube","u19"},
    {"great rhombicube","u21"},
    {"snub icosidodecahedron","u29"},
    {"snub disicosidodecahedron","u32"},
    {"dodecadodecahedron","u36"},
    {"great truncated dodecahedron","u37"},
    {"ditrigonary dodecadodecahedron","u41"},
    {"great dodekified icosidodecahedron","u42"},
    {"small dodekified icosidodecahedron","u43"},
    {"icosified dodecadodecahedron","u44"},
    {"great ditrigonary icosidodecahedron","u47"},
    {"great icosified icosidodecahedron","u48"},
    {"great truncated icosahedron","u55"},
    {"small stellatruncated dodecahedron","u58"},
    {"stellatruncated dodecadodecahedron","u59"},
    {"vertisnub dodecadodecahedron","u60"},
    {"great dodekicosahedron","u63"},
    {"great snub icosidisdodecahedron","u64"},
    {"great stellatruncated dodecahedron","u66"},
    {"stellatruncated icosidodecahedron","u68"},
    {"great vertisnub icosidodecahedron","u69"},
    {"great disnub disicosidisdodecahedron","u75"},
    // G Hart https://www.georgehart.com/virtual-polyhedra/uniform-index.html
    {"great dodecadodecahedron", "u36"},
    {"pentagrammatic prism","u78"},
    {"pentagrammatic antiprism","u79"},
    {"pentagrammatic crossed antiprism","u80"},
    // other alternate names found in Wikipedia
    {"triangular pyramid","u1"},
    {"cantic cube","u2"},
    {"allelotetratetrahedron","u3"},
    {"hemicuboctahedron","u4"},
    {"hemihexahedron","u4"},
    {"truncated hexahedron","u9"},
    {"small rhombicuboctahedron","u10"},
    {"convex great rhombicuboctahedron","u11"},
    {"rhombitruncated cuboctahedron","u11"},
    {"omnitruncated cube","u11"},
    {"cantitruncated cube","u11"},
    {"beveled cube","u11"},
    {"cubus simus","u12"},
    {"cuboctatruncated cuboctahedron","u16"},
    {"nonconvex great rhombicuboctahedron","u17"},
    {"quasirhombicuboctahedron","u17"},
    {"quasitruncated hexahedron","u19"},
    {"quasitruncated cuboctahedron","u20"},
    {"stellatruncated cuboctahedron","u20"},
    {"convex great rhombicosidodecahedron","u28"},
    {"rhombitruncated icosidodecahedron","u28"},
    {"great rhombicosidodecahedron","u28"},
    {"omnitruncated dodecahedron","u28"},
    {"omnitruncated icosahedron","u28"},
    {"small ditrigonary icosidodecahedron","u30"},
    {"small triambic icosidodecahedron","u30"},
    {"small icosified icosidodecahedron","u31"},
    {"holosnub icosahedron","u32"},
    {"small dodekicosidodecahedron","u33"},
    {"triambic dodecadodecahedron","u41"},
    {"icosidodecatruncated icosidodecahedron","u45"},
    {"great triambic icosidodecahedron","u47"},
    {"small icosahemidodecahedron","u49"},
    {"truncated great icosahedron","u55"}, // lists as title name!
    {"quasitruncated small stellated dodecahedron","u58"},
    {"quasitruncated dodecadodecahedron","u59"},
    {"great dodecahemiicosahedron","u62"},
    {"great snub dodekicosidodecahedron","u64"},
    {"small dodecahemiicosahedron","u65"},
    {"quasitruncated great stellated dodecahedron","u66"},
    {"nonconvex great rhombicosidodecahedron","u67"},
    {"quasirhombicosidodecahedron","u67"},
    {"great quasitruncated icosidodecahedron","u68"},
    {"great icosahemidodecahedron","u71"},
    {"retrosnub disicosidodecahedron","u72"},
    {"retrosnub ditrigonary icosidodecahedron","u72"},
    {"retroholosnub icosahedron","u72"},
    {"small inverted retrosnub icosicosidodecahedron","u72"},
    {"great inverted retrosnub icosidodecahedron","u74"},
    {"millers monster","u75"},
    // George Olshevsky Monster Names
    {"hastur","u32"},
    {"chaugnar faugn","u40"},
    {"dagon","u46"},
    {"tsathoggua","u57"},
    {"nyarlathotep","u60"},
    {"shub niggurath","u64"},
    {"cthulhu","u69"},
    {"yog sothoth","u72"},
    {"azathoth","u74"},
};

// alternative names for uniform polyhedra duals (for names following ud_)
// alternative name is lowercase, no abbreviations, words separated by spaces
std::unordered_map<string, string> ud_alt_names = {
    // duals
    // G Hart https://www.georgehart.com/virtual-polyhedra/uniform-index.html
    {"triakistetrahedron","u2"},
    {"tetrakishexahedron","u8"},
    {"triakisoctahedron","u9"},
    {"trapezoidal icositetrahedron","u10"},
    {"disdyakisdodecahedron","u11"},
    {"tetradyakishexahedron","u16"},
    {"great trapezoidal icositetrahedron","u17"},
    {"great triakisoctahedron","u19"},
    {"great disdyakisdodecahedron","u20"},
    {"pentakisdodecahedron","u25"},
    {"triakisicosahedron","u26"},
    {"trapezoidal hexecontahedron","u27"},
    {"disdyakistriacontahedron","u28"},
    {"hexagonal hexecontahedron","u32"},
    {"small stellapentakisdodecahedron","u37"},
    {"medial trapezoidal hexecontahedron","u38"},
    {"great dodecacronic hexecontahedron","u42"},
    {"small dodecacronic hexecontahedron","u43"},
    {"tridyakisicosahedron","u45"},
    {"great stellapentakisdodecahedron","u55"},
    {"great pentakisdodekahedron","u58"}, // legacy mispelling
    {"great pentakisdodecahedron","u58"},
    {"medial disdyakistriacontahedron","u59"},
    {"great dodecahemiicosacron","u62"},
    {"small dodecahemiicosacron","u65"},
    {"great triakisicosahedron","u66"},
    {"great trapezoidal hexecontahedron","u67"},
    {"great disdyakistriacontahedron","u68"},
    {"small hexagrammatic hexecontahedron","u72"},
    {"great pentagrammatic hexecontahedron","u74"},
    {"pentagonal trapezohedron","u77"},
    {"pentagrammatic dipyramid","u78"},
    {"pentagrammatic trapezohedron","u79"},
    {"pentagrammatic concave trapezohedron","u80"},
    // Wenninger names via Wikipedia
    {"hexahedron","u5"},
    // other alternate names found in Wikipedia
    {"kistetrahedron","u2"},
    {"tetrahexahedron","u8"},
    {"hextetrahedron","u8"},
    {"tetrakis cube","u8"},
    {"kiscube","u8"},
    {"disdyakis hexahedron","u8"},
    {"hexakis tetrahedron","u8"},
    {"kisoctahedron","u9"},
    {"trisoctahedron","u9"},
    {"trigonal trisoctahedron","u9"},
    {"tetragonal icosikaitetrahedron","u10"},
    {"tetragonal trisoctahedron","u10"},
    {"strombic icositetrahedron","u10"},
    {"lanceolar icositetrahedron","u10"},
    {"lanceolar disdodecahedron","u10"},
    {"hexoctahedron","u11"},
    {"hexakis octahedron","u11"},
    {"octakis cube","u11"},
    {"octakis hexahedron","u11"},
    {"kisrhombic dodecahedron","u11"},
    {"pentagonal icosikaitetrahedron","u12"},
    {"great sagittal disdodecahedron","u17"},
    {"great strombic icositetrahedron","u17"},
    {"small dipteral disdodecahedron","u18"},
    {"great dipteral disdodecahedron","u21"},
    {"triacontahedron","u24"},
    {"kisdodecahedron","u25"},
    {"kisicosahedron","u26"},
    {"strombic hexecontahedron","u27"},
    {"tetragonal hexecontahedron","u27"},
    {"hexakis icosahedron","u28"},
    {"decakis dodecahedron","u28"},
    {"kisrhombic triacontahedron","u28"},
    {"small lanceal trisicosahedron","u31"},
    {"small sagittal ditriacontahedron","u33"},
    {"midly rhombic triacontahedron","u36"},
    {"small stellated triacontahedron","u36"},
    {"midly triambic icosahedron","u41"},
    {"great lanceal trisicosahedron","u42"},
    {"fat star","u43"},
    {"midly sagittal ditriacontahedron","u44"},
    {"midly dentoid ditriacontahedron","u46"},
    {"medial triambic icosahedron","u47"},
    {"great sagittal trisicosahedron","u48"},
    {"small dipteral trisicosahedron","u50"},
    {"great stellated triacontahedron", "u54"},
    {"great astropentakis dodecahedron","u55"},
    {"midly dipteral ditriacontahedron","u56"},
    {"great petaloid ditriacontahedron","u57"},
    {"midly petaloid ditriacontahedron","u60"},
    {"great lanceal ditriacontahedron","u61"},
    {"great dipteral trisicosahedron","u63"},
    {"great astroid ditriacontahedron","u64"},
    {"great sagittal ditriacontahedron","u67"},
    {"great sagittal hexecontahedron","u67"},
    {"great strombic hexecontahedron","u67"},
    {"trisdyakis icosahedron","u68"},
    {"petaloidal trisicosahedron","u69"},
    {"great dipteral ditriacontahedron","u73"},
    {"great dentoid ditriacontahedron","u74"},
    {"pentagrammic trapezohedron","u79"},
    {"pentagrammic concave trapezohedron","u80"},
};

// alternative names for johnson solids (for names following j_)
// alternative name is lowercase, no abbreviations, words separated by spaces
std::unordered_map<string, string> j_alt_names = {
    {"lesser dome","j4"},
    {"diminished icosahedron","j11"},
    {"triangular dipyramid","j12"},
    {"hexadeltahedron","j12"},
    {"pentagonal dipyramid","j13"},
    {"decadeltahedron","j13"},
    {"elongated triangular dipyramid","j14"},
    {"triakis triangular prism","j14"},
    {"elongated octahedron","j15"},
    {"pencil cube","j15"},
    {"12 faced pencil cube","j15"},
    {"pentakis pentagonal prism","j16"},
    {"heccaidecadeltahedron","j17"},
    {"tetrakis square antiprism","j17"},
    {"hexakaidecadeltahedron","j17"},
    {"anticuboctahedron","j27"},
    {"twisted cuboctahedron","j27"},
    {"disheptahedron","j27"},
    {"cantellated triangular prism","j35"},
    {"millers solid","j37"},
    {"pseudo rhombicuboctahedron","j37"},
    {"gyrate rhombicuboctahedron","j37"},
    {"cantellated pentagonal prism","j38"},
    {"monolaterotruncated triangular bipyramid","j49"},
    {"tetracaidecadeltahedron,","j51"},
    {"tetrakis triangular prism","j51"},
    {"siamese dodecahedron","j84"},
    {"triangular dodecahedron","j84"},
    {"trigonal dodecahedron","j84"},
    {"dodecadeltahedron","j84"},
    {"pentakis elongated gyrobifastigium","j90"}
};

// alternative names for johnson solids duals (for names following jd_)
// alternative name is lowercase, no abbreviations, words separated by spaces
std::unordered_map<string, string> jd_alt_names = {
    {"triangular prism","j12"},
    {"pentagonal prism","j13"},
    {"triangular bifrustum","j14"},
    {"square bifrustum","j15"},
    {"pentagonal bifrustum","j16"},
    {"elongated tetragonal disphenoid","j26"},
    {"schmitt conway danzer biprism","j26"},
    {"trapezo rhombic dodecahedron","j27"},
    {"elongated square trapezohedron","j29"},
    {"elongated pentagonal trapezohedron","j31"},
    {"trapezo rhombic triacontahedron","j34"},
    {"pseudo deltoidal icositetrahedron","j37"},
    {"truncated triangular bipyramid","j51"},
    {"associahedron k5","j51"},
    {"monolaterotruncated pentagonal bipyramid","j52"},
    {"parabilaterotruncated pentagonal bipyramid","j53"},
    {"monolaterotruncated hexagonal bipyramid","j54"},
    {"parabilaterotruncated hexagonal bipyramid","j55"},
    {"alternate truncated hexagonal bipyramid","j57"},
    {"gyroelongated pentagonal bifrustum","j58"},
    {"elongated gyrobifastigium","j84"},
    {"truncated snub disphenoid","j90"}
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
   {"hexac",   "hexecontahedron"},
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
// clang-format on

bool process_letter_dual_prefix(string &name, const char *letters,
                                bool &is_dual)
{
  // Prefix format
  // letter from list optionally followed by d/D, and no mixed case
  // on success
  //    if second character in name is d/D, remove it
  //    is_dual is set to indicate if second character was d/D

  if (name.size() < 2) // minimum length for a model name
    return false;

  if (!strchr(letters, name[0])) // does name start with a char from letters
    return false;

  is_dual = strchr("dD", name[1]); // is the second character d/D
  if (is_dual) {
    // is dual model name too short, or has mixed case prefix
    if (name.size() < 3 || isupper(name[0]) != isupper(name[1]))
      return false;
    name.erase(1, 1); // remove second character (d/D) from name
  }

  return true;
}

int make_resource_uniform(Geometry &geom, string name, bool is_std,
                          string *error_msg = nullptr)
{
  // Model name formats:
  //   uNumber,  UNumber,  u_Name
  //   udNumber, UDNumber, ud_Name
  //   kNumber,  KNumber,  k_Name
  //   kdNumber, KDNumber, kd_Name
  //   wNumber,  WNumber,  w_name
  //   wdNumber, WDNumber, wd_name

  bool is_dual;
  if (!process_letter_dual_prefix(name, "uUkKwW", is_dual))
    return -1; // unhandled, not uniform name

  Uniform uni;
  int sym_no;
  if (read_int(name.c_str() + 1, &sym_no)) {
    // bypass for wenninger stellations, found in next section
    // still need the dD
    if (strchr("wW", name[0]) &&
        (sym_no == 19 || (sym_no >= 23 && sym_no <= 40) ||
         (sym_no >= 42 && sym_no <= 66)))
      return -1;

    // Uniform, Kaleido, Coxeter, or Wenninger number
    sym_no = uni.lookup_sym_no(name, is_dual);
    if (sym_no == -1) {
      if (error_msg)
        *error_msg = msg_str("invalid uniform polyhedron %sname",
                             (is_dual) ? "dual " : "");
      return 1; // fail
    }
  }
  else if (strchr(RES_SEPARATOR, name[1])) {
    if (isupper(name[0]))
      return -1; // unhandled, upper case prefix, not UC full name

    string expanded;

    if (is_dual) {
      expanded = expand_abbrevs(name, ud_abbrevs,
                                sizeof(ud_abbrevs) / sizeof(ud_abbrevs[0]));

      // check alternative names
      auto it = ud_alt_names.find(to_resource_name(expanded));
      if (it != ud_alt_names.end())
        expanded = it->second; // set name to the usual name for the model
    }
    else {
      expanded = expand_abbrevs(name, u_abbrevs,
                                sizeof(u_abbrevs) / sizeof(u_abbrevs[0]));

      // check alternative names
      auto it = u_alt_names.find(to_resource_name(expanded));
      if (it != u_alt_names.end())
        expanded = it->second; // set name to the usual name for the model
    }
    sym_no = uni.lookup_sym_no(expanded.c_str(), is_dual);
    if (sym_no == -1) {
      if (error_msg)
        *error_msg = msg_str("invalid uniform polyhedron %sname",
                             (is_dual) ? "dual " : "");
      return 1; // fail
    }
  }
  else
    return -1; // not uniform name

  uni.get_poly(geom, sym_no);

  if (!is_std) {
    double e_len = geom.edge_vec(geom.faces(0, 0), geom.faces(0, 1)).len();
    geom.transform(Trans3d::scale(1 / e_len));
    set_resource_polygon_color(geom);
  }

  // this is here if 'd' was the second character
  if (is_dual)
    make_resource_dual(geom, is_std);

  return 0; // name found
}

int make_resource_wenninger(Geometry &geom, string name, bool is_std,
                            string *error_msg = nullptr)
{
  // Model name formats:
  //   wNumber,  WNumber
  //   wdNumber, WDNumber

  bool is_dual;
  if (!process_letter_dual_prefix(name, "wW", is_dual))
    return -1; // unhandled, not wenninger name

  Wenninger wenn;
  int sym_no;
  if (read_int(name.c_str() + 1, &sym_no)) {
    sym_no--;
    if (sym_no < 0 || !wenn.test_valid(sym_no)) {
      if (error_msg)
        *error_msg = "wenninger stellation number out of range";
      return 1; // fail
    }
  }
  else
    return -1; // not integer

  wenn.get_poly(geom, sym_no);

  // if is_std, have to strip built in color
  if (is_std) {
    geom.colors(VERTS).clear();
    geom.colors(EDGES).clear();
    geom.colors(FACES).clear();
  }
  // color by symmetry or compound looks better
  else {
    Coloring clrng(&geom);
    ColorMap *cmap = colormap_from_name("compound");
    clrng.add_cmap(cmap);

    GeometryInfo info(geom);
    if (info.num_parts() > 1)
      // color by compound
      clrng.f_parts(true);
    else {
      // color by symmetry
      Symmetry sym;
      vector<vector<set<int>>> sym_equivs;
      sym.init(geom, &sym_equivs);
      clrng.f_sets(sym_equivs[2], true);
    }

    // edges can be confusing. set them to invisible
    Coloring(&geom).vef_one_col(Color::invisible, Color::invisible, Color());
  }

  // this is here if 'd' was the second character
  if (is_dual)
    make_resource_dual(geom, is_std);

  return 0; // name found
}

int make_resource_uniform_compound(Geometry &geom, string name, bool is_std,
                                   string *error_msg = nullptr)
{
  // Model name formats:
  //   ucNumber,  UCNumber,  uc_Name
  // any model that takes parameters can be followed by _Parameters

  if (name.size() < 3 ||
      (name.compare(0, 2, "uc") != 0 && name.compare(0, 2, "UC") != 0))
    return -1; // not uniform compound name

  double angle = NAN;
  int n = -1;
  int d = -1;
  int k = -1;

  UniformCompound uniform_compounds;

  // to restore name if no parameters found
  string name_copy = name;

  // if complex parms, parse them
  if (strpbrk(name.c_str(), RES_SEPARATOR)) {
    int ret = uniform_compounds.parse_uc_args(name, angle, n, d, k, error_msg);
    if (ret > 0)
      return 1; // fail
  }

  // restore unmangled name if parse_uc_args found no parameters
  if (std::isnan(angle) && (n == -1) && (d == -1) && (k == -1))
    name = name_copy;

  int sym_no;
  if (read_int(name.c_str() + 2, &sym_no)) {
    sym_no--;
    if (sym_no < 0 || sym_no >= uniform_compounds.get_last_uc()) {
      if (error_msg)
        *error_msg = "uniform compound number out of range";
      return 1; // fail
    }
  }
  else if (strchr(RES_SEPARATOR, name[2])) {
    if (isupper(name[0]))
      return -1; // unhandled, upper case prefix, not UC full name

    string expanded = expand_abbrevs(
        name, uc_abbrevs, sizeof(uc_abbrevs) / sizeof(uc_abbrevs[0]));
    sym_no = uniform_compounds.lookup_sym_no(expanded.c_str());
    if (sym_no == -1)
      return 1; // fail
  }
  else
    return -1; // not uniform compound name

  int ret =
      uniform_compounds.set_uc_args(sym_no + 1, angle, n, d, k, error_msg);
  if (ret > 0)
    return 1; // fail

  uniform_compounds.get_poly(geom, sym_no, angle, n, d, k, is_std);

  if (!is_std) {
    Coloring clrng(&geom);
    ColorMap *cmap = colormap_from_name("compound");
    clrng.add_cmap(cmap);

    clrng.v_apply_cmap();
    clrng.e_apply_cmap();
    clrng.f_apply_cmap();
  }

  return 0; // name found
}

int make_resource_johnson(Geometry &geom, string name, bool is_std,
                          string *error_msg = nullptr)
{
  // Model name formats:
  //   jNumber,  JNumber,  j_Name
  //   jdNumber, JDNumber, jd_Name
  // any of these can be followed by _raw

  bool is_dual;
  if (!process_letter_dual_prefix(name, "jJ", is_dual))
    return -1; // unhandled, not johnson name

  // If model name ends in _raw then don't symmetry align
  bool sym_align = true;
  if (name.size() > 4 && name.substr(name.size() - 4) == "_raw") {
    name.resize(name.size() - 4);
    sym_align = false;
  }

  Johnson json;
  int sym_no;
  if (read_int(name.c_str() + 1, &sym_no)) {
    sym_no--;
    if (sym_no < 0 || sym_no >= json.get_last_J()) {
      if (error_msg)
        *error_msg = "Johnson polyhedron number out of range";
      return 1; // fail
    }
  }
  else if (strchr(RES_SEPARATOR, name[1])) {
    if (isupper(name[0]))
      return -1; // unhandled, upper case prefix, not johnson full name

    string expanded = expand_abbrevs(name, j_abbrevs,
                                     sizeof(j_abbrevs) / sizeof(j_abbrevs[0]));

    // check alternate names
    if (is_dual) {
      auto it = jd_alt_names.find(to_resource_name(expanded));
      if (it != jd_alt_names.end())
        expanded = it->second; // set name to the usual name for the model
    }
    else {
      auto it = j_alt_names.find(to_resource_name(expanded));
      if (it != j_alt_names.end())
        expanded = it->second; // set name to the usual name for the model
    }

    sym_no = json.lookup_sym_no(expanded.c_str());

    // if expanded name turned into a j number (alternate names), must add 1
    if (strchr("jJ", expanded[0]))
      sym_no++;

    if (sym_no == -1) {
      if (error_msg)
        *error_msg = "invalid Johnson polyhedron name";
      return 1; // fail
    }
  }
  else
    return -1; // not johnson name

  json.get_poly(geom, sym_no);

  // this is here if 'd' was the second character
  if (is_dual)
    make_resource_dual(geom, is_std);

  if (!is_std)
    set_resource_polygon_color(geom);

  if (sym_align)
    geom.sym_align();

  return 0; // name found
}

class res_Coloring : public Coloring {
public:
  res_Coloring(Geometry *geo = nullptr) : Coloring(geo) {}
  void f_all_angles(bool as_values);
  void f_max_angle(bool as_values);
};

void res_Coloring::f_all_angles(bool as_values)
{
  int faces_sz = get_geom()->faces().size();
  for (int i = 0; i < faces_sz; i++) {
    vector<double> f_angs;
    get_geom()->face_angles_lengths(i, &f_angs);
    sort(f_angs.rbegin(), f_angs.rend());
    double ang_sum = 0;
    for (unsigned int j = 0; j < f_angs.size(); j++)
      ang_sum += j * f_angs[j];
    int idx = (int)(rad2deg(50 * ang_sum /
                            (f_angs.size() * (f_angs.size() + 1) / 2)) +
                    0.5);
    if (as_values)
      get_geom()->colors(FACES).set(i, get_col(idx));
    else
      get_geom()->colors(FACES).set(i, idx);
  }
}

int make_resource_geodesic(Geometry &geom, string name, bool is_std,
                           string *error_msg)
{
  if (name.compare(0, 4, "geo_") != 0)
    return -1; // unhandled, not geodesic name

  Geometry base;
  int offset = 4;
  if (name[offset] == 't') {
    base.read_resource("std_tet");
    offset++;
  }
  else if (name[offset] == 'o') {
    base.read_resource("std_oct");
    offset++;
  }
  else if (name[offset] == 'i') {
    base.read_resource("std_ico");
    offset++;
  }
  else if (name[offset] == 'T') {
    double X = 0.25;
    double Y = sqrt(3.0) / 12;
    double Z = 0.8;
    base.add_vert(Vec3d(-X, -Y, Z));   // 0
    base.add_vert(Vec3d(X, -Y, Z));    // 1
    base.add_vert(Vec3d(0, 2 * Y, Z)); // 3

    vector<int> face(3);
    for (int i = 0; i < 3; i++)
      face[i] = i;
    base.add_face(face);
    offset++;
  }
  else if (isdigit(name[offset])) {
    base.read_resource("std_ico");
  }
  else {
    if (error_msg)
      *error_msg = "geodesic base polyhedron not i, o, t or T";
    return 1; // fail
  }

  Split num_strs(name.c_str() + offset, "_");
  int nums[2] = {1, 0};
  if (num_strs.size() > 2) {
    if (error_msg)
      *error_msg = "geodesic division: must be given by 1 or 2 integers";
    return 1; // fail
  }

  if (num_strs.size() > 0 && !read_int(num_strs[0], &nums[0])) {
    if (error_msg)
      *error_msg =
          msg_str("geodesic division:%s division number not an integer",
                  (num_strs.size() > 1) ? " first" : "");
    return 1; // fail
  }

  if (num_strs.size() > 1 && !read_int(num_strs[1], &nums[1])) {
    if (error_msg)
      *error_msg = "geodesic division: second division number not an integer";
    return 1; // fail
  }

  if (!make_geodesic_sphere(geom, base, nums[0], nums[1])) {
    if (error_msg)
      *error_msg = "geodesic division: invalid division";
    return 1; // fail
  }

  if (!is_std) {
    res_Coloring clrng(&geom);
    ColorMap *cmap = colormap_from_name(
        "rng256_R1:0.7:1:0.5:1G1:0.4:1:0.3:1B1:0.2:1:0.7:1%");
    clrng.add_cmap(cmap);
    clrng.f_all_angles(true);
  }

  return 0; // name found
}

// if name ends in "_d" then make dual of resource
void process_for_dual(string &name, bool &dual_flag)
{
  int name_sz = name.size();
  if (name_sz > 2 && name[name_sz - 2] == '_' && name[name_sz - 1] == 'd') {
    name.resize(name_sz - 2);
    dual_flag = !dual_flag;
  }
}

static void get_arrow(Geometry &geom, const Symmetry &sym)
{
  geom.clear_all();
  geom.add_vert(Vec3d(0.0, 0.0, 0.0));
  geom.add_vert(Vec3d(3.0, 0.0, 0.0));
  geom.add_vert(Vec3d(2.0, 1.0, 0.0));
  geom.add_vert(Vec3d(2.0, 0.5, 0.0));
  geom.add_vert(Vec3d(0.0, 0.5, 0.0));
  geom.add_vert(Vec3d(0.0, 0.0, 0.5));
  geom.add_vert(Vec3d(3.0, 0.0, 0.5));
  geom.add_vert(Vec3d(2.0, 0.5, 0.25));
  geom.add_vert(Vec3d(0.0, 0.5, 0.25));

  geom.add_face({1, 2, 3, 4, 0}, Color(0.9, 0.9, 0.7));
  geom.add_face({6, 2, 7, 8, 5}, Color(0.0, 0.2, 0.6));
  geom.add_face({0, 1, 6, 5}, Color(0.9, 0.9, 0.7));
  geom.add_face({1, 2, 6}, Color(0.9, 0.9, 0.7));
  geom.add_face({2, 3, 7}, Color(0.9, 0.9, 0.7));
  geom.add_face({3, 4, 8, 7}, Color(0.9, 0.9, 0.7));
  geom.add_face({4, 0, 5, 8}, Color(0.9, 0.9, 0.7));

  Trans3d trans = Trans3d::translate(Vec3d(0.1, 0.2, 0.4));
  int fold = sym.get_nfold();

  switch (sym.get_sym_type()) {
  case Symmetry::C1:
  case Symmetry::Cs:
  case Symmetry::Ci:
    break;

  case Symmetry::C:
  case Symmetry::Ch:
  case Symmetry::S:
    trans = Trans3d::rotate(Vec3d::Z, -M_PI / 2) *
            Trans3d::translate(Vec3d(-1.6, -0.2, 0.0)) * trans;
    if (fold > 1)
      trans = Trans3d::translate(Vec3d(1.6 / sin(M_PI / fold), 0, 0)) * trans;
    break;

  case Symmetry::Cv:
  case Symmetry::D:
  case Symmetry::Dh:
    trans = Trans3d::rotate(Vec3d::Z, -M_PI / 2) * trans;
    if (fold > 1)
      trans = Trans3d::translate(Vec3d(3.2 / tan(M_PI / fold), 0, 0)) * trans;
    break;

  case Symmetry::Dv:
    if (fold > 1)
      trans = Trans3d::rotate(Vec3d::Z, 0.5 * M_PI / fold) *
              Trans3d::translate(Vec3d(3.2 / tan(M_PI / fold), 0, 0)) *
              Trans3d::rotate(Vec3d::Z, -M_PI / 2) * trans;
    break;

  case Symmetry::T:
  case Symmetry::Td:
  case Symmetry::Th:
    trans =
        Trans3d::align(Vec3d::Z, Vec3d::X, Vec3d(1, 1, 1), Vec3d(0, -1, 1)) *
        Trans3d::translate(Vec3d(0, 0.5, 3)) *
        Trans3d::rotate(Vec3d::Z, M_PI / 2) * trans;
    break;

  case Symmetry::O:
  case Symmetry::Oh:
    // trans = Trans3d::transl(Vec3d(0, 3, 4)) * trans;
    trans = Trans3d::translate(Vec3d(1, 0, 4)) * trans;
    break;

  case Symmetry::I:
  case Symmetry::Ih:
    trans = Trans3d::rotate(Vec3d::Z, Vec3d(0, 1, 1.618)) *
            Trans3d::translate(Vec3d(0, 1, 5.5)) *
            Trans3d::rotate(Vec3d::Z, M_PI / 2) * trans;
    break;
  }

  geom.transform(trans);
}

int make_resource_sym(Geometry &geom, string name, bool is_std,
                      string *error_msg = nullptr)
{
  if (name.compare(0, 4, "sym_") != 0)
    return -1; // unhandles, not sym_ name

  Symmetry sym;
  Status stat = sym.init(name.substr(4), Trans3d());
  if (stat.is_error()) {
    if (error_msg)
      *error_msg = stat.msg();
    return 1;
  }

  Geometry arrow;
  get_arrow(arrow, sym);
  sym_repeat(geom, arrow, sym.get_trans());

  if (is_std)
    geom.clear_cols();

  return 0;
}

int make_resource_pgon(Geometry &geom, string name, bool is_std,
                       string *error_msg)
{
  if (name.size() < 4)
    return -1; // unhandles, not a polygon res name

  int num_sides;
  int step = 1;
  char *pnam = &name[0]; // don't use name as string again
  char *pnum = pnam + 3;
  char *p = strchr(pnum, '/');
  if (p != nullptr) {
    *p++ = '\0';
    if (!read_int(p, &step))
      return -1;
  }
  if (!read_int(pnum, &num_sides))
    return -1;
  if (num_sides < 2)
    return -1;
  if (step < 1)
    return -1;
  if (step % num_sides == 0)
    return -1;

  *pnum = '\0';

  Polygon pgon(num_sides, step);
  int type;
  int subtype = Polygon::sub_default;
  if (strcmp("pri", pnam) == 0)
    type = Polygon::prism;
  else if (strcmp("ant", pnam) == 0)
    type = Polygon::antiprism;
  else if (strcmp("pyr", pnam) == 0)
    type = Polygon::pyramid;
  else if (strcmp("dip", pnam) == 0)
    type = Polygon::dipyramid;
  else if (strcmp("cup", pnam) == 0)
    type = Polygon::cupola;
  else if (strcmp("ort", pnam) == 0)
    type = Polygon::orthobicupola;
  else if (strcmp("gyr", pnam) == 0)
    type = Polygon::gyrobicupola;
  else if (strcmp("snu", pnam) == 0)
    type = Polygon::snub_antiprism;
  else if (strcmp("pol", pnam) == 0) {
    type = Polygon::dihedron;
    subtype = Polygon::sub_dihedron_polygon;
  }
  else
    return -1; // unhandled

  pgon.set_type(type);
  pgon.set_subtype(subtype);

  if (!uni_pgon(geom, pgon)) {
    if (error_msg)
      *error_msg = "polyhedron cannot have unit edges";
    return 1; // failure
  }

  if (!is_std)
    set_resource_polygon_color(geom);

  return 0; // success
}

int make_resource_schwarz(Geometry &geom, string name, bool is_std,
                          string *error_msg = nullptr)
{
  if (name.size() < (7 + 5) || name.compare(0, 7, "schwarz") != 0)
    return -1; // unhandles, not schwarz name

  if (name.find(':') != string::npos || name.find('|') != string::npos) {
    if (error_msg)
      *error_msg = "schwarz name: cannot contain bar character";
    return 1; // fail
  }

  string symbol_str = name.substr(7);
  for (char &i : symbol_str) {
    if (i == '_') // allow _ in place of space
      i = ' ';
  }

  bool as_poly = false;
  if (symbol_str.size() && symbol_str[symbol_str.size() - 1] == 'p') {
    as_poly = true;
    symbol_str.resize(symbol_str.size() - 1);
  }

  symbol_str += "|"; // add to make sure symbol will be valid for wythoff
  Status stat;
  Wythoff wyt(symbol_str.c_str(), &stat);
  if (stat.is_error()) {
    if (error_msg)
      *error_msg = "schwarz name: invalid symbol: " + stat.msg();
    return 1; // fail
  }

  if (as_poly)
    wyt.make_tri_poly(geom);
  else
    wyt.make_tri(geom);

  if (!is_std) {
    res_Coloring clrng(&geom);
    ColorMap *cmap = colormap_from_name("map_grey20:ivory");
    clrng.add_cmap(cmap);
    clrng.f_apply_cmap();
  }

  return 0; // name found
}

static int make_resource_wythoff(Geometry &geom, string name, bool is_std,
                                 string *error_msg = nullptr)
{
  if (name.size() < (7 + 5) || name.compare(0, 7, "wythoff") != 0)
    return -1; // unhandled, not wythoff name

  string symbol_str = name.substr(7);
  for (char &i : symbol_str) {
    if (i == '_') // allow _ in place of space
      i = ' ';
    else if (i == ':') // allow : in place of |
      i = '|';
  }

  Status stat;
  Wythoff wyt(symbol_str.c_str(), &stat);
  if (stat.is_error()) {
    if (error_msg)
      *error_msg = "wythoff name: invalid symbol: " + stat.msg();
    return 1; // fail
  }

  if (!(stat = wyt.make_poly(geom))) {
    if (error_msg)
      *error_msg = "wythoff name: invalid symbol: " + stat.msg();
    return 1; // fail
  }
  if (stat.is_warning() && error_msg)
    *error_msg = "wythoff name: invalid symbol: " + stat.msg();
  // fprintf(stderr, "warning: wythoff name: %s\n", stat.c_msg());

  if (!is_std) {
    double e_len = geom.edge_vec(geom.faces(0, 0), geom.faces(0, 1)).len();
    geom.transform(Trans3d::scale(1 / e_len));
    res_Coloring clrng(&geom);
    ColorMap *cmap = colormap_from_name("map_red:blue:yellow:ivory");
    clrng.add_cmap(cmap);
    clrng.f_apply_cmap();
    clrng.e_apply_cmap();
  }
  else {
    geom.clear_cols();
    geom.clear(EDGES);
  }

  return 0; // name found
}

}; // anonymous namespace

static void rh_dodecahedron(Geometry &geom, bool is_std = false)
{
  geom.clear_all();
  if (is_std) {
    Geometry geom2;
    geom2.read_resource("std_cube");
    geom.add_verts(geom2.verts());
    geom2.read_resource("std_oct");
    geom2.transform(Trans3d::scale(2));
    geom.add_verts(geom2.verts());
    geom.add_hull();
    normalized_face_list(geom);
  }
  else {
    make_resource_uniform(geom, "U7", is_std);
    make_resource_dual(geom, is_std);
  }
}

static void rh_triacontahedron(Geometry &geom, bool is_std = false)
{
  if (is_std) {
    geom.clear_all();
    Geometry geom2;
    geom2.read_resource("std_ico");
    geom.add_verts(geom2.verts());
    geom2.read_resource("std_dod");
    geom.add_verts(geom2.verts());
    geom.add_hull();
    normalized_face_list(geom);
  }
  else {
    make_resource_uniform(geom, "U24", is_std);
    make_resource_dual(geom, is_std);
  }
}

static void rh_enneacontahedron(Geometry &geom, bool is_std = false)
{
  geom.clear_all();
  Geometry geom2;
  geom2.read_resource("std_dod");
  if (is_std)
    geom2.transform(Trans3d::scale(1 / (2 * phi * phi)));
  else
    geom2.transform(Trans3d::scale(1 / geom2.verts(0).len()));
  make_zonohedron(geom, geom2.verts());
  normalized_face_list(geom);
  if (!is_std)
    set_resource_polygon_color(geom);
}

static void rh_hexecontahedron(Geometry &geom, bool is_std = false)
{
  geom.clear_all();
  Geometry geom_face;
  geom_face.add_vert(Vec3d(0, phi + 1, 1));
  geom_face.add_vert(Vec3d(1 / phi, phi, 1));
  geom_face.add_vert(Vec3d(0, 1 / phi, 1));
  geom_face.add_vert(Vec3d(-1 / phi, phi, 1));
  geom_face.add_face({0, 1, 2, 3});
  if (!is_std) // Make unit edges
    geom_face.transform(Trans3d::scale(1 / geom_face.edge_vec(0, 1).len()));
  sym_repeat(geom, geom_face, Symmetry(Symmetry::I));
  merge_coincident_elements(geom, "vef", anti::epsilon);
  normalized_face_list(geom);
  if (!is_std)
    set_resource_polygon_color(geom);
}

// clang-format off
// http://www.math.unm.edu/~vageli/papers/FLEX/Szilassi.pdf
static void csaszar(Geometry &geom, bool is_std = false)
{
  double crds[2][7][3] = {// Coordintes attributed to Csaszar
                          {{3, -3, 0},
                           {3, 3, 1},
                           {1, 2, 3},
                           {-1, -2, 3},
                           {-3, -3, 1},
                           {-3, 3, 0},
                           {0, 0, 15}},
                          // Coordinates C1
                          {{4 * sqrt(15), 0, 0},
                           {0, 8, 4},
                           {-1, 2, 11},
                           {1, -2, 11},
                           {0, -8, 4},
                           {-4 * sqrt(15), 0, 0},
                           {0, 0, 20}}};
  geom.clear_all();
  for (int i = 0; i < 7; i++)
    geom.add_vert(Vec3d(crds[is_std][i]));

  geom.add_face({0, 1, 5});
  geom.add_face({0, 3, 1});
  geom.add_face({4, 2, 1});
  geom.add_face({3, 0, 2});
  geom.add_face({1, 6, 5});
  geom.add_face({2, 6, 1});
  geom.add_face({0, 6, 2});
  geom.add_face({5, 4, 0});
  geom.add_face({5, 2, 4});
  geom.add_face({1, 3, 4});
  geom.add_face({2, 5, 3});
  geom.add_face({4, 6, 0});
  geom.add_face({3, 6, 4});
  geom.add_face({5, 6, 3});

  if (!is_std)
    set_resource_unique_color(geom);
}

static void szilassi(Geometry &geom, bool is_std = false)
{
  double crds[14][3] = {
      // Coordinates, Tom Ace http://www.minortriad.com/szilassi.html
      {0.0, -12.6, 12.0},
      {0.0, 12.6, 12.0},
      {-7.0, 0.0, -2.0},
      {-2.0, 5.0, 8.0},
      {-12.0, 0.0, -12.0},
      {-7.0, -2.5, -2.0},
      {2.0, -5.0, 8.0},
      {7.0, 0.0, -2.0},
      {12.0, 0.0, -12.0},
      {7.0, 2.5, -2.0},
      {-3.75, -3.75, 3.0},
      {3.75, 3.75, 3.0},
      {-4.5, 2.5, -2.0},
      {4.5, -2.5, -2.0},
  };

  geom.clear_all();
  for (auto &crd : crds)
    geom.add_vert(Vec3d(crd));

  geom.add_face({0, 2, 5, 4, 3, 1});
  geom.add_face({1, 7, 9, 8, 6, 0});
  geom.add_face({0, 6, 10, 11, 12, 2});
  geom.add_face({1, 3, 11, 10, 13, 7});
  geom.add_face({4, 5, 13, 10, 6, 8});
  geom.add_face({8, 9, 12, 11, 3, 4});
  geom.add_face({5, 2, 12, 9, 7, 13});

  if (!is_std)
    set_resource_unique_color(geom);
}

int make_resource_misc_poly(Geometry &geom, string name, bool is_std,
                            string * /*error_msg = nullptr*/)
{
   map<string, std_model_func> models;
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
   models["rhombic_hexecontahedron"]     = rh_hexecontahedron;
   models["rh_hexecontahedron"]          = rh_hexecontahedron;
   models["rh_hex"]                      = rh_hexecontahedron;
   models["szilassi"]                    = szilassi;
   models["csaszar"]                     = csaszar;

  auto mi = models.find(name);
  if (mi != models.end()) {
    mi->second(geom, is_std);
    return 0; // name found
  }

  return -1; // not handled
}
// clang-format on

Status make_resource_geom(Geometry &geom, string name)
{
  geom.clear_all();

  if (!name.size())
    return Status::error("no name given");

  // if name ends in underscore, don't let it pass
  if (name.substr(name.size() - 1) == "_")
    return Status::error("not found (or invalid)");

  // common spelling difference for hexe(a)contahedron
  name = std::regex_replace(name, std::regex("hexacon"), "hexecon");

  bool make_dual = false;
  process_for_dual(name, make_dual);

  string error_msg;
  bool geom_ok = false;

  // RK: standard keyword "std_" now used for all built-in polyhedra
  // truncate "std_" but keep track of it
  bool is_std = (name.size() > 3 && name.substr(0, 4) == "std_");
  if (is_std)
    name = name.substr(4);

  // Look for an internal alternative name
  auto it = alt_names.find(to_resource_name(name));
  if (it != alt_names.end()) {
    name = it->second; // set name to the usual name for the model
    process_for_dual(name, make_dual);
  }

  // Only check model name if it has two or more chars, and no '.'
  if (name.size() > 1 && name.find('.') == string::npos) {
    vector<int (*)(Geometry &, string, bool, string *)> make_funcs = {
        make_resource_misc_poly,        make_resource_pgon,
        make_resource_uniform,          make_resource_wenninger,
        make_resource_uniform_compound, make_resource_johnson,
        make_resource_geodesic,         make_resource_wythoff,
        make_resource_schwarz,          make_resource_sym};

    for (auto make_func : make_funcs) {
      int ret = make_func(geom, name, is_std, &error_msg);
      if (ret == 0) {
        geom_ok = true;
        break;
      }
      else if (ret > 0)
        return Status::error(error_msg);
    }
  }

  // Catch any failure to make a valid geometry
  if (!geom.is_set())
    geom_ok = false;

  if (!geom_ok) {
    if (error_msg.empty())
      error_msg = "not found (or invalid)";
  }

  if (geom_ok && make_dual)
    make_resource_dual(geom, is_std);

  if (geom_ok)
    return (error_msg.empty()) ? Status::ok() : Status::warning(error_msg);
  else
    return Status::error(error_msg);
}
