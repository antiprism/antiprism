/*
   Copyright (c) 2012-2016, Adrian Rossiter

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
   Name: help.h
   Description: help text
   Project: Antiprism - http://www.antiprism.com
*/

#ifndef HELP_H
#define HELP_H

const char *help_help = R"(Help topics
===========
Type 'off_util -H topic_name' for help on each topic

colour:       selecting colours for models
   col_val:      colour value formats
   col_names:    colour names
   col_map:      colour map formats and resource colour maps
   expressions:  mathematical expressions for floating point arguments
   symmetry:     symmetry features and option arguments
   models:       built in models
   common_polys: common polyhedra
   uniform:      uniform polyhedra (including wythoff_ models)
   ud:           uniform dual polyhedra
   wenninger:    uniform polyhedra by Wenninger number and included stellations
   johnson:      johnson polyhedra
   uc:           uniform compounds
   bowers:       supported Bowers short name notations
   polygon:      prisms, pyramids, antiprisms, etc
   geodesic:     geodesic spheres
   sym_models:   symmetry example models
   schwarz:      Schwarz triangles
   std_polys:    information on usual coordinates
   uniform_syms: alternate names for uniform polyhedra
   ud_syms:      alternate names for uniform dual polyhedra
   johnson_syms: alternate names for johnson solids and their duals)";

const char *help_models = R"(Models
======
Antiprism includes a number of built in models. These may be used by 
any program that expects an OFF file, by passing the model name instead
of a file name e.g. antiview u1. Adding _d to a resource model name will
give the dual of the model. Prefixing std_ to a resource name will give
an uncoloured standard or usual representation of the polyhedron.

See the help topic for each model type for more details
   common_polys: common polyhedra
   uniform:      uniform polyhedra (including wythoff_ models)
   ud:           uniform dual polyhedra
   johnson:      johnson polyhedra
   uc:           uniform compounds
   bowers:       supported Bowers short name notations
   polygon:      prisms, pyramids, antiprisms, etc
   geodesic:     geodesic spheres
   sym_models:   symmetry example models
   schwarz:      Schwarz triangles
   std_polys:    information on usual coordinates)";

const char *help_common_polys = R"(Common Polyhedra
================
Common polyhedra are given by name only. Some names have their own 
construction, while others are just a synonym for another resource 
model name

   Constructed:
     rhombic_dodecahedron, rd
     rhombic_triacontahedron, rt
     rhombic_enneacontahedron, re
     rhombic_hexecontahedron, rh_hex
     szilassi
     csaszar

   Synonyms:
     tetrahedron, tet                          u1
     truncated_tetrahedron, tr_tet             u2
     octahedron, oct                           u5
     cube                                      u6
     cuboctahedron, cubo                       u7
     truncated_octahedron, tr_oct              u8
     truncated_cube, tr_cube                   u9
     rhombicuboctahedron, rh_cubo              u10
     truncated_cuboctahedron, tr_cubo          u11
     snub_cube, sn_cube                        u12
     icosahedron, icosa, ico                   u22
     dodecahedron, dod                         u23
     icosidodecahedron, icosid                 u24
     truncated_icosahedron, tr_ico             u25
     truncated_dodecahedron, tr_dod            u26
     rhombicosidodecahedron, rh_icosid         u27
     truncated_icosidodecahedron, tr_icosid    u28
     snub_dodecahedron, sn_dod                 u29
     small_stellated_dodecahedron, sm_st_dod   u34
     great_dodecahedron, gr_dod                u35
     great_stellated_dodecahedron, gr_st_dod   u52
     great_icosahedron, gr_ico                 u53
     triakis_tetrahedron, tri_tet              u2_d
     triakis_octahedron, tri_oct               u9_d
     tetrakis_hexahedron, tetr_hex             u8_d
     disdyakis_hexahedron, disd_hex            u8_d
     deltoidal_icositetrahedron, delt_icosit   u10_d
     trapezoidal_icositetrahedron, trap_icosit u10_d
     hexakis_octahedron, hex_oct               u11_d
     disdyakis_dodecahedron, disd_dod          u11_d
     pentagonal_icositetrahedron, pen_icosit   u12_d
     triakis_icosahedron, tri_ico              u26_d
     pentakis_dodecahedron, pent_dod           u25_d
     deltoidal_hexecontahedron, delt_hexec     u27_d
     trapezoidal_hexecontahedron, trap_hexec   u27_d
     hexakis_icosahedron, hex_ico              u28_d
     disdyakis_triacontahedron, disd_tri       u28_d
     pentagonal_hexecontahedron, pen_hexec     u29_d

     e.g. tet, great_icosahedron)";

const char *help_uniform = R"(Uniform Polyhedra
=================
Uniform polyhedra can be specified by

o   A U number e.g. u8

o   u_ followed by the specific Wythoff Symbol in the list below,
    e.g. 'u_2 4|3', which will be used as a lookup for the
    corresponding uniform 

o   wythoff_ followed by a general Wythoff Symbol, where '_' and ':'
    maybe used instead of ' ', and ':', e.g. wythoff_2_4:3. The model
    is constructed directly from the symbol. Note that some degenerate
    snub models may not be constructed correctly. The faces around each
    triangle vertex are coloured, in order, red, blue, yellow, and snub
    triangles are coloured white.

o   u_ followed by a name (see the list below) Use '_' instead of
    a space to avoid having to quote the model name.

    When giving a name the following abbreviations can be used
       tr:     truncated
       sm:     small
       gr:     great
       st:     stellated
       sn:     snub
       tet:    tetrahedron
       ico:    icosahedron
       icosa:  icosahedron
       dod:    dodecahedron
       oct:    octahedron
       cubo:   cuboctahedron
       icosid: icosidodecahedron

    e.g. u_truncated_octahedron, u_tr_octahedron, u_tr_oct

o   Common polyhedra can be given by name only (see help for
    'common_polys' for list) e.g. tet, cube

o   Uniform prisms are named 'pri' followed by the polygon
    fraction e.g pri5, pri5/2

o   Uniform antiprisms are named 'ant' followed by the polygon fraction
    e.g ant5, ant5/2, ant5/3

Uniform List:

U No.  Wythoff Sym    Name
-----  -----------    ------------------------
U1           3|2 3    tetrahedron
U2           2 3|3    truncated tetrahedron
U3         3/2 3|3    octahemioctahedron
U4         3/2 3|2    tetrahemihexahedron
U5           4|2 3    octahedron
U6           3|2 4    cube
U7           2|3 4    cuboctahedron
U8           2 4|3    truncated octahedron
U9           2 3|4    truncated cube
U10          3 4|2    rhombicuboctahedron
U11         2 3 4|    truncated cuboctahedron
U12         |2 3 4    snub cube
U13        3/2 4|4    small cubicuboctahedron
U14        3 4|4/3    great cubicuboctahedron
U15        4/3 4|3    cubohemioctahedron
U16       4/3 3 4|    cubitruncated cuboctahedron
U17        3/2 4|2    great rhombicuboctahedron
U18       3/2 2 4|    small rhombihexahedron
U19        2 3|4/3    stellated truncated hexahedron
U20       4/3 2 3|    great truncated cuboctahedron
U21     4/3 3/2 2|    great rhombihexahedron
U22          5|2 3    icosahedron
U23          3|2 5    dodecahedron
U24          2|3 5    icosidodecahedron
U25          2 5|3    truncated icosahedron
U26          2 3|5    truncated dodecahedron
U27          3 5|2    rhombicosidodecahedron
U28         2 3 5|    truncated icosidodecahedron
U29         |2 3 5    snub dodecahedron
U30        3|5/2 3    small ditrigonal icosidodecahedron
U31        5/2 3|3    small icosicosidodecahedron
U32       |5/2 3 3    small snub icosicosidodecahedron
U33        3/2 5|5    small dodecicosidodecahedron
U34        5|2 5/2    small stellated dodecahedron
U35        5/2|2 5    great dodecahedron
U36        2|5/2 5    dodecadodecahedron
U37        2 5/2|5    truncated great dodecahedron
U38        5/2 5|2    rhombidodecadodecahedron
U39       2 5/2 5|    small rhombidodecahedron
U40       |2 5/2 5    snub dodecadodecahedron
U41        3|5/3 5    ditrigonal dodecadodecahedron
U42        3 5|5/3    great ditrigonal dodecicosidodecahedron
U43        5/3 3|5    small ditrigonal dodecicosidodecahedron
U44        5/3 5|3    icosidodecadodecahedron
U45       5/3 3 5|    icositruncated dodecadodecahedron
U46       |5/3 3 5    snub icosidodecadodecahedron
U47        3/2|3 5    great ditrigonal icosidodecahedron
U48        3/2 5|3    great icosicosidodecahedron
U49        3/2 3|5    small icosihemidodecahedron
U50       3/2 3 5|    small dodecicosahedron
U51        5/4 5|5    small dodecahemidodecahedron
U52        3|2 5/2    great stellated dodecahedron
U53        5/2|2 3    great icosahedron
U54        2|5/2 3    great icosidodecahedron
U55        2 5/2|3    great truncated icosahedron
U56       2 5/2 3|    rhombicosahedron
U57       |2 5/2 3    great snub icosidodecahedron
U58        2 5|5/3    small stellated truncated dodecahedron
U59       5/3 2 5|    truncated dodecadodecahedron
U60       |5/3 2 5    inverted snub dodecadodecahedron
U61      5/2 3|5/3    great dodecicosidodecahedron
U62      5/3 5/2|3    small dodecahemicosahedron
U63     5/3 5/2 3|    great dodecicosahedron
U64     |5/3 5/2 3    great snub dodecicosidodecahedron
U65        5/4 5|3    great dodecahemicosahedron
U66        2 3|5/3    great stellated truncated dodecahedron
U67        5/3 3|2    great rhombicosidodecahedron
U68       5/3 2 3|    great truncated icosidodecahedron
U69       |5/3 2 3    great inverted snub icosidodecahedron
U70    5/3 5/2|5/3    great dodecahemidodecahedron
U71      3/2 3|5/3    great icosihemidodecahedron
U72   |3/2 3/2 5/2    small retrosnub icosicosidodecahedron
U73     3/2 5/3 2|    great rhombidodecahedron
U74     |3/2 5/3 2    great retrosnub icosidodecahedron
U75  3/2 5/3 3 5/2    great dirhombicosidodecahedron
U76          2 5|2    pentagonal prism
U77         |2 2 5    pentagonal antiprism
U78        2 5/2|2    pentagrammic prism
U79       |2 2 5/2    pentagrammic antiprism
U80       |2 2 5/3    pentagrammic crossed antiprism)";

const char *help_wenninger = R"(Wenninger Number
================
Polyhedra can be listed by Wenninger Number and includes Wenninger Stellations

o   A W number e.g. w8

Wenninger List:

W No.  Name
-----  ------------------------
W1     tetrahedron
W2     octahedron
W3     cube
W4     icosahedron
W5     dodecahedron
W6     truncated tetrahedron
W7     truncated octahedron
W8     truncated cube
W9     truncated icosahedron
W10    truncated dodecahedron
W11    cuboctahedron
W12    icosidodecahedron
W13    rhombicuboctahedron
W14    rhombicosidodecahedron
W15    truncated cuboctahedron
W16    truncated icosidodecahedron
W17    snub cube
W18    snub dodecahedron
W19    compound of two tetrahedra, also uc4
W20    small stellated dodecahedron
W21    great dodecahedron
W22    great stellated dodecahedron
W23    compound of five octahedra, also uc17
W24    compound of five tetrahedra, also uc5
W25    compound of ten tetrahedra, also uc6
W26    small triambic icosahedron, also u30_d
W27    second stellation of icosahedron
W28    third stellation of icosahedron, excavated dodecahedron
W29    fourth stellation of icosahedron
W30    fifth stellation of icosahedron
W31    sixth stellation of icosahedron
W32    seventh stellation of icosahedron
W33    eighth stellafion of icosahedron, great triambic icosahedron, also u47_d
W34    ninth stellation of icosahedron
W35    tenth stellation of icosahedron
W36    eleventh stellation of icosahedron
W37    twelfth stellation of icosahedron
W38    thirteenth stellation of icosahedron
W39    fourteenth stellation of icosahedron
W40    fifteenth stellation of icosahedron
W41    great icosahedron
W42    final stellation of the icosahedron
W43    compound of cube and octahedron
W44    second stellation of cuboctahedron
W45    third stellation of cuboctahedron
W46    fourth stellation of cuboctahedron
W47    compound of dodecahedron and icosahedron
W48    second stellation of icosidodecahedron
W49    third stellation of icosidodecahedron
W50    compound of small stellated dodecahedron and triakis icosahedron
W51    compound of small stellated dodecahedron and five octahedra
W52    sixth stellation of icosidodecahedron
W53    seventh stellation of icosidodecahedron
W54    compound of five tetrahedra and great dodecahedron
W55    ninth stellation of icosidodecahedron
W56    tenth stellation of icosidodecahedron
W57    eleventh stellation of icosidodecahedron
W58    twelfth stellation of icosidodecahedron
W59    thirteenth stellation of icosidodecahedron
W60    fourteenth stellation of icosidodecahedron
W61    compound of great stellated dodecahedron and great icosahedron
W62    fifteenth stellation of icosidodecahedron
W63    sixteenth stellation of icosidodecahedron
W64    seventeenth stellation of icosidodecahedron
W65    eighteenth stellation of icosidodecahedron
W66    nineteenth stellation of icosidodecahedron
W67    tetrahemihexahedron
W68    octahemioctahedron
W69    small cubicuboctahedron
W70    small ditrigonal icosidodecahedron
W71    small icosicosidodecahedron
W72    small dodecicosidodecahedron
W73    dodecadodecahedron
W74    small rhombidodecahedron
W75    truncated great dodecahedron
W76    rhombidodecadodecahedron
W77    great cubicuboctahedron
W78    cubohemioctahedron
W79    cubitruncated cuboctahedron
W80    ditrigonal dodecadodecahedron
W81    great ditrigonal dodecicosidodecahedron
W82    small ditrigonal dodecicosidodecahedron
W83    icosidodecadodecahedron
W84    icositruncated dodecadodecahedron
W85    great rhombicuboctahedron
W86    small rhombihexahedron
W87    great ditrigonal icosidodecahedron
W88    great icosicosidodecahedron
W89    small icosihemidodecahedron
W90    small dodecicosahedron
W91    small dodecahemidodecahedron
W92    stellated truncated hexahedron
W93    great truncated cuboctahedron
W94    great icosidodecahedron
W95    great truncated icosahedron
W96    rhombicosahedron
W97    small stellated truncated dodecahedron
W98    truncated dodecadodecahedron
W99    great dodecicosidodecahedron
W100   small dodecahemicosahedron
W101   great dodecicosahedron
W102   great dodecahemicosahedron
W103   great rhombihexahedron
W104   great stellated truncated dodecahedron
W105   great rhombicosidodecahedron
W106   great icosihemidodecahedron
W107   great dodecahemidodecahedron
W108   great truncated icosidodecahedron
W109   great rhombidodecahedron
W110   small snub icosicosidodecahedron
W111   snub dodecadodecahedron
W112   snub icosidodecadodecahedron
W113   great inverted snub icosidodecahedron
W114   inverted snub dodecadodecahedron
W115   great snub dodecicosidodecahedron
W116   great snub icosidodecahedron
W117   great retrosnub icosidodecahedron
W118   small retrosnub icosicosidodecahedron
W119   great dirhombicosidodecahedron)";

const char *help_uniform_syms = R"(Uniform Synonyms
================
Uniform polyhedra have aquired many alternate names over the years. This is a
list of many names that have been seen and will work in Antiprism.

o   u_ followed by the name (see the list below) Use '_' instead of
    a space to avoid having to quote the model name.

    Synonyms:
      triangular pyramid                               u1
      cantic cube                                      u2
      allelotetratetrahedron                           u3
      hemicuboctahedron                                u4
      hemihexahedron                                   u4
      truncated hexahedron                             u9
      small rhombicuboctahedron                        u10
      beveled cube                                     u11
      cantitruncated cube                              u11
      convex great rhombicuboctahedron                 u11
      omnitruncated cube                               u11
      rhombitruncated cuboctahedron                    u11
      cubus simus                                      u12
      snub cuboctahedron                               u12
      cuboctatruncated cuboctahedron                   u16
      nonconvex great rhombicuboctahedron              u17
      quasirhombicuboctahedron                         u17
      small rhombicube                                 u18
      quasitruncated hexahedron                        u19
      stellatruncated cube                             u19
      quasitruncated cuboctahedron                     u20
      stellatruncated cuboctahedron                    u20
      great rhombicube                                 u21
      convex great rhombicosidodecahedron              u28
      great rhombicosidodecahedron                     u28
      omnitruncated dodecahedron                       u28
      omnitruncated icosahedron                        u28
      rhombitruncated icosidodecahedron                u28
      snub icosidodecahedron                           u29
      small ditrigonary icosidodecahedron              u30
      small triambic icosidodecahedron                 u30
      small icosified icosidodecahedron                u31
      hastur                                           u32
      holosnub icosahedron                             u32
      snub disicosidodecahedron                        u32
      small dodekicosidodecahedron                     u33
      dodecadodecahedron                               u36
      great dodecadodecahedron                         u36
      great truncated dodecahedron                     u37
      chaugnar faugn                                   u40
      ditrigonary dodecadodecahedron                   u41
      triambic dodecadodecahedron                      u41
      great dodekified icosidodecahedron               u42
      small dodekified icosidodecahedron               u43
      icosified dodecadodecahedron                     u44
      icosidodecatruncated icosidodecahedron           u45
      dagon                                            u46
      great ditrigonary icosidodecahedron              u47
      great triambic icosidodecahedron                 u47
      great icosified icosidodecahedron                u48
      small icosahemidodecahedron                      u49
      great truncated icosahedron                      u55
      truncated great icosahedron                      u55
      tsathoggua                                       u57
      quasitruncated small stellated dodecahedron      u58
      small stellatruncated dodecahedron               u58
      quasitruncated dodecadodecahedron                u59
      stellatruncated dodecadodecahedron               u59
      nyarlathotep                                     u60
      vertisnub dodecadodecahedron                     u60
      great dodecahemiicosahedron                      u62
      great dodekicosahedron                           u63
      great snub dodekicosidodecahedron                u64
      great snub icosidisdodecahedron                  u64
      shub niggurath                                   u64
      small dodecahemiicosahedron                      u65
      great stellatruncated dodecahedron               u66
      quasitruncated great stellated dodecahedron      u66
      nonconvex great rhombicosidodecahedron           u67
      quasirhombicosidodecahedron                      u67
      great quasitruncated icosidodecahedron           u68
      stellatruncated icosidodecahedron                u68
      cthulhu                                          u69
      great vertisnub icosidodecahedron                u69
      great icosahemidodecahedron                      u71
      retroholosnub icosahedron                        u72
      retrosnub disicosidodecahedron                   u72
      retrosnub ditrigonary icosidodecahedron          u72
      small inverted retrosnub icosicosidodecahedron   u72
      yog sothoth                                      u72
      azathoth                                         u74
      great inverted retrosnub icosidodecahedron       u74
      great disnub disicosidisdodecahedron             u75
      millers monster                                  u75
      pentagrammatic prism                             u78
      pentagrammatic antiprism                         u79
      pentagrammatic crossed antiprism                 u80)";

const char *help_uniform_duals = R"(Uniform Dual Polyhedra
======================
Uniform dual polyhedra can be specified by

o   ud followed by a U number e.g. ud8

o   ud_ followed by a name (see the list below) Use '_' instead of
    a space to avoid having to quote the model name.

    When giving a name the following abbreviations can be used
       sm:       small
       gr:       great
       st:       stellated
       inv:      inverted
       delt:     deltoidal
       pen:      pentagonal
       med:      medial
       triam:    triambic
       ditrig:   ditrigonal
       tri:      triakis
       tetr:     tetrakis
       pent:     pentakis
       hex:      hexakis
       disd:     disdyakis
       tet:      tetrahedron
       ico:      icosahedron
       icosa:    icosahedron
       dod:      dodecahedron
       oct:      octahedron
       cubo:     cuboctahedron
       hexa:     hexahedron
       hexec:    hexecontahedron
       hexac:    hexacontahedron
       icositet: icositetrahedron

    e.g. ud_triakis_tetrahedron, ud_tri_tetrahedron, ud_tri_tet

o   Common polyhedra can be given by name only (see help for
    'common_polys' for list) e.g. triakis_tetrahedron, triakis_octahedron

o   The duals of uniform prisms are named 'pri' followed by the polygon
    fraction and ending with '_d' e.g pri5_d, pri5/2_d

o   The duals of uniform antiprisms are named 'ant' followed by the polygon
    fraction and ending with '_d' e.g ant5_d, ant5/2_d, ant5/3_d

Uniform Dual List:

UD No. Name
------ ------------------------
UD1    tetrahedron
UD2    triakis tetrahedron
UD3    octahemioctacron
UD4    tetrahemihexacron
UD5    cube
UD6    octahedron
UD7    rhombic dodecahedron
UD8    tetrakis hexahedron
UD9    triakis octahedron
UD10   deltoidal icositetrahedron
UD11   disdyakis dodecahedron
UD12   pentagonal icositetrahedron
UD13   small hexacronic icositetrahedron
UD14   great hexacronic icositetrahedron
UD15   hexahemioctacron
UD16   tetradyakis hexahedron
UD17   great deltoidal icositetrahedron
UD18   small rhombihexacron
UD19   great triakis octahedron
UD20   great disdyakis dodecahedron
UD21   great rhombihexacron
UD22   dodecahedron
UD23   icosahedron
UD24   rhombic triacontahedron
UD25   pentakis dodecahedron
UD26   triakis icosahedron
UD27   deltoidal hexecontahedron
UD28   disdyakis triacontahedron
UD29   pentagonal hexecontahedron
UD30   small triambic icosahedron
UD31   small icosacronic hexecontahedron
UD32   small hexagonal hexecontahedron
UD33   small dodecacronic hexecontahedron
UD34   great dodecahedron
UD35   small stellated dodecahedron
UD36   medial rhombic triacontahedron
UD37   small stellapentakis dodecahedron
UD38   medial deltoidal hexecontahedron
UD39   small rhombidodecacron
UD40   medial pentagonal hexecontahedron
UD41   medial triambic icosahedron
UD42   great ditrigonal dodecacronic hexecontahedron
UD43   small ditrigonal dodecacronic hexecontahedron
UD44   medial icosacronic hexecontahedron
UD45   tridyakis icosahedron
UD46   medial hexagonal hexecontahedron
UD47   great triambic icosahedron
UD48   great icosacronic hexecontahedron
UD49   small icosihemidodecacron
UD50   small dodecicosacron
UD51   small dodecahemidodecacron
UD52   great icosahedron
UD53   great stellated dodecahedron
UD54   great rhombic triacontahedron
UD55   great stellapentakis dodecahedron
UD56   rhombicosacron
UD57   great pentagonal hexecontahedron
UD58   great pentakis dodecahedron
UD59   medial disdyakis triacontahedron
UD60   medial inverted pentagonal hexecontahedron
UD61   great dodecacronic hexecontahedron
UD62   small dodecahemicosacron
UD63   great dodecicosacron
UD64   great hexagonal hexecontahedron
UD65   great dodecahemicosacron
UD66   great triakis icosahedron
UD67   great deltoidal hexecontahedron
UD68   great disdyakis triacontahedron
UD69   great inverted pentagonal hexecontahedron
UD70   great dodecahemidodecacron
UD71   great icosihemidodecacron
UD72   small hexagrammic hexecontahedron
UD73   great rhombidodecacron
UD74   great pentagrammic hexecontahedron
UD75   great dirhombicosidodecacron
UD76   pentagonal dipyramid
UD77   pentagonal deltohedron
UD78   pentagrammic dipyramid
UD79   pentagrammic deltohedron
UD80   pentagrammic concave deltohedron)";

const char *help_uniform_dual_syms = R"(Uniform Dual Synonyms
================
Uniform dual polyhedra have aquired many alternate names over the years. This
is a list of many names that have been seen and will work in Antiprism.

o   ud_ followed by the name (see the list below) Use '_' instead of
    a space to avoid having to quote the model name.

    Synonyms:
      kistetrahedron                                     u2_d
      triakistetrahedron                                 u2_d
      hexahedron                                         u5_d
      disdyakis hexahedron                               u8_d
      hexakis tetrahedron                                u8_d
      hextetrahedron                                     u8_d
      kiscube                                            u8_d
      tetrahexahedron                                    u8_d
      tetrakis cube                                      u8_d
      tetrakishexahedron                                 u8_d
      kisoctahedron                                      u9_d
      triakisoctahedron                                  u9_d
      trigonal trisoctahedron                            u9_d
      trisoctahedron                                     u9_d
      lanceolar disdodecahedron                          u10_d
      lanceolar icositetrahedron                         u10_d
      strombic icositetrahedron                          u10_d
      tetragonal icosikaitetrahedron                     u10_d
      tetragonal trisoctahedron                          u10_d
      trapezoidal icositetrahedron                       u10_d
      disdyakisdodecahedron                              u11_d
      hexakis octahedron                                 u11_d
      hexoctahedron                                      u11_d
      kisrhombic dodecahedron                            u11_d
      octakis cube                                       u11_d
      octakis hexahedron                                 u11_d
      pentagonal icosikaitetrahedron                     u12_d
      tetradyakishexahedron                              u16_d
      great sagittal disdodecahedron                     u17_d
      great strombic icositetrahedron                    u17_d
      great trapezoidal icositetrahedron                 u17_d
      small dipteral disdodecahedron                     u18_d
      great triakisoctahedron                            u19_d
      great disdyakisdodecahedron                        u20_d
      great dipteral disdodecahedron                     u21_d
      triacontahedron                                    u24_d
      kisdodecahedron                                    u25_d
      pentakisdodecahedron                               u25_d
      kisicosahedron                                     u26_d
      triakisicosahedron                                 u26_d
      strombic hexecontahedron                           u27_d
      tetragonal hexecontahedron                         u27_d
      trapezoidal hexecontahedron                        u27_d
      decakis dodecahedron                               u28_d
      disdyakistriacontahedron                           u28_d
      hexakis icosahedron                                u28_d
      kisrhombic triacontahedron                         u28_d
      small lanceal trisicosahedron                      u31_d
      hexagonal hexecontahedron                          u32_d
      small sagittal ditriacontahedron                   u33_d
      midly rhombic triacontahedron                      u36_d
      small stellated triacontahedron                    u36_d
      small stellapentakisdodecahedron                   u37_d
      medial trapezoidal hexecontahedron                 u38_d
      midly triambic icosahedron                         u41_d
      great dodecacronic hexecontahedron                 u42_d
      great lanceal trisicosahedron                      u42_d
      fat star                                           u43_d
      small dodecacronic hexecontahedron                 u43_d
      midly sagittal ditriacontahedron                   u44_d
      tridyakisicosahedron                               u45_d
      midly dentoid ditriacontahedron                    u46_d
      medial triambic icosahedron                        u47_d
      great sagittal trisicosahedron                     u48_d
      small dipteral trisicosahedron                     u50_d
      great stellated triacontahedron                    u54_d
      great astropentakis dodecahedron                   u55_d
      great stellapentakisdodecahedron                   u55_d
      midly dipteral ditriacontahedron                   u56_d
      great petaloid ditriacontahedron                   u57_d
      great pentakisdodecahedron                         u58_d
      great pentakisdodekahedron                         u58_d
      medial disdyakistriacontahedron                    u59_d
      midly petaloid ditriacontahedron                   u60_d
      great lanceal ditriacontahedron                    u61_d
      great dodecahemiicosacron                          u62_d
      great dipteral trisicosahedron                     u63_d
      great astroid ditriacontahedron                    u64_d
      small dodecahemiicosacron                          u65_d
      great triakisicosahedron                           u66_d
      great sagittal ditriacontahedron                   u67_d
      great sagittal hexecontahedron                     u67_d
      great strombic hexecontahedron                     u67_d
      great trapezoidal hexecontahedron                  u67_d
      great disdyakistriacontahedron                     u68_d
      trisdyakis icosahedron                             u68_d
      petaloidal trisicosahedron                         u69_d
      small hexagrammatic hexecontahedron                u72_d
      great dipteral ditriacontahedron                   u73_d
      great dentoid ditriacontahedron                    u74_d
      great pentagrammatic hexecontahedron               u74_d
      pentagonal trapezohedron                           u77_d
      pentagrammatic dipyramid                           u78_d
      pentagrammatic trapezohedron                       u79_d
      pentagrammic trapezohedron                         u79_d
      pentagrammatic concave trapezohedron               u80_d
      pentagrammic concave trapezohedron                 u80_d)";

const char *help_johnson = R"(Johnson Polyhedra
=================
Norman Johnson says there are 17 elementary Johnson solids. They cannot be
made by augmenting another solid.

(1-6, 63, 80, 83-86, 88-92)
square pyramid (J1)
pentagonal pyramid (J2)
triangular cupola (J3)
square cupola (J4)
pentagonal cupola (J5)
pentagonal rotunda (J6)
tridiminished icosahedron (J63)
parabidiminished rhombicosidodecahedron (J80)
tridiminished rhombicosidodecahedron (J83)
snub disphenoid (J84)
snub square antiprism (J85)
sphenocorona (J86)
sphenomegacorona (J88)
hebesphenomegacorona (J89)
disphenocingulum (J90)
bilunabirotunda (J91)
triangular hebesphenorotunda (J92)

Johnson polyhedra can be specified by a

o   A J number e.g. j8

o   j_ followed by a name (see the list below) Use '_' instead of
    a space to avoid having to quote the model name.

    When giving a name the following abbreviations can be used
       tri:   triangular
       sq:    square
       squ:   square
       pe:    pentagonal
       pen:   pentagonal
       el:    elongated
       ge:    gyroelongated
       tr:    truncated
       au:    augmented
       ba:    biaugmenbted
       ta:    triaugmented

    e.g. j_elongated_triangular_pyramid, j_el_tri_pyramid

The Johnson polyhedron will be aligned with its symmetry group, but
appending '_raw' to the name will preserve the original construction
alignment.

Johnson List:

J No.  Name
-----  -------------------------------
J1     square pyramid
J2     pentagonal pyramid
J3     triangular cupola
J4     square cupola
J5     pentagonal cupola
J6     pentagonal rotunda
J7     elongated triangular pyramid
J8     elongated square pyramid
J9     elongated pentagonal pyramid
J10    gyroelongated square pyramid
J11    gyroelongated pentagonal pyramid
J12    triangular dipyramid
J13    pentagonal dipyramid
J14    elongated triangular dipyramid
J15    elongated square dipyramid
J16    elongated pentagonal dipyramid
J17    gyroelongated square dipyramid
J18    elongated triangular cupola
J19    elongated square cupola
J20    elongated pentagonal cupola
J21    elongated pentagonal rotunda
J22    gyroelongated triangular cupola
J23    gyroelongated square cupola
J24    gyroelongated pentagonal cupola
J25    gyroelongated pentagonal rotunda
J26    gyrobifastigium
J27    triangular orthobicupola
J28    square orthobicupola
J29    square gyrobicupola
J30    pentagonal orthobicupola
J31    pentagonal gyrobicupola
J32    pentagonal orthocupolarotunda
J33    pentagonal gyrocupolarotunda
J34    pentagonal orthobirotunda
J35    elongated triangular orthobicupola
J36    elongated triangular gyrobicupola
J37    elongated square gyrobicupola
J38    elongated pentagonal orthobicupola
J39    elongated pentagonal gyrobicupola
J40    elongated pentagonal orthocupolarotunda
J41    elongated pentagonal gyrocupolarotunda
J42    elongated pentagonal orthobirotunda
J43    elongated pentagonal gyrobirotunda
J44    gyroelongated triangular bicupola
J45    gyroelongated square bicupola
J46    gyroelongated pentagonal bicupola
J47    gyroelongated pentagonal cupolarotunda
J48    gyroelongated pentagonal birotunda
J49    augmented triangular prism
J50    biaugmented triangular prism
J51    triaugmented triangular prism
J52    augmented pentagonal prism
J53    biaugmented pentagonal prism
J54    augmented hexagonal prism
J55    parabiaugmented hexagonal prism
J56    metabiaugmented hexagonal prism
J57    triaugmented hexagonal prism
J58    augmented dodecahedron
J59    parabiaugmented dodecahedron
J60    metabiaugmented dodecahedron
J61    triaugmented dodecahedron
J62    metabidiminished icosahedron
J63    tridiminished icosahedron
J64    augmented tridiminished icosahedron
J65    augmented truncated tetrahedron
J66    augmented truncated cube
J67    biaugmented truncated cube
J68    augmented truncated dodecahedron
J69    parabiaugmented truncated dodecahedron
J70    metabiaugmented truncated dodecahedron
J71    triaugmented truncated dodecahedron
J72    gyrate rhombicosidodecahedron
J73    parabigyrate rhombicosidodecahedron
J74    metabigyrate rhombicosidodecahedron
J75    trigyrate rhombicosidodecahedron
J76    diminished rhombicosidodecahedron
J77    paragyrate diminished rhombicosidodecahedron
J78    metagyrate diminished rhombicosidodecahedron
J79    bigyrate diminished rhombicosidodecahedron
J80    parabidiminished rhombicosidodecahedron
J81    metabidiminished rhombicosidodecahedron
J82    gyrate bidiminished rhombicosidodecahedron
J83    tridiminished rhombicosidodecahedron
J84    snub disphenoid
J85    snub square antiprism
J86    sphenocorona
J87    augmented sphenocorona
J88    sphenomegacorona
J89    hebesphenomegacorona
J90    disphenocingulum
J91    bilunabirotunda
J92    triangular hebesphenorotunda)";

const char *help_johnson_syms = R"(Johnson Solid Synonyms
======================
Johnson Solids have aquired many alternate names over the years. This
is a list of many names that have been seen and will work in Antiprism.

o   j_ followed by the name (see the list below) Use '_' instead of
    a space to avoid having to quote the model name.

    Synonyms:
      lesser dome                                      j4
      diminished icosahedron                           j11
      hexadeltahedron                                  j12
      triangular dipyramid                             j12
      decadeltahedron                                  j13
      pentagonal dipyramid                             j13
      elongated triangular dipyramid                   j14
      triakis triangular prism                         j14
      elongated octahedron                             j15
      pencil cube                                      j15
      12 faced pencil cube                             j15
      pentakis pentagonal prism                        j16
      heccaidecadeltahedron                            j17
      hexakaidecadeltahedron                           j17
      tetrakis square antiprism                        j17
      anticuboctahedron                                j27
      disheptahedron                                   j27
      twisted cuboctahedron                            j27
      cantellated triangular prism                     j35
      gyrate rhombicuboctahedron                       j37
      millers solid,                                   j37
      pseudo rhombicuboctahedron                       j37
      cantellated pentagonal prism                     j38
      monolaterotruncated triangular bipyramid         j49
      tetracaidecadeltahedron,                         j51
      tetrakis triangular prism                        j51
      dodecadeltahedron                                j84
      siamese dodecahedron                             j84
      triangular dodecahedron                          j84
      trigonal dodecahedron                            j84
      pentakis elongated gyrobifastigium               j90

o   jd_ followed by the name (see the list below) Use '_' instead of
    a space to avoid having to quote the model name.
    
    Synonyms for duals:    
      triangular prism                                 j12_d
      pentagonal prism                                 j13_d
      triangular bifrustum                             j14_d
      square bifrustum                                 j15_d
      pentagonal bifrustum                             j16_d
      elongated tetragonal disphenoid                  j26_d
      schmitt conway danzer biprism                    j26_d
      trapezo rhombic dodecahedron                     j27_d
      elongated square trapezohedron                   j29_d
      elongated pentagonal trapezohedron               j31_d
      trapezo rhombic triacontahedron                  j34_d
      pseudo deltoidal icositetrahedron                j37_d
      associahedron k5                                 j51_d
      truncated triangular bipyramid                   j51_d
      monolaterotruncated pentagonal bipyramid         j52_d
      parabilaterotruncated pentagonal bipyramid       j53_d
      monolaterotruncated hexagonal bipyramid          j54_d
      parabilaterotruncated hexagonal bipyramid        j55_d
      alternate truncated hexagonal bipyramid          j57_d
      gyroelongated pentagonal bifrustum               j58_d
      elongated gyrobifastigium                        j84_d
      truncated snub disphenoid                        j90_d)";

const char *help_polygon = R"(Polygon-based Polyhedra
=======================
Polygon-based polyhedra have unit edges and are specified by their
abbreviated class name followed by the base polygon, given as a fraction.
(Note that some models, e.g. pyr7, cannot be made to have unit edges.)

The class names are
   pol:     polygon
   pri:     prism
   ant:     antiprism
   pyr:     pyramid
   dip:     dipyramid
   cup:     cupola
   ort:     orthobicupola
   gyr:     gyrobicupola
   snu:     snub antiprism

   e.g. pri5, pyr7/2)";

const char *help_std_polys = R"(Standard Polyhedra
==================
Any built-in polyhedron can be generated in a standard form by prefixing
std_ before the name. Any polyhedron generated this way will have no colour
added. For uniform polyhedra, uniform duals, and uniform compounds, their
coordinates are left as calculated and not given unit edge lengths

Some notable models and their coordinates are
   tetrahedron, tet:                 0, 1
   truncated_tetrahedron, tr_tet:    1, 3
   cube:                             1
   truncated_cube, tr_cube:          1, sqrt(2)-1
   snub_cube, sn_cube:               expressions using cube roots
   octahedron, oct:                  0, 1
   truncated_octahedron, tr_oct:     0, 1, 2
   dodecahedron, dod:                0, phi, 1/phi
   snub_dodecahedron, sn_dod:        expressions using phi and cube roots
   icosahedron, ico:                 0, 1, phi
   cuboctahedron, cubo:              0, 1
   truncated_cuboctahedron, tr_cubo: 1, 1+sqrt(2), 1+2*sqrt(2)
   truncated_icosidodecahedron, tr_icosid:
                                     phi and integer based
   rhombicuboctahedron, rhombicubo, rh_cubo:
                                     1, 1+sqrt(2)
   rhombicosidodecahedron, rhombicosid, rh_icosid:
                                     1, phi, phi^2, phi^3, 2phi, 2+phi

   rhombic_dodecahedron, rd           0, 1, 2
   rhombic_triacontahedron, rt        0, 1, phi, 1/phi
   rhombic_enneacontahedron, re       0, 1, expresions using phi

   e.g. std_truncated_tetrahedron, std_tr_tet, std_rh_cubo)";

const char *help_uniform_compounds = R"(Uniform Compounds
=================
Uniform compounds can be specified by

o   A UC number e.g. uc5

o   uc_ followed by a name (see the list below) Use '_' instead of
    a space to avoid having to quote the model name.

    When giving a name the following abbreviations can be used
       tr:     truncated
       sm:     small
       gr:     great
       st:     stellated
       sn:     snub
       tet:    tetrahedra
       ico:    icosahedra
       icosa:  icosahedra
       dod:    dodecahedra
       oct:    octahedra
       cubo:   cuboctahedra
       icosid: icosidodecahedra
       pri:    prisms
       ant:    antiprisms
       rot:    rotational

    e.g. uc_5_tetrahedra, uc_5_tet, uc_2_tr_tet, uc_20_octahedra_rot

    When a compound is listed as rotational, an angle can be supplied after
    an underscore. If no angle is supplied, a random angle is generated

    e.g. uc2_a30, uc7_a22.5, uc11_a11.315, uc28_a0

    For uc20 to uc25 n/d and k can be supplied after an underscore.
    If no n/d or k is supplied then random values are generated

    e.g. uc21_n5/2k4, uc22_n7/3, uc23_n7/3k3, uc25_n7/4k3

    Because uc20, uc22, and uc24 an angle can also be supplied.
    Note that the order of a, n/d and k does not matter

    e.g. uc20_n9/2, uc20_k5, uc22_n7/3k2a4, uc22_a7n5/3k3, uc24_a3k5n7/4


Uniform Compound List:

UC No  Name
-----  -------------------------------
UC1    6 tetrahedra rotational
UC2    12 tetrahedra rotational
UC3    6 tetrahedra
UC4    2 tetrahedra
UC5    5 tetrahedra
UC6    10 tetrahedra
UC7    6 cubes rotational
UC8    3 cubes
UC9    5 cubes
UC10   4 octahedra rotational
UC11   8 octahedra rotational
UC12   4 octahedra
UC13   20 octahedra rotational
UC14   20 octahedra
UC15   10 octahedra 1
UC16   10 octahedra 2
UC17   5 octahedra
UC18   5 tetrahemihexahedra
UC19   20 tetrahemihexahedra
UC20   2k n d gonal prisms rotational
UC21   k n d gonal prisms
UC22   2k n odd d gonal antiprisms rotational
UC23   k n odd d gonal antiprisms
UC24   2k n even d gonal antiprisms rotational
UC25   k n even d gonal antiprisms
UC26   12 pentagonal antiprisms rotational
UC27   6 pentagonal antiprisms
UC28   12 pentagrammic crossed antiprisms rotational
UC29   6 pentagrammic crossed antiprisms
UC30   4 triangular prisms
UC31   8 triangular prisms
UC32   10 triangular prisms
UC33   20 triangular prisms
UC34   6 pentagonal prisms
UC35   12 pentagonal prisms
UC36   6 pentagrammic prisms
UC37   12 pentagrammic prisms
UC38   4 hexagonal prisms
UC39   10 hexagonal prisms
UC40   6 decagonal prisms
UC41   6 decagrammic prisms
UC42   3 square antiprisms
UC43   6 square antiprisms
UC44   6 pentagrammic antiprisms
UC45   12 pentagrammic antiprisms
UC46   2 icosahedra
UC47   5 icosahedra
UC48   2 great dodecahedra
UC49   5 great dodecahedra
UC50   2 small stellated dodecahedra
UC51   5 small stellated dodecahedra
UC52   2 great icosahedra
UC53   5 great icosahedra
UC54   2 truncated tetrahedra
UC55   5 truncated tetrahedra
UC56   10 truncated tetrahedra
UC57   5 truncated cubes
UC58   5 stellated truncated hexahedra
UC59   5 cuboctahedra
UC60   5 cubohemioctahedra
UC61   5 octahemioctahedra
UC62   5 rhombicuboctahedra
UC63   5 small rhombihexahedra
UC64   5 small cubicuboctahedra
UC65   5 great cubicuboctahedra
UC66   5 great rhombihexahedra
UC67   5 great rhombicuboctahedra
UC68   2 snub cubes
UC69   2 snub dodecahedra
UC70   2 great snub icosidodecahedra
UC71   2 great inverted snub icosidodecahedra
UC72   2 great retrosnub icosidodecahedra
UC73   2 snub dodecadodecahedra
UC74   2 inverted snub dodecadodecahedra
UC75   2 snub icosidodecadodecahedra)";

const char *help_geodesic = R"(Geodesic Spheres
================
A geodesic sphere is specified by a base polyhedron and a pattern
given by two integers. The pattern numbers determined the class of
geodesic sphere, and more specifically the 'stepping' for moving
between original base vertices in the geodesic model.

o   geo_
    optionally followed by a letter to indicate the base polyhedron
       - i: icosahedron (default)
       - o: octahedron
       - t: tetrahdron
    followed by the first pattern number
    optionally followed by a second pattern number (default 0)

Class I patterns are of the form m,0 or 0,n (the pattern frequency
   is m or n, and this figure is also the 'frequency')
   e.g. geo_4, geo_o3, geo_o3_0, geo_t0_3
Class II patterns are of the form m,m (the pattern frequency is m,i
   but the 'frequency' is often described as 2m)
   e.g. geo_4_4, geo_o3_3
Class III patterns are of the form m,n with m>0, n>0 and m=/=n (the
   the pattern frequency is the greatest common divisor of m and n)
   e.g. geo_4_1, geo_o3_2, geo_t4_2)";

const char *help_sym_models = R"(Symmetry Example Models
=======================
These models are made as a symmetric arrangement of arrows. The
number of arrows in the model is the size of the given symmetry
group. To display the symmetry elements view the models with, for
example, 'antiview -s a sym_D3h'

The models are given by:
o   sym_ followed by a Schoenflies symbol from
            Cs  - mirror
            Ci  - inversion
            Cn  - cyclic rotational
            Cnv - cyclic rotational with vertical mirror
            Cnh - cyclic rotational with horizontal mirror
            Dn  - dihedral rotational
            Dnv - dihedral rotational with vertical mirror
            Dnh - dihedral rotational with horizontal mirror
            Sn  - cyclic rotational (n/2-fold) with inversion
            T   - tetrahedral rotational
            Td  - tetrahedral rotational with mirror
            Th  - tetrahedral rotational with inversion
            O   - octahedral rotational
            Oh  - octahedral rotational with mirror
            I   - icosahedral rotational
            Ih  - icosahedral rotational with mirror

   e.g. sym_Cs, sym_C3, sym_D3h, sym_S6, sym_I)";

const char *help_color = R"(Colour
======
Colouring a model usually involves specifying individual colours. A colour
can be specified directly in a number of formats, or can be the name of a
colour that Antiprism recognises, or can be an index number that will
generally be turned into a final colour by looking it up in a colour map.
Antiprism includes some colour maps, and facilities for creating new ones.

See these help topics for more details
   col_val :     colour value formats
   col_names :   colour names
   col_map :     colour map formats and resource colour maps)";

const char *help_color_val = R"(Colour Values
=============
Colour values may be given in a number of formats.

Red, Green, Blue and Alpha components (RGBA):
   3 or 4 decimal numbers separated by commas in the range 0.0-1.0
       (e.g. orange is '1.0,0.5,0.0')
   3 or 4 integers separated by commas in the range 0 - 255
       (e.g. orange is '255,128,0')
   x (or X or #) followed by 3 or 4 pairs of hexadecimal integers
       (e.g. orange is 'xFF8000', 'x' alone, is transparent black and is
        a shortcut for 'invisible') 
Hue, Saturation, Value and Alpha components (HSVA):
   H or h, immediately followed by 3 or 4 decimal numbers separated by
   commas in the range 0.0-1.0. For initial h the hue has range
   0.0 to 360.0
       (e.g. orange is 'H0.07,1.0,1.0', or 'h25,1.0,1.0)
Colour name:
   Common colour names. See 'off_util -H col_names' for the full list.
       (e.g. 'orange')
Colour index:
   A positive integer that can be used as index into a colour map
       (e.g. '5'))";

const char *help_color_names = R"(Named Colours
=============

Colour names can be used to specify a colour given on the command
line or in a colour map file.

The names are those used in X11's rgb.txt. The list of colours is also
given in the Antiprism X11 colour map file included with the
resources. Other names: 'invisible' corresponds to transparent black,
and is used to indicate that display elements should be excluded from
display (a colour value of 'x' can also be used for 'invisible'); 'none'
is a colour with no colour data, setting an element to this colour removes
any colour data for the element

The colour names are (with their integer RGB description):

snow (255,250,250)               ghostwhite (248,248,255)
whitesmoke (245,245,245)         gainsboro (220,220,220)
floralwhite (255,250,240)        oldlace (253,245,230)
linen (250,240,230)              antiquewhite (250,235,215)
papayawhip (255,239,213)         blanchedalmond (255,235,205)
bisque (255,228,196)             peachpuff (255,218,185)
navajowhite (255,222,173)        moccasin (255,228,181)
cornsilk (255,248,220)           ivory (255,255,240)
lemonchiffon (255,250,205)       seashell (255,245,238)
honeydew (240,255,240)           mintcream (245,255,250)
azure (240,255,255)              aliceblue (240,248,255)
lavender (230,230,250)           lavenderblush (255,240,245)
mistyrose (255,228,225)          white (255,255,255)
black (0,0,0)                    darkslategray (47,79,79)
dimgray (105,105,105)            slategray (112,128,144)
lightslategray (119,136,153)     gray (190,190,190)
lightgray (211,211,211)          midnightblue (25,25,112)
navy (0,0,128)                   navyblue (0,0,128)
cornflowerblue (100,149,237)     darkslateblue (72,61,139)
slateblue (106,90,205)           mediumslateblue (123,104,238)
lightslateblue (132,112,255)     mediumblue (0,0,205)
royalblue (65,105,225)           blue (0,0,255)
dodgerblue (30,144,255)          deepskyblue (0,191,255)
skyblue (135,206,235)            lightskyblue (135,206,250)
steelblue (70,130,180)           lightsteelblue (176,196,222)
lightblue (173,216,230)          powderblue (176,224,230)
paleturquoise (175,238,238)      darkturquoise (0,206,209)
mediumturquoise (72,209,204)     turquoise (64,224,208)
cyan (0,255,255)                 lightcyan (224,255,255)
cadetblue (95,158,160)           mediumaquamarine (102,205,170)
aquamarine (127,255,212)         darkgreen (0,100,0)
darkolivegreen (85,107,47)       darkseagreen (143,188,143)
seagreen (46,139,87)             mediumseagreen (60,179,113)
lightseagreen (32,178,170)       palegreen (152,251,152)
springgreen (0,255,127)          lawngreen (124,252,0)
green (0,255,0)                  chartreuse (127,255,0)
mediumspringgreen (0,250,154)    greenyellow (173,255,47)
limegreen (50,205,50)            yellowgreen (154,205,50)
forestgreen (34,139,34)          olivedrab (107,142,35)
darkkhaki (189,183,107)          khaki (240,230,140)
palegoldenrod (238,232,170)      lightgoldenrodyellow (250,250,210)
lightyellow (255,255,224)        yellow (255,255,0)
gold (255,215,0)                 lightgoldenrod (238,221,130)
goldenrod (218,165,32)           darkgoldenrod (184,134,11)
rosybrown (188,143,143)          indianred (205,92,92)
saddlebrown (139,69,19)          sienna (160,82,45)
peru (205,133,63)                burlywood (222,184,135)
beige (245,245,220)              wheat (245,222,179)
sandybrown (244,164,96)          tan (210,180,140)
chocolate (210,105,30)           firebrick (178,34,34)
brown (165,42,42)                darksalmon (233,150,122)
salmon (250,128,114)             lightsalmon (255,160,122)
orange (255,165,0)               darkorange (255,140,0)
coral (255,127,80)               lightcoral (240,128,128)
tomato (255,99,71)               orangered (255,69,0)
red (255,0,0)                    hotpink (255,105,180)
deeppink (255,20,147)            pink (255,192,203)
lightpink (255,182,193)          palevioletred (219,112,147)
maroon (176,48,96)               mediumvioletred (199,21,133)
violetred (208,32,144)           magenta (255,0,255)
violet (238,130,238)             plum (221,160,221)
orchid (218,112,214)             mediumorchid (186,85,211)
darkorchid (153,50,204)          darkviolet (148,0,211)
blueviolet (138,43,226)          purple (160,32,240)
mediumpurple (147,112,219)       thistle (216,191,216)
snow1 (255,250,250)              snow2 (238,233,233)
snow3 (205,201,201)              snow4 (139,137,137)
seashell1 (255,245,238)          seashell2 (238,229,222)
seashell3 (205,197,191)          seashell4 (139,134,130)
antiquewhite1 (255,239,219)      antiquewhite2 (238,223,204)
antiquewhite3 (205,192,176)      antiquewhite4 (139,131,120)
bisque1 (255,228,196)            bisque2 (238,213,183)
bisque3 (205,183,158)            bisque4 (139,125,107)
peachpuff1 (255,218,185)         peachpuff2 (238,203,173)
peachpuff3 (205,175,149)         peachpuff4 (139,119,101)
navajowhite1 (255,222,173)       navajowhite2 (238,207,161)
navajowhite3 (205,179,139)       navajowhite4 (139,121,94)
lemonchiffon1 (255,250,205)      lemonchiffon2 (238,233,191)
lemonchiffon3 (205,201,165)      lemonchiffon4 (139,137,112)
cornsilk1 (255,248,220)          cornsilk2 (238,232,205)
cornsilk3 (205,200,177)          cornsilk4 (139,136,120)
ivory1 (255,255,240)             ivory2 (238,238,224)
ivory3 (205,205,193)             ivory4 (139,139,131)
honeydew1 (240,255,240)          honeydew2 (224,238,224)
honeydew3 (193,205,193)          honeydew4 (131,139,131)
lavenderblush1 (255,240,245)     lavenderblush2 (238,224,229)
lavenderblush3 (205,193,197)     lavenderblush4 (139,131,134)
mistyrose1 (255,228,225)         mistyrose2 (238,213,210)
mistyrose3 (205,183,181)         mistyrose4 (139,125,123)
azure1 (240,255,255)             azure2 (224,238,238)
azure3 (193,205,205)             azure4 (131,139,139)
slateblue1 (131,111,255)         slateblue2 (122,103,238)
slateblue3 (105,89,205)          slateblue4 (71,60,139)
royalblue1 (72,118,255)          royalblue2 (67,110,238)
royalblue3 (58,95,205)           royalblue4 (39,64,139)
blue1 (0,0,255)                  blue2 (0,0,238)
blue3 (0,0,205)                  blue4 (0,0,139)
dodgerblue1 (30,144,255)         dodgerblue2 (28,134,238)
dodgerblue3 (24,116,205)         dodgerblue4 (16,78,139)
steelblue1 (99,184,255)          steelblue2 (92,172,238)
steelblue3 (79,148,205)          steelblue4 (54,100,139)
deepskyblue1 (0,191,255)         deepskyblue2 (0,178,238)
deepskyblue3 (0,154,205)         deepskyblue4 (0,104,139)
skyblue1 (135,206,255)           skyblue2 (126,192,238)
skyblue3 (108,166,205)           skyblue4 (74,112,139)
lightskyblue1 (176,226,255)      lightskyblue2 (164,211,238)
lightskyblue3 (141,182,205)      lightskyblue4 (96,123,139)
slategray1 (198,226,255)         slategray2 (185,211,238)
slategray3 (159,182,205)         slategray4 (108,123,139)
lightsteelblue1 (202,225,255)    lightsteelblue2 (188,210,238)
lightsteelblue3 (162,181,205)    lightsteelblue4 (110,123,139)
lightblue1 (191,239,255)         lightblue2 (178,223,238)
lightblue3 (154,192,205)         lightblue4 (104,131,139)
lightcyan1 (224,255,255)         lightcyan2 (209,238,238)
lightcyan3 (180,205,205)         lightcyan4 (122,139,139)
paleturquoise1 (187,255,255)     paleturquoise2 (174,238,238)
paleturquoise3 (150,205,205)     paleturquoise4 (102,139,139)
cadetblue1 (152,245,255)         cadetblue2 (142,229,238)
cadetblue3 (122,197,205)         cadetblue4 (83,134,139)
turquoise1 (0,245,255)           turquoise2 (0,229,238)
turquoise3 (0,197,205)           turquoise4 (0,134,139)
cyan1 (0,255,255)                cyan2 (0,238,238)
cyan3 (0,205,205)                cyan4 (0,139,139)
darkslategray1 (151,255,255)     darkslategray2 (141,238,238)
darkslategray3 (121,205,205)     darkslategray4 (82,139,139)
aquamarine1 (127,255,212)        aquamarine2 (118,238,198)
aquamarine3 (102,205,170)        aquamarine4 (69,139,116)
darkseagreen1 (193,255,193)      darkseagreen2 (180,238,180)
darkseagreen3 (155,205,155)      darkseagreen4 (105,139,105)
seagreen1 (84,255,159)           seagreen2 (78,238,148)
seagreen3 (67,205,128)           seagreen4 (46,139,87)
palegreen1 (154,255,154)         palegreen2 (144,238,144)
palegreen3 (124,205,124)         palegreen4 (84,139,84)
springgreen1 (0,255,127)         springgreen2 (0,238,118)
springgreen3 (0,205,102)         springgreen4 (0,139,69)
green1 (0,255,0)                 green2 (0,238,0)
green3 (0,205,0)                 green4 (0,139,0)
chartreuse1 (127,255,0)          chartreuse2 (118,238,0)
chartreuse3 (102,205,0)          chartreuse4 (69,139,0)
olivedrab1 (192,255,62)          olivedrab2 (179,238,58)
olivedrab3 (154,205,50)          olivedrab4 (105,139,34)
darkolivegreen1 (202,255,112)    darkolivegreen2 (188,238,104)
darkolivegreen3 (162,205,90)     darkolivegreen4 (110,139,61)
khaki1 (255,246,143)             khaki2 (238,230,133)
khaki3 (205,198,115)             khaki4 (139,134,78)
lightgoldenrod1 (255,236,139)    lightgoldenrod2 (238,220,130)
lightgoldenrod3 (205,190,112)    lightgoldenrod4 (139,129,76)
lightyellow1 (255,255,224)       lightyellow2 (238,238,209)
lightyellow3 (205,205,180)       lightyellow4 (139,139,122)
yellow1 (255,255,0)              yellow2 (238,238,0)
yellow3 (205,205,0)              yellow4 (139,139,0)
gold1 (255,215,0)                gold2 (238,201,0)
gold3 (205,173,0)                gold4 (139,117,0)
goldenrod1 (255,193,37)          goldenrod2 (238,180,34)
goldenrod3 (205,155,29)          goldenrod4 (139,105,20)
darkgoldenrod1 (255,185,15)      darkgoldenrod2 (238,173,14)
darkgoldenrod3 (205,149,12)      darkgoldenrod4 (139,101,8)
rosybrown1 (255,193,193)         rosybrown2 (238,180,180)
rosybrown3 (205,155,155)         rosybrown4 (139,105,105)
indianred1 (255,106,106)         indianred2 (238,99,99)
indianred3 (205,85,85)           indianred4 (139,58,58)
sienna1 (255,130,71)             sienna2 (238,121,66)
sienna3 (205,104,57)             sienna4 (139,71,38)
burlywood1 (255,211,155)         burlywood2 (238,197,145)
burlywood3 (205,170,125)         burlywood4 (139,115,85)
wheat1 (255,231,186)             wheat2 (238,216,174)
wheat3 (205,186,150)             wheat4 (139,126,102)
tan1 (255,165,79)                tan2 (238,154,73)
tan3 (205,133,63)                tan4 (139,90,43)
chocolate1 (255,127,36)          chocolate2 (238,118,33)
chocolate3 (205,102,29)          chocolate4 (139,69,19)
firebrick1 (255,48,48)           firebrick2 (238,44,44)
firebrick3 (205,38,38)           firebrick4 (139,26,26)
brown1 (255,64,64)               brown2 (238,59,59)
brown3 (205,51,51)               brown4 (139,35,35)
salmon1 (255,140,105)            salmon2 (238,130,98)
salmon3 (205,112,84)             salmon4 (139,76,57)
lightsalmon1 (255,160,122)       lightsalmon2 (238,149,114)
lightsalmon3 (205,129,98)        lightsalmon4 (139,87,66)
orange1 (255,165,0)              orange2 (238,154,0)
orange3 (205,133,0)              orange4 (139,90,0)
darkorange1 (255,127,0)          darkorange2 (238,118,0)
darkorange3 (205,102,0)          darkorange4 (139,69,0)
coral1 (255,114,86)              coral2 (238,106,80)
coral3 (205,91,69)               coral4 (139,62,47)
tomato1 (255,99,71)              tomato2 (238,92,66)
tomato3 (205,79,57)              tomato4 (139,54,38)
orangered1 (255,69,0)            orangered2 (238,64,0)
orangered3 (205,55,0)            orangered4 (139,37,0)
red1 (255,0,0)                   red2 (238,0,0)
red3 (205,0,0)                   red4 (139,0,0)
deeppink1 (255,20,147)           deeppink2 (238,18,137)
deeppink3 (205,16,118)           deeppink4 (139,10,80)
hotpink1 (255,110,180)           hotpink2 (238,106,167)
hotpink3 (205,96,144)            hotpink4 (139,58,98)
pink1 (255,181,197)              pink2 (238,169,184)
pink3 (205,145,158)              pink4 (139,99,108)
lightpink1 (255,174,185)         lightpink2 (238,162,173)
lightpink3 (205,140,149)         lightpink4 (139,95,101)
palevioletred1 (255,130,171)     palevioletred2 (238,121,159)
palevioletred3 (205,104,137)     palevioletred4 (139,71,93)
maroon1 (255,52,179)             maroon2 (238,48,167)
maroon3 (205,41,144)             maroon4 (139,28,98)
violetred1 (255,62,150)          violetred2 (238,58,140)
violetred3 (205,50,120)          violetred4 (139,34,82)
magenta1 (255,0,255)             magenta2 (238,0,238)
magenta3 (205,0,205)             magenta4 (139,0,139)
orchid1 (255,131,250)            orchid2 (238,122,233)
orchid3 (205,105,201)            orchid4 (139,71,137)
plum1 (255,187,255)              plum2 (238,174,238)
plum3 (205,150,205)              plum4 (139,102,139)
mediumorchid1 (224,102,255)      mediumorchid2 (209,95,238)
mediumorchid3 (180,82,205)       mediumorchid4 (122,55,139)
darkorchid1 (191,62,255)         darkorchid2 (178,58,238)
darkorchid3 (154,50,205)         darkorchid4 (104,34,139)
purple1 (155,48,255)             purple2 (145,44,238)
purple3 (125,38,205)             purple4 (85,26,139)
mediumpurple1 (171,130,255)      mediumpurple2 (159,121,238)
mediumpurple3 (137,104,205)      mediumpurple4 (93,71,139)
thistle1 (255,225,255)           thistle2 (238,210,238)
thistle3 (205,181,205)           thistle4 (139,123,139)
gray0 (0,0,0)                    gray1 (3,3,3)
gray2 (5,5,5)                    gray3 (8,8,8)
gray4 (10,10,10)                 gray5 (13,13,13)
gray6 (15,15,15)                 gray7 (18,18,18)
gray8 (20,20,20)                 gray9 (23,23,23)
gray10 (26,26,26)                gray11 (28,28,28)
gray12 (31,31,31)                gray13 (33,33,33)
gray14 (36,36,36)                gray15 (38,38,38)
gray16 (41,41,41)                gray17 (43,43,43)
gray18 (46,46,46)                gray19 (48,48,48)
gray20 (51,51,51)                gray21 (54,54,54)
gray22 (56,56,56)                gray23 (59,59,59)
gray24 (61,61,61)                gray25 (64,64,64)
gray26 (66,66,66)                gray27 (69,69,69)
gray28 (71,71,71)                gray29 (74,74,74)
gray30 (77,77,77)                gray31 (79,79,79)
gray32 (82,82,82)                gray33 (84,84,84)
gray34 (87,87,87)                gray35 (89,89,89)
gray36 (92,92,92)                gray37 (94,94,94)
gray38 (97,97,97)                gray39 (99,99,99)
gray40 (102,102,102)             gray41 (105,105,105)
gray42 (107,107,107)             gray43 (110,110,110)
gray44 (112,112,112)             gray45 (115,115,115)
gray46 (117,117,117)             gray47 (120,120,120)
gray48 (122,122,122)             gray49 (125,125,125)
gray50 (127,127,127)             gray51 (130,130,130)
gray52 (133,133,133)             gray53 (135,135,135)
gray54 (138,138,138)             gray55 (140,140,140)
gray56 (143,143,143)             gray57 (145,145,145)
gray58 (148,148,148)             gray59 (150,150,150)
gray60 (153,153,153)             gray61 (156,156,156)
gray62 (158,158,158)             gray63 (161,161,161)
gray64 (163,163,163)             gray65 (166,166,166)
gray66 (168,168,168)             gray67 (171,171,171)
gray68 (173,173,173)             gray69 (176,176,176)
gray70 (179,179,179)             gray71 (181,181,181)
gray72 (184,184,184)             gray73 (186,186,186)
gray74 (189,189,189)             gray75 (191,191,191)
gray76 (194,194,194)             gray77 (196,196,196)
gray78 (199,199,199)             gray79 (201,201,201)
gray80 (204,204,204)             gray81 (207,207,207)
gray82 (209,209,209)             gray83 (212,212,212)
gray84 (214,214,214)             gray85 (217,217,217)
gray86 (219,219,219)             gray87 (222,222,222)
gray88 (224,224,224)             gray89 (227,227,227)
gray90 (229,229,229)             gray91 (232,232,232)
gray92 (235,235,235)             gray93 (237,237,237)
gray94 (240,240,240)             gray95 (242,242,242)
gray96 (245,245,245)             gray97 (247,247,247)
gray98 (250,250,250)             gray99 (252,252,252)
gray100 (255,255,255)            darkgray (169,169,169)
darkblue (0,0,139)               darkcyan (0,139,139)
darkmagenta (139,0,139)          darkred (139,0,0)
lightgreen (144,238,144))";

const char *help_ColorMap = R"(Colour Maps
===========
OFF elements may be coloured by index numbers. Colour maps provide a
way to convert these index numbers into colour values and are used
in command options such as off_color -m, antiview -m and n_icons -m.

The maps may use the Antiprism colour map format, Gimp Palette format
or Fractint format. The Antiprism format has lines of the form
index_number = color_value # comment_text, e.g. '2 = 1.0,0.0,0.0 # red'
anything after # is a comment and ignored. Blank lines are ignored

A map may be given by several maps separated by ',' e.g. 'map1,map2'.
The maps are tried in order until a conversion is found for an index
number

Map modifiers
-------------
A map may be modified by remapping its own entries. The modifiers
can be given in any order but the general form is

   map_name+shift*step%wrap

Any index, idx, is mapped to the colour value with index
   wrap is 0:      shift + idx*step
   wrap is not 0: (shift + idx*step) % wrap   [where % is modulus]
The defaults are shift=0, step=1, wrap=0. If a bare % is given then the
the wrap value will be the largest index number in the map plus one
e.g. cmap+1  : for index 10 get the colour value for cmap index 11
     cmap*2  : for index 10 get the colour value for cmap index 20
     cmap%6 : for index 10 get the colour value for cmap index 4

Resource Maps
-------------
Internal (see below for format):
   spread
      Gives a range of colours each differing from the last few. Useful
      to colour elements whose index numbers have been set sequentially.
   map
      makes a colour map on the command line. Map entries are separated
      by ':' and each entry corresponds to a line in an Antiprism format
      map. Colours given by components may also have the components
      separated by '/' (as spaces will require quoting or escaping).
   rnd, rand, random
      A random map, with colours selected within certain ranges
      (default: component ranges H0:1S0.7:1V0.7:1).
   rng, range
      A map made by ranging between component values
      (default: size 256, component ranges H0:1S0.9V0.9).
   remap
      A map of index numbers to themselves. Use with the map modifiers
      to remap index numbers.
   null
      An empty map.
   deal
      A map (default: size 256) containing a random shuffle of the
      values 0 to packsize-1, packsize is the same as size by default,
      but can be changed by adding _packsize (sequential deals are used
      if this is less than size), e.g. deal100, deal_3 
   grey, greyw
      greyscales (default: size 256), grey runs from black to white
      and greyw is wrappable and runs from black to white to black again.
   rainbow
      A seven color rainbow map (default: size 256), which include a brightness
      _value range from -1 (all black) to 1 (all white), e.g. rainbow_0.5
   uniform
      used to colour the uniform, Johnson and polygon-based resource
      models (applied with off_color -f A -m uniform)
   compound
      used to colour the uniform compound resource models (applied
      with off_color -f K -v F -e F -m compound)

   These maps contain a mapping for every index number (except default for
   rng is 256 entries). Follow the map name immediately by a number to set
   a particular size for the map (e.g. rnd64). The random and range maps
   (with the optional size specification) can be followed by '_' and then
   a letter from HSVA or RGBA (upper or lower case) to specify a component.
   This is followed by floating point numbers numbers separated by ':'.
   These numbers must be in the range 0.0 to 1.0. H may be larger than
   1.0 to allow hue ranges that cross the 1.0/0.0 boundary. Lower case 'h'
   gives the hue in the range 0-360. The random map allows one number for
   a fixed value, or two numbers for a range. The range map places the
   components at equal steps throughout the map and interpolates between
   the values
   e.g. rnd - default random map
        rnd64 - random map with 64 entries
        rnd256_S0V0:1 - random map with 256 grey entries
        rng - default range map, like a rainbow map with 256 entries
        rng16 - rainbow map with 16 entries
        rng_S0V0:1 - greyscale with 256 entries
        rng_R0:1:1G1:1:0B1:0:1 - runs from cyan to yellow to magenta.

External (in resource directory 'col_maps'):
   x11 (549 colours)
       Broad range of named colours from X11's rgb.txt. The names may
       also be used to specify colours on the command line.
   vga (16 colours)
       Original VGA named colours.
   html (140 colours)
       Broad range of colours. These are the colours that can be used
       by name in HTML.
   ms (48 colours)
       A colour map based on the Microsoft colour dialog.
   iscc (267 colours)
       Colour centroids (http://tx4.us/nbs-iscc.htm)
   rainbowc (192 colours)
       A rainbow map, with cyan but not green
   rainbowg (192 colours)
       A rainbow map, with green but not cyan
   spectrum (401 colours)
       An approximate visible spectrum)";

const char *help_symmetry = R"(Symmetry
========
Several programs have features involving symmetry:
   poly_kscope - repeats a model symmetrically, used to make compounds
   off_align - repeats a model symmetrically, used to augment polyhedra
   off_trans - aligns a model according to symmetry
   off_report - prints symmetry information for a model
   off_color - colours by symmetry orbit
   antiview, off2vrml, off2pov - displays symmetry elements

The program option parameters are organised around the following ideas

  Symmetry group
    A set (group) of (Euclidean) transformations that carry a
    polyhedron onto itself, described in general form using
    Schoenflies notation (see below) e.g. Oh, D3v.

  Full symmetry group of a polyhedron
    The set of all (Euclidean) transformations that carry a
    polyhedron onto itself.

  Symmetry subgroup, or subsymmetry
    A set of transformations from a symmetry group which, considered
    alone, also form a symmetry group, e.g. a cube has Oh symmetry
    and has a 3-fold axis corresponding to a C3v subgroup.

  Symmetry orbit of an element
    A set of equivalent elements, those elements that this element is
    carried on to by a symmetry or subsymmetry of the model.

  Standard alignment of a symmetry
    A symmetry group could be aligned anywhere in the coordinate
    system, but there are particular alignments that fit nicely
    with the coordinate axes, and these are used as the 'standard'
    alignments in Antiprism. They are a way of associating a symbol
    like D3v with a fixed set of transformations.

  Conjugation subtype of a subsymmetry
    This is an integer used to distinguish subgroups which are
    not carried onto each other (by conjugation) by the
    transformations of the parent symmetry group. For example a
    cube has a 2 fold axis through mid-edge and a 2-fold axis
    through a face centre. There is no symmetry of the cube that
    carries one onto the other and so they will have different
    subtype numbers. Geometrically, they look different in the cube.

  Symmetry realignment
    If you align a polyhedron with the standard set of symmetries
    for its full symmetry group there is often more than one
    distinct way to achieve this (a transformation not in the
    symmetry group that transforms the symmetry group onto
    itself). For example, if you align a cereal box-like cuboid
    naturally with the coordinate axes there are 6 possibilities
    i.e. the centres of the three rectangle types can lie on any
    of the axes, with 3x2x1 = 6. Possibilities for some polyhedra
    are infinite e.g. the symmetry group of a pyramid does not change
    when it is translated along its principal axis. The realignment
    is given by a series of colon separated numbers, the first
    number selects from a finite set of realignments, and the
    following numbers are decimals to control rotations and
    translations as follows
      axial rotation:    1 number  - degrees around principle axis
      full rotation:     3 numbers - degrees around x, y and z axes
      axial translation: 1 number  - distance to translate along principal
                                     axis
      plane translation: 2 numbers - distance to translate along two
                                     orthogonal directions in (mirror) plane
      full translation:  3 numbers - distance to translate along x, y and z
                                     axes

  Schoenflies notation
    Used to specify symmetry groups. The standard alignments have,
    preferentially, a centre (fixed point) on the origin a principal
    rotational axes on the z-axis, a dihedral axis on the x-axis, a
    mirror normal on the y-axis (except Cs has a mirror normal on the
    z-axis). The polyhedral symmetry types have a 3-fold axis on (1,1,1).
    In the following list of symbols, when a type contains 'n' this must
    be replaced by an integer (giving an n-fold axis), and for S this
    integer must be even.
         C1  - identity
         Cs  - mirror
         Ci  - inversion
         Cn  - cyclic rotational
         Cnv - cyclic rotational with vertical mirror
         Cnh - cyclic rotational with horizontal mirror
         Dn  - dihedral rotational
         Dnv - dihedral rotational with vertical mirror
         Dnh - dihedral rotational with horizontal mirror
         Sn  - cyclic rotational (n/2-fold) with inversion
         T   - tetrahedral rotational
         Td  - tetrahedral rotational with mirror
         Th  - tetrahedral rotational with inversion
         O   - octahedral rotational
         Oh  - octahedral rotational with mirror
         I   - icosahedral rotational
         Ih  - icosahedral rotational with mirror)";

const char *help_bowers = R"(Bowers Short Name Notations
===========================
Jonathan Bowers has created short names for the Uniforms, Uniform Compounds,
Johnson polyhedra, and prisms. These notations are supported by Antiprism.

Short Name   Symbol and Name
----------   ---------------
tet          U1   tetrahedron
tut          U2   truncated tetrahedron
oho          U3   octahemioctahedron
thah         U4   tetrahemihexahedron
oct          U5   octahedron
cube         U6   cube
co           U7   cuboctahedron
toe          U8   truncated octahedron
tic          U9   truncated cube
sirco        U10  rhombicuboctahedron
girco        U11  truncated cuboctahedron
snic         U12  snub cube
socco        U13  small cubicuboctahedron
gocco        U14  great cubicuboctahedron
cho          U15  cubohemioctahedron
cotco        U16  cubitruncated cuboctahedron
querco       U17  great rhombicuboctahedron
sroh         U18  small rhombihexahedron
quith        U19  stellated truncated hexahedron
quitco       U20  great truncated cuboctahedron
groh         U21  great rhombihexahedron
ike          U22  icosahedron
doe          U23  dodecahedron
id           U24  icosidodecahedron
ti           U25  truncated icosahedron
tid          U26  truncated dodecahedron
srid         U27  rhombicosidodecahedron
grid         U28  truncated icosidodecahedron
snid         U29  snub dodecahedron
sidtid       U30  small ditrigonal icosidodecahedron
siid         U31  small icosicosidodecahedron
seside       U32  small snub icosicosidodecahedron
saddid       U33  small dodecicosidodecahedron
sissid       U34  small stellated dodecahedron
gad          U35  great dodecahedron
did          U36  dodecadodecahedron
tigid        U37  truncated great dodecahedron
raded        U38  rhombidodecadodecahedron
sird         U39  small rhombidodecahedron
siddid       U40  snub dodecadodecahedron
ditdid       U41  ditrigonal dodecadodecahedron
gidditdid    U42  great ditrigonal dodecicosidodecahedron
sidditdid    U43  small ditrigonal dodecicosidodecahedron
ided         U44  icosidodecadodecahedron
idtid        U45  icositruncated dodecadodecahedron
sided        U46  snub icosidodecadodecahedron
gidtid       U47  great ditrigonal icosidodecahedron
giid         U48  great icosicosidodecahedron
seihid       U49  small icosihemidodecahedron
siddy        U50  small dodecicosahedron
sidhid       U51  small dodecahemidodecahedron
gissid       U52  great stellated dodecahedron
gike         U53  great icosahedron
gid          U54  great icosidodecahedron
tiggy        U55  great truncated icosahedron
ri           U56  rhombicosahedron
gosid        U57  great snub icosidodecahedron
quitsissid   U58  small stellated truncated dodecahedron
quitdid      U59  truncated dodecadodecahedron
isdid        U60  inverted snub dodecadodecahedron
gaddid       U61  great dodecicosidodecahedron
sidhei       U62  small dodecahemicosahedron
giddy        U63  great dodecicosahedron
gisdid       U64  great snub dodecicosidodecahedron
gidhei       U65  great dodecahemicosahedron
quitgissid   U66  great stellated truncated dodecahedron
qrid         U67  great rhombicosidodecahedron
gaquatid     U68  great truncated icosidodecahedron
gisid        U69  great inverted snub icosidodecahedron
gidhid       U70  great dodecahemidodecahedron
geihid       U71  great icosihemidodecahedron
sirsid       U72  small retrosnub icosicosidodecahedron
gird         U73  great rhombidodecahedron
girsid       U74  great retrosnub icosidodecahedron
gidrid       U75  great dirhombicosidodecahedron
pip          U76  pentagonal prism
pap          U77  pentagonal antiprism
stip         U78  pentagrammic prism
stap         U79  pentagrammic antiprism
starp        U80  pentagrammic crossed antiprism

sis          UC1   6 tetrahedra rotational
dis          UC2   12 tetrahedra rotational
snu          UC3   6 tetrahedra
so           UC4   2 tetrahedra
ki           UC5   5 tetrahedra
e            UC6   10 tetrahedra
risdoh       UC7   6 cubes rotational
rah          UC8   3 cubes
rhom         UC9   5 cubes
dissit       UC10  4 octahedra rotational
doso         UC11  8 octahedra rotational
sno          UC12  4 octahedra
addasi       UC13  20 octahedra rotational
dasi         UC14  20 octahedra
gissi        UC15  10 octahedra 1
si           UC16  10 octahedra 2
se           UC17  5 octahedra
hirki        UC18  5 tetrahemihexahedra
sapisseri    UC19  20 tetrahemihexahedra
gadsid       UC26  12 pentagonal antiprisms rotational
gassid       UC27  6 pentagonal antiprisms
gidasid      UC28  12 pentagrammic crossed antiprisms rotational
gissed       UC29  6 pentagrammic crossed antiprisms
ro           UC30  4 triangular prisms
dro          UC31  8 triangular prisms
kri          UC32  10 triangular prisms
dri          UC33  20 triangular prisms
kred         UC34  6 pentagonal prisms
dird         UC35  12 pentagonal prisms
gikrid       UC36  6 pentagrammic prisms
giddird      UC37  12 pentagrammic prisms
griso        UC38  4 hexagonal prisms
rosi         UC39  10 hexagonal prisms
rassid       UC40  6 decagonal prisms
grassid      UC41  6 decagrammic prisms
gassic       UC42  3 square antiprisms
gidsac       UC43  6 square antiprisms
sassid       UC44  6 pentagrammic antiprisms
sadsid       UC45  12 pentagrammic antiprisms
siddo        UC46  2 icosahedra
sne          UC47  5 icosahedra
presipsido   UC48  2 great dodecahedra
presipsi     UC49  5 great dodecahedra
passipsido   UC50  2 small stellated dodecahedra
passipsi     UC51  5 small stellated dodecahedra
sirsido      UC52  2 great icosahedra
sirsei       UC53  5 great icosahedra
tisso        UC54  2 truncated tetrahedra
taki         UC55  5 truncated tetrahedra
te           UC56  10 truncated tetrahedra
harie        UC57  5 truncated cubes
quahri       UC58  5 stellated truncated hexahedra
arie         UC59  5 cuboctahedra
gari         UC60  5 cubohemioctahedra
iddei        UC61  5 octahemioctahedra
rasseri      UC62  5 rhombicuboctahedra
rasher       UC63  5 small rhombihexahedra
rahrie       UC64  5 small cubicuboctahedra
raquahri     UC65  5 great cubicuboctahedra
rasquahr     UC66  5 great rhombihexahedra
rasquahpri   UC67  5 great rhombicuboctahedra
disco        UC68  2 snub cubes
dissid       UC69  2 snub dodecahedra
giddasid     UC70  2 great snub icosidodecahedra
gidsid       UC71  2 great inverted snub icosidodecahedra
gidrissid    UC72  2 great retrosnub icosidodecahedra
disdid       UC73  2 snub dodecadodecahedra
idisdid      UC74  2 inverted snub dodecadodecahedra
desided      UC75  2 snub icosidodecadodecahedra

squippy      J1   square pyramid
peppy        J2   pentagonal pyramid
tricu        J3   triangular cupola
squicu       J4   square cupola
pecu         J5   pentagonal cupola
pero         J6   pentagonal rotunda
etripy       J7   elongated triangular pyramid
esquippy     J8   elongated square pyramid
epeppy       J9   elongated pentagonal pyramid
gyesp        J10  gyroelongated square pyramid
gyepip       J11  gyroelongated pentagonal pyramid
tridpy       J12  triangular dipyramid
pedpy        J13  pentagonal dipyramid
etidpy       J14  elongated triangular dipyramid
esquidpy     J15  elongated square dipyramid
epedpy       J16  elongated pentagonal dipyramid
gyesqidpy    J17  gyroelongated square dipyramid
etcu         J18  elongated triangular cupola
escu         J19  elongated square cupola
epcu         J20  elongated pentagonal cupola
epro         J21  elongated pentagonal rotunda
gyetcu       J22  gyroelongated triangular cupola
gyescu       J23  gyroelongated square cupola
gyepcu       J24  gyroelongated pentagonal cupola
gyepro       J25  gyroelongated pentagonal rotunda
gybef        J26  gyrobifastigium
tobcu        J27  triangular orthobicupola
squobcu      J28  square orthobicupola
squigybcu    J29  square gyrobicupola
pobcu        J30  pentagonal orthobicupola
pegybcu      J31  pentagonal gyrobicupola
pocuro       J32  pentagonal orthocupolarotunda
pegycuro     J33  pentagonal gyrocupolarotunda
pobro        J34  pentagonal orthobirotunda
etobcu       J35  elongated triangular orthobicupola
etigybcu     J36  elongated triangular gyrobicupola
esquigybcu   J37  elongated square gyrobicupola
epobcu       J38  elongated pentagonal orthobicupola
epigybcu     J39  elongated pentagonal gyrobicupola
epocuro      J40  elongated pentagonal orthocupolarotunda
epgycuro     J41  elongated pentagonal gyrocupolarotunda
epobro       J42  elongated pentagonal orthobirotunda
epgybro      J43  elongated pentagonal gyrobirotunda
gyetibcu     J44  gyroelongated triangular bicupola
gyesquibcu   J45  gyroelongated square bicupola
gyepibcu     J46  gyroelongated pentagonal bicupola
gyepcuro     J47  gyroelongated pentagonal cupolarotunda
gyepabro     J48  gyroelongated pentagonal birotunda
autip        J49  augmented triangular prism
bautip       J50  biaugmented triangular prism
tautip       J51  triaugmented triangular prism
aupip        J52  augmented pentagonal prism
baupip       J53  biaugmented pentagonal prism
auhip        J54  augmented hexagonal prism
pabauhip     J55  parabiaugmented hexagonal prism
mabauhip     J56  metabiaugmented hexagonal prism
tauhip       J57  triaugmented hexagonal prism
aud          J58  augmented dodecahedron
pabaud       J59  parabiaugmented dodecahedron
mabaud       J60  metabiaugmented dodecahedron
taud         J61  triaugmented dodecahedron
mibdi        J62  metabidiminished icosahedron
teddi        J63  tridiminished icosahedron
auteddi      J64  augmented tridiminished icosahedron
autut        J65  augmented truncated tetrahedron
autic        J66  augmented truncated cube
bautic       J67  biaugmented truncated cube
autid        J68  augmented truncated dodecahedron
pabautid     J69  parabiaugmented truncated dodecahedron
mabautid     J70  metabiaugmented truncated dodecahedron
tautid       J71  triaugmented truncated dodecahedron
gyrid        J72  gyrate rhombicosidodecahedron
pabgyrid     J73  parabigyrate rhombicosidodecahedron
mabgyrid     J74  metabigyrate rhombicosidodecahedron
tagyrid      J75  trigyrate rhombicosidodecahedron
dirid        J76  diminished rhombicosidodecahedron
pagydrid     J77  paragyrate diminished rhombicosidodecahedron
magydrid     J78  metagyrate diminished rhombicosidodecahedron
bagydrid     J79  bigyrate diminished rhombicosidodecahedron
pabidrid     J80  parabidiminished rhombicosidodecahedron
mabidrid     J81  metabidiminished rhombicosidodecahedron
gybadrid     J82  gyrate bidiminished rhombicosidodecahedron
tedrid       J83  tridiminished rhombicosidodecahedron
snadow       J84  snub disphenoid
snisquap     J85  snub square antiprism
waco         J86  sphenocorona
auwaco       J87  augmented sphenocorona
wamco        J88  sphenomegacorona
hawmco       J89  hebesphenomegacorona
dawci        J90  disphenocingulum
bilbiro      J91  bilunabirotunda
thawro       J92  triangular hebesphenorotunda

trip         pri3     triangular prism
pip          pri5     pentagonal prism
stip         pri5/2   pentagrammic prism
hip          pri6     hexagonal prism
hep          pri7     heptagonal prism
ship         pri7/2   heptagrammic prism
giship       pri7/3   heptagrammic prism
op           pri8     octagonal prism
stop         pri8/3   octagrammic prism
ep           pri9     nonagonal prism
step         pri9/2   nonagrammic prism
gistep       pri9/4   nongrammic prism
dip          pri10    decagonal prism
stiddip      pri10/3  decagrammic prism

squap        ant4     square antiprism
pap          ant5     pentagonal antiprism
stap         ant5/2   pentagrammic antiprism
starp        ant5/3   pentagrammic crossed antiprism
hap          ant6     hexagonal antiprism
heap         ant7     heptagonal antiprism
sthap        ant7/2   heptagrammic antiprism
gisthap      ant7/3   heptagrammic antiprism
gisthirp     ant7/4   heptagrammic crossed antiprism
oap          ant8     octagonal antiprism
stoap        ant8/3   octagrammic antiprism
storp        ant8/5   octagrammic crossed antiprism
eap          ant9     nonagonal antiprism
steap        ant9/2   nonagonal antiprism
gisteap      ant9/4   nonagonal antiprism
gisterp      ant9/5   nonagonal crossed antiprism
dap          ant10    decagonal antiprism
stiddap      ant10/3  decagrammic antiprism)";

const char *help_schwarz = R"(Schwarz triangles
=================
A Schwarz triangle, aligned with with the corresponding minimal symmetry
group, can be specified as follows:

schwarz_ followed by the three fractions separated by spaces or '_'. The
fractions are in the form used in the Wythoff symbol. If a triangle has
an angle of PIn/d, its fraction is d/n, e.g. the Mobius triangle for
octahedral symmetry has angles PI/2, PI/3 and PI/4 and is given by
schwarz_2_3_4, a spherical tetrahedron face has angles of 2PI/3 and is
given by schwarz_3/2_3/2_3/2. Add 'p' to repeat the triangle to make a
polyhedron, e.g. schwarz_5_3_2p.)";

const char *help_expreval = R"(Expression Evaluation
=====================
Any program that accepts a floating point argument will also accept the
argument as a mathematical expression. This expression will be used to
calculate a numeric value, that will then passed to the program.

Several of the characters that can appear in an expression may be
intercepted by the shell for its own use, and the expression will
often need to be enclosed in quotes, e.g.

   off_trans -T '1+sqrt(2)',0,0 cube

Arithmetic operators, + - * / ^, may be used, and various functions,
including: sqrt(), sin(), cos(), tan(), asin(), acos(), atan(), deg().
The full list of functions is included further below.

The following constants are defined: convenience square roots of the
form rt2, rt3, rt5, etc, pi = 3.14159..., phi = (sqrt(5)+1)/2 = 1.61803...

A limited set of variables, var0, var1, ...var9, may be used in
expressions (default value 0.0). Assign a value to a variable with
'=', and separate sub-expressions with ';'. The result of the last
sub-expression is returned as the result of the whole expression.

Examples:
   Expresion                  Result
   1                            1.0
   1+1                          2.0
   1+2*3                        7.0
   (1+2)*3                      9.0
   sin(30)                      0.5
   deg(pi/3)                   60.0
   sum(1;2;3)                   6.0
   var0=2;var1=3;var0^var1      8.0


The code that handles the expression evaluation is muParse, written
by Ingo Berg: http://muparser.beltoforion.de/ . The following is from
the muParser project documentation (with slight changes):

Built-in functions
The following table gives an overview of the functions supported by the
default implementation. It lists the function names, the number of
arguments and a brief description. Separate function arguments with ';'.

   Name    Argc.  Explanation
   sin     1      sine function (degrees)
   cos     1      cosine function (degrees)
   tan     1      tangens function (degrees)
   asin    1      arcus sine function (degrees)
   acos    1      arcus cosine function (degrees)
   atan    1      arcus tangens function (degrees)
   sinh    1      hyperbolic sine function
   cosh    1      hyperbolic cosine
   tanh    1      hyperbolic tangens function
   asinh   1      hyperbolic arcus sine function
   acosh   1      hyperbolic arcus tangens function
   atanh   1      hyperbolic arcur tangens function
   log2    1      logarithm to the base 2
   log10   1      logarithm to the base 10
   log     1      logarithm to the base 10
   ln      1      logarithm to base e (2.71828...)
   exp     1      e raised to the power of x
   sqrt    1      square root of a value
   sign    1      sign function -1 if x<0; 1 if x>0
   rint    1      round to nearest integer
   deg     1      convert from radians to degrees
   rad     1      convert from degrees to radians
   abs     1      absolute value
   min     var.   min of all arguments
   max     var.   max of all arguments
   sum     var.   sum of all arguments
   avg     var.   mean value of all arguments

Built-in binary operators
The following table lists the default binary operators supported by the
parser.

   Operator  Meaning                   Priority
   =         assignment (*see below)   -1
   &&        logical and                1
   ||        logical or                 2
   <=        less or equal              4
   >=        greater or equal           4
   !=        not equal                  4
   ==        equal                      4
   >         greater than               4
   <         less than                  4
   +         addition                   5
   -         subtraction                5
   *         multiplication             6
   */        division                   6
   ^         raise x to the power of y  7
*The assignment operator is special since it changes one of its "
arguments
and can only by applied to variables.

Other operators
muParser has built in support for the if then else operator. It uses
lazy evaluation in order to make sure only the necessary branch of the
expression is evaluated.

   Operator    Meaning                 Remarks
   ?:          if then else operator   C++ style syntax)";

#endif // HELP_H
