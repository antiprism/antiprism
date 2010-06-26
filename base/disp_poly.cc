/*
   Copyright (c) 2003-2008, Adrian Rossiter

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

/* \file disp_poly.cc
   \brief display a polyhedron as plane faces, edge rods and vertex balls and
   with element number labels
*/


#include <string.h>

#include "utils.h"
#include "math_utils.h"
#include "scene.h"
#include "symmetry.h"
#include "disp_poly.h"
#include "vrml_writer.h"
#include "pov_writer.h"
#include "gl_writer.h"
      

disp_poly::disp_poly(): triangulate(true), face_alpha(-1), use_lines(false)  
{ 
}

col_val disp_poly::get_def_v_col()
{
   col_val col = v().get_col();
   if(col.is_idx())
      col = get_v_clrng().get_col(col.get_idx());

   return col.is_val() ? col : col_val(1.0,0.5,0.0);
}

col_val disp_poly::get_def_e_col()
{
   col_val col = e().get_col();
   if(col.is_idx())
      col = get_e_clrng().get_col(col.get_idx());

   return col.is_val() ? col : col_val(0.8,0.6,0.8);
}

col_val disp_poly::get_def_f_col()
{
   col_val col = f().get_col();
   if(col.is_idx())
      col = get_f_clrng().get_col(col.get_idx());

   return col.is_val() ? col : col_val(0.8,0.9,0.9);
}

void disp_poly::geom_changed()
{ 
   disp_geom = sc_geom->get_geom();
   vector<int> face_map;
   if(triangulate)
      disp_geom.triangulate(col_val::invisible, &face_map);
   else {
      face_map.resize(sc_geom->get_geom().faces().size());
      for(unsigned int i=0; i<face_map.size(); i++)
         face_map[i] = i;
   }
   face_map.push_back(disp_geom.faces().size()); // add end marker

   if(face_alpha>0) {
      for(unsigned int i=0; i<face_map.size()-1; i++) {
         col_val col = sc_geom->get_geom().get_f_col(i);
         col = col_val(col[0], col[1], col[2], face_alpha<256?face_alpha:255);
         for(int f_idx=face_map[i]; f_idx<face_map[i+1]; f_idx++)
            disp_geom.set_f_col(f_idx, col);
      }
   }
}


void disp_poly::set_triangulate(bool tri)
{
   if(triangulate!=tri) {
      triangulate=tri;
      geom_changed();
   }

}

void disp_poly::set_face_alpha(int alpha)
{ 
   if(face_alpha!=alpha) {
      face_alpha = alpha;
      geom_changed();
   }
}


void disp_poly::animate()
{
   for(int i=0; i<3; i++) {
      unsigned int msecs = clrngs[i].get_cycle_msecs();
      if(msecs && cmap_tmrs[i].finished()) {
         clrngs[i].cycle_map_cols();
         cmap_tmrs[i].inc_timer(msecs/1000.0);
      }
   }   
}


// --------------------------------------------------------------
// disp_poly - vrml 

void disp_poly::vrml_protos(FILE *ofile)
{
   col_val vcol = get_def_v_col();
   fprintf(ofile, 
"\n"
"PROTO V_%s [\n"
"   field SFVec3f C 0 0 0    # centre\n"
"   field SFColor clr %s     # colour\n"
"   field SFFloat trn %.4f     # transparency\n"
"]\n"
"{\n"
"   Transform {\n"
"      translation IS C\n"
"      children [\n"
"         Shape {\n"
"            appearance Appearance {\n"
"               material Material {\n"
"                  diffuseColor IS clr\n"
"                  transparency IS trn\n"
"               }\n"
"            }\n"
"            geometry Sphere {\n"
"               radius %g\n"
"            }\n"
"         }\n"
"      ]\n"
"   }\n"
"}\n", dots2underscores(sc_geom->get_name()).c_str(),
       vrml_col(vcol).c_str(), vcol.get_transd(),
       get_vert_rad());

   col_val ecol = get_def_e_col();
   fprintf(ofile, 
"\n"
"PROTO E_%s [\n"
"   field SFVec3f C 0 0 0     # centre\n"
"   field SFRotation R 1 0 0 0     # rotation\n"
"   field SFColor clr %s    # colour\n"
"   field SFFloat trn %.4f     # transparency\n"
"   field SFFloat rad %g    # radius\n"
"   field SFFloat ht 1   # height\n"
"]\n"
"{\n"
"   Transform {\n"
"      translation IS C\n"
"      rotation IS R\n"
"      children [\n"
"         Shape {\n"
"            appearance Appearance {\n"
"               material Material {\n"
"                  diffuseColor IS clr\n"
"                  transparency IS trn\n"
"               }\n"
"            }\n"
"            geometry Cylinder {\n"
"               radius IS rad\n"
"               height IS ht\n"
"            }\n"
"         }\n"
"      ]\n"
"   }\n"
"}\n", dots2underscores(sc_geom->get_name()).c_str(),
       vrml_col(ecol).c_str(), ecol.get_transd(),
       get_edge_rad());


   col_val fcol = get_def_f_col();
   fprintf(ofile, 
"\n"
"PROTO F0_%s [\n"
"   field MFInt32 ci [0 0 0 -1]  # coordinate index node\n"
"   field SFNode vc NULL         # coords\n"
"   field SFColor clr %s  # colour\n"
"   field SFFloat trn %.4f  # transparency\n"
"]\n"
"{\n"
"   Shape {\n"
"      appearance Appearance {\n"
"         material Material {\n"
"            diffuseColor IS clr\n"
"            transparency IS trn\n"
"         }\n"
"      }\n"
"      geometry IndexedFaceSet  {\n"
"         colorPerVertex FALSE\n"
"         solid FALSE\n"
"         coord IS vc\n"
"         coordIndex IS ci\n"
"      }\n"
"  }\n"
"}\n", dots2underscores(sc_geom->get_name()).c_str(),
   vrml_col(fcol).c_str(), fcol.get_transd());

   fprintf(ofile, 
"\n"
"PROTO F_%s [\n"
"   field MFInt32 ci [0 0 0 -1]  # coordinate index node\n"
"   field SFNode vc NULL         # coords\n"
"   field MFColor clrs [0 0 0]  # colours\n"
//"   field SFColor clr %s  # colour\n"
"   field SFFloat trn %.4f  # transparency\n"
"]\n"
"{\n"
"   Shape {\n"
"      appearance Appearance {\n"
"         material Material {\n"
//"            diffuseColor IS clr\n"
"            transparency IS trn\n"
"         }\n"
"      }\n"
"      geometry IndexedFaceSet  {\n"
"         colorPerVertex FALSE\n"
"         solid FALSE\n"
"         coord IS vc\n"
"         coordIndex IS ci\n"
"         color Color { color IS clrs }\n"
"      }\n"
"  }\n"
"}\n", dots2underscores(sc_geom->get_name()).c_str(),
       /*vrml_col(fcol).c_str(),*/
       fcol.get_transd());
}

void disp_poly::vrml_trans_begin(FILE *ofile, const scene &scen)
{
   fprintf(ofile,
         "# scene transformations\n"
         "Transform {\n"
         "   translation %s\n"
         "   children [\n\n"
         "# forget indentation and carry on...\n\n",
   vrml_vec(-scen.cur_camera().get_lookat()).c_str());
}

 
void disp_poly::vrml_coords(FILE *ofile, int sig_digits)
{
   fprintf(ofile, "# Vertex Coordinates\n");
   fprintf(ofile,
"Shape {\n"
"   geometry IndexedFaceSet {\n"
"      coord DEF CRDS Coordinate {\n"
"         point [\n");

   const vector<vec3d> &vs = disp_geom.verts();
   if(vs.size()) {
      for(unsigned int i=0; i<vs.size()-1; i++)
         fprintf(ofile, "\t%s,\n", vrml_vec(vs[i], sig_digits).c_str());
   fprintf(ofile, "\t%s\n", vrml_vec(vs[vs.size()-1], sig_digits).c_str());
   }
       
   fprintf(ofile,
"         ]\n"\
"      }\n"\
"   }\n"\
"}\n\n");
}


void disp_poly::vrml_verts(FILE *ofile, int sig_digits)
{
   fprintf(ofile, "# Vertex elements\n");
   
   const vector<vec3d> &vs = disp_geom.verts();
   for(unsigned int i=0; i<sc_geom->get_geom().verts().size(); i++) {
      if(disp_geom.get_v_col(i).is_inv())
         continue;
      fprintf(ofile, "V_%s { C %s ",
            dots2underscores(sc_geom->get_name()).c_str(),
            vrml_vec(vs[i], sig_digits).c_str());
      col_val col = disp_geom.get_v_col(i);
      if(col.is_idx())
         col = get_f_clrng().get_col(col.get_idx());
      if(col.is_val())
         fprintf(ofile, "clr %s trn %.4f", vrml_col(col).c_str(),
               col.get_transd());
      fprintf(ofile, "}\n");
   }

   fprintf(ofile, "\n\n\n");
}

void disp_poly::vrml_verts_l(FILE *ofile)
{
   fprintf(ofile,
"# Vertex elements\n"
"\n"
"Shape {\n"
"   appearance Appearance {\n"
"      material Material {\n"
"         emissiveColor %s\n"
"      }\n"
"   }\n"
"   geometry PointSet {\n"
"      coord USE CRDS\n"
"   }\n"
"}\n"
"\n", vrml_col(get_def_v_col()).c_str());
   
   fprintf(ofile, "\n\n\n");
}



void disp_poly::vrml_edges(FILE *ofile)
{
   
   fprintf(ofile, "# Edge elements\n");
   
   const vector<vec3d> &vs = disp_geom.verts();
   const vector<vector<int> > &es = disp_geom.edges();
   for(unsigned int i=0; i<es.size(); i++) {
      if(disp_geom.get_e_col((int)i).is_inv())
         continue;
      vec3d mid = (vs[es[i][0]] +  vs[es[i][1]])/2.0 ;
      vec3d dir = vs[es[i][0]] -  vs[es[i][1]];
      double ht =  dir.mag() ;
      dir /= ht; // to unit
      double ang = -acos(safe_for_trig(dir[1])); // angle between dir and y-axis
      vec3d axis = vcross(dir, vec3d(0,1,0)).unit();  // axis
      fprintf(ofile, "E_%s { C %s R %s %g\n\t  ht %g ",
            dots2underscores(sc_geom->get_name()).c_str(),
            vrml_vec(mid, 8 /*sig_digits*/).c_str(),
            vrml_vec(axis,8 /*sig_digits*/).c_str(), ang, ht);
      
      col_val col = disp_geom.get_e_col((int)i);
      if(col.is_idx())
         col = get_e_clrng().get_col(col.get_idx());
      if(col.is_val())
         fprintf(ofile, "clr %s trn %.4f", vrml_col(col).c_str(),
               col.get_transd());
      fprintf(ofile, "}\n");
   }

   fprintf(ofile, "\n\n\n");
}


void disp_poly::vrml_edges_l(FILE *ofile)
{
   fprintf(ofile, "# Edge elements\n");
   fprintf(ofile,
"Shape {\n"
"   appearance Appearance {\n"
"      material Material {\n"
"         emissiveColor %s\n"
"      }\n"
"   }\n", vrml_col(get_def_e_col()).c_str());
   fprintf(ofile,
"   geometry IndexedLineSet {\n"
"      colorPerVertex FALSE\n"
"      coord USE CRDS\n"
"      coordIndex [\n");
 
   const vector<vector<int> > &es = disp_geom.edges();
   for(unsigned int i=0; i<es.size(); i++) {
      fprintf(ofile, "%d %d -1  ", es[i][0], es[i][1]);
      if(!(i%6))
         fprintf(ofile, "\n");
   }
   fprintf(ofile,
"      ]\n"
"   }\n"
"}\n");

   fprintf(ofile, "\n\n\n");
}


void disp_poly::vrml_faces(FILE *ofile)
{
   fprintf(ofile, "# Face elements\n");
   
   map<int, vector<int> > f_alpha;
   const vector<vector<int> > &fs = disp_geom.faces();
   for(unsigned int i=0; i<fs.size(); i++) {
      if(fs[i].size()<3)  // skip degenerate polygons
         continue;
      int alpha = -1;

      col_val col = disp_geom.get_f_col((int)i);
      if(col.is_idx())
         col = get_f_clrng().get_col(col.get_idx());
      if(col.is_inv())
         continue;
      if(col.is_val())
         alpha = col[3];
      f_alpha[alpha].push_back(i);
   }
   
   map<int, vector<int> >::iterator mi;
   for(mi=f_alpha.begin(); mi!=f_alpha.end(); mi++) {
      fprintf(ofile, "F%s_%s { vc USE CRDS ci [ ", mi->first<0?"0":"",
            dots2underscores(sc_geom->get_name()).c_str() );
      int f_cnt=0;
      for(unsigned int i=0; i<mi->second.size(); i++) {
         int idx = mi->second[i];
         if(fs[idx].size()<3)  // skip digons
            continue;
         for(unsigned int j=0; j<fs[idx].size(); j++)
            fprintf(ofile, "%d ", fs[idx][j]);
         fprintf(ofile, "-1  ");
         if(!((++f_cnt)%4))
            fprintf(ofile, "\n");
      }
      fprintf(ofile, "   ]\n\t");
      if(mi->first >= 0) {
         fprintf(ofile,"clrs [ ");
         f_cnt = 0;
         for(unsigned int i=0; i<mi->second.size(); i++) {
            int idx = mi->second[i];
            if(fs[idx].size()<3)  // skip degenerate polygons
               continue;
            col_val col = disp_geom.get_f_col((int)idx);
            if(col.is_idx())
               col = get_f_clrng().get_col(col.get_idx());
            fprintf(ofile, "%s, ", vrml_col(col).c_str());
            if(!((++f_cnt)%3))
               fprintf(ofile, "\n\t");
         }
         fprintf(ofile, " %s", "0 0 0 ]"); // dummy color for the last ','
         fprintf(ofile, " trn %.4f", 1 - col_val::i2f(mi->first));
      }
   fprintf(ofile, " }\n");
   }

   fprintf(ofile, "\n\n\n");
}

void disp_poly::vrml_trans_end(FILE *ofile)
{

   fprintf(ofile, "# close scene transformations\n"
"   ]\n"
"}\n");

}


void disp_poly::vrml_geom(FILE *ofile, const scene &scen, int sig_digits)
{
   if(disp_geom.verts().size()==0)   // Don't write out empty geometries
      return;
   
   vrml_protos(ofile);
   vrml_trans_begin(ofile, scen);

   if(f().get_show() || use_lines)
      vrml_coords(ofile, sig_digits);
   if(v().get_show()) {
      if(use_lines)
         vrml_verts_l(ofile);
      else
         vrml_verts(ofile, sig_digits);
   }
   if(e().get_show()) {
      if(use_lines)
         vrml_edges_l(ofile);
      else
         vrml_edges(ofile);
   }
   if(f().get_show())
      vrml_faces(ofile);
   
   vrml_trans_end(ofile);
}      

// --------------------------------------------------------------
// disp_poly - pov 


void disp_poly::pov_default_vals(FILE *ofile)
{
   fprintf(ofile,
         "#declare PtsCentre = %s;\n"
         "#declare PtsWidth = %g;\n"
         "#declare PtsBallRad = %g;\n"
         "\n",
         pov_vec(sc_geom->get_centre()).c_str(),
         sc_geom->get_width(),
         sc_geom->get_v_ball_rad() );

   fprintf(ofile,
         "// Display flags\n"
         "#declare show = 1; // Show object, may be 1 - show, 0 hide\n"
         "\n"
         "   // Show elements of a type values may be 1 - show, 0 - hide\n"
         "   #declare verts_show = %d;\n"
         "   #declare edges_show = %d;\n"
         "   #declare faces_show = %d;\n"
         "\n",
         v().get_show(),
         e().get_show(),
         f().get_show() );
  
   fprintf(ofile, 
         "// Display values\n"
         "   // Size (or radius) of elements\n"
         "   #declare vert_sz = %g; // %g\n"
         "   #declare edge_sz = %g; // %g\n"
         "   #declare face_sz = %g; // %g\n"
         "\n",
         get_vert_rad(), get_vert_rad(),
         get_edge_rad(), get_edge_rad(),
         f().get_size(), f().get_size() );

   fprintf(ofile, 
         "   // Colour of elements (used to set up default textures\n"
         "   #declare vert_col = %s; // %s\n"
         "   #declare edge_col = %s; // %s\n"
         "   #declare face_col = %s; // %s\n"
         "\n",
         pov_col(get_def_v_col()).c_str(), pov_col(get_def_v_col()).c_str(),
         pov_col(get_def_e_col()).c_str(), pov_col(get_def_e_col()).c_str(),
         pov_col(get_def_f_col()).c_str(), pov_col(get_def_f_col()).c_str() );

   fprintf(ofile, 
         "   // Texture of elements\n"
         "   #declare vert_tex=texture{ pigment{ rgbt vert_col}}\n"
         "   #declare edge_tex=texture{ pigment{ rgbt edge_col}}\n"
         "   #declare face_tex=texture{ pigment{ rgbt face_col}}\n"
         "\n"
         "#declare col_map = array[1]; // Default colourmap\n"
         "#declare tex_map = array[1]; // Default texmap\n"
         "\n" );
}

void disp_poly::pov_disp_macros(FILE *ofile)
{
      fprintf(ofile, 
            "#macro disp_vertex(vertex, col)\n"
            "   default_disp_vertex(vertex, col)\n"
            "#end\n"
            "\n"
            "#macro disp_edge(edge, col)\n"
            "   default_disp_edge(edge, col)\n"
            "#end\n"
            "\n"
            "#macro disp_face(face_no, idx, col)\n"
            "   default_disp_face(face_no, idx, col)\n"
            "#end\n"
            "\n"
            "#macro disp_extra()\n"
            "   default_disp_extra()\n"
            "#end\n");
}


void disp_poly::pov_include_files(FILE *ofile)
{
   fprintf(ofile, "#if(file_exists(\"default_off_i.inc\")) #include \"default_off_i.inc\" #end\n");
   string name(dots2underscores(sc_geom->get_name()));
   if(name != "default_off")
   fprintf(ofile, "#if(file_exists(\"%s_i.inc\")) #include \"%s_i.inc\" #end\n",
                  name.c_str(), name.c_str());
   for(unsigned int i=0; i<includes.size(); i++)
      fprintf(ofile, "#include \"%s\"\n", includes[i].c_str());
   
   fprintf(ofile, "\n");
}


void disp_poly::pov_vert_arrays(FILE *ofile, int sig_digits)
{
   const vector<vec3d> &vs = disp_geom.verts();
   fprintf(ofile, 
         "// Array of vertex coordinates\n"
         "#declare num_verts = %lu;\n",
         (unsigned long)vs.size());

   if(!vs.size())
      return;
   
   fprintf(ofile, "#declare verts = array [num_verts] {\n");
   for(unsigned int i=0; i<vs.size(); i++)
      fprintf(ofile, "   %s%s", pov_vec(vs[i], sig_digits).c_str(),
            (i<vs.size()-1) ? ",\n" : "");
   fprintf(ofile, "\n}\n\n");

   fprintf(ofile, "// Array of vertex colours\n"
                  "#declare v_cols = array [num_verts]\n");
   for(unsigned int i=0; i<vs.size(); i++) {
      col_val col = disp_geom.get_v_col((int)i);
      if(col.is_idx())
         col = get_v_clrng().get_col(col.get_idx());
      if(col.is_set())
         fprintf(ofile, "#declare v_cols[%d]=%s;\n", i,pov_col(col).c_str());
   }

   fprintf(ofile, "\n\n\n");
}

void disp_poly::pov_edge_arrays(FILE *ofile)
{
   const vector<vector<int> > &es = disp_geom.edges();
   fprintf(ofile, "// Array of edge indexes\n"
                  "#declare num_edges = %lu;\n", (unsigned long)es.size());
   if(!es.size())
      return;
   
   fprintf(ofile, "#declare edges = array [num_edges][2] {\n");
   for(unsigned int i=0; i<es.size(); i++)
      fprintf(ofile, "   {%d, %d}%s", es[i][0], es[i][1],
            (i<es.size()-1) ? ",\n" : "");
   fprintf(ofile, "\n}\n\n");

   fprintf(ofile, "// Array of edge colours\n"
                  "#declare e_cols = array [num_edges]\n");
   for(unsigned int i=0; i<es.size(); i++) {
      col_val col = disp_geom.get_e_col((int)i);
      if(col.is_idx())
         col = get_e_clrng().get_col(col.get_idx());
      if(col.is_set())
         fprintf(ofile, "#declare e_cols[%d]=%s;\n", i,pov_col(col).c_str());
   }

   fprintf(ofile, "\n\n\n");
}

void disp_poly::pov_face_arrays(FILE *ofile)
{
   const vector<vector<int> > &fs = disp_geom.faces();
   
   int num_face_items = 0;
   for(unsigned int i=0; i<fs.size(); i++)
      num_face_items += fs[i].size()+1;   // add 1 for the numper of points
   
   fprintf(ofile, "// Array of face vertex counts and indexes\n"
                  "#declare num_faces = %lu;\n"
                  "#declare num_face_items = %d;\n",
                  (unsigned long)fs.size(),
                  num_face_items);
   if(!fs.size())
      return;
   fprintf(ofile, "#declare faces = array [num_face_items] {");
   for(unsigned int i=0; i<fs.size(); i++) {
      fprintf(ofile, "\n   %lu, ", (unsigned long)fs[i].size());
      for(unsigned int j=0; j<fs[i].size(); j++)
         fprintf(ofile, "%d%s", fs[i][j],
               ((i==fs.size()-1) && (j==fs[i].size()-1)) ? "" : ", ");
   }
   
   fprintf(ofile, "\n}\n\n");

   fprintf(ofile, "// Array of face colours\n"
                  "#declare f_cols = array [num_faces]\n");
   for(unsigned int i=0; i<fs.size(); i++) {
      col_val col = disp_geom.get_f_col((int)i);
      if(col.is_idx())
         col = get_f_clrng().get_col(col.get_idx());
      if(col.is_set())
         fprintf(ofile, "#declare f_cols[%d]=%s;\n", i,pov_col(col).c_str());
   }

   fprintf(ofile, "\n\n\n");
}


void disp_poly::pov_elements(FILE *ofile, int sig_digits)
{
   pov_vert_arrays(ofile, sig_digits);
   pov_edge_arrays(ofile);
   pov_face_arrays(ofile);
}      


void disp_poly::pov_col_maps(FILE *ofile)
{
   fprintf(ofile,
"// Colour Maps - redefine these, normally in an include file\n"
"   #declare col_map = array[1];\n"
"   #declare tex_map = array[1];\n"
"   #declare vert_col_map = col_map;\n"
"   #declare vert_tex_map = tex_map;\n"
"   #declare edge_col_map = col_map;\n"
"   #declare edge_tex_map = tex_map;\n"
"   #declare face_col_map = col_map;\n"
"   #declare face_tex_map = tex_map;\n\n");
}

void disp_poly::pov_object(FILE *ofile)
{
   fprintf(ofile,
"#if (show)\n"
//"union {\n"
"#declare NoColour = <-1, -1, -1, 0>; // Indicates no colour has been set" 
"// Display vertex elements\n"
"#if(verts_show)\n"
"   #declare i=0;\n"
"   #while (i<num_verts)\n"
"      #declare col = NoColour;\n"
"      #ifdef (v_cols[i]) #declare col=v_cols[i]+<0,0,0,0>; #end\n"
"         #if (col.x!=0 | col.y!=0 | col.z!=0 | col.t!=1)\n"
"            disp_vertex(i, col)\n"
"         #end\n"
"      #declare i=i+1;\n"
"      #end\n"
"   #end // (verts_show)\n"
"\n"
"// Display edge elements\n"
"#if (edges_show)\n"
"   #declare i=0;\n"
"   #while (i<num_edges)\n"
"      #declare col = NoColour;\n"
"      #ifdef (e_cols[i]) #declare col=e_cols[i]+<0,0,0,0>; #end\n"
"         #if (col.x!=0 | col.y!=0 | col.z!=0 | col.t!=1)\n"
"            disp_edge(i, col)\n"
"         #end\n"
"      #declare i=i+1;\n"
"      #end\n"
"   #end // (edges_show)\n"
"\n"
"// Display face elements\n"
"#if (faces_show)\n"
"   #declare face_no=0;"
"   #declare idx=0;\n"
"   #while (face_no<num_faces)\n"
"      #declare col = NoColour;\n"
"      #ifdef (f_cols[face_no]) #declare col=f_cols[face_no]+<0,0,0,0>; #end\n"
"         #if (col.x!=0 | col.y!=0 | col.z!=0 | col.t!=1)\n"
"            disp_face(face_no, idx, col)\n"
"         #end\n"
"      #declare idx = idx + faces[idx] + 1;\n"
"      #declare face_no=face_no+1;\n"
"      #end\n"
"   #end // (faces_show)\n"
"\n"
"// Extra object\n"
"disp_extra()\n"
"\n"
//"}\n\n"
"#end // (show)\n");
}


void disp_poly::pov_geom(FILE *ofile, const scene &, int sig_digits)
{
   if(disp_geom.verts().size()==0)   // Don't write out empty geometries
      return;
   pov_default_vals(ofile);
   pov_disp_macros(ofile);
   pov_elements(ofile, sig_digits);
   pov_col_maps(ofile);
   pov_include_files(ofile);
   pov_object(ofile);
}

#ifdef HAVE_CONFIG_H
   #include "../config.h"
#endif


#if FOUND_GLU == 1

#include <GL/gl.h>
#include <GL/glu.h>


extern unsigned char stippleMask[17][128];

void gl_set_material(col_val col = vec3d(.5,.5,.5), bool trans=true,
      bool both_sides=true);
 
void gl_set_material(col_val col, bool trans, bool both_sides)
{
   vec4d cv = col.get_vec4d();
   int sides = both_sides ? GL_FRONT_AND_BACK : GL_FRONT;
   GLfloat f_specular[] = { cv[0]/2, cv[1]/2, cv[2]/2, 1.0 };
   GLfloat f_diffuse[] = { cv[0], cv[1], cv[2], 1.0 };
   GLfloat f_shininess[] = { 100.0 };
   glMaterialfv(sides, GL_SPECULAR, f_specular);
   glMaterialfv(sides, GL_DIFFUSE, f_diffuse);
   glMaterialfv(sides, GL_SHININESS, f_shininess);
   if(trans)
      glPolygonStipple(stippleMask[int(cv[3]*16+0.5)]);
}

 

void disp_poly::gl_verts(const scene &scen)
{
   double v_rad = get_vert_rad();
   double extra = 10*v_rad/scen.get_width();
   int long_div = int(11*(1+extra));
   int lat_div = int(7*(1+extra));
   GLUquadric *quad = gluNewQuadric();
   GLuint sph = glGenLists(1);
   glNewList(sph, GL_COMPILE);
   gluSphere(quad, v_rad, long_div, lat_div);
   glEndList();
   
   const vector<vec3d> &verts = disp_geom.verts();
   for(unsigned int i=0; i<verts.size(); i++) { 
      col_val col = disp_geom.get_v_col((int)i);
      if(col.is_idx())
         col = get_v_clrng().get_col(col.get_idx());
      if(!col.is_val())
         col = get_def_v_col(); // use default
      if(col.is_inv())
         continue;

      gl_set_material(col, get_elem_trans(), false);
      glPushMatrix();
      glTranslated(verts[i][0], verts[i][1], verts[i][2]);
      glCallList(sph);
      glPopMatrix();
   }
   glDeleteLists(sph, 1);
   gluDeleteQuadric(quad);
}

static void edge_cyl_trans(vec3d p1, vec3d p2)
{
   vec3d p1to2 = p2 - p1;
   glTranslated(p1[0], p1[1], p1[2]);
   mat3d rot = mat3d::rot(vec3d(0,0,1), p1to2);
   glMultMatrixd(rot.transpose().get_m());
   glScaled(1, 1, p1to2.mag());
} 


void disp_poly::gl_edges(const scene &scen)
{
   double e_rad = get_edge_rad();
   double extra = 10*e_rad/scen.get_width();
   int long_div = int(11*(1+extra));
   int lat_div = 1;
   
   GLUquadric *quad = gluNewQuadric();
   GLuint cyl = glGenLists(1);
   glNewList(cyl, GL_COMPILE);
   gluCylinder(quad, e_rad, e_rad, 1, long_div, lat_div);
   glEndList();
   
   const vector<vec3d> &verts = disp_geom.verts();
   const vector<vector<int> > &edges = disp_geom.edges();
   for(unsigned int i=0; i<edges.size(); i++) { 
      col_val col = disp_geom.get_e_col((int)i);
      if(col.is_idx())
         col = get_e_clrng().get_col(col.get_idx());
      if(!col.is_val())
         col = get_def_e_col(); // use default
      if(col.is_inv())
         continue;

      gl_set_material(col, get_elem_trans(), false);
      glPushMatrix();
      edge_cyl_trans(verts[edges[i][0]], verts[edges[i][1]]);
      glCallList(cyl);
      glPopMatrix();
   }
   glDeleteLists(cyl, 1);
   gluDeleteQuadric(quad);
}


void disp_poly::gl_faces(const scene &)
{
   const vector<vec3d> &verts = disp_geom.verts();
   const vector<vector<int> > &faces = disp_geom.faces();
   for(unsigned int i=0; i<faces.size(); i++) { 
      if(faces[i].size()<3)
         continue;
      col_val col = disp_geom.get_f_col((int)i);
      if(col.is_idx())
         col = get_f_clrng().get_col(col.get_idx());
      if(!col.is_val())
         col = get_def_f_col(); // use default
      if(col.is_inv())
         continue;

      gl_set_material(col, get_elem_trans());
      glBegin(GL_POLYGON);
      vec3d norm = face_norm(verts,faces[i]);
      glNormal3dv(norm.get_v());
      for(unsigned int j=0; j<faces[i].size(); j++)
         glVertex3dv(verts[faces[i][j]].get_v());
      glEnd();
   }
}

#else

void disp_poly::gl_verts(const scene &) {}
void disp_poly::gl_edges(const scene &) {}
void disp_poly::gl_faces(const scene &) {}

#endif // FOUND_GLU == 1


void disp_poly::gl_geom(const scene &scen)
{
   if(v().get_show())
      gl_verts(scen);
   if(f().get_show())
      gl_faces(scen);
   if(e().get_show())
      gl_edges(scen);
}











// --------------------------------------------------------------
// disp_num_labels 

disp_num_labels::disp_num_labels()
{
   v().set_col(col_val(0.5,0.0,0.0));
   e().set_col(col_val(0.0,0.5,0.0));
   f().set_col(col_val(0.0,0.0,0.5));
}


void disp_num_labels::pov_geom(FILE *file, const scene &scen, int sig_dgts)
{
   if(file || &scen || sig_dgts)
      fprintf(stderr, "\n");
}


void disp_num_labels::vrml_protos(FILE *ofile, const scene &scen)
{
   col_val bg = scen.get_bg_col();
   bool bg_dark = (bg[0]+bg[1]+bg[2])<1.5;
   vec3d txt_col = vec3d(bg_dark, bg_dark, bg_dark);
   double txt_sz = sc_geom->get_width()/30;

   fprintf(ofile, 
"\n"
"PROTO LAB [\n"
"   field SFColor clr %s"
"   field MFString txt \"\"\n"
"   field SFVec3f pos 0 0 0\n"
"]\n"
"{\n"
"   Transform {\n"
"      translation IS pos\n"
"      children [\n"
"         Billboard {\n"
"            axisOfRotation 0 0 0\n"
"            children [\n"
"               Shape {\n"
"                  geometry Text { string IS txt fontStyle FontStyle { size %g justify \"MIDDLE\"} }\n"
"                  appearance Appearance {\n"
"                     material Material {\n"
"                        diffuseColor IS clr\n"
"                     }\n"
"                  }\n"
"               }\n"
"            ]\n"
"         }\n"
"      ]\n"
"   }\n"
"}\n",  vrml_col(txt_col).c_str(), txt_sz);

}


void disp_num_labels::vrml_verts(FILE *ofile)
{
   const col_geom_v &geom = sc_geom->get_geom();
   fprintf(ofile, "# Vertex number labels\n");
   int v_sz = geom.verts().size();
   for(int i=0; i<v_sz; i++) {
      if(geom.get_v_col((int)i).is_inv())
         continue;
      fprintf(ofile, "LAB { txt \"%d\" pos %s }\n",
            i, vrml_vec(sc_geom->get_v_label_pos(i), 4).c_str());
   }
   fprintf(ofile, "\n\n\n");
}

void disp_num_labels::vrml_edges(FILE *ofile)
{
   fprintf(ofile, "\n\n\n");
}

void disp_num_labels::vrml_faces(FILE *ofile)
{
   fprintf(ofile, "\n\n\n");
}

void disp_num_labels::vrml_geom(FILE *ofile, const scene &scen, int)
{
   vrml_protos(ofile, scen);
   if(v().get_show())
      vrml_verts(ofile);
   if(e().get_show())
      vrml_edges(ofile);
   if(f().get_show())
      vrml_faces(ofile);
}      


#if FOUND_GLU == 1

#include <GL/gl.h>
#include <GL/glu.h>


static void draw_text(char *str, double font_sz, vec3d pos, int halign=1,
      int valign=1, bool bill=true,
      const anti_StrokeFont*font=ANTI_STROKE_ROMAN)
{  
   const float drop = 33.33;
   const float height = 152.38;
   float scale = font_sz/height;
   float width = antiStrokeLength(font, (unsigned char *)str);
  
   float off_x=0, off_y=0;
   switch(halign) {
      case 0: // left
         off_x=0;
         break;
      case 1: // centre
         off_x -= 0.5 * width;
         break;
      case 2: // right
         off_x -= width;
         break;
   }
   
   switch(valign) {
      case 0: // bottom
         off_y -= -drop;
         break;
      case 1: // centre
         off_y -= height/2 - drop;
         break;
      case 2: // right
         off_y -= height - drop;
         break;
      case 3: // right
         off_y = 0;
         break;
   }
   
   glLineWidth(1.5);
   glPushMatrix();
   glTranslated(pos[0], pos[1], pos[2]);
   
   if(bill) {
      float m[16];
      glGetFloatv(GL_MODELVIEW_MATRIX , m);
      for(int i=0; i<3; i++ ) 
         for(int j=0; j<3; j++ ) {
            if ( i==j )
               m[i*4+j] = 1;
            else
               m[i*4+j] = 0;
         }
      glLoadMatrixf(m);
   }
  
   glScalef(scale, scale, scale);
   glTranslatef(off_x, off_y, 0);
   glNormal3f(0,0,1);
   for (char *c=str; *c ; c++)
      antiStrokeCharacter(font, *c);
   glPopMatrix();
}



static void gl_write_label(char *label, vec3d pos, const camera &cam)
{
   draw_text(label, cam.get_text_sz(pos), pos);
}


void disp_num_labels::gl_verts(const scene &scen)
{
   col_geom_v &geom = sc_geom->get_geom();
   gl_set_material(get_label_col(v().get_col()), get_elem_trans());
   char label[64];
   const vector<vec3d> &verts = geom.verts();
   for(unsigned int i=0; i<verts.size(); i++) {
      if(geom.get_v_col((int)i).is_inv())
         continue;
      sprintf(label, "%u", i);
      gl_write_label(label, sc_geom->get_v_label_pos(i), scen.cur_camera());
   }
}

void disp_num_labels::gl_edges(const scene &scen)
{
   col_geom_v &geom = sc_geom->get_geom();
   gl_set_material(get_label_col(e().get_col()), get_elem_trans());
   char label[64];
   const vector<vector<int> > &edges = geom.edges();
   for(unsigned int i=0; i<edges.size(); i++) { 
      if(geom.get_e_col((int)i).is_inv())
         continue;
      sprintf(label, "%u", i);
      gl_write_label(label, sc_geom->get_e_label_pos(i), scen.cur_camera());
   }
}


void disp_num_labels::gl_faces(const scene &scen)
{
   col_geom_v &geom = sc_geom->get_geom();
   gl_set_material(get_label_col(f().get_col()), get_elem_trans());
   char label[64];
   const vector<vector<int> > &faces = geom.faces();
   for(unsigned int i=0; i<faces.size(); i++) { 
      if(geom.get_f_col((int)i).is_inv())
         continue;
      sprintf(label, "%u", i);
      gl_write_label(label, sc_geom->get_f_label_pos(i), scen.cur_camera());
   }
}

#else

void disp_num_labels::gl_verts(const scene &) {}
void disp_num_labels::gl_edges(const scene &) {}
void disp_num_labels::gl_faces(const scene &) {}

#endif // FOUND_GLU == 1


void disp_num_labels::gl_geom(const scene &scen)
{
   //camera->set_label_rot();
   if(v().get_show())
      gl_verts(scen);
   if(f().get_show())
      gl_faces(scen);
   if(e().get_show())
      gl_edges(scen);
}


// --------------------------------------------------------------
// disp_sym

disp_sym::disp_sym(): show_axes(false), show_mirrors(false),
   show_rotrefls(false), sym()
{
   v().set_show(false);
   e().set_show(false);
}
      
void initialise_unitialised_sym();

void add_ring(col_geom_v &geom, float i_rad, float o_rad, int steps,
      col_val col, const mat3d &trans)
{
   int curr_idx = geom.verts().size();
   for(int i=0; i<steps; i++) {
      for(int h=0; h<2; h++) {
         double ang = 2*M_PI*i/steps;
         geom.add_vert(trans * vec3d(i_rad*cos(ang), i_rad*sin(ang), 0));
         geom.add_vert(trans * vec3d(o_rad*cos(ang), o_rad*sin(ang), 0));
         vector<int> face(4);
         int offsets[] = {0,1,5,4};
         for(int j=0; j<4; j++)
            face[j] = curr_idx + (4*i+2*h+offsets[j]) % (4*steps);
         geom.add_col_face(face, col);
      }
   }
   
}

void add_mirror_elem(col_geom_v &geom, float rad, const mat3d &trans)
{
   int steps = 30;
   float out_rad = 1.01*rad;
   float height = .005*rad;
   col_val mirror_col = col_val(0.6,0.6,0.6);

   int curr_idx = geom.verts().size();
   for(int i=0; i<steps; i++) {
      for(int h=0; h<2; h++) {
         float ht = (1-2*h)*height;
         double ang = 2*M_PI*i/steps;
         geom.add_vert(trans * vec3d(rad*cos(ang), rad*sin(ang), ht));
         geom.add_vert(trans * vec3d(out_rad*cos(ang), out_rad*sin(ang), 0));
         vector<int> face(4);
         int offsets[] = {0,1,5,4};
         for(int j=0; j<4; j++)
            face[j] = curr_idx + (4*i+2*h+offsets[j]) % (4*steps);
         geom.add_col_face(face, mirror_col);
      }
   }
   
}

void add_rotrefl_elem(col_geom_v &geom, float rad, const mat3d &trans)
{
   int steps = 60;
   float out_rad = 1.01*rad;
   float height = .01*rad;
   col_val mirror_col = col_val(1.0,1.0,1.0);

   int curr_idx = geom.verts().size();
   for(int i=0; i<steps; i++) {
      for(int h=0; h<2; h++) {
         float ht = (1-2*h)*height;
         double ang = 2*M_PI*i/steps;
         geom.add_vert(trans * vec3d(rad*cos(ang), rad*sin(ang), ht));
         geom.add_vert(trans * vec3d(out_rad*cos(ang), out_rad*sin(ang), 0));
         if((i+h)%2) {
            vector<int> face(4);
            int offsets[] = {0,1,5,4};
            for(int j=0; j<4; j++)
               face[j] = curr_idx + (4*i+2*h+offsets[j]) % (4*steps);
            geom.add_col_face(face, mirror_col);
         }
      }
   }
   
}

void axis_cap(col_geom_v &geom, const sch_axis &sym, col_val col)
{
   geom.clear_all();
   col_val alt_col = col_val(0.4,0.4,0.4);
   int fold = sym.get_nfold();
   int typ = sym.get_sym_type();
   if(typ==sch_sym::S)
      fold /= 2;

   float rad2 = 0.6;
   double ang_inc = 2*M_PI/fold;
   double extra_inc = 0;
   if(typ==sch_sym::Cv || typ==sch_sym::Dv || typ==sch_sym::Dh ) {
      rad2 = 1.0;
      extra_inc = ang_inc/4;
   }

     
   vector<int> face;
   geom.add_vert(vec3d(0,0,0));
   for(int i=0; i<fold; i++) {
      double ang = i*ang_inc + extra_inc;
      vector<int> face(3);
      face[0] = 0;
      face[1] = geom.verts().size();
      geom.add_vert(vec3d(cos(ang), sin(ang), 0));
      face[2] = geom.verts().size();
      geom.add_vert(vec3d(rad2*cos(ang+ang_inc/2), rad2*sin(ang+ang_inc/2), 0));
      geom.add_col_face(face, col);
      switch(typ) {
         case sch_sym::C:
         case sch_sym::Cv:
         case sch_sym::Dh:
            break;
         case sch_sym::S:
            face[1] = geom.verts().size();
            geom.add_vert(vec3d(cos(ang+ang_inc/2), sin(ang+ang_inc/2), 0));
            face[2] = geom.verts().size();
            geom.add_vert(vec3d(rad2*cos(ang+ang_inc),rad2*sin(ang+ang_inc),0));
            geom.add_col_face(face, alt_col);
            break;
         case sch_sym::D:
         case sch_sym::Dv:
            face[1] = geom.verts().size();
            geom.add_vert(vec3d(rad2*cos(ang+ang_inc/2),rad2*sin(ang+ang_inc/2),0));
            face[2] = geom.verts().size();
            geom.add_vert(vec3d(cos(ang+ang_inc), sin(ang+ang_inc), 0));
            geom.add_col_face(face, alt_col);
            break;
      }
   }

   if(typ==sch_sym::Ch || typ==sch_sym::Dh)
      add_ring(geom, 0.8, 1.1, 12, alt_col, mat3d::transl(vec3d(0,0,-0.1)));
}


void add_axis_elem(col_geom_v &geom, const sch_axis &sym, float rad,
      const mat3d &trans)
{
   int fold = sym.get_nfold();
   if(fold<2)
      return;

   int typ = sym.get_sym_type();
   if(typ==sch_sym::S)
      fold /= 2;

   col_val cols[6] = { col_val(0.6,0.3,0.0), col_val(), col_val(0.8,0.8,0.2),
      col_val(0.3,0.8,0.3), col_val(0.6,0.0,0.0), col_val(0.0,0.0,0.6) };
   col_val col = (fold<=5) ? cols[fold] : cols[0];

   float radius = 0.007*rad; // axis radius
   float ht = 1.03*rad;      // half the axis height
   int steps = 6;
   int curr_idx = geom.verts().size();
   for(int i=0; i<steps; i++) {
      double ang = 2*M_PI*i/steps;
      geom.add_vert(trans * vec3d(radius*cos(ang), radius*sin(ang), ht));
      geom.add_vert(trans * vec3d(radius*cos(ang), radius*sin(ang), -ht));
      vector<int> face(4);
      for(int j=0; j<4; j++)
         face[j] = curr_idx + (2*i + ((j<2)?j:(5-j))) % (2*steps);
      geom.add_col_face(face, col);
   }

   col_geom_v cap, cap2;
   axis_cap(cap, sym, col);
   cap.transform(mat3d::scale(rad*0.06));
   cap2 = cap;
   cap.transform(trans * mat3d::transl(vec3d(0,0,ht)));
   //cap2.transform(trans * mat3d::transl(vec3d(0,0,-ht)));
   cap2.transform(trans * mat3d::inversion() * mat3d::transl(vec3d(0,0,ht)));
   geom.append(cap);
   geom.append(cap2);
   
}

void disp_sym::disp_changed()
{
   disp_geom.clear_all();
   if(!sc_geom)
      return;

   if(sym.get_sym_type()==sch_sym::unknown &&
         (show_axes || show_mirrors || show_rotrefls))
      sym.init(sc_geom->get_geom());
   
   double rad = 1.05*sc_geom->get_width()/2;
   vec3d cent = sc_geom->get_centre();
   const set<sch_axis> &axes = sym.get_axes(); 
   set<sch_axis>::const_iterator ax;
   for(ax=axes.begin(); ax!=axes.end(); ++ax) {
      mat3d trans = mat3d::transl(cent);
      if(ax->get_nfold()==2 && (sym.get_sym_type()==sch_sym::D ||
            sym.get_sym_type()==sch_sym::Dv || sym.get_sym_type()==sch_sym::Dh))
         trans *= mat3d::alignment(vec3d::z, vec3d::x,
               ax->get_axis(), mat3d::inverse(sym.get_to_std())*vec3d::z-cent);
      else if(ax->get_perp().is_set())
         trans *= mat3d::alignment(vec3d::z, vec3d::x,
               ax->get_axis(), ax->get_perp());
      else
         trans *=  mat3d::rot(vec3d::z, ax->get_axis());
      int sym_type = ax->get_sym_type();
      if(show_rotrefls && (sym_type==sch_sym::S  || sym_type==sch_sym::Dv ) ) {
         add_rotrefl_elem(disp_geom, rad, trans);
      }
      if(show_axes) {
         add_axis_elem(disp_geom, *ax, rad, trans);
      }
   }
   
   if(show_mirrors) {
      const set<vec3d> &mirrors = sym.get_mirrors();
      mat3d turn;
      set<vec3d>::const_iterator mir;
      for(mir=mirrors.begin(); mir!=mirrors.end(); ++mir) {
         mat3d trans = mat3d::transl(cent) * turn * mat3d::rot(vec3d::z, *mir);
         add_mirror_elem(disp_geom, 1.05*sc_geom->get_width()/2, trans);
      }
   }

}



void disp_sym::geom_changed()
{
   disp_changed();
}

void disp_sym::vrml_geom(FILE *ofile, const scene &scen, int /*sig_dgts*/)
{
   disp_poly::vrml_geom(ofile, scen, 4);
}

void disp_sym::pov_geom(FILE *ofile, const scene &scen, int /*sig_dgts*/)
{
   disp_poly::pov_geom(ofile, scen, 4);
}

void disp_sym::gl_geom(const scene &scen)
{
   disp_poly::gl_geom(scen);
}



// --------------------------------------------------------------
// Other functions

void view_opts::set_view_vals(scene &scen)
{
   scen = scen_defs;
   scen.add_camera(cam_defs);
 
   char errmsg[MSG_SZ];
   for(unsigned int i=0; i<ifiles.size(); i++) {
      col_geom_v geom;
      if(!geom.read(ifiles[i], errmsg))
         error(errmsg);
      if(*errmsg)
         warning(errmsg);
   
      if((geom_defs.e().get_col()!=col_val(0,0,0,0)))
         geom.add_missing_impl_edges();
      
      scene_geom sc_geom;
      sc_geom.set_scene(&scen);
      sc_geom.add_disp(geom_defs);
      sc_geom.set_label(lab_defs);
      sc_geom.set_sym(sym_defs);
      //fprintf(stderr, "before geom = %p\n", sc_geom.get_sym
      sc_geom.set_geom(geom);


      if(ifiles[i]!="")
         sc_geom.set_name(basename2(ifiles[i].c_str()));
      else
         sc_geom.set_name("stdin");

      scen.add_geom(sc_geom);
      
      // Use element sizes from first geometry
      scen.get_geoms()[i].get_disps()[0]->
         e().set_size(scen.get_geoms()[0].get_disps()[0]->get_edge_rad());
      scen.get_geoms()[i].get_disps()[0]->
         v().set_size(scen.get_geoms()[0].get_disps()[0]->get_vert_rad());
   }
}
   
const char *view_opts::help_view_text =
"  -v <rad>  radius of vertex spheres, or 'b' to have radius of balls\n"
"            of the maximum size without overlap (default: ball_rad/15)\n"
"  -e <rad>  radius of edge cylinders (default: vertex_rad/1.5)\n"
"  -V <col>  default vertex colour, in form 'R,G,B,A' (3 or 4 values\n"
"            0.0-1.0, or 0-255) or hex 'xFFFFFF' (default: 1.0,0.5,0.0)\n"
"  -E <col>  default edge colour, in form 'R,G,B,A' (3 or 4 values\n"
"            0.0-1.0, or 0-255) or hex 'xFFFFFF', 'x' to hide implicit edges\n"
"            (default: 0.8,0.6,0.8)\n"
"  -F <col>  default face colour, in form 'R,G,B,A' (3 or 4 values\n"
"            0.0-1.0, or 0-255) or hex 'xFFFFFF' (default: 0.8,0.9,0.9)\n"
"  -x <elms> hide elements. The element string can include v, e and f\n"
"            to hide vertices, edges and faces\n"
"  -n <elms> show element index number labels. The element string can\n"
"            include v, e and f to label vertices, edges and faces\n"
"  -s <syms> show symmetry elements. The element string can include\n"
"               x - rotation axes\n"
"               m - mirror planes\n"
"               r - rotation-reflection planes\n"
"               a - all elements (same as xmr)\n"
"  -m <maps> a comma separated list of colour maps used to transform colour\n"
"            indexes, a part consisting of letters from v, e, f, selects \n"
"            the element types to apply the map list to (default 'vef').\n"
;

const char *view_opts::help_scene_text =
"  -D <dist> distance to camera\n"
"  -C <cent> centre of points, in form 'X,Y,Z'\n"
"  -L <look> point to look at, in form 'X,Y,Z'\n"
"            (default, points centre)\n"
"  -R <rot>  rotate about axes through centre of points, in\n"
"            form 'X-ang,Y-ang,Z-ang' (degrees)\n"
"  -B <col>  background colour, in form 'R,G,B,A' (3 or 4 values\n"
"            0.0-1.0, or 0-255) or hex 'xFFFFFF'\n";

const char *view_opts::help_prec_text =
"  -d <dgts> number of significant digits (default 17) or if negative\n"
"            then the number of digits after the decimal point\n"
"  -t <type> display type for faces 0 - native polygons, 1 - triangulate\n"
"            polygons (default)\n";

bool view_opts::read_disp_option(char opt, char *optarg, char *errmsg,
            vector<string> &warnings)
{
   *errmsg = '\0';
   char errmsg2[MSG_SZ];
   double val;
   int num;
   vec3d vec;
   col_val col;

   bool handled = true;
   switch(opt) {
      case 'v':
         if(strcmp(optarg, "b")==0) {
            geom_defs.v().set_size(geom_disp::rad_ball);
         }
         else {
            if(!read_double(optarg, &val, errmsg2))
               strcpy(errmsg, errmsg2);
            else if(val < 0)
               strcpy(errmsg, "vertex sphere radius cannot be negative");
            else
               geom_defs.v().set_size(val);
         }
         break;

         case 'e':
            if(!read_double(optarg, &val, errmsg2))
               strcpy(errmsg, errmsg2);
            else if(val < 0)
               strcpy(errmsg, "edge cylinder radius cannot be negative");
            else
               geom_defs.e().set_size(val);
            break;

         case 'V':
            if(!col.read(optarg, errmsg2))
               strcpy(errmsg, errmsg2);
            else
               geom_defs.v().set_col(col);
            break;
            
         case 'E':
            if(!col.read(optarg, errmsg2))
               strcpy(errmsg, errmsg2);
            else
               geom_defs.e().set_col(col);
            break;
         
         case 'F':
            if(!col.read(optarg, errmsg2))
               strcpy(errmsg, errmsg2);
            else
               geom_defs.f().set_col(col);
            break;

         case 'm':
            if(!read_colorings(geom_defs.get_clrngs(), optarg, errmsg2))
               strcpy(errmsg, errmsg2);
            else if(*errmsg2)
               warnings.push_back(errmsg2);
            break;

         case 'x':
            if(strspn(optarg, "vef") != strlen(optarg))
               snprintf(errmsg2, MSG_SZ, "elements to hide are '%s' must be "
                     "from v, e, or f", optarg);
            else {
               if(strchr(optarg, 'v'))
                  geom_defs.v().set_show(false);
               if(strchr(optarg, 'e'))
                  geom_defs.e().set_show(false);
               if(strchr(optarg, 'f'))
                  geom_defs.f().set_show(false);
            }
            break;

         case 'n':
            if(strspn(optarg, "vef") != strlen(optarg))
               snprintf(errmsg2, MSG_SZ, "elements to label are '%s' must be "
                     "from v, e, f", optarg);
            else {
               if(strchr(optarg, 'v'))
                  lab_defs.v().set_show(true);
               if(strchr(optarg, 'e'))
                  lab_defs.e().set_show(true);
               if(strchr(optarg, 'f'))
                  lab_defs.f().set_show(true);
            }
            break;

         case 's':
            if(strspn(optarg, "axmr") != strlen(optarg))
               snprintf(errmsg2, MSG_SZ, "symmetry elements to show are"
                     "'%s' must be from a (all), x, m, r", optarg);
            else {
               sym_defs.set_show_axes(
                     strchr(optarg,'x') || strchr(optarg,'a'));
               sym_defs.set_show_mirrors(
                     strchr(optarg,'m') || strchr(optarg,'a'));
               sym_defs.set_show_rotrefls(
                     strchr(optarg,'r') || strchr(optarg,'a'));
            }
            break;

         case 't':
            if(!read_int(optarg, &num, errmsg2))
               strcpy(errmsg, errmsg2);
            else if(num < 0 || num >1)
               strcpy(errmsg, "display type for faces must be 0 or 1");
            else
               geom_defs.set_triangulate(num);
            break;


         case 'w':
            if(!read_double(optarg, &val, errmsg))
               strcpy(errmsg, errmsg2);
            else if(val <= 0.0)
               strcpy(errmsg, "width must be a positive number");
            else
               cam_defs.set_width(val);
            break;

         case 'I':
            if(!read_double(optarg, &val, errmsg))
               strcpy(errmsg, errmsg2);
            else if(val <= 0.0)
               strcpy(errmsg, "infinity distance must be a positive number");
            else
               scen_defs.set_inf_dist(val);
            break;

         case 'D':
            if(!read_double(optarg, &val, errmsg))
               strcpy(errmsg, errmsg2);
            else if(val <= 0.0)
               strcpy(errmsg, "distance must be a positive number");
            else
               cam_defs.set_distance(val);
            break;

         case 'C':
            if(!vec.read(optarg, errmsg2))
               strcpy(errmsg, errmsg2);
            else
               cam_defs.set_centre(vec);
            break;
            
         case 'L':
            if(!vec.read(optarg, errmsg2))
               strcpy(errmsg, errmsg2);
            else
               cam_defs.set_lookat(vec);
            break;
            
         case 'R':
            if(!vec.read(optarg, errmsg2))
               strcpy(errmsg, errmsg2);
            else
               cam_defs.set_rotation(mat3d::rot((vec)*deg2rad()));
            break;
         
         case 'P':
            if(!read_double(optarg, &val, errmsg2))
               strcpy(errmsg, errmsg2);
            if(val <= 0)
               strcpy(errmsg, "perspective factor must be a positive number");
            else
               cam_defs.set_persp(val);

            break;
 
         case 'B':
            if(!col.read(optarg, errmsg2))
               strcpy(errmsg, errmsg2);
            else
               scen_defs.set_bg_col(col);
            break;
            
         default:
            handled = false;
   }

   return handled;
}
 
