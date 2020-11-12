/*
   Copyright (c) 2003-2016, Adrian Rossiter

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

/* \file PovWriter.cc
   Export to POV-Ray format
*/

#include "povwriter.h"
#include "displaypoly.h"
#include "mathutils.h"

#include <string>
#include <vector>

using std::string;
using std::vector;

namespace anti {

// ----------------------------------------------------------------
// functions

string pov_vec(double x, double y, double z, int sig_digits)
{
  char buf[256];
  if (sig_digits > 0)
    snprintf(buf, 256, "<%.*g, %.*g, %.*g>", sig_digits, x, sig_digits, y,
             sig_digits, z);
  else
    snprintf(buf, 256, "<%.*f, %.*f, %.*f>", -sig_digits, x, -sig_digits, y,
             -sig_digits, z);
  return buf;
}

string pov_col(const Color &col)
{
  char buf[MSG_SZ];
  *buf = '\0';
  if (col.is_value()) {
    Vec4d cv = col.get_vec4d();
    snprintf(buf, MSG_SZ, "<%g, %g, %g, %g>", cv[0], cv[1], cv[2], 1 - cv[3]);
  }
  if (col.is_index())
    snprintf(buf, MSG_SZ, "<%d, -1, 0, 0>", col.get_index());
  return buf;
}

// ------------------------------------------------------------------
// DisplayPoly_pov

// ---------------------------------------------------------------------
// PovWriter

void PovWriter::write(FILE *ofile, const Scene &scen, int sig_dgts)
{
  if (o_type != 'o') {
    scene_header(ofile, scen);
    cameras(ofile, scen);
    common_macros(ofile);
    include_files(ofile, scen);
    scene_setup(ofile);
    stereo_setup(ofile);
    geometry_objects(ofile, scen, sig_dgts);
    stereo_clipping(ofile);
    camera_lights(ofile);
  }
  else
    geometry_objects(ofile, scen, sig_dgts);
}

void PovWriter::set_o_type(char otype) { o_type = otype; }

void PovWriter::scene_header(FILE *ofile, const Scene &scen)
{
  fprintf(ofile, "#version 3.6;\n\n");
  fprintf(ofile, "#include \"colors.inc\"\n\n");
  fprintf(ofile,
          "// Scene Width - maximum distance between any 2 points\n"
          "#declare SceneWidth = %g;\n\n",
          scen.get_bound_sph().get_width());
  fprintf(ofile,
          "// Scene Centre\n"
          "#declare SceneCentre = %s;\n\n",
          pov_vec(scen.get_bound_sph().get_centre()).c_str());
  fprintf(ofile,
          "// The values below may be changed (original values in brackets):\n"
          "\n"
          "// Stereo type may be 0 (mono) or 1 (single image stereo)\n"
          "// or 2 (double image file stereo). If 2 then you must use a\n"
          "// povray option like +KFF2 which will produce the left and right\n"
          "// stereo views in separate image files.\n");
  fprintf(ofile, "#declare StereoType = %d;\n\n", get_stereo_type());
  fprintf(ofile,
          "// Shadow, 0 - shadows mono no shadows stero,"
          "1 - no shadows, 2 - shadows\n"
          "#declare Shadow = %d;\n\n",
          get_shadow());
  fprintf(ofile,
          "// Background colour (%s)\n"
          "#declare BgColour = %s; \n\n",
          pov_col(scen.get_bg_col()).c_str(),
          pov_col(scen.get_bg_col()).c_str());
  fprintf(ofile, "// Max Trace Level (10), increase if black areas appear when "
                 "using transparency\n"
                 "#declare MaxTraceLevel = 10; \n\n");
  fprintf(ofile,
          "// Aspect Ratio (1.3333) set on command line with "
          "Declare=AspectRatio=1.3333\n"
          "#ifndef(AspectRatio) #declare AspectRatio = 1.33333; #end\n\n");

  fprintf(ofile,
          "// Vertex numbering\n"
          "#declare TextSize = %g; // if 0 use value calc'd from geom\n"
          "#declare TextColour = %s; //\n"
          "#declare FontFile = \"cyrvetic.ttf\"\n\n",
          scen.get_width() / 20, pov_col(text_colour).c_str());
}

void PovWriter::common_macros(FILE *ofile)
{
  fprintf(ofile,
          "\n"
          "\n"
          "#macro v_equal(v1,v2) ((v1.x=v2.x)&(v1.y=v2.y)&(v1.z=v2.z)) #end\n"
          "\n"
          "#macro col_to_tex(col, elem_tex_map, elem_col_map, def_tex)\n"
          "   #local typ=0;\n"
          "   #if(col.y>=0) #local typ=1; #end\n"
          "   #if(!typ & col.y=-1)\n"
          "      #if(col.x<dimension_size(elem_tex_map,1))\n"
          "         #ifdef(elem_tex_map[col.x]) #local typ=2; #end\n"
          "      #end\n"
          "      #if(!typ & col.x<dimension_size(elem_col_map,1))\n"
          "         #ifdef(elem_col_map[col.x]) #local typ=3; #end\n"
          "      #end\n"
          "   #end\n"
          "   #switch(typ)\n"
          "      #case(0) texture{ def_tex } #break;\n"
          "      #case(1) texture{ pigment{ rgbt col}} #break;\n"
          "      #case(2) texture{ elem_tex_map[-col.x]} #break;\n"
          "      #case(3) texture{ pigment{ color elem_col_map[-col.x]}} "
          "#break;\n"
          "   #end\n"
          "#end\n"
          "\n"
          "\n"
          "#macro default_disp_vertex(vertex, col)\n"
          "   sphere{ verts[vertex] vert_sz"
          " col_to_tex(col, vert_tex_map, vert_col_map, vert_tex) }\n"
          "#end\n"
          "\n"
          "\n"
          "#macro default_disp_edge(edge, col)\n"
          "   #if(!v_equal(verts[edges[edge][0]], verts[edges[edge][1]]) )\n"
          "      cylinder{verts[edges[edge][0]] verts[edges[edge][1]] edge_sz\n"
          "         col_to_tex(col, edge_tex_map, edge_col_map, edge_tex)\n"
          "      }\n"
          "   #end\n"
          "#end\n"
          "\n"
          "\n"
          "#macro default_disp_face(face_no, idx, col)\n"
          "   disp_face_triangles(face_no, idx, col)\n"
          "#end\n"
          "\n"
          "\n"
          "#macro default_disp_extra()\n"
          "#end\n"
          "\n"
          "\n"
          "#macro disp_face_polygon(face_no, idx, col)\n"
          "   #if (faces[idx]>2)\n"
          "      polygon{ faces[idx]+1 \n"
          "      #local i=0;\n"
          "      #while (i< faces[idx])\n"
          "         verts[faces[idx+i+1]]\n"
          "         #local i = i+1;\n"
          "         #end\n"
          "      verts[faces[idx+1]]\n"
          "      col_to_tex(col, face_tex_map, face_col_map, face_tex)\n"
          "      }\n"
          "   #end\n"
          "#end\n"
          "\n"
          "\n"
          "#macro disp_face_flattened_polygon(face_no, idx, col)\n"
          "   #local p1 = verts[faces[idx+1]];\n"
          "   #local i=1;\n"
          "   #while (i<faces[idx]-1)\n"
          "      #local norm = vcross(verts[faces[idx+i+1]]-p1,"
          "verts[faces[idx+i+2]]-p1);\n"
          "      #local sin_a = vlength(norm)/(vlength(verts[faces[idx+2]]-p1)*"
          "vlength(verts[faces[idx+3]]-p1));\n"
          "      #if(abs(sin_a)>1e-6)\n"
          "         #local norm=vnormalize(norm);\n"
          "         #local i=faces[idx];\n"
          "      #else\n"
          "         #local sin_a=0;\n"
          "      #end\n"
          "      #local i=i+1;\n"
          "   #end\n"
          "   #if(sin_a!=0)\n"
          "      polygon{ faces[idx]+1 p1\n"
          "      #local i=1;\n"
          "      #while (i<faces[idx])\n"
          "         #local vec = verts[faces[idx+i+1]]-p1;\n"
          "         p1 + vec - vdot(vec, norm)*norm\n"
          "         #local i = i+1;\n"
          "         #end\n"
          "      p1\n"
          "      col_to_tex(col, face_tex_map, face_col_map, face_tex)\n"
          "      }\n"
          "   #end\n"
          "#end\n"
          "\n"
          "#macro disp_face_triangles(face_no, idx, col)\n"
          "   #local centroid=0;\n"
          "   #local i=0;\n"
          "   #while (i< faces[idx])\n"
          "      #local centroid=centroid+verts[faces[idx+i+1]];\n"
          "      #local i = i+1;\n"
          "      #end\n"
          "   #local centroid=centroid/faces[idx];\n"
          "   #local i=0;\n"
          "   union {\n"
          "   #while (i< faces[idx])\n"
          "      triangle { centroid verts[faces[idx+i+1]] "
          "verts[faces[idx+mod(i+1, faces[idx])+1]] }\n"

          "      #local i = i+1;\n"
          "      #end\n"
          "   col_to_tex(col, face_tex_map, face_col_map, face_tex)\n"
          "   }\n"
          "#end\n"
          "\n"
          "\n"
          "#macro disp_elem_label(pos, txt, col)\n"
          "   #if(TextSize=0)\n"
          "      #declare TextSize = SceneWidth/40;\n"
          "      #end\n"
          "   text {\n"
          "      ttf FontFile txt 1 0\n"
          "      pigment { rgbt col }\n"
          "      scale <TextSize, TextSize, TextSize/20>\n"
          "      translate -<TextSize/2, TextSize/2, TextSize/40>\n"
          "      #declare vec = vnormalize(pos - PtsCentre);\n"
          "      #declare loc =  pos + (1.2*vert_sz+TextSize)*vec;\n"
          "      translate SceneCentre + vrotate(loc - SceneCentre, Rotation)\n"
          "      transform {\n"
          "         translate -SceneCentre\n"
          "         rotate Rotation\n"
          "         translate SceneCentre\n"
          "         inverse\n"
          "      }\n"
          "   }\n"
          "\n"
          "#end\n"
          "\n"
          "\n\n");
}

void PovWriter::cameras(FILE *ofile, const Scene &scen)
{
  int cur_cam_idx = 0;
  for (unsigned int i = 0; i < scen.get_cameras().size(); i++) {
    if (&scen.get_cameras()[i] == &scen.cur_camera())
      cur_cam_idx = i;
  }
  fprintf(ofile, "#ifndef(CameraNumber) #declare CameraNumber = %d; #end\n\n",
          cur_cam_idx);
  fprintf(ofile, "#switch(CameraNumber)\n");

  for (unsigned int i = 0; i < scen.get_cameras().size(); i++) {
    Camera cam = scen.get_cameras()[i];
    fprintf(ofile, "#case (%u) // %s\n", i, scen.get_camera_name(i).c_str());

    string dist_txt;
    if (cam.get_distance())
      dist_txt = dtostr(cam.get_distance());
    else
      dist_txt = "   1.2 * SceneWidth";
    fprintf(ofile,
            "   // Distance from viewer to LookAt point (%s)\n"
            "   #declare Distance = 0.9*%s;\n\n",
            dist_txt.c_str(), dist_txt.c_str());

    fprintf(ofile,
            "   // Rotation Centre\n"
            "   #declare RotCentre = %s;\n\n",
            pov_vec(cam.get_centre()).c_str());
    fprintf(ofile,
            "   // Width of perspective (%g) if 0 use default\n"
            "   #declare PerspFactor = %g;\n\n",
            cam.get_persp(), cam.get_persp());
    fprintf(ofile,
            "   // View point, where camera looks (%s)\n"
            "   #declare LookAt = %s;\n\n",
            pov_vec(cam.get_lookat()).c_str(),
            pov_vec(cam.get_lookat()).c_str());

    Vec3d angs = -cam.get_rotation().get_euler() * rad2deg();
    fprintf(ofile,
            "   // Rotation about points centre (%s)\n"
            "   #declare Rotation = %s;\n\n",
            pov_vec(angs).c_str(), pov_vec(angs).c_str());
    fprintf(ofile,
            "   // Shadow, 0 - shadows mono no shadows stero,"
            " 1 - no shadows, 2 - shadows\n"
            "   #declare Shadow = %d;\n\n",
            shadow);
    fprintf(ofile, "   #break\n\n");
  }
  fprintf(ofile, "#else\n"
                 "   #debug \"CameraNumber out of range\"\n\n"
                 "#end        // switch(CameraNumber)\n\n");
}

void PovWriter::include_files(FILE *ofile, const Scene &scen)
{
  fprintf(ofile, "#if(file_exists(\"default_pov.inc\")) #include "
                 "\"default_pov.inc\" #end\n");
  string name = dots2underscores(scen.get_name());
  if (name != "default_pov")
    fprintf(ofile, "#if(file_exists(\"%s.inc\")) #include \"%s.inc\" #end\n",
            name.c_str(), name.c_str());
  for (auto &include : includes)
    fprintf(ofile, "#include \"%s\"\n", include.c_str());
}

void PovWriter::scene_setup(FILE *ofile)
{
  fprintf(
      ofile,
      "\n"
      "// ###########  non-configurable section ############\n"
      "\n"
      "\n"
      "#if (StereoType=1 | StereoType=3)\n"
      "   #declare Distance = Distance*1.4;\n"
      "#end\n"
      "\n"
      "#if (PerspFactor=0)\n"
      "   #if (StereoType=1 | StereoType=3)\n"
      "      #declare PerspFactor=4;\n"
      "   #else\n"
      "      #declare PerspFactor=2;\n"
      "   #end\n"
      "#end\n"
      "\n"
      "// Stereo offset (offset of images or camera from mono position)\n"
      "#declare StereoOffset = Distance/3.3;\n"
      "\n"
      "#declare CamLocOff  = <0, 0, Distance*PerspFactor>;\n"
      "\n"
      "#if (StereoType=2)  // Stereo using separate image file for each view\n"
      "   #if (clock=0)\n"
      "      #declare CamLocOff = CamLocOff + < StereoOffset, 0, 0>;\n"
      "   #else\n"
      "      #declare CamLocOff = CamLocOff + <-StereoOffset, 0, 0>;\n"
      "   #end\n"
      "#end\n"
      "#if (StereoType=3) // Four views, so move the camera back a bit more\n"
      "   #declare CamLocOff = CamLocOff*1.2;\n"
      "#end\n"
      "\n");
}

void PovWriter::stereo_setup(FILE *ofile)
{
  fprintf(ofile, "#if (StereoType=0|StereoType=2) // Mono, or Stereo using two "
                 "image files\n"
                 "   #declare ArrSize = 1;\n"
                 "   #declare Offsets = array[ArrSize] {<0, 0, 0>}\n"
                 "#end\n"
                 "#if (StereoType=1)          // Single Image Stereo\n"
                 "   #declare ArrSize = 2;\n"
                 "   #declare Offsets = array[ArrSize] {<+StereoOffset, 0, 0>, "
                 "<-StereoOffset, 0, 0>}\n"
                 "#end\n"
                 "#if (StereoType=3)          // Tetra View\n"
                 "   #declare ArrSize = 4;\n"
                 "   #declare Offsets = array[ArrSize] {<-StereoOffset, "
                 "+StereoOffset, 0>, <+StereoOffset, +StereoOffset, 0>, "
                 "<+StereoOffset, -StereoOffset, 0>, <-StereoOffset, "
                 "-StereoOffset, 0>}\n"
                 "   #declare OrigRot = Rotation;\n"
                 "   #declare Rots = array [ArrSize] {<0, 0, 0>, <109, 0, 0>, "
                 "<109, 120, 0>, <109, 240, 0>}\n"
                 "#end\n"
                 "\n"
                 "#declare Off=0;\n"
                 "#while(Off<ArrSize)\n"
                 "   #if (StereoType=3)\n"
                 "      #declare Rotation = OrigRot + Rots[Off];\n"
                 "   #end\n"
                 "\n"
                 "   object {\n");
}

void PovWriter::geometry_objects(FILE *ofile, const Scene &scen, int sig_dgts)
{
  fprintf(ofile, "      // Objects\n"
                 "      union {\n");
  const vector<SceneGeometry> &sgeoms = scen.get_geoms();
  for (const auto &sgeom : sgeoms) {
    if (o_type == 's' || o_type == 't') {
      string fname = dots2underscores(sgeom.get_name()) + ".inc";
      fprintf(ofile, "         #include \"%s\"\n", fname.c_str());
      if (o_type == 's') { // write geometries into their own files
        FILE *gfile = fopen(fname.c_str(), "w");
        if (gfile == nullptr) {
          fprintf(stderr, "could not open geometry output file \'%s\'",
                  fname.c_str());
          exit(1);
        }
        const vector<GeometryDisplay *> &disps = sgeom.get_disps();
        for (auto disp : disps)
          disp->pov_geom(gfile, scen, sig_dgts);
        if (sgeom.get_label())
          sgeom.get_label()->pov_geom(gfile, scen, sig_dgts);
        if (sgeom.get_sym())
          sgeom.get_sym()->pov_geom(gfile, scen, sig_dgts);
        fclose(gfile);
      }
    }
    else {
      const vector<GeometryDisplay *> &disps = sgeom.get_disps();
      for (auto disp : disps)
        disp->pov_geom(ofile, scen, sig_dgts);
      if (sgeom.get_label())
        sgeom.get_label()->pov_geom(ofile, scen, sig_dgts);
      if (sgeom.get_sym())
        sgeom.get_sym()->pov_geom(ofile, scen, sig_dgts);
    }
  }

  for (auto &obj_include : obj_includes)
    fprintf(ofile, "         #include \"%s\"\n", obj_include.c_str());

  if (o_type != 'o')
    fprintf(ofile, "\n"
                   "         translate - SceneCentre\n"
                   "         rotate Rotation\n"
                   "         translate SceneCentre + Offsets[Off]\n");

  fprintf(ofile, "      }\n\n");
}

void PovWriter::stereo_clipping(FILE *ofile)
{
  fprintf(ofile,
          "      // Don't do slow clipping unless there is chance of overlap\n"
          "      #if ((StereoType=1|StereoType=3) & ((Distance < "
          "1.2*SceneWidth*1.4) |"
          "(vlength(LookAt-SceneCentre))) )\n");
  fprintf(ofile,
          "         bounded_by {\n"
          "            box { SceneWidth*(<1, 1, 1>), -SceneWidth*<1, 1, 1>\n"
          "               translate SceneWidth*<select(Offsets[Off].x, -1, 0, "
          "1), select(Offsets[Off].y, -1, 0, 1), 0>\n"
          "               translate Distance*PerspFactor*z\n"
          "               rotate "
          "<-3*atan(Offsets[Off].y/(Distance*PerspFactor)), "
          "3*atan(Offsets[Off].x/(Distance*PerspFactor)),0>\n"
          "               translate -Distance*PerspFactor*z+LookAt\n"
          "            }\n"
          "         }\n"
          "         clipped_by { bounded_by }\n"
          "      #end\n"
          "   } // close of object\n"
          "   #declare Off = Off+1;\n"
          "#end\n"
          "\n");
}

void PovWriter::camera_lights(FILE *ofile)
{
  fprintf(
      ofile,
      "\n"
      "background {color BgColour}\n"
      "global_settings{max_trace_level MaxTraceLevel}\n"
      "\n"
      "#ifndef(ExcludeDefCamera) #declare ExcludeDefCamera=0; #end\n"
      "#if(!ExcludeDefCamera)\n"
      "  camera { location LookAt + CamLocOff\n"
      "           look_at LookAt\n"
      "           direction <0 , 0, PerspFactor>\n"
      "           right <-1*AspectRatio , 0, 0>\n"
      "         }\n"
      "#end\n"
      "\n"
      "\n"
      "#ifndef(ExcludeDefLights) #declare ExcludeDefLights=0; #end\n"
      "#if(!ExcludeDefLights)\n"
      "   light_source {LookAt + <0, 0, SceneWidth*10> color <0.7, 0.7, 0.7> \n"
      "      #if (Shadow=1|(!Shadow & StereoType))\n"
      "         shadowless\n"
      "      #end\n"
      "   }\n"
      "\n"
      "   light_source {LookAt + <0 ,2*SceneWidth, SceneWidth*10> color White "
      "\n"
      "      #if (Shadow=1|(!Shadow & StereoType))\n"
      "         shadowless\n"
      "      #end\n"
      "   }\n"
      "#end\n"
      "\n");
}

} // namespace anti
