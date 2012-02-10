// Antiprism Resource File - http://www.antiprism.com
// This file may be copied, modified and redistributed
//
// Converts faces to a ring with radial bars to the edge mid-points,
// and edges to a ring that joins the bars meeting there (Adrian Rossiter)
// Original idea proposed by Dick Fischbeck http://www.freewebtown.com/randome/


#include "metals.inc"

//======== configurable variables ==========================
//Set these variables to overide the defaults

//Strut radius is edge radius * strut factor;
//#declare strut_factor = 1; // default 1

//Ring radius is edge radius * ring factor;
//#declare ring_factor = 2; // default 2

//Strut and ring material
//#declare strut_tex = texture { T_Brass_4A } // default texture { T_Brass_4A }

//======== end of configurable variables ====================

#ifndef(strut_factor) #declare strut_factor = 1; #end
#ifndef(ring_factor) #declare ring_factor = 2; #end
#ifndef(strut_tex) #declare strut_tex = texture { T_Brass_4A } #end


#declare strut_rad = edge_sz*strut_factor;
#declare ring_rad = edge_sz*ring_factor;


#macro rot_y_to_vec(vec)
   #if(vec.z=0) #local vec = vec + <1e-6*(abs(vec.y)+abs(vec.z)), 0, 0>; #end
   #if(vec.y=0) #local vec = vec + <0, 1e-6*(abs(vec.x)+abs(vec.z)), 0>; #end
   #if(vec.z=0) #local vec = vec + <0, 0, 1e-6*(abs(vec.x)+abs(vec.y))>; #end
   #local xz_len = abs(vec.z)/abs(vec.z)*sqrt(vec.x*vec.x + vec.z*vec.z);
   //#debug concat("vec = ", vstr(3, vec, ", ", 3, 3), "\n")
   <degrees(atan2(xz_len, vec.y)),
   degrees(atan2(vec.x, vec.z)),
   0>
#end

#macro face_normal(idx)
   #local p1 = verts[faces[idx+1]];
   #local i=1;
   #while (i<faces[idx]-1)
      #local norm = vcross(verts[faces[idx+i+1]]-p1,verts[faces[idx+i+2]]-p1);
      #local sin_a = vlength(norm)/(vlength(verts[faces[idx+2]]-p1)*vlength(verts[faces[idx+3]]-p1));
      #if(abs(sin_a)>1e-6) #local norm=vnormalize(norm); #local i=faces[idx];
      #else #local sin_a=0;
         #end
      #local i=i+1;
      #end
   #if(sin_a!=0)
      #local ret_norm = norm;
   #else
      #local ret_norm = <0, 0, 0>;
      #end
   ret_norm
#end

#macro disp_face(face_no, idx, col)
   #local centroid=0;
   #local i=0;
   #while (i< faces[idx])
      #local centroid=centroid+verts[faces[idx+i+1]];
      #local i = i+1;
      #end
   #local centroid=centroid/faces[idx];

   #local norm = face_normal(idx);
   torus { ring_rad strut_rad
      rotate rot_y_to_vec(norm)
      translate centroid
      texture{ strut_tex }
      }

   #local i=0;
   #while (i < faces[idx])
      #local p0 = verts[faces[idx+i+1]];
      #local p1 = verts[faces[idx + mod(i+1, faces[idx])+1]];
      #local mid_edge = (p0+p1)/2;
      #local strut_v = vnormalize(centroid-mid_edge);
      cylinder{ centroid-strut_v*ring_rad
                mid_edge+strut_v*ring_rad
                strut_rad
         texture{ strut_tex }
      }
      #local i = i+1;
      #end
#end


#macro disp_edge(edge, col)
   object {
      torus { ring_rad strut_rad }
      rotate rot_y_to_vec(verts[edges[edge][0]] - verts[edges[edge][1]])
      translate (verts[edges[edge][0]] + verts[edges[edge][1]])/2
      texture{ strut_tex }
   }
#end

//#macro disp_face(face, idx, col) #end
//#macro disp_edge(edge, col) #end
#macro disp_vertex(vert, col) #end

