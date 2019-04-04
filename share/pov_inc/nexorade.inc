// Antiprism Resource File - http://www.antiprism.com
// This file may be copied, modified and redistributed
//
// Render nexorade units (Adrian Rossiter)
// off2pov -v controls the radius and -e the strut length
// Example command (must use use 'rotegrity -O f' to ensure
// that faces are included in the POV file, and must use
// 'off2pov -t no_tri' to ensure the faces are not triangulated)
// off_color -e S geo_3_1 | rotegrity -a n -O f -o nex_geo_3_1.off
// off2pov -t no_tri -v 0.00789 -e 0.45 -i nexorade.inc -o nex_geo_3_1.pov nex_geo_3_1.off
// povray +a +p +H600 +W800 nex_geo_3_1.pov

#include "textures.inc"
#include "transforms.inc"

#declare Rand = seed(0);

#macro local_col_to_tex(col, elem_tex_map, elem_col_map, def_tex)
   #local typ=0;
   #if(col.y>=0) #local typ=1; #end
   #if(!typ & col.y=-1)
      #if(col.x<dimension_size(elem_tex_map,1))
         #ifdef(elem_tex_map[col.x]) #local typ=2; #end
      #end
      #if(!typ & col.x<dimension_size(elem_col_map,1))
         #ifdef(elem_col_map[col.x]) #local typ=3; #end
      #end
   #end
   #switch(typ)
      #case(0) texture{ def_tex } #break;
      #case(1) texture{
        pigment{ wrinkles
          scale vert_sz/2
          rotate <90*rand(Rand),90*rand(Rand),90*rand(Rand)>
          translate <rand(Rand),rand(Rand),rand(Rand)>
          color_map{
            [0.0 rgbt col*0.8]
            [1.0 rgbt col]
          }
        }
      }
      #break;
      #case(2) texture{ elem_tex_map[-col.x]} #break;
      #case(3) texture{ pigment{ color elem_col_map[-col.x]}} #break;
   #end
#end


#ifndef(strut_len) #declare strut_len = edge_sz; #end

#macro disp_face(face_no, idx, col)
  #local v0 = verts[faces[idx+1]];
  #local v1 = verts[faces[idx+2]];
  #local e_vec = v1-v0;
  #local e_len = vlength(e_vec);
  #local e_dir = e_vec/e_len;
  #local e_mid = (v0 + v1)/2;
  #local end0 = e_mid - e_dir*strut_len/2;
  #local end1 = e_mid + e_dir*strut_len/2;
  
  #local col_tex = local_col_to_tex(f_cols[face_no], face_tex_map, face_col_map, face_tex)
  cylinder { end0, end1, vert_sz texture { col_tex } }
#end


// extra object
//sphere{ <0, 0, 0> 1-2*edge_sz texture { pigment {rgbt White} } }
//plane { -<0, 0, 1> 0.1 texture { pigment {rgbt White} } }

#macro disp_vertex(vert, col) #end
#macro disp_edge(edge, col) #end

