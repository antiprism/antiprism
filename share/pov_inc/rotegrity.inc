// Antiprism Resource File - http://www.antiprism.com
// This file may be copied, modified and redistributed
//
// Render rotegrity units (Adrian Rossiter)
// Example command (must use use 'rortegrit -O f' to ensure
// that the faces are included in the POV file, and must use
// 'off2pov -t no_tri' to ensure the faces are not triangulated)
// off_color -e S geo_3_1 | rotegrity -O f | off2pov -t no_tri -v 0.02 -i rotegrity.inc > rot_geo_3_1.pov


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

#declare incr = vert_sz/2;

#macro disp_face(face_no, idx, col)
  #local norm = vnormalize(face_normal(idx));
  #local v0 = verts[faces[idx+1]];
  #local v1 = verts[faces[idx+2]];
  #local edge_vec = v1 - v0;
  #local edge_mid_dir = vnormalize((v0 + v1)/2);
  #local edge_dir = vnormalize(edge_vec);
  #local axis = vcross(norm, edge_dir);
  #local cos_a = sqrt(1 - pow(vlength(edge_vec/2), 2)); // a=half central angle
  #local d = 1 - (1+incr)*cos_a;        // arc height of shallow arc
  #local e = (1+incr)*sqrt(1-cos_a*cos_a);  // half chord length of shallow arc
  #local rt = (e*e + d*d)/(2*d);
  #if (rt>0)
    #local R = sqrt(rt);   // radius of shallow arc
  #else
    #local R = 10000;      // far away for flat unit
  #end
  #local col_tex = local_col_to_tex(f_cols[face_no], face_tex_map, face_col_map, face_tex)
  intersection {
    difference {
      sphere { <0, 0, 0> R+incr/2 }
      sphere { <0, 0, 0> R }
      translate -edge_mid_dir*(R-1)
    }
    union {
      intersection {
        plane { norm vert_sz }
        plane { -norm vert_sz }
        plane { vcross(v0, norm) 0 inverse}
        plane { -vcross(v1, norm) 0 inverse}
      }
      cylinder{ <0, 0, 0> (1+2*incr)*v0 vert_sz }
      cylinder{ <0, 0, 0> (1+2*incr)*v1 vert_sz }
    }
    texture { col_tex }
  }
#end

#macro disp_vertex(vert, col)
  #local dir = vnormalize(verts[vert]);
  #local rad = vert_sz*0.5;
  #local ht = incr;
  union {
    cylinder { (1+ht)*dir dir, rad }
    sphere { <0, 0, 0>
      rad
      Axial_Scale_Trans(dir, 0.5)
      translate (1+ht)*dir
    }
    sphere { <0, 0, 0>
      rad
      Axial_Scale_Trans(dir, 0.5)
      translate dir
    }
    texture {Silver2}
  }
#end

// extra object
sphere{ <0, 0, 0> 1 texture{ pigment{ color <1, 1, 1>} } }

#macro disp_edge(edge, col) #end

