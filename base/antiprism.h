#ifndef ANTIPRISM_H
#define ANTIPRISM_H

/**\mainpage Antiprism Library Documentation
 * \section features Features
 *
 * The Antiprism library provides a framework for working
 * with polyhedra. It includes:
 * - import of OFF, LG3D and simple text coordinate files
 * - export of OFF, VRML, POV, LG3D and text coordinate files
 * - vector and matrix operations
 * - geometric utilities
 * - analysis of polyhedra
 * - polyhedron operations
 * - colouring operations
 * - symmetry operations
 * - models of commonly used polyhedra
 *
 */

/*!\file antiprism.h
 * \brief Includes all the headers need to use the Antiprism library
 */

#include "const.h"
#include "vec3d.h"
#include "mat3d.h"
#include "vec_utils.h"
#include "vec4d.h"
#include "mat4d.h"
#include "col_val.h"
#include "col_geom.h"
#include "geom.h"
#include "col_map.h"
#include "coloring.h"
#include "symmetry.h"
#include "transforms.h"
#include "geom_utils.h"
#include "math_utils.h"
#include "ultragetopt.h"
#include "utils.h"
#include "info.h"
#include "bbox.h"
#include "timing.h"
#include "scene.h"
#include "disp_poly.h"
#include "vrml_writer.h"
#include "pov_writer.h"
#include "polygons.h"

#endif // ANTIPRISM_H
