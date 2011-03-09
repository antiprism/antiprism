/*
** License Applicability. Except to the extent portions of this file are
** made subject to an alternative license as permitted in the SGI Free
** Software License B, Version 1.1 (the "License"), the contents of this
** file are subject only to the provisions of the License. You may not use
** this file except in compliance with the License. You may obtain a copy
** of the License at Silicon Graphics, Inc., attn: Legal Services, 1600
** Amphitheatre Parkway, Mountain View, CA 94043-1351, or at:
** 
** http://oss.sgi.com/projects/FreeB
** 
** Note that, as provided in the License, the Software is distributed on an
** "AS IS" basis, with ALL EXPRESS AND IMPLIED WARRANTIES AND CONDITIONS
** DISCLAIMED, INCLUDING, WITHOUT LIMITATION, ANY IMPLIED WARRANTIES AND
** CONDITIONS OF MERCHANTABILITY, SATISFACTORY QUALITY, FITNESS FOR A
** PARTICULAR PURPOSE, AND NON-INFRINGEMENT.
** 
** Original Code. The Original Code is: OpenGL Sample Implementation,
** Version 1.2.1, released January 26, 2000, developed by Silicon Graphics,
** Inc. The Original Code is Copyright (c) 1991-2000 Silicon Graphics, Inc.
** Copyright in any portions created by third parties is as indicated
** elsewhere herein. All Rights Reserved.
** 
** Additional Notice Provisions: The application programming interfaces
** established by SGI in conjunction with the Original Code are The
** OpenGL(R) Graphics System: A Specification (Version 1.2.1), released
** April 1, 1999; The OpenGL(R) Graphics System Utility Library (Version
** 1.3), released November 4, 1998; and OpenGL(R) Graphics with the X
** Window System(R) (Version 1.3), released October 19, 1998. This software
** was created using the OpenGL(R) version 1.2.1 Sample Implementation
** published by SGI, but has not been independently verified as being
** compliant with the OpenGL(R) version 1.2.1 Specification.
**
*/
/*
** Author: Eric Veach, July 1994.
**
** $Date$ $Revision$
** $Header: //depot/main/gfx/lib/localglu/libtess/mesh.h#5 $
*/

#ifndef __mesh_h_
#define __mesh_h_

#include "glu.h"

typedef struct localGLUmesh localGLUmesh; 

typedef struct localGLUvertex localGLUvertex;
typedef struct localGLUface localGLUface;
typedef struct localGLUhalfEdge localGLUhalfEdge;

typedef struct ActiveRegion ActiveRegion;	/* Internal data */

/* The mesh structure is similar in spirit, notation, and operations
 * to the "quad-edge" structure (see L. Guibas and J. Stolfi, Primitives
 * for the manipulation of general subdivisions and the computation of
 * Voronoi diagrams, ACM Transactions on Graphics, 4(2):74-123, April 1985).
 * For a simplified description, see the course notes for CS348a,
 * "Mathematical Foundations of Computer Graphics", available at the
 * Stanford bookstore (and taught during the fall quarter).
 * The implementation also borrows a tiny subset of the graph-based approach
 * use in Mantyla's Geometric Work Bench (see M. Mantyla, An Introduction
 * to Sold Modeling, Computer Science Press, Rockville, Maryland, 1988).
 *
 * The fundamental data structure is the "half-edge".  Two half-edges
 * go together to make an edge, but they point in opposite directions.
 * Each half-edge has a pointer to its mate (the "symmetric" half-edge Sym),
 * its origin vertex (Org), the face on its left side (Lface), and the
 * adjacent half-edges in the CCW direction around the origin vertex
 * (Onext) and around the left face (Lnext).  There is also a "next"
 * pointer for the global edge list (see below).
 *
 * The notation used for mesh navigation:
 *	Sym   = the mate of a half-edge (same edge, but opposite direction)
 *	Onext = edge CCW around origin vertex (keep same origin)
 *	Dnext = edge CCW around destination vertex (keep same dest)
 *	Lnext = edge CCW around left face (dest becomes new origin)
 *	Rnext = edge CCW around right face (origin becomes new dest)
 *
 * "prev" means to substitute CW for CCW in the definitions above.
 *
 * The mesh keeps global lists of all vertices, faces, and edges,
 * stored as doubly-linked circular lists with a dummy header node.
 * The mesh stores pointers to these dummy headers (vHead, fHead, eHead).
 *
 * The circular edge list is special; since half-edges always occur
 * in pairs (e and e->Sym), each half-edge stores a pointer in only
 * one direction.  Starting at eHead and following the e->next pointers
 * will visit each *edge* once (ie. e or e->Sym, but not both).
 * e->Sym stores a pointer in the opposite direction, thus it is
 * always true that e->Sym->next->Sym->next == e.
 *
 * Each vertex has a pointer to next and previous vertices in the
 * circular list, and a pointer to a half-edge with this vertex as
 * the origin (NULL if this is the dummy header).  There is also a
 * field "data" for client data.
 *
 * Each face has a pointer to the next and previous faces in the
 * circular list, and a pointer to a half-edge with this face as
 * the left face (NULL if this is the dummy header).  There is also
 * a field "data" for client data.
 *
 * Note that what we call a "face" is really a loop; faces may consist
 * of more than one loop (ie. not simply connected), but there is no
 * record of this in the data structure.  The mesh may consist of
 * several disconnected regions, so it may not be possible to visit
 * the entire mesh by starting at a half-edge and traversing the edge
 * structure.
 *
 * The mesh does NOT support isolated vertices; a vertex is deleted along
 * with its last edge.  Similarly when two faces are merged, one of the
 * faces is deleted (see __localgl_meshDelete below).  For mesh operations,
 * all face (loop) and vertex pointers must not be NULL.  However, once
 * mesh manipulation is finished, __localgl_MeshZapFace can be used to delete
 * faces of the mesh, one at a time.  All external faces can be "zapped"
 * before the mesh is returned to the client; then a NULL face indicates
 * a region which is not part of the output polygon.
 */

struct localGLUvertex {
  localGLUvertex	*next;		/* next vertex (never NULL) */
  localGLUvertex	*prev;		/* previous vertex (never NULL) */
  localGLUhalfEdge	*anEdge;	/* a half-edge with this origin */
  void		*data;		/* client's data */

  /* Internal data (keep hidden) */
  localGLdouble	coords[3];	/* vertex location in 3D */
  localGLdouble	s, t;		/* projection onto the sweep plane */
  long		pqHandle;	/* to allow deletion from priority queue */
};

struct localGLUface {
  localGLUface	*next;		/* next face (never NULL) */
  localGLUface	*prev;		/* previous face (never NULL) */
  localGLUhalfEdge	*anEdge;	/* a half edge with this left face */
  void		*data;		/* room for client's data */

  /* Internal data (keep hidden) */
  localGLUface	*trail;		/* "stack" for conversion to strips */
  localGLboolean	marked;		/* flag for conversion to strips */
  localGLboolean	inside;		/* this face is in the polygon interior */
};

struct localGLUhalfEdge {
  localGLUhalfEdge	*next;		/* doubly-linked list (prev==Sym->next) */
  localGLUhalfEdge	*Sym;		/* same edge, opposite direction */
  localGLUhalfEdge	*Onext;		/* next edge CCW around origin */
  localGLUhalfEdge	*Lnext;		/* next edge CCW around left face */
  localGLUvertex	*Org;		/* origin vertex (Overtex too long) */
  localGLUface	*Lface;		/* left face */

  /* Internal data (keep hidden) */
  ActiveRegion	*activeRegion;	/* a region with this upper edge (sweep.c) */
  int		winding;	/* change in winding number when crossing
                                   from the right face to the left face */
};

#define	Rface	Sym->Lface
#define Dst	Sym->Org

#define Oprev	Sym->Lnext
#define Lprev   Onext->Sym
#define Dprev	Lnext->Sym
#define Rprev	Sym->Onext
#define Dnext	Rprev->Sym	/* 3 pointers */
#define Rnext	Oprev->Sym	/* 3 pointers */


struct localGLUmesh {
  localGLUvertex	vHead;		/* dummy header for vertex list */
  localGLUface	fHead;		/* dummy header for face list */
  localGLUhalfEdge	eHead;		/* dummy header for edge list */
  localGLUhalfEdge	eHeadSym;	/* and its symmetric counterpart */
};

/* The mesh operations below have three motivations: completeness,
 * convenience, and efficiency.  The basic mesh operations are MakeEdge,
 * Splice, and Delete.  All the other edge operations can be implemented
 * in terms of these.  The other operations are provided for convenience
 * and/or efficiency.
 *
 * When a face is split or a vertex is added, they are inserted into the
 * global list *before* the existing vertex or face (ie. e->Org or e->Lface).
 * This makes it easier to process all vertices or faces in the global lists
 * without worrying about processing the same data twice.  As a convenience,
 * when a face is split, the "inside" flag is copied from the old face.
 * Other internal data (v->data, v->activeRegion, f->data, f->marked,
 * f->trail, e->winding) is set to zero.
 *
 * ********************** Basic Edge Operations **************************
 *
 * __localgl_meshMakeEdge( mesh ) creates one edge, two vertices, and a loop.
 * The loop (face) consists of the two new half-edges.
 *
 * __localgl_meshSplice( eOrg, eDst ) is the basic operation for changing the
 * mesh connectivity and topology.  It changes the mesh so that
 *	eOrg->Onext <- OLD( eDst->Onext )
 *	eDst->Onext <- OLD( eOrg->Onext )
 * where OLD(...) means the value before the meshSplice operation.
 *
 * This can have two effects on the vertex structure:
 *  - if eOrg->Org != eDst->Org, the two vertices are merged together
 *  - if eOrg->Org == eDst->Org, the origin is split into two vertices
 * In both cases, eDst->Org is changed and eOrg->Org is untouched.
 *
 * Similarly (and independently) for the face structure,
 *  - if eOrg->Lface == eDst->Lface, one loop is split into two
 *  - if eOrg->Lface != eDst->Lface, two distinct loops are joined into one
 * In both cases, eDst->Lface is changed and eOrg->Lface is unaffected.
 *
 * __localgl_meshDelete( eDel ) removes the edge eDel.  There are several cases:
 * if (eDel->Lface != eDel->Rface), we join two loops into one; the loop
 * eDel->Lface is deleted.  Otherwise, we are splitting one loop into two;
 * the newly created loop will contain eDel->Dst.  If the deletion of eDel
 * would create isolated vertices, those are deleted as well.
 *
 * ********************** Other Edge Operations **************************
 *
 * __localgl_meshAddEdgeVertex( eOrg ) creates a new edge eNew such that
 * eNew == eOrg->Lnext, and eNew->Dst is a newly created vertex.
 * eOrg and eNew will have the same left face.
 *
 * __localgl_meshSplitEdge( eOrg ) splits eOrg into two edges eOrg and eNew,
 * such that eNew == eOrg->Lnext.  The new vertex is eOrg->Dst == eNew->Org.
 * eOrg and eNew will have the same left face.
 *
 * __localgl_meshConnect( eOrg, eDst ) creates a new edge from eOrg->Dst
 * to eDst->Org, and returns the corresponding half-edge eNew.
 * If eOrg->Lface == eDst->Lface, this splits one loop into two,
 * and the newly created loop is eNew->Lface.  Otherwise, two disjoint
 * loops are merged into one, and the loop eDst->Lface is destroyed.
 *
 * ************************ Other Operations *****************************
 *
 * __localgl_meshNewMesh() creates a new mesh with no edges, no vertices,
 * and no loops (what we usually call a "face").
 *
 * __localgl_meshUnion( mesh1, mesh2 ) forms the union of all structures in
 * both meshes, and returns the new mesh (the old meshes are destroyed).
 *
 * __localgl_meshDeleteMesh( mesh ) will free all storage for any valid mesh.
 *
 * __localgl_meshZapFace( fZap ) destroys a face and removes it from the
 * global face list.  All edges of fZap will have a NULL pointer as their
 * left face.  Any edges which also have a NULL pointer as their right face
 * are deleted entirely (along with any isolated vertices this produces).
 * An entire mesh can be deleted by zapping its faces, one at a time,
 * in any order.  Zapped faces cannot be used in further mesh operations!
 *
 * __localgl_meshCheckMesh( mesh ) checks a mesh for self-consistency.
 */

localGLUhalfEdge	*__localgl_meshMakeEdge( localGLUmesh *mesh );
int		__localgl_meshSplice( localGLUhalfEdge *eOrg, localGLUhalfEdge *eDst );
int		__localgl_meshDelete( localGLUhalfEdge *eDel );

localGLUhalfEdge	*__localgl_meshAddEdgeVertex( localGLUhalfEdge *eOrg );
localGLUhalfEdge	*__localgl_meshSplitEdge( localGLUhalfEdge *eOrg );
localGLUhalfEdge	*__localgl_meshConnect( localGLUhalfEdge *eOrg, localGLUhalfEdge *eDst );

localGLUmesh		*__localgl_meshNewMesh( void );
localGLUmesh		*__localgl_meshUnion( localGLUmesh *mesh1, localGLUmesh *mesh2 );
void		__localgl_meshDeleteMesh( localGLUmesh *mesh );
void		__localgl_meshZapFace( localGLUface *fZap );

#ifdef NDEBUG
#define		__localgl_meshCheckMesh( mesh )
#else
void		__localgl_meshCheckMesh( localGLUmesh *mesh );
#endif

#endif
