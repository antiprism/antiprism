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
** $Header: //depot/main/gfx/lib/localglu/libtess/tess.c#7 $
*/

#include "gluos.h"
#include <stddef.h>
#include <assert.h>
#include <setjmp.h>
#include "memalloc.h"
#include "tess.h"
#include "mesh.h"
#include "normal.h"
#include "sweep.h"
#include "tessmono.h"
#include "render.h"

#define localGLU_TESS_DEFAULT_TOLERANCE 0.0
#define localGLU_TESS_MESH		100112	/* void (*)(localGLUmesh *mesh)	    */

#define TRUE 1
#define FALSE 0

/*ARGSUSED*/ static void localGLAPI noBegin( localGLenum /*type*/ ) {}
/*ARGSUSED*/ static void localGLAPI noEdgeFlag( localGLboolean /*boundaryEdge*/ ) {}
/*ARGSUSED*/ static void localGLAPI noVertex( void * /*data*/ ) {}
/*ARGSUSED*/ static void localGLAPI noEnd( void ) {}
/*ARGSUSED*/ static void localGLAPI noError( localGLenum /*errnum*/ ) {}
/*ARGSUSED*/ static void localGLAPI noCombine( localGLdouble /*coords*/[3], void * /*data*/[4],
                            localGLfloat /*weight*/[4], void ** /*dataOut*/ ) {}
/*ARGSUSED*/ static void localGLAPI noMesh( localGLUmesh * /*mesh*/ ) {}


/*ARGSUSED*/ void localGLAPI __localgl_noBeginData( localGLenum /*type*/,
					     void * /*polygonData*/ ) {}
/*ARGSUSED*/ void localGLAPI __localgl_noEdgeFlagData(
                   localGLboolean /*boundaryEdge*/, 
				       void * /*polygonData*/ ) {}
/*ARGSUSED*/ void localGLAPI __localgl_noVertexData( void * /*data*/,
					      void * /*polygonData*/ ) {}
/*ARGSUSED*/ void localGLAPI __localgl_noEndData( void * /*polygonData*/ ) {}
/*ARGSUSED*/ void localGLAPI __localgl_noErrorData( localGLenum /*errnum*/,
					     void * /*polygonData*/ ) {}
/*ARGSUSED*/ void localGLAPI __localgl_noCombineData(
                      localGLdouble /*coords*/[3],
					       void * /*data*/[4],
					       localGLfloat /*weight*/[4],
					       void ** /*outData*/,
					       void * /*polygonData*/ ) {}

/* Half-edges are allocated in pairs (see mesh.c) */
typedef struct { localGLUhalfEdge e, eSym; } EdgePair;

#define MAX(a,b)	((a) > (b) ? (a) : (b))
#define MAX_FAST_ALLOC	(MAX(sizeof(EdgePair), \
			 MAX(sizeof(localGLUvertex),sizeof(localGLUface))))


localGLUtesselator * localGLAPI
localgluNewTess( void )
{
  localGLUtesselator *tess;

  /* Only initialize fields which can be changed by the api.  Other fields
   * are initialized where they are used.
   */

  if (memInit( MAX_FAST_ALLOC ) == 0) {
     return 0;			/* out of memory */
  }
  tess = (localGLUtesselator *)memAlloc( sizeof( localGLUtesselator ));
  if (tess == NULL) {
     return 0;			/* out of memory */
  }

  tess->state = T_DORMANT;

  tess->normal[0] = 0;
  tess->normal[1] = 0;
  tess->normal[2] = 0;

  tess->relTolerance = localGLU_TESS_DEFAULT_TOLERANCE;
  tess->windingRule = localGLU_TESS_WINDING_ODD;
  tess->flagBoundary = FALSE;
  tess->boundaryOnly = FALSE;

  tess->callBegin = &noBegin;
  tess->callEdgeFlag = &noEdgeFlag;
  tess->callVertex = &noVertex;
  tess->callEnd = &noEnd;

  tess->callError = &noError;
  tess->callCombine = &noCombine;
  tess->callMesh = &noMesh;

  tess->callBeginData= &__localgl_noBeginData;
  tess->callEdgeFlagData= &__localgl_noEdgeFlagData;
  tess->callVertexData= &__localgl_noVertexData;
  tess->callEndData= &__localgl_noEndData;
  tess->callErrorData= &__localgl_noErrorData;
  tess->callCombineData= &__localgl_noCombineData;

  tess->polygonData= NULL;

  return tess;
}

static void MakeDormant( localGLUtesselator *tess )
{
  /* Return the tessellator to its original dormant state. */

  if( tess->mesh != NULL ) {
    __localgl_meshDeleteMesh( tess->mesh );
  }
  tess->state = T_DORMANT;
  tess->lastEdge = NULL;
  tess->mesh = NULL;
}

#define RequireState( tess, s )   if( tess->state != s ) GotoState(tess,s)

static void GotoState( localGLUtesselator *tess, enum TessState newState )
{
  while( tess->state != newState ) {
    /* We change the current state one level at a time, to get to
     * the desired state.
     */
    if( tess->state < newState ) {
      switch( tess->state ) {
      case T_DORMANT:
	CALL_ERROR_OR_ERROR_DATA( localGLU_TESS_MISSING_BEGIN_POLYGON );
	localgluTessBeginPolygon( tess, NULL );
	break;
      case T_IN_POLYGON:
	CALL_ERROR_OR_ERROR_DATA( localGLU_TESS_MISSING_BEGIN_CONTOUR );
	localgluTessBeginContour( tess );
	break;
      case T_IN_CONTOUR:  // No action, to remove warnings
	break;
      }
    } else {
      switch( tess->state ) {
      case T_IN_CONTOUR:
	CALL_ERROR_OR_ERROR_DATA( localGLU_TESS_MISSING_END_CONTOUR );
	localgluTessEndContour( tess );
	break;
      case T_IN_POLYGON:
	CALL_ERROR_OR_ERROR_DATA( localGLU_TESS_MISSING_END_POLYGON );
	/* localgluTessEndPolygon( tess ) is too much work! */
	MakeDormant( tess );
	break;
      case T_DORMANT:  // No action, to remove warnings
	break;
      }
    }
  }
}


void localGLAPI
localgluDeleteTess( localGLUtesselator *tess )
{
  RequireState( tess, T_DORMANT );
  memFree( tess );
}


void localGLAPI
localgluTessProperty( localGLUtesselator *tess, localGLenum which, localGLdouble value )
{
  localGLenum windingRule;

  switch( which ) {
  case localGLU_TESS_TOLERANCE:
    if( value < 0.0 || value > 1.0 ) break;
    tess->relTolerance = value;
    return;

  case localGLU_TESS_WINDING_RULE:
    windingRule = (localGLenum) value;
    if( windingRule != value ) break;	/* not an integer */

    switch (windingRule) {
    case localGLU_TESS_WINDING_ODD:
    case localGLU_TESS_WINDING_NONZERO:
    case localGLU_TESS_WINDING_POSITIVE:
    case localGLU_TESS_WINDING_NEGATIVE:
    case localGLU_TESS_WINDING_ABS_GEQ_TWO:
      tess->windingRule = windingRule;
      return;
    default:
      break;
    }
    break;

  case localGLU_TESS_BOUNDARY_ONLY:
    tess->boundaryOnly = (value != 0);
    return;

  default:
    CALL_ERROR_OR_ERROR_DATA( localGLU_INVALID_ENUM );
    return;
  }
  CALL_ERROR_OR_ERROR_DATA( localGLU_INVALID_VALUE );
}

/* Returns tessellator property */
void localGLAPI
localgluGetTessProperty( localGLUtesselator *tess, localGLenum which, localGLdouble *value )
{
   switch (which) {
   case localGLU_TESS_TOLERANCE:
      /* tolerance should be in range [0..1] */
      assert(0.0 <= tess->relTolerance && tess->relTolerance <= 1.0);
      *value= tess->relTolerance;
      break;    
   case localGLU_TESS_WINDING_RULE:
      assert(tess->windingRule == localGLU_TESS_WINDING_ODD ||
	     tess->windingRule == localGLU_TESS_WINDING_NONZERO ||
	     tess->windingRule == localGLU_TESS_WINDING_POSITIVE ||
	     tess->windingRule == localGLU_TESS_WINDING_NEGATIVE ||
	     tess->windingRule == localGLU_TESS_WINDING_ABS_GEQ_TWO);
      *value= tess->windingRule;
      break;
   case localGLU_TESS_BOUNDARY_ONLY:
      assert(tess->boundaryOnly == TRUE || tess->boundaryOnly == FALSE);
      *value= tess->boundaryOnly;
      break;
   default:
      *value= 0.0;
      CALL_ERROR_OR_ERROR_DATA( localGLU_INVALID_ENUM );
      break;
   }
} /* localgluGetTessProperty() */

void localGLAPI
localgluTessNormal( localGLUtesselator *tess, localGLdouble x, localGLdouble y, localGLdouble z )
{
  tess->normal[0] = x;
  tess->normal[1] = y;
  tess->normal[2] = z;
}

void localGLAPI
localgluTessCallback( localGLUtesselator *tess, localGLenum which, void (localGLAPI *fn)())
{
  switch( which ) {
  case localGLU_TESS_BEGIN:
    tess->callBegin = (fn == NULL) ? &noBegin : (void (localGLAPI *)(localGLenum)) fn;
    return;
  case localGLU_TESS_BEGIN_DATA:
    tess->callBeginData = (fn == NULL) ?
	&__localgl_noBeginData : (void (localGLAPI *)(localGLenum, void *)) fn;
    return;
  case localGLU_TESS_EDGE_FLAG:
    tess->callEdgeFlag = (fn == NULL) ? &noEdgeFlag :
					(void (localGLAPI *)(localGLboolean)) fn;
    /* If the client wants boundary edges to be flagged,
     * we render everything as separate triangles (no strips or fans).
     */
    tess->flagBoundary = (fn != NULL);
    return;
  case localGLU_TESS_EDGE_FLAG_DATA:
    tess->callEdgeFlagData= (fn == NULL) ?
	&__localgl_noEdgeFlagData : (void (localGLAPI *)(localGLboolean, void *)) fn; 
    /* If the client wants boundary edges to be flagged,
     * we render everything as separate triangles (no strips or fans).
     */
    tess->flagBoundary = (fn != NULL);
    return;
  case localGLU_TESS_VERTEX:
    tess->callVertex = (fn == NULL) ? &noVertex :
				      (void (localGLAPI *)(void *)) fn;
    return;
  case localGLU_TESS_VERTEX_DATA:
    tess->callVertexData = (fn == NULL) ?
	&__localgl_noVertexData : (void (localGLAPI *)(void *, void *)) fn;
    return;
  case localGLU_TESS_END:
    tess->callEnd = (fn == NULL) ? &noEnd : (void (localGLAPI *)(void)) fn;
    return;
  case localGLU_TESS_END_DATA:
    tess->callEndData = (fn == NULL) ? &__localgl_noEndData : 
                                       (void (localGLAPI *)(void *)) fn;
    return;
  case localGLU_TESS_ERROR:
    tess->callError = (fn == NULL) ? &noError : (void (localGLAPI *)(localGLenum)) fn;
    return;
  case localGLU_TESS_ERROR_DATA:
    tess->callErrorData = (fn == NULL) ?
	&__localgl_noErrorData : (void (localGLAPI *)(localGLenum, void *)) fn;
    return;
  case localGLU_TESS_COMBINE:
    tess->callCombine = (fn == NULL) ? &noCombine :
	(void (localGLAPI *)(localGLdouble [3],void *[4], localGLfloat [4], void ** )) fn;
    return;
  case localGLU_TESS_COMBINE_DATA:
    tess->callCombineData = (fn == NULL) ? &__localgl_noCombineData :
                                           (void (localGLAPI *)(localGLdouble [3],
						     void *[4], 
						     localGLfloat [4], 
						     void **,
						     void *)) fn;
    return;
  case localGLU_TESS_MESH:
    tess->callMesh = (fn == NULL) ? &noMesh : (void (localGLAPI *)(localGLUmesh *)) fn;
    return;
  default:
    CALL_ERROR_OR_ERROR_DATA( localGLU_INVALID_ENUM );
    return;
  }
}

static int AddVertex( localGLUtesselator *tess, localGLdouble coords[3], void *data )
{
  localGLUhalfEdge *e;

  e = tess->lastEdge;
  if( e == NULL ) {
    /* Make a self-loop (one vertex, one edge). */

    e = __localgl_meshMakeEdge( tess->mesh );
    if (e == NULL) return 0;
    if ( !__localgl_meshSplice( e, e->Sym ) ) return 0;
  } else {
    /* Create a new vertex and edge which immediately follow e
     * in the ordering around the left face.
     */
    if (__localgl_meshSplitEdge( e ) == NULL) return 0;
    e = e->Lnext;
  }

  /* The new vertex is now e->Org. */
  e->Org->data = data;
  e->Org->coords[0] = coords[0];
  e->Org->coords[1] = coords[1];
  e->Org->coords[2] = coords[2];
  
  /* The winding of an edge says how the winding number changes as we
   * cross from the edge''s right face to its left face.  We add the
   * vertices in such an order that a CCW contour will add +1 to
   * the winding number of the region inside the contour.
   */
  e->winding = 1;
  e->Sym->winding = -1;

  tess->lastEdge = e;

  return 1;
}


static void CacheVertex( localGLUtesselator *tess, localGLdouble coords[3], void *data )
{
  CachedVertex *v = &tess->cache[tess->cacheCount];

  v->data = data;
  v->coords[0] = coords[0];
  v->coords[1] = coords[1];
  v->coords[2] = coords[2];
  ++tess->cacheCount;
}


static int EmptyCache( localGLUtesselator *tess )
{
  CachedVertex *v = tess->cache;
  CachedVertex *vLast;

  tess->mesh = __localgl_meshNewMesh();
  if (tess->mesh == NULL) return 0;

  for( vLast = v + tess->cacheCount; v < vLast; ++v ) {
    if ( !AddVertex( tess, v->coords, v->data ) ) return 0;
  }
  tess->cacheCount = 0;
  tess->emptyCache = FALSE;

  return 1;
}


void localGLAPI
localgluTessVertex( localGLUtesselator *tess, localGLdouble coords[3], void *data )
{
  int i, tooLarge = FALSE;
  localGLdouble x, clamped[3];

  RequireState( tess, T_IN_CONTOUR );

  if( tess->emptyCache ) {
    if ( !EmptyCache( tess ) ) {
       CALL_ERROR_OR_ERROR_DATA( localGLU_OUT_OF_MEMORY );
       return;
    }
    tess->lastEdge = NULL;
  }
  for( i = 0; i < 3; ++i ) {
    x = coords[i];
    if( x < - localGLU_TESS_MAX_COORD ) {
      x = - localGLU_TESS_MAX_COORD;
      tooLarge = TRUE;
    }
    if( x > localGLU_TESS_MAX_COORD ) {
      x = localGLU_TESS_MAX_COORD;
      tooLarge = TRUE;
    }
    clamped[i] = x;
  }
  if( tooLarge ) {
    CALL_ERROR_OR_ERROR_DATA( localGLU_TESS_COORD_TOO_LARGE );
  }

  if( tess->mesh == NULL ) {
    if( tess->cacheCount < TESS_MAX_CACHE ) {
      CacheVertex( tess, clamped, data );
      return;
    }
    if ( !EmptyCache( tess ) ) {
       CALL_ERROR_OR_ERROR_DATA( localGLU_OUT_OF_MEMORY );
       return;
    }
  }
  if ( !AddVertex( tess, clamped, data ) ) {
       CALL_ERROR_OR_ERROR_DATA( localGLU_OUT_OF_MEMORY );
  }
}


void localGLAPI
localgluTessBeginPolygon( localGLUtesselator *tess, void *data )
{
  RequireState( tess, T_DORMANT );

  tess->state = T_IN_POLYGON;
  tess->cacheCount = 0;
  tess->emptyCache = FALSE;
  tess->mesh = NULL;

  tess->polygonData= data;
}


void localGLAPI
localgluTessBeginContour( localGLUtesselator *tess )
{
  RequireState( tess, T_IN_POLYGON );

  tess->state = T_IN_CONTOUR;
  tess->lastEdge = NULL;
  if( tess->cacheCount > 0 ) {
    /* Just set a flag so we don't get confused by empty contours
     * -- these can be generated accidentally with the obsolete
     * NextContour() interface.
     */
    tess->emptyCache = TRUE;
  }
}


void localGLAPI
localgluTessEndContour( localGLUtesselator *tess )
{
  RequireState( tess, T_IN_CONTOUR );
  tess->state = T_IN_POLYGON;
}

void localGLAPI
localgluTessEndPolygon( localGLUtesselator *tess )
{
  localGLUmesh *mesh;

  if (setjmp(tess->env) != 0) {	
     /* come back here if out of memory */
     CALL_ERROR_OR_ERROR_DATA( localGLU_OUT_OF_MEMORY );
     return;
  }

  RequireState( tess, T_IN_POLYGON );
  tess->state = T_DORMANT;

  if( tess->mesh == NULL ) {
    if( ! tess->flagBoundary && tess->callMesh == &noMesh ) {

      /* Try some special code to make the easy cases go quickly
       * (eg. convex polygons).  This code does NOT handle multiple contours,
       * intersections, edge flags, and of course it does not generate
       * an explicit mesh either.
       */
      if( __localgl_renderCache( tess )) {
	tess->polygonData= NULL; 
	return;
      }
    }
    if ( !EmptyCache( tess ) ) longjmp(tess->env,1); /* could've used a label*/
  }

  /* Determine the polygon normal and project vertices onto the plane
   * of the polygon.
   */
  __localgl_projectPolygon( tess );

  /* __localgl_computeInterior( tess ) computes the planar arrangement specified
   * by the given contours, and further subdivides this arrangement
   * into regions.  Each region is marked "inside" if it belongs
   * to the polygon, according to the rule given by tess->windingRule.
   * Each interior region is guaranteed be monotone.
   */
  if ( !__localgl_computeInterior( tess ) ) {
     longjmp(tess->env,1);	/* could've used a label */
  }

  mesh = tess->mesh;
  if( ! tess->fatalError ) {
    int rc = 1;

    /* If the user wants only the boundary contours, we throw away all edges
     * except those which separate the interior from the exterior.
     * Otherwise we tessellate all the regions marked "inside".
     */
    if( tess->boundaryOnly ) {
      rc = __localgl_meshSetWindingNumber( mesh, 1, TRUE );
    } else {
      rc = __localgl_meshTessellateInterior( mesh ); 
    }
    if (rc == 0) longjmp(tess->env,1);	/* could've used a label */

    __localgl_meshCheckMesh( mesh );

    if( tess->callBegin != &noBegin || tess->callEnd != &noEnd
       || tess->callVertex != &noVertex || tess->callEdgeFlag != &noEdgeFlag 
       || tess->callBeginData != &__localgl_noBeginData 
       || tess->callEndData != &__localgl_noEndData
       || tess->callVertexData != &__localgl_noVertexData
       || tess->callEdgeFlagData != &__localgl_noEdgeFlagData )
    {
      if( tess->boundaryOnly ) {
	__localgl_renderBoundary( tess, mesh );  /* output boundary contours */
      } else {
	__localgl_renderMesh( tess, mesh );	   /* output strips and fans */
      }
    }
    if( tess->callMesh != &noMesh ) {

      /* Throw away the exterior faces, so that all faces are interior.
       * This way the user doesn't have to check the "inside" flag,
       * and we don't need to even reveal its existence.  It also leaves
       * the freedom for an implementation to not generate the exterior
       * faces in the first place.
       */
      __localgl_meshDiscardExterior( mesh );
      (*tess->callMesh)( mesh );		/* user wants the mesh itself */
      tess->mesh = NULL;
      tess->polygonData= NULL;
      return;
    }
  }
  __localgl_meshDeleteMesh( mesh );
  tess->polygonData= NULL;
  tess->mesh = NULL;
}


/*XXXblythe unused function*/
#if 0
void localGLAPI
localgluDeleteMesh( localGLUmesh *mesh )
{
  __localgl_meshDeleteMesh( mesh );
}
#endif



/*******************************************************/

/* Obsolete calls -- for backward compatibility */

void localGLAPI
localgluBeginPolygon( localGLUtesselator *tess )
{
  localgluTessBeginPolygon( tess, NULL );
  localgluTessBeginContour( tess );
}


/*ARGSUSED*/
void localGLAPI
localgluNextContour( localGLUtesselator *tess, localGLenum /*type*/ )
{
  localgluTessEndContour( tess );
  localgluTessBeginContour( tess );
}


void localGLAPI
localgluEndPolygon( localGLUtesselator *tess )
{
  localgluTessEndContour( tess );
  localgluTessEndPolygon( tess );
}
