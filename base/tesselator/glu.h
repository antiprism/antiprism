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
** Additional Notice Provisions: This software was created using the
** OpenGL(R) version 1.2.1 Sample Implementation published by SGI, but has
** not been independently verified as being compliant with the OpenGL(R)
** version 1.2.1 Specification.
*/

#ifndef __glu_h__
#define __glu_h__

#include "../const.h"
#include "gl.h"

#ifdef __cplusplus
extern "C" {
#endif

/*************************************************************/

/* ErrorCode */
#define localGLU_INVALID_ENUM                   100900
#define localGLU_INVALID_VALUE                  100901
#define localGLU_OUT_OF_MEMORY                  100902
#define localGLU_INVALID_OPERATION              100904

/* TessCallback */
#define localGLU_TESS_BEGIN                     100100
#define localGLU_BEGIN                          100100
#define localGLU_TESS_VERTEX                    100101
#define localGLU_VERTEX                         100101
#define localGLU_TESS_END                       100102
#define localGLU_END                            100102
#define localGLU_TESS_ERROR                     100103
#define localGLU_TESS_EDGE_FLAG                 100104
#define localGLU_EDGE_FLAG                      100104
#define localGLU_TESS_COMBINE                   100105
#define localGLU_TESS_BEGIN_DATA                100106
#define localGLU_TESS_VERTEX_DATA               100107
#define localGLU_TESS_END_DATA                  100108
#define localGLU_TESS_ERROR_DATA                100109
#define localGLU_TESS_EDGE_FLAG_DATA            100110
#define localGLU_TESS_COMBINE_DATA              100111

/* TessContour */
#define localGLU_CW                             100120
#define localGLU_CCW                            100121
#define localGLU_INTERIOR                       100122
#define localGLU_EXTERIOR                       100123
#define localGLU_UNKNOWN                        100124

/* TessProperty */
#define localGLU_TESS_WINDING_RULE              100140
#define localGLU_TESS_BOUNDARY_ONLY             100141
#define localGLU_TESS_TOLERANCE                 100142

/* TessError */
#define localGLU_TESS_ERROR1                    100151
#define localGLU_TESS_ERROR2                    100152
#define localGLU_TESS_ERROR3                    100153
#define localGLU_TESS_ERROR4                    100154
#define localGLU_TESS_ERROR5                    100155
#define localGLU_TESS_ERROR6                    100156
#define localGLU_TESS_ERROR7                    100157
#define localGLU_TESS_ERROR8                    100158
#define localGLU_TESS_MISSING_BEGIN_POLYGON     100151
#define localGLU_TESS_MISSING_BEGIN_CONTOUR     100152
#define localGLU_TESS_MISSING_END_POLYGON       100153
#define localGLU_TESS_MISSING_END_CONTOUR       100154
#define localGLU_TESS_COORD_TOO_LARGE           100155
#define localGLU_TESS_NEED_COMBINE_CALLBACK     100156

/* TessWinding */
#define localGLU_TESS_WINDING_ODD               TESS_WINDING_ODD
#define localGLU_TESS_WINDING_NONZERO           TESS_WINDING_NONZERO
#define localGLU_TESS_WINDING_POSITIVE          TESS_WINDING_POSITIVE
#define localGLU_TESS_WINDING_NEGATIVE          TESS_WINDING_NEGATIVE
#define localGLU_TESS_WINDING_ABS_GEQ_TWO       TESS_WINDING_ABS_GEQ_TWO
/*
#define localGLU_TESS_WINDING_ODD               100130
#define localGLU_TESS_WINDING_NONZERO           100131
#define localGLU_TESS_WINDING_POSITIVE          100132
#define localGLU_TESS_WINDING_NEGATIVE          100133
#define localGLU_TESS_WINDING_ABS_GEQ_TWO       100134
*/

/*************************************************************/


#ifdef __cplusplus
struct localGLUtesselator;
#else
typedef struct localGLUtesselator localGLUtesselator;
#endif

typedef struct localGLUtesselator localGLUtesselatorObj;
typedef struct localGLUtesselator localGLUtriangulatorObj;

#define localGLU_TESS_MAX_COORD 1.0e150

extern void localgluBeginPolygon (localGLUtesselator* tess);
extern void localgluDeleteTess (localGLUtesselator* tess);
extern void localgluEndPolygon (localGLUtesselator* tess);
extern const localGLubyte * localgluErrorString (localGLenum error);
extern const localGLubyte * localgluGetString (localGLenum name);
extern void localgluGetTessProperty (localGLUtesselator* tess, localGLenum which, localGLdouble* data);
extern localGLUtesselator* localgluNewTess (void);
extern void localgluNextContour (localGLUtesselator* tess, localGLenum type);
extern void localgluTessBeginContour (localGLUtesselator* tess);
extern void localgluTessBeginPolygon (localGLUtesselator* tess, localGLvoid* data);
extern void localgluTessCallback (localGLUtesselator* tess, localGLenum which, _localGLfuncptr CallBackFunc);
extern void localgluTessEndContour (localGLUtesselator* tess);
extern void localgluTessEndPolygon (localGLUtesselator* tess);
extern void localgluTessNormal (localGLUtesselator* tess, localGLdouble valueX, localGLdouble valueY, localGLdouble valueZ);
extern void localgluTessProperty (localGLUtesselator* tess, localGLenum which, localGLdouble data);
extern void localgluTessVertex (localGLUtesselator* tess, localGLdouble *location, localGLvoid* data);

#ifdef __cplusplus
}
#endif

#endif /* __localglu_h__ */
