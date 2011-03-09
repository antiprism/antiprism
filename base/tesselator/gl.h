#ifndef __gl_h_
#define __gl_h_

#ifdef __cplusplus
extern "C" {
#endif

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

typedef unsigned int localGLenum;
typedef unsigned char localGLboolean;
typedef unsigned char localGLubyte;
typedef float localGLfloat;
typedef double localGLdouble;
typedef void localGLvoid;
typedef void (*_localGLfuncptr)();

/*************************************************************/

/* BeginMode */
#define localGL_POINTS                         0x0000
#define localGL_LINES                          0x0001
#define localGL_LINE_LOOP                      0x0002
#define localGL_LINE_STRIP                     0x0003
#define localGL_TRIANGLES                      0x0004
#define localGL_TRIANGLE_STRIP                 0x0005
#define localGL_TRIANGLE_FAN                   0x0006
#define localGL_QUADS                          0x0007
#define localGL_QUAD_STRIP                     0x0008
#define localGL_POLYGON                        0x0009

#ifdef __cplusplus
}
#endif

#endif /* __gl_h_ */
