/* $Id: kaleido.h,v 3.27 2002-01-06 16:21:30+02 rl Exp $ */
/*
 *****************************************************************************
 * kaleido
 *
 *	Kaleidoscopic construction of uniform polyhedra
 *	Copyright © 1991-2002 Dr. Zvi Har'El <rl@math.technion.ac.il>
 *
 *	Redistribution and use in source and binary forms, with or without
 *	modification, are permitted provided that the following conditions
 *	are met:
 *
 *	1. Redistributions of source code must retain the above copyright
 *	   notice, this list of conditions and the following disclaimer.
 *
 *	2. Redistributions in binary form must reproduce the above copyright
 *	   notice, this list of conditions and the following disclaimer in
 *	   the documentation and/or other materials provided with the
 *	   distribution.
 *
 *	3. The end-user documentation included with the redistribution,
 *	   if any, must include the following acknowledgment:
 *		"This product includes software developed by
 *		 Dr. Zvi Har'El (http://www.math.technion.ac.il/~rl/)."
 *	   Alternately, this acknowledgment may appear in the software itself,
 *	   if and wherever such third-party acknowledgments normally appear.
 *
 *	This software is provided 'as-is', without any express or implied
 *	warranty.  In no event will the author be held liable for any
 *	damages arising from the use of this software.
 *
 *	Author:
 *		Dr. Zvi Har'El,
 *		Technion, Israel Institue of Technology,
 *		Haifa 32000, Israel.
 *		E-Mail: rl@math.technion.ac.il
 *****************************************************************************
 */
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <errno.h>


#ifndef MAXLONG
#define MAXLONG 0x7FFFFFFF
#endif
#ifndef MAXDIGITS
#define MAXDIGITS 10 /* (int)log10((double)MAXLONG) + 1 */
#endif
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#ifndef DBL_EPSILON
#define DBL_EPSILON 2.2204460492503131e-16
#endif
#define BIG_EPSILON 3e-2
#define DEG (180/M_PI)
#define AZ M_PI/7 /* axis azimuth */
#define EL M_PI/17 /* axis elevation */


#define VAL (argp = s[1] ? s + 1 : (argc--, *++argv), \
		s = "-", !argp ? usage() : argp)
#define Err(x) {\
	error.message = x;\
	error.line = __LINE__;\
	error.file = __FILE__;\
	return 0;\
    }
#define Free(lvalue) {\
	if (lvalue) {\
	    free((char*) lvalue);\
	    lvalue=0;\
	}\
    }
#define Matfree(lvalue,n) {\
	if (lvalue) \
	    matfree((char*) lvalue, n);\
	lvalue=0;\
    }
#define Malloc(lvalue,n,type) {\
	if (!(lvalue = (type*) malloc((n) * sizeof(type)))) \
	    Err(0)\
    }
#define Realloc(lvalue,n,type) {\
	if (!(lvalue = (type*) realloc(lvalue, (n) * sizeof(type)))) \
	    Err(0)\
    }
#define Calloc(lvalue,n,type) {\
	if (!(lvalue = (type*) calloc(n, sizeof(type))))\
	    Err(0)\
    }
#define Matalloc(lvalue,n,m,type) {\
	if (!(lvalue = (type**) matalloc(n, (m) * sizeof(type))))\
	    Err(0)\
    }
#define Sprintfrac(lvalue,x) {\
	if (!(lvalue=sprintfrac(x)))\
	    return 0;\
    }
#define numerator(x) (frac(x), frax.n)
#define denominator(x) (frac(x), frax.d)
#define compli(x) (frac(x), (double) frax.n / (frax.n-frax.d))

typedef struct {
    double x, y, z;
} Vector;

typedef struct {
    /* NOTE: some of the int's can be replaced by short's, char's,
	    or even bit fields, at the expense of readability!!!*/
    int index; /* index to the standard list, the array uniform[] */
    int N; /* number of faces types (atmost 5)*/
    int M; /* vertex valency  (may be big for dihedral polyhedra) */
    int V; /* vertex count */
    int E; /* edge count */
    int F; /* face count */
    int D; /* density */
    int chi; /* Euler characteristic */
    int g; /* order of symmetry group */
    int K; /* symmetry type: D=2, T=3, O=4, I=5 */
    int hemi;/* flag hemi polyhedron */
    int onesided;/* flag onesided polyhedron */
    int even; /* removed face in pqr| */
    int *Fi; /* face counts by type (array N)*/
    int *rot; /* vertex configuration (array M of 0..N-1) */
    int *snub; /* snub triangle configuration (array M of 0..1) */
    int *firstrot; /* temporary for vertex generation (array V) */
    int *anti; /* temporary for direction of ideal vertices (array E) */
    int *ftype; /* face types (array F) */
    int **e; /* edges (matrix 2 x E of 0..V-1)*/
    int **dual_e; /* dual edges (matrix 2 x E of 0..F-1)*/
    int **incid; /* vertex-face incidence (matrix M x V of 0..F-1)*/
    int **adj; /* vertex-vertex adjacency (matrix M x V of 0..V-1)*/
    double p[4]; /* p, q and r; |=0 */
    double minr; /* smallest nonzero inradius */
    double gon; /* basis type for dihedral polyhedra */
    double *n; /* number of side of a face of each type (array N) */
    double *m; /* number of faces at a vertex of each type (array N) */
    double *gamma; /* fundamental angles in radians (array N) */
    char *polyform; /* printable Wythoff symbol */
    char *config; /* printable vertex configuration */
    char *name; /* name, standard or manifuctured */
    char *dual_name; /* dual name, standard or manifuctured */
    Vector *v; /* vertex coordinates (array V) */
    Vector *f; /* face coordinates (array F)*/
} Polyhedron;

typedef struct { /* See uniform.h for explanation of the fields */
	const char *Wythoff, *name, *dual;
	short Coxeter, Wenninger;
} Uniform;

typedef struct {
    long n,d;
} Fraction;

typedef struct {
    const char *message;
    int line;
    const char *file;
} Error;

extern Polyhedron *kaleido(char *sym, Uniform *uniform, int last_uniform);
extern Polyhedron *polyalloc(void);
extern Vector rotate(Vector vertex, Vector axis, double angle);
extern Vector sum3(Vector a, Vector b, Vector c);
extern Vector scale(double k, Vector a);
extern Vector sum(Vector a, Vector b);
extern Vector diff(Vector a, Vector b);
extern Vector pole (double r, Vector a, Vector b, Vector c);
extern Vector cross(Vector a, Vector b);
extern double dot(Vector a, Vector b);
extern int same(Vector a, Vector b, double epsilon);
extern char *picfile(Vector *v, int V, int **e, int E, int *anti, int index,
	char *star, char *subtitle, char *subsub, char *prefix,
	int need_rotation);
extern char *sprintfrac(double x);
extern char *nextsym(int input, int index, int *pargc, char ***pargv);
extern void frac(double x);
extern void matfree(void *mat, int rows);
extern void *matalloc(int rows, int row_size);
extern void more(void);
extern void printvec(FILE *fp, Vector v, int digits);
extern void rgbcolor(FILE *fp, int n);

int get_uniform_list(Uniform **uniform);

//extern Uniform uniform[];
extern Fraction frax;
extern Error error;
//extern int last_uniform;
extern int gcont;
