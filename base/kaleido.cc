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
 *		Deptartment of Mathematics,
 *		Technion, Israel Institue of Technology,
 *		Haifa 32000, Israel.
 *		E-Mail: rl@math.technion.ac.il
 *****************************************************************************
 */

/* Changes 16-01-05: Adrian Rossiter <adrian_r@terra.es>
   Added function prototypes and braces to avoid warning messages
 */

#include "kaleido.h"

#define GR_DIRH_IDX 74

int unpacksym(char *sym, Polyhedron *P, Uniform *uniform, int last_uniform);
int moebius(Polyhedron *P);
int decompose(Polyhedron *P);
int newton(Polyhedron *P);
int exceptions(Polyhedron *P);
int count(Polyhedron *P);
int vertices(Polyhedron *P);
int faces(Polyhedron *P);
int mod (int i, int j);

Error error;

Polyhedron *kaleido(char *sym, Uniform *uniform, int last_uniform)
{
   Polyhedron *P;
   /*
    * Allocate a Polyhedron structure P.
    */
   if (!(P = polyalloc()))
      return 0;

   /*
    * Unpack input symbol into P.
    */
   if (!unpacksym(sym, P, uniform, last_uniform))
      return 0;
   /*
    * Find Mebius triangle, its density and Euler characteristic.
    */
   if (!moebius(P))
      return 0;
   
   /*
    * Decompose Schwarz triangle.
    */
   if (!decompose(P))
      return 0;
   /*
    * Solve Fundamental triangles, optionally printing approximations.
    */
   if (!newton(P))
      return 0;
   /*
    * Deal with exceptional polyhedra.
    */
   if (!exceptions(P))
      return 0;
   /*
    * Count edges and faces, update density and characteristic if needed.
    */
   if (!count(P))
      return 0;
   /*
    * Compute coordinates.
    */
   if (!vertices(P))
      return 0;
   if (!faces (P))
      return 0;
   
   return P;
}

/*
 * Allocate a blank Polyhedron structure and initialize some of its nonblank
 * fields.
 *
 * Array and matrix field are allocated when needed.
 */
Polyhedron *
polyalloc()
{
   Polyhedron *P;
   Calloc(P, 1, Polyhedron)
   P->index = -1;
   P->even = -1;
   P->K = 2;
   return P;
}

/*
 * Free the struture allocated by polyalloc(), as well as all the array and
 * matrix fields.
 */
void
polyfree(Polyhedron *P)
{
   Free(P->Fi)
   Free(P->n)
   Free(P->m)
   Free(P->gamma)
   Free(P->rot)
   Free(P->snub)
   Free(P->firstrot)
   Free(P->anti)
   Free(P->ftype)
   Free(P->polyform)
   Free(P->config)
   if (P->index < 0) {
      Free(P->name)
      Free(P->dual_name)
   }
   
   Free(P->v)
   Free(P->f)
   Matfree(P->e, 2)
   Matfree(P->dual_e, 2)
   Matfree(P->incid, P->M)
   Matfree(P->adj, P->M)
   free(P);
}

void *
matalloc(int rows, int row_size)
{
   void **mat;
   int i = 0;
   if (!(mat = (void **)malloc(rows * sizeof (void *))))
      return 0;
   while ((mat[i] = malloc(row_size)) && ++i < rows)
      ;
   if (i == rows)
      return (void *)mat;
   while (--i >= 0)
      free(mat[i]);
   free(mat);
   return 0;
}

void
matfree(void *mat, int rows)
{
   while (--rows >= 0)
      free(((void **)mat)[rows]);
   free(mat);
}

/*
 * Unpack input symbol: Wythoff symbol or an index to uniform[]. The symbol is
 * a # followed by a number, or a three fractions and a bar in some order. We
 * allow no bars only if it result from the input symbol #80.
 */
int
unpacksym(char *sym, Polyhedron *P, Uniform *uniform, int last_uniform)
{
   int i = 0, n, d, bars = 0;
   char c;
   while ((c = *sym++) && isspace(c))
      ;
   if (!c) Err("no data")
   if (c == '#') {
      while ((c = *sym++) && isspace(c))
         ;
      if (!c)
         Err("no digit after #")
      if (!isdigit(c))
         Err("not a digit")
      n = c - '0';
      while ((c = *sym++) && isdigit(c))
         n = n * 10 + c - '0';
      if (!n)
	      Err("zero index")
	   if (n > last_uniform)
         Err("index too big")
	   sym--;
	   while ((c = *sym++) && isspace(c))
	      ;
      if (c)
	      Err("data exceeded")
      char wyth[256];
      strncpy(wyth, uniform[P->index = n - 1].Wythoff, 256);
      sym = wyth;
   } else
	sym--;
    for (;;) {
	while ((c = *sym++) && isspace(c))
	    ;
	if (!c) {
	    if (i == 4 && (bars || P->index == GR_DIRH_IDX))
		return 1;
	    if (!bars)
		Err("no bars")
	    Err("not enough fractions")
	}
	if (i == 4)
	    Err("data exceeded")
	if (c == '|'){
	    if (++bars > 1)
		Err("too many bars")
	    P->p[i++] = 0;
	    continue;
	}
	if (!isdigit(c))
	    Err("not a digit")
	n = c - '0';
	while ((c = *sym++) && isdigit(c))
	    n = n * 10 + c - '0';
	if (c && isspace (c))
	    while ((c = *sym++) && isspace(c))
		;
	if (c != '/') {
	    sym--;
	    if ((P->p[i++] = n) <= 1)
		Err("fraction<=1")
	    continue;
	}
	while ((c = *sym++) && isspace(c))
	    ;
	if (!c || !isdigit(c))
	    return 0;
	d = c - '0';
	while ((c = *sym++) && isdigit(c))
	    d = d * 10 + c - '0';
	if (!d)
	    Err("zero denominator")
	sym--;
	if ((P->p[i++] = (double) n / d) <= 1)
	    Err("fraction<=1")
    }
}

/*
 * Using Wythoff symbol (p|qr, pq|r, pqr| or |pqr), find the Moebius triangle
 * (2 3 K) (or (2 2 n)) of the Schwarz triangle (pqr), the order g of its
 * symmetry group, its Euler characteristic chi, and its covering density D.
 * g is the number of copies of (2 3 K) covering the sphere, i.e.,
 *
 *		g * pi * (1/2 + 1/3 + 1/K - 1) = 4 * pi
 *
 * D is the number of times g copies of (pqr) cover the sphere, i.e.
 *
 *		D * 4 * pi = g * pi * (1/p + 1/q + 1/r - 1)
 *
 * chi is V - E + F, where F = g is the number of triangles, E = 3*g/2 is the
 * number of triangle edges, and V = Vp+ Vq+ Vr, with Vp = g/(2*np) being the
 * number of vertices with angle pi/p (np is the numerator of p).
 */
int
moebius(Polyhedron *P)
{
    int twos = 0, j, len = 1;
/*
 * Arrange Wythoff symbol in a presentable form. In the same time check the
 * restrictions on the three fractions: They all have to be greater then one,
 * and the numerators 4 or 5 cannot occur together.  We count the ocurrences of
 * 2 in `two', and save the largest numerator in `P->K', since they reflect on
 * the symmetry group.
 */
    P->K = 2;
    if (P->index == GR_DIRH_IDX) {
	Malloc(P->polyform, ++len, char)
	strcpy(P->polyform, "|");
    } else
	Calloc(P->polyform, len, char)
    for (j = 0; j < 4; j++) {
	if (P->p[j]) {
	    char *s;
	    Sprintfrac(s, P->p[j])
	    if (j && P->p[j-1]) {
		Realloc(P->polyform, len += strlen (s) + 1, char)
		strcat(P->polyform, " ");
	    } else
		Realloc (P->polyform, len += strlen (s), char)
	    strcat(P->polyform, s);
	    free(s);
	    if (P->p[j] != 2) {
		int k;
		if ((k = numerator (P->p[j])) > P->K) {
		    if (P->K == 4)
			break;
		    P->K = k;
		} else if (k < P->K && k == 4)
		    break;
	    } else
		twos++;
	} else  {
	    Realloc(P->polyform, ++len, char)
	    strcat(P->polyform, "|");
	}
    }
/*
 * Find the symmetry group P->K (where 2, 3, 4, 5 represent the dihedral,
 * tetrahedral, octahedral and icosahedral groups, respectively), and its order
 * P->g.
 */
    if (twos >= 2) {/* dihedral */
	P->g = 4 * P->K;
	P->K = 2;
    } else {
	if (P->K > 5)
	    Err("numerator too large")
	P->g = 24 * P->K / (6 - P->K);
    }
/*
 * Compute the nominal density P->D and Euler characteristic P->chi.
 * In few exceptional cases, these values will be modified later.
 */
    if (P->index != GR_DIRH_IDX) {
	int i;
	P->D = P->chi = - P->g;
	for (j = 0; j < 4; j++) if (P->p[j]) {
	    P->chi += i = P->g / numerator(P->p[j]);
	    P->D += i * denominator(P->p[j]);
	}
	P->chi /= 2;
	P->D /= 4;
	if (P->D <= 0)
	    Err("nonpositive density")
    }
    return 1;
}

/*
 * Decompose Schwarz triangle into N right triangles and compute the vertex
 * count V and the vertex valency M.  V is computed from the number g of
 * Schwarz triangles in the cover, divided by the number of triangles which
 * share a vertex. It is halved for one-sided polyhedra, because the
 * kaleidoscopic construction really produces a double orientable covering of
 * such polyhedra. All q' q|r are of the "hemi" type, i.e. have equatorial {2r}
 * faces, and therefore are (except 3/2 3|3 and the dihedra 2 2|r) one-sided. A
 * well known example is 3/2 3|4, the "one-sided heptahedron". Also, all p q r|
 * with one even denominator have a crossed parallelogram as a vertex figure,
 * and thus are one-sided as well.
 */
int
decompose(Polyhedron *P)
{
    int j, J, *s, *t;
    if (!P->p[1]) { /* p|q r */
	P->N = 2;
	P->M = 2 * numerator(P->p[0]);
	P->V = P->g / P->M;
	Malloc(P->n, P->N, double)
	Malloc(P->m, P->N, double)
	Malloc(P->rot, P->M, int)
	s = P->rot;
	for (j = 0; j < 2; j++) {
	    P->n[j] = P->p[j+2];
	    P->m[j] = P->p[0];
	}
	for (j = P->M / 2; j--;) {
	    *s++ = 0;
	    *s++ = 1;
	}
    } else if (!P->p[2]) { /* p q|r */
	P->N = 3;
	P->M = 4;
	P->V = P->g / 2;
	Malloc(P->n, P->N, double)
	Malloc(P->m, P->N, double)
	Malloc(P->rot, P->M, int)
	s = P->rot;
	P->n[0] = 2 * P->p[3];
	P->m[0] = 2;
	for (j = 1; j < 3; j++) {
	    P->n[j] = P->p[j-1];
	    P->m[j] = 1;
	    *s++ = 0;
	    *s++ = j;
	}
	if (fabs(P->p[0] - compli (P->p[1])) < DBL_EPSILON) {/* p = q' */
		    /* P->p[0]==compli(P->p[1]) should work.  However, MSDOS
		     * yeilds a 7e-17 difference! Reported by Jim Buddenhagen
		     * <jb1556@daditz.sbc.com> */
	    P->hemi = 1;
	    P->D = 0;
	    if (P->p[0] != 2 && !(P->p[3] == 3 && (P->p[0] == 3 ||
			    P->p[1] == 3))) {
		P->onesided = 1;
		P->V /= 2;
		P->chi /= 2;
	    }
	}
    } else if (!P->p[3]) { /* p q r| */
	P->M = P->N = 3;
	P->V = P->g;
	Malloc(P->n, P->N, double)
	Malloc(P->m, P->N, double)
	Malloc(P->rot, P->M, int)
	s = P->rot;
	for (j = 0; j < 3; j++) {
	    if (!(denominator(P->p[j]) % 2)) {
		/* what happens if there is more then one even denominator? */
		if (P->p[(j+1)%3] != P->p[(j+2)%3]) { /* needs postprocessing */
		    P->even = j;/* memorize the removed face */
		    P->chi -= P->g / numerator(P->p[j]) / 2;
		    P->onesided = 1;
		    P->D = 0;
		} else {/* for p = q we get a double 2 2r|p */
		/* noted by Roman Maeder <maeder@inf.ethz.ch> for 4 4 3/2| */
		/* Euler characteristic is still wrong */
		    P->D /= 2;
		}
		P->V /= 2;
	    }
	    P->n[j] = 2 * P->p[j];
	    P->m[j] = 1;
	    *s++ = j;
	}
    } else { /* |p q r - snub polyhedron */
	P->N = 4;
	P->M = 6;
	P->V = P->g / 2;/* Only "white" triangles carry a vertex */
	Malloc(P->n, P->N, double)
	Malloc(P->m, P->N, double)
	Malloc(P->rot, P->M, int)
	Malloc(P->snub, P->M, int)
	s = P->rot;
	t = P->snub;
	P->m[0] = P->n[0] = 3;
	for (j = 1; j < 4; j++) {
	    P->n[j] = P->p[j];
	    P->m[j] = 1;
	    *s++ = 0;
	    *s++ = j;
	    *t++ = 1;
	    *t++ = 0;
	}
    }
/*
 * Sort the fundamental triangles (using bubble sort) according to decreasing
 * n[i], while pushing the trivial triangles (n[i] = 2) to the end.
 */
    J = P->N - 1;
    while (J) {
	int last;
	last = J;
	J = 0;
	for (j = 0; j < last; j++) {
	    if ((P->n[j] < P->n[j+1] || P->n[j] == 2) && P->n[j+1] != 2) {
		int i;
		double temp;
		temp = P->n[j];
		P->n[j] = P->n[j+1];
		P->n[j+1] = temp;
		temp = P->m[j];
		P->m[j] = P->m[j+1];
		P->m[j+1] = temp;
		for (i = 0; i < P->M; i++) {
		    if (P->rot[i] == j)
			P->rot[i] = j+1;
		    else if (P->rot[i] == j+1)
			P->rot[i] = j;
		}
		if (P->even != -1) {
		    if (P->even == j)
			P->even = j+1;
		    else if (P->even == j+1)
			P->even = j;
		}
		J = j;
	    }
	}
    }
/*
 *  Get rid of repeated triangles.
 */
    for (J = 0; J < P->N && P->n[J] != 2;J++) {
	int k, i;
	for (j = J+1; j < P->N && P->n[j]==P->n[J]; j++)
	    P->m[J] += P->m[j];
	k = j - J - 1;
	if (k) {
	    for (i = j; i < P->N; i++) {
		P->n[i - k] = P->n[i];
		P->m[i - k] = P->m[i];
	    }
	    P->N -= k;
	    for (i = 0; i < P->M; i++) {
		if (P->rot[i] >= j)
		    P->rot[i] -= k;
		else if (P->rot[i] > J)
		    P->rot[i] = J;
	    }
	    if (P->even >= j)
		P->even -= k;
	}
    }
/*
 * Get rid of trivial triangles.
 */
    if (!J)
	J = 1; /* hosohedron */
    if (J < P->N) {
	int i;
	P->N = J;
	for (i = 0; i < P->M; i++) {
	    if (P->rot[i] >= P->N) {
		for (j = i + 1; j < P->M; j++) {
		    P->rot[j-1] = P->rot[j];
		    if (P->snub)
			P->snub[j-1] = P->snub[j];
		}
		P->M--;
	    }
	}
    }
/*
 * Truncate arrays
 */
    Realloc(P->n, P->N, double)
    Realloc(P->m, P->N, double)
    Realloc(P->rot, P->M, int)
    if (P->snub)
	Realloc(P->snub, P->M, int)
    return 1;
}


/*
 * Solve the fundamental right spherical triangles.
 * If need_approx is set, print iterations on standard error.
 */
int
newton(Polyhedron *P)
{
/*
 * First, we find initial approximations.
 */
    int j;
    double cosa;
    Malloc(P->gamma, P->N, double)
    if (P->N == 1) {
	P->gamma[0] = M_PI / P->m[0];
	return 1;
    }
    for (j = 0; j < P->N; j++)
	P->gamma[j] = M_PI / 2 - M_PI / P->n[j];
    errno = 0; /* may be non-zero from some reason */
/*
 * Next, iteratively find closer approximations for gamma[0] and compute
 * other gamma[j]'s from Napier's equations.
 */
    for (;;) {
	double delta = M_PI, sigma = 0;
	for (j = 0; j < P->N; j++) {
	    delta -= P->m[j] * P->gamma[j];
	}
	if (fabs(delta) < 11 * DBL_EPSILON)
	    return 1;
	/* On a RS/6000, fabs(delta)/DBL_EPSILON may occilate between 8 and
	 * 10. Reported by David W. Sanderson <dws@ssec.wisc.edu> */
	for (j = 0; j < P->N; j++)
	    sigma += P->m[j] * tan(P->gamma[j]);
	P->gamma[0] += delta * tan(P->gamma[0]) / sigma;
	if (P->gamma[0] < 0 || P->gamma[0] > M_PI)
	    Err("gamma out of bounds")
	cosa = cos(M_PI / P->n[0]) / sin(P->gamma[0]);
	for (j = 1; j < P->N; j++)
		P->gamma[j] = asin(cos(M_PI / P->n[j]) / cosa);
	if (errno)
	    Err(0)
    }
}

/*
 * Postprocess pqr| where r has an even denominator (cf. Coxeter &al. Sec.9).
 * Remove the {2r} and add a retrograde {2p} and retrograde {2q}.
 */
int
exceptions(Polyhedron *P)
{
    int j;
    if (P->even != -1) {
	P->M = P->N = 4;
	Realloc(P->n, P->N, double)
	Realloc(P->m, P->N, double)
	Realloc(P->gamma, P->N, double)
	Realloc(P->rot, P->M, int)
	for (j = P->even + 1; j < 3; j++) {
	    P->n[j-1] = P->n[j];
	    P->gamma[j-1] = P->gamma[j];
	}
	P->n[2] = compli(P->n[1]);
	P->gamma[2] = - P->gamma[1];
	P->n[3] = compli(P->n[0]);
	P->m[3] = 1;
	P->gamma[3] = - P->gamma[0];
	P->rot[0] = 0;
	P->rot[1] = 1;
	P->rot[2] = 3;
	P->rot[3] = 2;
    }

/*
 * Postprocess the last polyhedron |3/2 5/3 3 5/2 by taking a |5/3 3 5/2,
 * replacing the three snub triangles by four equatorial squares and adding
 * the missing {3/2} (retrograde triangle, cf. Coxeter &al. Sec. 11).
 */
    if (P->index == GR_DIRH_IDX) {
	P->N = 5;
	P->M = 8;
	Realloc(P->n, P->N, double)
	Realloc(P->m, P->N, double)
	Realloc(P->gamma, P->N, double)
	Realloc(P->rot, P->M, int)
	Realloc(P->snub, P->M, int)
	P->hemi = 1;
	P->D = 0;
	for (j = 3; j; j--) {
	    P->m[j] = 1;
	    P->n[j] = P->n[j-1];
	    P->gamma[j] = P->gamma[j-1];
	}
	P->m[0] = P->n[0] = 4;
	P->gamma[0] = M_PI / 2;
	P->m[4] = 1;
	P->n[4] = compli(P->n[1]);
	P->gamma[4] = - P->gamma[1];
	for (j = 1; j < 6; j += 2) P->rot[j]++;
	P->rot[6] = 0;
	P->rot[7] = 4;
	P->snub[6] = 1;
	P->snub[7] = 0;
    }
    return 1;
}

/*
 * Compute edge and face counts, and update D and chi.  Update D in the few
 * cases the density of the polyhedron is meaningful but different than the
 * density of the corresponding Schwarz triangle (cf. Coxeter &al., p. 418 and
 * p. 425).
 * In these cases, spherical faces of one type are concave (bigger than a
 * hemisphere), and the actual density is the number of these faces less the
 * computed density.  Note that if j != 0, the assignment gamma[j] = asin(...)
 * implies gamma[j] cannot be obtuse.  Also, compute chi for the only
 * non-Wythoffian polyhedron.
 */
int
count(Polyhedron *P)
{
    int j, temp;
    Malloc(P->Fi, P->N, int)
    for (j = 0; j < P->N; j++) {
	P->E += temp = P->V * numerator(P->m[j]);
	P->F += P->Fi[j] = temp / numerator(P->n[j]);
    }
    P->E /= 2;
    if (P->D && P->gamma[0] > M_PI / 2)
	P->D = P->Fi[0] - P->D;
    if (P->index == GR_DIRH_IDX)
	P->chi = P->V - P->E + P->F;
    return 1;
}

/*
 * Compute polyhedron vertices and vertex adjecency lists.
 * The vertices adjacent to v[i] are v[adj[0][i], v[adj[1][i], ...
 * v[adj[M-1][i], ordered counterclockwise.  The algorith is a BFS on the
 * vertices, in such a way that the vetices adjacent to a givem vertex are
 * obtained from its BFS parent by a cyclic sequence of rotations. firstrot[i]
 * points to the first  rotaion in the sequence when applied to v[i]. Note that
 * for non-snub polyhedra, the rotations at a child are opposite in sense when
 * compared to the rotations at the parent. Thus, we fill adj[*][i] from the
 * end to signify clockwise rotations. The firstrot[] array is not needed for
 * display thus it is freed after being used for face computations below.
 */
int
vertices(Polyhedron *P)
{
    int i, newV = 2;
    double cosa;
    Malloc(P->v, P->V, Vector)
    Matalloc(P->adj, P->M, P->V, int)
    Malloc(P->firstrot, P->V, int) /* temporary , put in Polyhedron
		    structure so that may be freed on error */
    cosa = cos(M_PI / P->n[0]) / sin(P->gamma[0]);
    P->v[0].x = 0;
    P->v[0].y = 0;
    P->v[0].z = 1;
    P->firstrot[0] = 0;
    P->adj[0][0] = 1;
    P->v[1].x = 2 * cosa * sqrt(1 - cosa * cosa);
    P->v[1].y = 0;
    P->v[1].z = 2 * cosa * cosa - 1;
    if (!P->snub) {
	P->firstrot[1] = 0;
	P->adj[0][1] = -1;/* start the other side */
	P->adj[P->M-1][1] = 0;
    } else {
	P->firstrot[1] = P->snub[P->M-1] ? 0 : P->M-1 ;
	P->adj[0][1] = 0;
    }
    for (i = 0; i < newV; i++) {
	int j, k;
	int last, one, start, limit;
	if (P->adj[0][i] == -1) {
	    one = -1; start = P->M-2; limit = -1;
	} else {
	    one = 1; start = 1; limit = P->M;
	}
	k = P->firstrot[i];
	for (j = start; j != limit; j += one) {
	    Vector temp;
	    int J;
	    temp = rotate (P->v[P->adj[j-one][i]], P->v[i],
		one * 2 * P->gamma[P->rot[k]]);
	    for (J=0; J<newV && !same(P->v[J],temp,BIG_EPSILON); J++)
		;/* noop */
	    P->adj[j][i] = J;
	    last = k;
	    if (++k == P->M)
		k = 0;
	    if (J == newV) { /* new vertex */
		if (newV == P->V) Err ("too many vertices")
		P->v[newV++] = temp;
		if (!P->snub) {
		    P->firstrot[J] = k;
		    if (one > 0) {
			P->adj[0][J] = -1;
			P->adj[P->M-1][J] = i;
		    } else {
			P->adj[0][J] = i;
		    }
		} else {
		    P->firstrot[J] = !P->snub[last] ? last :
			!P->snub[k] ? (k+1)%P->M : k ;
		    P->adj[0][J] = i;
		}
	    }
	}
    }
    return 1;
}

/*
 * Compute polyhedron faces (dual vertices) and incidence matrices.
 * For orientable polyhedra, we can distinguish between the two faces meeting
 * at a given directed edge and identify the face on the left and the face on
 * the right, as seen from the outside.  For one-sided polyhedra, the vertex
 * figure is a papillon (in Coxeter &al.  terminology, a crossed parallelogram)
 * and the two faces meeting at an edge can be identified as the side face
 * (n[1] or n[2]) and the diagonal face (n[0] or n[3]).
 */
int
faces(Polyhedron *P)
{
    int i, newF = 0;
    Malloc (P->f, P->F, Vector)
    Malloc (P->ftype, P->F, int)
    Matalloc (P->incid, P->M, P->V, int)
    P->minr = 1 / fabs (tan (M_PI / P->n[P->hemi]) * tan (P->gamma[P->hemi]));
    for (i = P->M; --i>=0;) {
	int j;
	for (j = P->V; --j>=0;)
	    P->incid[i][j] = -1;
    }
    for (i = 0; i < P->V; i++) {
	int j;
	for (j = 0; j < P->M; j++) {
	    int i0, J;
	    int pap=0;/* papillon edge type */
	    if (P->incid[j][i] != -1)
		continue;
	    P->incid[j][i] = newF;
	    if (newF == P->F)
		Err("too many faces")
	    P->f[newF] = pole(P->minr, P->v[i], P->v[P->adj[j][i]],
		    P->v[P->adj[mod(j + 1, P->M)][i]]);
	    P->ftype[newF] = P->rot[mod(P->firstrot[i] + (P->adj[0][i] <
			P->adj[P->M - 1][i] ?  j : -j - 2), P->M)];
	    if (P->onesided)
		pap = (P->firstrot[i] + j) % 2;
	    i0 = i;
	    J = j;
	    for (;;) {
		int k;
		k = i0;
		if ((i0 = P->adj[J][k]) == i) break;
		for (J = 0; J < P->M && P->adj[J][i0] != k; J++)
		    ;/* noop */
		if (J == P->M)
		    Err("too many faces")
		if (P->onesided && (J + P->firstrot[i0]) % 2 == pap) {
		    P->incid [J][i0] = newF;
		    if (++J >= P->M)
			J = 0;
		} else {
		    if (--J < 0)
			J = P->M - 1;
		    P->incid [J][i0] = newF;
		}
	    }
	newF++;
	}
    }
    Free(P->firstrot)
    Free(P->rot)
    Free(P->snub);
    return 1;
}

/*
 * compute the mathematical modulus function.
 */
int
mod (int i, int j)
{
    return (i%=j)>=0?i:j<0?i-j:i+j;
}

Fraction frax;

/*
 * Find the numerator and the denominator using the Euclidean algorithm.
 */
void
frac(double x)
{
    static Fraction zero = {0,1}, inf = {1,0};
    Fraction r0, r;
    long f;
    double s = x;
    r = zero;
    frax = inf;
    for (;;) {
	if (fabs(s) > (double) MAXLONG)
	    return;
	f = (long) floor (s);
	r0 = r;
	r = frax;
	frax.n = frax.n * f + r0.n;
	frax.d = frax.d * f + r0.d;
	if (x == (double)frax.n/(double)frax.d)
	    return;
	s = 1 / (s - f);
    }
}

/*
 * Print a fraction.
 */
void
printfrac(double x)
{
    frac (x);
    if (!frax.d) {
	    printf ("infinity");
	    return;
    }
    printf ("%ld",frax.n);
    if (frax.d!=1) printf ("/%ld",frax.d);
}

char *
sprintfrac(double x)
{
    char *s;
    frac (x);
    if (!frax.d) {
	Malloc(s, sizeof ("infinity"), char)
	strcpy(s, "infinity");
    } else if (frax.d == 1) {
	char n[MAXDIGITS + 1];
	sprintf(n, "%ld", frax.n);
	Malloc(s, strlen (n) +  1, char)
	strcpy(s, n);
    } else {
	char n[MAXDIGITS + 1], d[MAXDIGITS + 1];
	sprintf(n, "%ld", frax.n);
	sprintf(d, "%ld", frax.d);
	Malloc(s, strlen (n) + strlen (d) + 2, char)
	sprintf(s, "%s/%s", n, d);
    }
    return s;
}

double
dot(Vector a, Vector b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

Vector
scale(double k, Vector a)
{
    a.x *= k;
    a.y *= k;
    a.z *= k;
    return a;
}

Vector
diff(Vector a, Vector b)
{
    a.x -= b.x;
    a.y -= b.y;
    a.z -= b.z;
    return a;
}

Vector
cross(Vector a, Vector b)
{
    Vector p;
    p.x = a.y * b.z - a.z * b.y;
    p.y = a.z * b.x - a.x * b.z;
    p.z = a.x * b.y - a.y * b.x;
    return p;
}

Vector
sum(Vector a, Vector b)
{
    a.x += b.x;
    a.y += b.y;
    a.z += b.z;
    return a;
}

Vector
sum3(Vector a, Vector b, Vector c)
{
    a.x += b.x + c.x;
    a.y += b.y + c.y;
    a.z += b.z + c.z;
    return a;
}

Vector
rotate(Vector vertex, Vector axis, double angle)
{
    Vector p;
    p = scale(dot (axis, vertex), axis);
    return sum3(p, scale(cos(angle), diff(vertex, p)),
	    scale(sin(angle), cross(axis, vertex)));
}

Vector x, y, z;


int
same(Vector a, Vector b, double epsilon)
{
    return fabs(a.x - b.x) < epsilon && fabs(a.y - b.y) < epsilon
	    && fabs(a.z - b.z) < epsilon;
}

/*
 * Compute the polar reciprocal of the plane containing a, b and c:
 *
 * If this plane does not contain the origin, return p such that
 *	dot(p,a) = dot(p,b) = dot(p,b) = r.
 *
 * Otherwise, return p such that
 *	dot(p,a) = dot(p,b) = dot(p,c) = 0
 * and
 *	dot(p,p) = 1.
 */
Vector
pole(double r, Vector a, Vector b, Vector c)
{
    Vector p;
    double k;
    p = cross(diff(b, a), diff(c, a));
    k = dot(p, a);
    if (fabs(k) < 1e-6)
	return scale(1 / sqrt(dot(p, p)), p);
    else
	return scale(r/ k , p);
}

Uniform uniform_list[] = {

/*******************************************************************************
*		Tetrahedral Schwarz Triangles
*******************************************************************************/

/* (2 3 3) (T1) */

/*1*/	{"3|2 3","tetrahedron",
			"tetrahedron",15,1},
/*2*/	{"2 3|3","truncated tetrahedron",
			"triakistetrahedron",16,6},

/* (3/2 3 3) (T2) */

/*3*/	{"3/2 3|3","octahemioctahedron",
			"octahemioctacron",37,68},

/* (3/2 2 3) (T3) */

/*4*/	{"3/2 3|2","tetrahemihexahedron",
			"tetrahemihexacron",36,67},

/******************************************************************************
*		Octahedral Schwarz Triangles
******************************************************************************/

/* (2 3 4) (O1) */

/*5*/	{"4|2 3","octahedron",
			"cube",17,2},
/*6*/	{"3|2 4","cube",
			"octahedron",18,3},
/*7*/	{"2|3 4","cuboctahedron",
			"rhombic dodecahedron",19,11},
/*8*/	{"2 4|3","truncated octahedron",
			"tetrakishexahedron",20,7},
/*9*/	{"2 3|4","truncated cube",
			"triakisoctahedron",21,8},
/*10*/	{"3 4|2","rhombicuboctahedron",
			"deltoidal icositetrahedron",22,13},
/*11*/	{"2 3 4|","truncated cuboctahedron",
			"disdyakisdodecahedron",23,15},
/*12*/	{"|2 3 4","snub cube",
			"pentagonal icositetrahedron",24,17},

/* (3/2 4 4) (O2b) */

/*13*/	{"3/2 4|4","small cubicuboctahedron",
			"small hexacronic icositetrahedron",38,69},

/* (4/3 3 4) (O4) */

/*14*/	{"3 4|4/3","great cubicuboctahedron",
			"great hexacronic icositetrahedron",50,77},
/*15*/	{"4/3 4|3","cubohemioctahedron",
			"hexahemioctacron",51,78},
/*16*/	{"4/3 3 4|","cubitruncated cuboctahedron",
			"tetradyakishexahedron",52,79},

/* (3/2 2 4) (O5) */

/*17*/	{"3/2 4|2","great rhombicuboctahedron",
			"great deltoidal icositetrahedron",59,85},
/*18*/	{"3/2 2 4|","small rhombihexahedron",
			"small rhombihexacron",60,86},

/* (4/3 2 3) (O7) */

/*19*/	{"2 3|4/3","stellated truncated hexahedron",
			"great triakisoctahedron",66,92},
/*20*/	{"4/3 2 3|","great truncated cuboctahedron",
			"great disdyakisdodecahedron",67,93},

/* (4/3 3/2 2) (O11) */

/*21*/	{"4/3 3/2 2|","great rhombihexahedron",
			"great rhombihexacron",82,103},

/******************************************************************************
*		Icosahedral Schwarz Triangles
******************************************************************************/

/* (2 3 5) (I1) */

/*22*/	{"5|2 3","icosahedron",
			"dodecahedron",25,4},
/*23*/	{"3|2 5","dodecahedron",
			"icosahedron",26,5},
/*24*/	{"2|3 5","icosidodecahedron",
			"rhombic triacontahedron",28,12},
/*25*/	{"2 5|3","truncated icosahedron",
			"pentakisdodecahedron",27,9},
/*26*/	{"2 3|5","truncated dodecahedron",
			"triakisicosahedron",29,10},
/*27*/	{"3 5|2","rhombicosidodecahedron",
			"deltoidal hexecontahedron",30,14},
/*28*/	{"2 3 5|","truncated icosidodecahedron",
			"disdyakistriacontahedron",31,16},
/*29*/	{"|2 3 5","snub dodecahedron",
			"pentagonal hexecontahedron",32,18},

/* (5/2 3 3) (I2a) */

/*30*/	{"3|5/2 3","small ditrigonal icosidodecahedron",
			"small triambic icosahedron",39,70},
/*31*/	{"5/2 3|3","small icosicosidodecahedron",
			"small icosacronic hexecontahedron",40,71},
/*32*/	{"|5/2 3 3","small snub icosicosidodecahedron",
			"small hexagonal hexecontahedron",41,110},

/* (3/2 5 5) (I2b) */

/*33*/	{"3/2 5|5","small dodecicosidodecahedron",
			"small dodecacronic hexecontahedron",42,72},

/* (2 5/2 5) (I3) */

/*34*/	{"5|2 5/2","small stellated dodecahedron",
			"great dodecahedron",43,20},
/*35*/	{"5/2|2 5","great dodecahedron",
			"small stellated dodecahedron",44,21},
/*36*/	{"2|5/2 5","great dodecadodecahedron",
			"medial rhombic triacontahedron",45,73},
/*37*/	{"2 5/2|5","truncated great dodecahedron",
			"small stellapentakisdodecahedron",47,75},
/*38*/	{"5/2 5|2","rhombidodecadodecahedron",
			"medial deltoidal hexecontahedron",48,76},
/*39*/	{"2 5/2 5|","small rhombidodecahedron",
			"small rhombidodecacron",46,74},
/*40*/	{"|2 5/2 5","snub dodecadodecahedron",
			"medial pentagonal hexecontahedron",49,111},

/* (5/3 3 5) (I4) */

/*41*/	{"3|5/3 5","ditrigonal dodecadodecahedron",
			"medial triambic icosahedron",53,80},
/*42*/	{"3 5|5/3","great ditrigonal dodecicosidodecahedron",
			"great ditrigonal dodecacronic hexecontahedron",54,81},
/*43*/	{"5/3 3|5","small ditrigonal dodecicosidodecahedron",
			"small ditrigonal dodecacronic hexecontahedron",55,82},
/*44*/	{"5/3 5|3","icosidodecadodecahedron",
			"medial icosacronic hexecontahedron",56,83},
/*45*/	{"5/3 3 5|","icositruncated dodecadodecahedron",
			"tridyakisicosahedron",57,84},
/*46*/	{"|5/3 3 5","snub icosidodecadodecahedron",
			"medial hexagonal hexecontahedron",58,112},

/* (3/2 3 5) (I6b) */

/*47*/	{"3/2|3 5","great ditrigonal icosidodecahedron",
			"great triambic icosahedron",61,87},
/*48*/	{"3/2 5|3","great icosicosidodecahedron",
			"great icosacronic hexecontahedron",62,88},
/*49*/	{"3/2 3|5","small icosihemidodecahedron",
			"small icosihemidodecacron",63,89},
/*50*/	{"3/2 3 5|","small dodecicosahedron",
			"small dodecicosacron",64,90},

/* (5/4 5 5) (I6c) */

/*51*/	{"5/4 5|5","small dodecahemidodecahedron",
			"small dodecahemidodecacron",65,91},

/* (2 5/2 3) (I7) */

/*52*/	{"3|2 5/2","great stellated dodecahedron",
			"great icosahedron",68,22},
/*53*/	{"5/2|2 3","great icosahedron",
			"great stellated dodecahedron",69,41},
/*54*/	{"2|5/2 3","great icosidodecahedron",
			"great rhombic triacontahedron",70,94},
/*55*/	{"2 5/2|3","great truncated icosahedron",
			"great stellapentakisdodecahedron",71,95},
/*56*/	{"2 5/2 3|","rhombicosahedron",
			"rhombicosacron",72,96},
/*57*/	{"|2 5/2 3","great snub icosidodecahedron",
			"great pentagonal hexecontahedron",73,113},

/* (5/3 2 5) (I9) */

/*58*/	{"2 5|5/3","small stellated truncated dodecahedron",
			"great pentakisdodekahedron",74,97},
/*59*/	{"5/3 2 5|","truncated dodecadodecahedron",
			"medial disdyakistriacontahedron",75,98},
/*60*/	{"|5/3 2 5","inverted snub dodecadodecahedron",
			"medial inverted pentagonal hexecontahedron",76,114},

/* (5/3 5/2 3) (I10a) */

/*61*/	{"5/2 3|5/3","great dodecicosidodecahedron",
			"great dodecacronic hexecontahedron",77,99},
/*62*/	{"5/3 5/2|3","small dodecahemicosahedron",
			"small dodecahemicosacron",78,100},
/*63*/	{"5/3 5/2 3|","great dodecicosahedron",
			"great dodecicosacron",79,101},
/*64*/	{"|5/3 5/2 3","great snub dodecicosidodecahedron",
			"great hexagonal hexecontahedron",80,115},

/* (5/4 3 5) (I10b) */

/*65*/	{"5/4 5|3","great dodecahemicosahedron",
			"great dodecahemicosacron",81,102},

/* (5/3 2 3) (I13) */

/*66*/	{"2 3|5/3","great stellated truncated dodecahedron",
			"great triakisicosahedron",83,104},
/*67*/	{"5/3 3|2","great rhombicosidodecahedron",
			"great deltoidal hexecontahedron",84,105},
/*68*/	{"5/3 2 3|","great truncated icosidodecahedron",
			"great disdyakistriacontahedron",87,108},
/*69*/	{"|5/3 2 3","great inverted snub icosidodecahedron",
			"great inverted pentagonal hexecontahedron",88,116},

/* (5/3 5/3 5/2) (I18a) */

/*70*/	{"5/3 5/2|5/3","great dodecahemidodecahedron",
			"great dodecahemidodecacron",86,107},

/* (3/2 5/3 3) (I18b) */

/*71*/	{"3/2 3|5/3","great icosihemidodecahedron",
			"great icosihemidodecacron",85,106},

/* (3/2 3/2 5/3) (I22) */

/*72*/	{"|3/2 3/2 5/2","small retrosnub icosicosidodecahedron",
			"small hexagrammic hexecontahedron",91,118},

/* (3/2 5/3 2) (I23) */

/*73*/	{"3/2 5/3 2|","great rhombidodecahedron",
			"great rhombidodecacron",89,109},
/*74*/	{"|3/2 5/3 2","great retrosnub icosidodecahedron",
			"great pentagrammic hexecontahedron",90,117},

/******************************************************************************
*		Last But Not Least
******************************************************************************/

/*75*/	{"3/2 5/3 3 5/2","great dirhombicosidodecahedron",
			"great dirhombicosidodecacron",92,119},

/******************************************************************************
*		Dihedral Schwarz Triangles (D5 only)
******************************************************************************/

/* (2 2 5) (D1/5) */

/*76*/	{"2 5|2","pentagonal prism",
			"pentagonal dipyramid",0,0},
/*77*/	{"|2 2 5","pentagonal antiprism",
			"pentagonal deltohedron",0,0},

/* (2 2 5/2) (D2/5) */

/*78*/	{"2 5/2|2","pentagrammic prism",
			"pentagrammic dipyramid",0,0},
/*79*/	{"|2 2 5/2","pentagrammic antiprism",
			"pentagrammic deltohedron",0,0},

/* (5/3 2 2) (D3/5) */

/*80*/	{"|2 2 5/3","pentagrammic crossed antiprism",
			"pentagrammic concave deltohedron",0,0}

};
 
int get_uniform_list(Uniform **uniform)
{
   *uniform = uniform_list;
   return sizeof (uniform_list) / sizeof (uniform_list[0]);
}



