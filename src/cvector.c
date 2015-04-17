/*
 * vector.c
 *
 *  Created on: 2015/03/14
 *      Author: utsugi
 */

#include <stdlib.h>
#include <math.h>

#include "cvector.h"
#include "private/util.h"

static cvector *
cvector_alloc (void)
{
	cvector	*cv = (cvector *) malloc (sizeof (cvector));
	cv->x = 0.;
	cv->y = 0.;
	cv->z = 0.;
	return cv;
}

cvector *
cvector_new (const double x, const double y, const double z)
{
	cvector	*cv = cvector_alloc ();
	cv->x = x;
	cv->y = y;
	cv->z = z;
	return cv;
}

cvector *
cvector_new_with_geodesic_poler (const double r, const double inc, const double dec)
{
	cvector	*cv = cvector_alloc ();
	double	t = deg2rad(inc);
	double	p = deg2rad (dec);
	cv->x = r * cos (t) * sin (p);
	cv->y = r * cos (t) * cos (p);
	cv->z = - r * sin (t);
	return cv;
}

void
cvector_set (cvector *cv, const double x, const double y, const double z)
{
	if (!cv) error_and_exit ("cvector_set", "cvector *cv is empty.", __FILE__, __LINE__);
	cv->x = x;
	cv->y = y;
	cv->z = z;
	return;
}

void
cvector_free (cvector *cv)
{
	if (cv) free (cv);
	return;
}

void
cvector_scale (cvector *x, const double alpha)
{
	x->x *= alpha;
	x->y *= alpha;
	x->z *= alpha;
	return;
}

cvector *
cvector_copy (const cvector *src)
{
	cvector	*dist = cvector_alloc ();
	if (!src) error_and_exit ("cvector_copy", "cvector *src is empty.", __FILE__, __LINE__);
	cvector_set (dist, src->x, src->y, src->z);
	return dist;
}

void
cvector_axpy (const double alpha, const cvector *x, cvector *y)
{
	if (!x) error_and_exit ("cvector_axpy", "cvector *x is empty.", __FILE__, __LINE__);
	if (!y) error_and_exit ("cvector_axpy", "cvector *y is empty.", __FILE__, __LINE__);
	y->x += alpha * x->x;
	y->y += alpha * x->y;
	y->z += alpha * x->z;
	return;
}
