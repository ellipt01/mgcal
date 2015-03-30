/*
 * vector.c
 *
 *  Created on: 2015/03/14
 *      Author: utsugi
 */

#include <stdlib.h>
#include <math.h>

#include "cvector.h"
#include "util.h"

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

cvector *
cvector_copy (const cvector *src)
{
	cvector	*dist = cvector_alloc ();
	if (!src) error_and_exit ("cvector_copy", "cvector *src is empty.", __FILE__, __LINE__);
	cvector_set (dist, src->x, src->y, src->z);
	return dist;
}

void
cvector_add (cvector *v1, const cvector *v2)
{
	if (!v1) error_and_exit ("cvector_add", "cvector *v1 is empty.", __FILE__, __LINE__);
	if (!v2) error_and_exit ("cvector_add", "cvector *v2 is empty.", __FILE__, __LINE__);
	v1->x += v2->x;
	v1->y += v2->y;
	v1->z += v2->z;
	return;
}

void
cvector_sub (cvector *v1, const cvector *v2)
{
	if (!v1) error_and_exit ("cvector_sub", "cvector *v1 is empty.", __FILE__, __LINE__);
	if (!v2) error_and_exit ("cvector_sub", "cvector *v2 is empty.", __FILE__, __LINE__);
	v1->x -= v2->x;
	v1->y -= v2->y;
	v1->z -= v2->z;
	return;
}
