/*
 * vector.c
 *
 *  Created on: 2015/03/14
 *      Author: utsugi
 */

#include <stdlib.h>
#include <math.h>

#include "vector3d.h"
#include "private/util.h"

static vector3d *
vector3d_alloc (void)
{
	vector3d	*v = (vector3d *) malloc (sizeof (vector3d));
	v->x = 0.;
	v->y = 0.;
	v->z = 0.;
	return v;
}

vector3d *
vector3d_new (const double x, const double y, const double z)
{
	vector3d	*v = vector3d_alloc ();
	v->x = x;
	v->y = y;
	v->z = z;
	return v;
}

vector3d *
vector3d_new_with_geodesic_poler (const double r, const double inc, const double dec)
{
	vector3d	*v = vector3d_alloc ();
	double	t = deg2rad(inc);
	double	p = deg2rad (dec);
	v->x = r * cos (t) * cos (p);
	v->y = r * cos (t) * sin (p);
	v->z = r * sin (t);
	return v;
}

void
vector3d_set (vector3d *v, const double x, const double y, const double z)
{
	if (!v) error_and_exit ("vector3d_set", "vector3d *v is empty.", __FILE__, __LINE__);
	v->x = x;
	v->y = y;
	v->z = z;
	return;
}

void
vector3d_free (vector3d *v)
{
	if (v) free (v);
	return;
}

void
vector3d_scale (vector3d *x, const double alpha)
{
	x->x *= alpha;
	x->y *= alpha;
	x->z *= alpha;
	return;
}

vector3d *
vector3d_copy (const vector3d *src)
{
	vector3d	*dist = vector3d_alloc ();
	if (!src) error_and_exit ("vector3d_copy", "vector3d *src is empty.", __FILE__, __LINE__);
	vector3d_set (dist, src->x, src->y, src->z);
	return dist;
}

void
vector3d_axpy (const double alpha, const vector3d *x, vector3d *y)
{
	if (!x) error_and_exit ("vector3d_axpy", "vector3d *x is empty.", __FILE__, __LINE__);
	if (!y) error_and_exit ("vector3d_axpy", "vector3d *y is empty.", __FILE__, __LINE__);
	y->x += alpha * x->x;
	y->y += alpha * x->y;
	y->z += alpha * x->z;
	return;
}

double
vector3d_dot (const vector3d *x, const vector3d *y)
{
	return x->x * y->x + x->y * y->y + x->z * y->z;
}

double
vector3d_nrm (const vector3d *c)
{
	return sqrt (vector3d_dot (c, c));
}
