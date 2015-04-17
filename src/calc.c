/*
 * calc.c
 *
 *  Created on: 2015/03/14
 *      Author: utsugi
 */

#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <float.h>

#include "cvector.h"
#include "source.h"
#include "calc.h"
#include "private/util.h"

double scale_factor = 1.;

static double
sign (const double val)
{
	return (val >= 0.) ? +1. : -1.;
}

static void
dipole_kernel (cvector *f, const double x, const double y, const double z, const cvector *mgz)
{
	double	fx, fy, fz;

	double	r  = sqrt (x * x + y * y + z * z);
	double	r3 = pow (r, 3.);
	double	r5 = pow (r, 5.);

	double	jx = mgz->x;
	double	jy = mgz->y;
	double	jz = mgz->z;

	fx =
		- jx * (1.0 / r3 - 3. * pow (x, 2.) / r5)
		+ jy * (3.0 * x * y / r5)
		+ jz * (3.0 * x * z / r5);

	fy =
		+ jx * (3.0 * y * x / r5)
		- jy * (1.0 / r3 - 3. * pow (y, 2.) / r5)
		+ jz * (3.0 * y * z / r5);

	fz =
		+ jx * (3.0 * z * x / r5)
		+ jy * (3.0 * z * y / r5)
		- jz * (1.0 / r3 - 3. * pow (z, 2.) / r5);

	cvector_set (f, fx, fy, fz);
	return;
}

static void
prism_kernel (cvector *f, const double x, const double y, const double z, const cvector *mgz)
{
	double	fx, fy, fz;
	double	lnx, lny, lnz;

	double	r = sqrt (x * x + y * y + z * z);

	if (fabs (y) <= DBL_EPSILON && fabs (z) <= DBL_EPSILON && x < 0.)
		lnx = - log (r - x);
	else
		lnx = log (r + x);

	if (fabs (x) <= DBL_EPSILON && fabs (z) <= DBL_EPSILON && y < 0.)
		lny = - log (r - y);
	else
		lny = log (r + y);

	if (fabs (x) <= DBL_EPSILON && fabs (y) <= DBL_EPSILON && z < 0.)
		lnz = - log (r - z);
	else
		lnz = log (r + z);

	{
		double	jx = mgz->x;
		double	jy = mgz->y;
		double	jz = mgz->z;

		fx = - jx * atan2 (y * z, x * r)
			+ jy * lnz
			+ jz * lny;

		fy = jx * lnz
			- jy * atan2 (x * z, y * r)
			+ jz * lnx;

		fz = jx * lny
			+ jy * lnx
			- jz * atan2 (x * y, z * r);
	}

	cvector_set (f, fx, fy, fz);
	return;
}

cvector *
dipole (const cvector *obs, const source *s)
{
	double		x, y, z;
	double		x0, y0, z0;
	cvector		*f;
	cvector		tmp;
	source_item	*cur;

	if (!obs) error_and_exit ("dipole", "cvector *obs is empty.", __FILE__, __LINE__);
	if (!s) error_and_exit ("dipole", "source *s is empty.", __FILE__, __LINE__);

	x0 = obs->x;
	y0 = obs->y;
	z0 = obs->z;

	f = cvector_new (0., 0., 0.);

	cur = s->begin;
	while (cur) {
		double	dv;

		if (!cur->pos) error_and_exit ("dipole", "position of source item is empty.", __FILE__, __LINE__);
		if (!cur->mgz) error_and_exit ("dipole", "magnetization of source item is empty.", __FILE__, __LINE__);

		dv = 1.0;
		if (cur->dim) {
			double	dx = cur->dim->x;
			double	dy = cur->dim->y;
			double	dz = cur->dim->z;
			dv = fabs (dx * dy * dz);
		}
		x = cur->pos->x;
		y = cur->pos->y;
		z = cur->pos->z;
		dipole_kernel (&tmp, x - x0, y - y0, z - z0, cur->mgz);
		cvector_axpy (dv, &tmp, f);	// f = f + dv * tmp
		cur = cur->next;
	}
	cvector_scale (f, scale_factor);
	return f;
}

cvector *
prism (const cvector *obs, const source *s)
{
	double		a[2], b[2], c[2];
	double		x, y, z;
	double		x0, y0, z0;
	cvector		*f;
	cvector		tmp[8];
	source_item	*cur;

	if (!obs) error_and_exit ("prism", "cvector *obs is empty.", __FILE__, __LINE__);
	if (!s) error_and_exit ("prism", "source *s is empty.", __FILE__, __LINE__);

	x0 = obs->x;
	y0 = obs->y;
	z0 = obs->z;

	f = cvector_new (0., 0., 0.);

	cur = s->begin;
	while (cur) {
		double	dx, dy, dz, dv;
		double	flag;

		if (!cur->pos) error_and_exit ("prism", "position of source item is empty.", __FILE__, __LINE__);
		if (!cur->dim) error_and_exit ("prism", "dimension of source item is empty.", __FILE__, __LINE__);
		if (!cur->mgz) error_and_exit ("prism", "magnetization of source item is empty.", __FILE__, __LINE__);

		dx = cur->dim->x;
		dy = cur->dim->y;
		dz = cur->dim->z;
		dv = dx * dy * dz;
		flag = sign (dv);

		x = cur->pos->x;
		y = cur->pos->y;
		z = cur->pos->z;

		a[0] = x - 0.5 * dx - x0;
		b[0] = y - 0.5 * dy - y0;
		c[0] = z - 0.5 * dz - z0;

		a[1] = a[0] + dx;
		b[1] = b[0] + dy;
		c[1] = c[0] + dz;

		prism_kernel (&tmp[0], a[1], b[1], c[1], cur->mgz);
		prism_kernel (&tmp[1], a[1], b[1], c[0], cur->mgz);
		prism_kernel (&tmp[2], a[1], b[0], c[1], cur->mgz);
		prism_kernel (&tmp[3], a[1], b[0], c[0], cur->mgz);
		prism_kernel (&tmp[4], a[0], b[1], c[1], cur->mgz);
		prism_kernel (&tmp[5], a[0], b[1], c[0], cur->mgz);
		prism_kernel (&tmp[6], a[0], b[0], c[1], cur->mgz);
		prism_kernel (&tmp[7], a[0], b[0], c[0], cur->mgz);

		f->x += flag * (tmp[0].x - tmp[1].x - tmp[2].x + tmp[3].x
			- tmp[4].x + tmp[5].x + tmp[6].x - tmp[7].x);

		f->y += flag * (tmp[0].y - tmp[1].y - tmp[2].y + tmp[3].y
			- tmp[4].y + tmp[5].y + tmp[6].y - tmp[7].y);

		f->z += flag * (tmp[0].z - tmp[1].z - tmp[2].z + tmp[3].z
			- tmp[4].z + tmp[5].z + tmp[6].z - tmp[7].z);

		cur = cur->next;
	}
	cvector_scale (f, scale_factor);
	return f;
}

static void
dipole_yz_kernel (cvector *f, const double y, const double z, const cvector *mgz)
{
	double	fy, fz;

	double	jy = mgz->y;
	double	jz = mgz->z;

	double	r2 = y * y + z * z;
	double	r4 = pow (r2, 2.);

	fy = jy * (1. / r2 - 2. * pow (y, 2.) / r4) - jz * 2. * y * z / r4;
	fz = - jy * 2. * y * z / r4 + jz * (1. / r2 - 2. * pow (z, 2.) / r4);

	cvector_set (f, 0., - 2. * fy, - 2. * fz);
	return;
}

static void
prism_yz_kernel (cvector *f, const double y, const double z, const cvector *mgz)
{
	double	fy, fz;

	double	jy = mgz->y;
	double	jz = mgz->z;

	double	r = sqrt (y * y + z * z);

	fy = jy * atan (z / y) + jz * log (r);
	fz = jy * log (r) + jz * atan (y / z);

	cvector_set (f, 0., - 2. * fy, - 2. * fz);
	return;
}

cvector *
dipole_yz (const cvector *obs, const source *s)
{
	double		y, z;
	double		y0, z0;
	cvector		*f;
	cvector		tmp;
	source_item	*cur;

	y0 = obs->y;
	z0 = obs->z;

	f = cvector_new (0., 0., 0.);

	cur = s->begin;
	while (cur) {
		double	ds;

		if (!cur->pos) error_and_exit ("dipole_yz", "position of source item is empty.", __FILE__, __LINE__);
		if (!cur->mgz) error_and_exit ("dipole_yz", "magnetization of source item is empty.", __FILE__, __LINE__);

		ds = 1.0;
		if (cur->dim) {
			double	dy = cur->dim->y;
			double	dz = cur->dim->z;
			ds = fabs (dy * dz);
		}
		y = cur->pos->y - y0;
		z = cur->pos->z - z0;

		dipole_yz_kernel (&tmp, y, z, cur->mgz);
		cvector_axpy (ds, &tmp, f);
		cur = cur->next;
	}
	cvector_scale (f, scale_factor);
	return f;
}

cvector *
prism_yz (const cvector *obs, const source *s)
{
	double		b[2], c[2];
	double		y, z;
	double		y0, z0;
	cvector		*f;
	cvector		tmp[4];
	source_item	*cur;

	y0 = obs->y;
	z0 = obs->z;

	f = cvector_new (0., 0., 0.);

	cur = s->begin;
	while (cur) {
		double	dy, dz, ds;
		double	flag;

		if (!cur->pos) error_and_exit ("prism_yz", "position of source item is empty.", __FILE__, __LINE__);
		if (!cur->dim) error_and_exit ("prism_yz", "dimension of source item is empty.", __FILE__, __LINE__);
		if (!cur->mgz) error_and_exit ("prism_yz", "magnetization of source item is empty.", __FILE__, __LINE__);

		dy = cur->dim->y;
		dz = cur->dim->z;
		ds = dy * dz;
		flag = sign (ds);

		y = cur->pos->y;
		z = cur->pos->z;

		b[0] = y - 0.5 * dy - y0;
		c[0] = z - 0.5 * dz - z0;

		b[1] = b[0] + dy;
		c[1] = c[0] + dz;

		prism_yz_kernel (&tmp[0], b[1], c[1], cur->mgz);
		prism_yz_kernel (&tmp[1], b[1], c[0], cur->mgz);
		prism_yz_kernel (&tmp[2], b[0], c[1], cur->mgz);
		prism_yz_kernel (&tmp[3], b[0], c[0], cur->mgz);

		f->x += flag * (tmp[0].x - tmp[1].x - tmp[2].x + tmp[3].x);

		f->y += flag * (tmp[0].y - tmp[1].y - tmp[2].y + tmp[3].y);

		f->z += flag * (tmp[0].z - tmp[1].z - tmp[2].z + tmp[3].z);

		cur = cur->next;
	}
	cvector_scale (f, scale_factor);
	return f;
}

typedef enum {
	MGCAL_X_COMPONENT,	// Hx
	MGCAL_Y_COMPONENT,	// Hy
	MGCAL_Z_COMPONENT,	// Hz
	MGCAL_TOTAL_FORCE	// f
} MgcalComponent;

double
total_force (const cvector *exf, const cvector *f)
{
	if (exf == NULL) error_and_exit ("total_force", "cvector *exf is empty.", __FILE__, __LINE__);
	if (f == NULL) error_and_exit ("total_force", "cvector *f is empty.", __FILE__, __LINE__);
	return exf->x * f->x + exf->y * f->y + exf->z * f->z;
}

static double
calc_component (const cvector *f, const source *src, MgcalComponent comp)
{
	double	val = 0.;
	switch (comp) {
	case MGCAL_X_COMPONENT:
		val = f->x;
		break;

	case MGCAL_Y_COMPONENT:
		val = f->y;
		break;

	case MGCAL_Z_COMPONENT:
		val = f->z;
		break;

	case MGCAL_TOTAL_FORCE:
		val = total_force (src->exf, f);
		break;

	}
	return val;
}

/*** dipole ***/
static double
component_dipole (const cvector *obs, const source *src, MgcalComponent comp)
{
	cvector	*f = dipole (obs, src);
	double	val = calc_component (f, src, comp);
	cvector_free (f);
	return val;
}

double
x_component_dipole (const cvector *obs, const source *src, void *data)
{
	return component_dipole (obs, src, MGCAL_X_COMPONENT);
}

double
y_component_dipole (const cvector *obs, const source *src, void *data)
{
	return component_dipole (obs, src, MGCAL_Y_COMPONENT);
}

double
z_component_dipole (const cvector *obs, const source *src, void *data)
{
	return component_dipole (obs, src, MGCAL_Z_COMPONENT);
}

double
total_force_dipole (const cvector *obs, const source *src, void *data)
{
	return component_dipole (obs, src, MGCAL_TOTAL_FORCE);
}

/*** prism ***/
static double
component_prism (const cvector *obs, const source *src, MgcalComponent comp)
{
	cvector	*f = prism (obs, src);
	double	val = calc_component (f, src, comp);
	cvector_free (f);
	return val;
}

double
x_component_prism (const cvector *obs, const source *src, void *data)
{
	return component_prism (obs, src, MGCAL_X_COMPONENT);
}

double
y_component_prism (const cvector *obs, const source *src, void *data)
{
	return component_prism (obs, src, MGCAL_Y_COMPONENT);
}

double
z_component_prism (const cvector *obs, const source *src, void *data)
{
	return component_prism (obs, src, MGCAL_Z_COMPONENT);
}

double
total_force_prism (const cvector *obs, const source *src, void *data)
{
	return component_prism (obs, src, MGCAL_TOTAL_FORCE);
}

/*** dipole yz ***/
static double
component_dipole_yz (const cvector *obs, const source *src, MgcalComponent comp)
{
	cvector	*f = dipole_yz (obs, src);
	double	val = calc_component (f, src, comp);
	cvector_free (f);
	return val;
}

double
y_component_dipole_yz (const cvector *obs, const source *src, void *data)
{
	return component_dipole_yz (obs, src, MGCAL_Y_COMPONENT);
}

double
z_component_dipole_yz (const cvector *obs, const source *src, void *data)
{
	return component_dipole_yz (obs, src, MGCAL_Z_COMPONENT);
}

double
total_force_dipole_yz (const cvector *obs, const source *src, void *data)
{
	return component_dipole_yz (obs, src, MGCAL_TOTAL_FORCE);
}

/*** prism yz ***/
static double
component_prism_yz (const cvector *obs, const source *src, MgcalComponent comp)
{
	cvector	*f = prism_yz (obs, src);
	double	val = calc_component (f, src, comp);
	cvector_free (f);
	return val;
}

double
y_component_prism_yz (const cvector *obs, const source *src, void *data)
{
	return component_prism_yz (obs, src, MGCAL_Y_COMPONENT);
}

double
z_component_prism_yz (const cvector *obs, const source *src, void *data)
{
	return component_prism_yz (obs, src, MGCAL_Z_COMPONENT);
}

double
total_force_prism_yz (const cvector *obs, const source *src, void *data)
{
	return component_prism_yz (obs, src, MGCAL_TOTAL_FORCE);
}
