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

#include "../include/vector3d.h"
#include "private/util.h"
#include "source.h"
#include "calc.h"

double scale_factor = 1.;

#define SIGN(a) ((a) < 0. ? -1. : +1.)

static void
transform (const double x1, const double y1, const double z1, double *x, double *y, double *z)
{
	if (x) *x = y1;
	if (y) *y = x1;
	if (z) *z = -z1;
	return;
}

static void
dipole_kernel (vector3d *f, const double x1, const double y1, const double z1, const vector3d *mgz)
{
	double	x, y, z;
	double	jx, jy, jz;

	double	fx, fy, fz;
	double	r, r3, r5;

	transform (x1, y1, z1, &x, &y, &z);
	transform (mgz->x, mgz->y, mgz->z, &jx, &jy, &jz);

	r  = sqrt (x * x + y * y + z * z);
	r3 = pow (r, 3.);
	r5 = pow (r, 5.);

	fx =
		- jx * (1. / r3 - 3. * pow (x, 2.) / r5)
		+ jy * (3. * x * y / r5)
		+ jz * (3. * x * z / r5);

	fy =
		+ jx * (3. * y * x / r5)
		- jy * (1. / r3 - 3. * pow (y, 2.) / r5)
		+ jz * (3. * y * z / r5);

	fz =
		+ jx * (3. * z * x / r5)
		+ jy * (3. * z * y / r5)
		- jz * (1. / r3 - 3. * pow (z, 2.) / r5);

	vector3d_set (f, fy, fx, - fz);
	return;
}

static void
prism_kernel (vector3d *f, const double x1, const double y1, const double z1, const vector3d *mgz)
{
	double	x, y, z;
	double	jx, jy, jz;

	double	fx, fy, fz;
	double	r, lnx, lny, lnz;

	transform (x1, y1, z1, &x, &y, &z);
	transform (mgz->x, mgz->y, mgz->z, &jx, &jy, &jz);

	r = sqrt (x * x + y * y + z * z);
	lnx = (fabs (r + x) > DBL_EPSILON) ? log (r + x) : - log (r - x);
	lny = (fabs (r + y) > DBL_EPSILON) ? log (r + y) : - log (r - y);
	lnz = (fabs (r + z) > DBL_EPSILON) ? log (r + z) : - log (r - z);

	fx = - jx * atan2 (y * z, x * r)
		+ jy * lnz
		+ jz * lny;

	fy = jx * lnz
		- jy * atan2 (x * z, y * r)
		+ jz * lnx;

	fz = jx * lny
		+ jy * lnx
		- jz * atan2 (x * y, z * r);

	vector3d_set (f, fy, fx, - fz);
	return;
}

vector3d *
dipole (const vector3d *obs, const source *s)
{
	double		x, y, z;
	double		x0, y0, z0;
	vector3d	*f;
	vector3d	tmp;
	source_item	*cur;

	if (!obs) error_and_exit ("dipole", "vector3d *obs is empty.", __FILE__, __LINE__);
	if (!s) error_and_exit ("dipole", "source *s is empty.", __FILE__, __LINE__);

	x0 = obs->x;
	y0 = obs->y;
	z0 = obs->z;

	f = vector3d_new (0., 0., 0.);

	cur = s->begin;
	while (cur) {
		double	dv;

		if (!cur->pos) error_and_exit ("dipole", "position of source item is empty.", __FILE__, __LINE__);
		if (!cur->mgz) error_and_exit ("dipole", "magnetization of source item is empty.", __FILE__, __LINE__);

		dv = 1.;
		if (cur->dim) {
			double	dx = cur->dim->x;
			double	dy = cur->dim->y;
			double	dz = cur->dim->z;
			double	dd = dx * dy;
			if (fabs (dz) > DBL_EPSILON) dd *= dz;
			dv = fabs (dd);
		}
		x = cur->pos->x;
		y = cur->pos->y;
		z = cur->pos->z;
		dipole_kernel (&tmp, x - x0, y - y0, z - z0, cur->mgz);
		vector3d_axpy (dv, &tmp, f);	// f = f + dv * tmp
		cur = cur->next;
	}
	vector3d_scale (f, scale_factor);
	return f;
}

vector3d *
prism (const vector3d *obs, const source *s)
{
	double		a[2], b[2], c[2];
	double		x, y, z;
	double		x0, y0, z0;
	vector3d	*f;
	vector3d	tmp[8];
	source_item	*cur;

	if (!obs) error_and_exit ("prism", "vector3d *obs is empty.", __FILE__, __LINE__);
	if (!s) error_and_exit ("prism", "source *s is empty.", __FILE__, __LINE__);

	x0 = obs->x;
	y0 = obs->y;
	z0 = obs->z;

	f = vector3d_new (0., 0., 0.);

	cur = s->begin;
	while (cur) {
		double	dx, dy, dz;
		double	flag;

		if (!cur->pos) error_and_exit ("prism", "position of source item is empty.", __FILE__, __LINE__);
		if (!cur->dim) error_and_exit ("prism", "dimension of source item is empty.", __FILE__, __LINE__);
		if (!cur->mgz) error_and_exit ("prism", "magnetization of source item is empty.", __FILE__, __LINE__);

		dx = cur->dim->x;
		dy = cur->dim->y;
		dz = cur->dim->z;
		flag = SIGN (dx) * SIGN (dy) * SIGN (dz);

		x = cur->pos->x;
		y = cur->pos->y;
		z = cur->pos->z;

		a[0] = x - 0.5 * dx - x0;
		b[0] = y - 0.5 * dy - y0;
		c[0] = z - 0.5 * dz - z0;

		a[1] = a[0] + dx;
		b[1] = b[0] + dy;
		c[1] = c[0] + dz;

		prism_kernel (&tmp[0], a[1], b[1], c[0], cur->mgz);
		prism_kernel (&tmp[2], a[1], b[0], c[0], cur->mgz);
		prism_kernel (&tmp[4], a[0], b[1], c[0], cur->mgz);
		prism_kernel (&tmp[6], a[0], b[0], c[0], cur->mgz);

		if (fabs (dz) < DBL_EPSILON) {
			vector3d_set (&tmp[1], 0., 0., 0.);
			vector3d_set (&tmp[3], 0., 0., 0.);
			vector3d_set (&tmp[5], 0., 0., 0.);
			vector3d_set (&tmp[7], 0., 0., 0.);
		} else {
			prism_kernel (&tmp[1], a[1], b[1], c[1], cur->mgz);
			prism_kernel (&tmp[3], a[1], b[0], c[1], cur->mgz);
			prism_kernel (&tmp[5], a[0], b[1], c[1], cur->mgz);
			prism_kernel (&tmp[7], a[0], b[0], c[1], cur->mgz);
		}

		f->x += flag * (tmp[0].x - tmp[1].x - tmp[2].x + tmp[3].x
			- tmp[4].x + tmp[5].x + tmp[6].x - tmp[7].x);

		f->y += flag * (tmp[0].y - tmp[1].y - tmp[2].y + tmp[3].y
			- tmp[4].y + tmp[5].y + tmp[6].y - tmp[7].y);

		f->z += flag * (tmp[0].z - tmp[1].z - tmp[2].z + tmp[3].z
			- tmp[4].z + tmp[5].z + tmp[6].z - tmp[7].z);

		cur = cur->next;
	}
	vector3d_scale (f, scale_factor);
	return f;
}

static void
dipole_xz_kernel (vector3d *f, const double x1, const double z1, const vector3d *mgz)
{
	double	y, z;
	double	jy, jz;

	double	fy, fz;
	double	r2, r4;

	transform (x1, 0., z1, NULL, &y, &z);
	transform (mgz->x, 0., mgz->z, NULL, &jy, &jz);

	r2 = y * y + z * z;
	r4 = pow (r2, 2.);

	fy = jy * (1. / r2 - 2. * pow (y, 2.) / r4) - jz * 2. * y * z / r4;
	fz = - jy * 2. * y * z / r4 + jz * (1. / r2 - 2. * pow (z, 2.) / r4);

	vector3d_set (f, - 2. * fy, 0., 2. * fz);
	return;
}

static void
prism_xz_kernel (vector3d *f, const double x1, const double z1, const vector3d *mgz)
{
	double	y, z;
	double	jy, jz;

	double	fy, fz;
	double	r;

	transform (x1, 0., z1, NULL, &y, &z);
	transform (mgz->x, 0., mgz->z, NULL, &jy, &jz);

	r = sqrt (y * y + z * z);

	fy = jy * atan (z / y) + jz * log (r);
	fz = jy * log (r) + jz * atan (y / z);

	vector3d_set (f, - 2. * fy, 0., 2. * fz);
	return;
}

vector3d *
dipole_xz (const vector3d *obs, const source *s)
{
	double		x, z;
	double		x0, z0;
	vector3d	*f;
	vector3d	tmp;
	source_item	*cur;

	x0 = obs->x;
	z0 = obs->z;

	f = vector3d_new (0., 0., 0.);

	cur = s->begin;
	while (cur) {
		double	ds;

		if (!cur->pos) error_and_exit ("dipole_xz", "position of source item is empty.", __FILE__, __LINE__);
		if (!cur->mgz) error_and_exit ("dipole_xz", "magnetization of source item is empty.", __FILE__, __LINE__);

		ds = 1.0;
		if (cur->dim) {
			double	dx = cur->dim->x;
			double	dz = cur->dim->z;
			ds = fabs (dx * dz);
		}
		x = cur->pos->x - x0;
		z = cur->pos->z - z0;

		dipole_xz_kernel (&tmp, x, z, cur->mgz);
		vector3d_axpy (ds, &tmp, f);
		cur = cur->next;
	}
	vector3d_scale (f, scale_factor);
	return f;
}

vector3d *
prism_xz (const vector3d *obs, const source *s)
{
	double		b[2], c[2];
	double		x, z;
	double		x0, z0;
	vector3d	*f;
	vector3d	tmp[4];
	source_item	*cur;

	x0 = obs->x;
	z0 = obs->z;

	f = vector3d_new (0., 0., 0.);

	cur = s->begin;
	while (cur) {
		double	dx, dz;
		double	flag;

		if (!cur->pos) error_and_exit ("prism_xz", "position of source item is empty.", __FILE__, __LINE__);
		if (!cur->dim) error_and_exit ("prism_xz", "dimension of source item is empty.", __FILE__, __LINE__);
		if (!cur->mgz) error_and_exit ("prism_xz", "magnetization of source item is empty.", __FILE__, __LINE__);

		dx = cur->dim->x;
		dz = cur->dim->z;
		flag = SIGN (dx) * SIGN (dz);

		x = cur->pos->x;
		z = cur->pos->z;

		b[0] = x - 0.5 * dx - x0;
		c[0] = z - 0.5 * dz - z0;

		b[1] = b[0] + dx;
		c[1] = c[0] + dz;

		prism_xz_kernel (&tmp[0], b[1], c[0], cur->mgz);
		prism_xz_kernel (&tmp[2], b[0], c[0], cur->mgz);

		if (fabs (dz) < DBL_EPSILON) {
			vector3d_set (&tmp[1], 0., 0., 0.);
			vector3d_set (&tmp[3], 0., 0., 0.);
		} else {
			prism_xz_kernel (&tmp[1], b[1], c[1], cur->mgz);
			prism_xz_kernel (&tmp[3], b[0], c[1], cur->mgz);
		}

		f->x += flag * (tmp[0].x - tmp[1].x - tmp[2].x + tmp[3].x);
		f->z += flag * (tmp[0].z - tmp[1].z - tmp[2].z + tmp[3].z);

		cur = cur->next;
	}
	vector3d_scale (f, scale_factor);
	return f;
}

typedef enum {
	MGCAL_X_COMPONENT,	// Hx
	MGCAL_Y_COMPONENT,	// Hy
	MGCAL_Z_COMPONENT,	// Hz
	MGCAL_TOTAL_FORCE	// f
} MgcalComponent;

double
total_force (const vector3d *exf, const vector3d *f)
{
	if (exf == NULL) error_and_exit ("total_force", "vector3d *exf is empty.", __FILE__, __LINE__);
	if (f == NULL) error_and_exit ("total_force", "vector3d *f is empty.", __FILE__, __LINE__);
	return exf->x * f->x + exf->y * f->y + exf->z * f->z;
}

static double
calc_component (const vector3d *f, const source *src, MgcalComponent comp)
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

		default:
			break;
	}
	return val;
}

/*** dipole ***/
static double
component_dipole (const vector3d *obs, const source *src, MgcalComponent comp)
{
	vector3d	*f = dipole (obs, src);
	double	val = calc_component (f, src, comp);
	vector3d_free (f);
	return val;
}

double
x_component_dipole (const vector3d *obs, const source *src, void *data)
{
	return component_dipole (obs, src, MGCAL_X_COMPONENT);
}

double
y_component_dipole (const vector3d *obs, const source *src, void *data)
{
	return component_dipole (obs, src, MGCAL_Y_COMPONENT);
}

double
z_component_dipole (const vector3d *obs, const source *src, void *data)
{
	return component_dipole (obs, src, MGCAL_Z_COMPONENT);
}

double
total_force_dipole (const vector3d *obs, const source *src, void *data)
{
	return component_dipole (obs, src, MGCAL_TOTAL_FORCE);
}

/*** prism ***/
static double
component_prism (const vector3d *obs, const source *src, MgcalComponent comp)
{
	vector3d	*f = prism (obs, src);
	double	val = calc_component (f, src, comp);
	vector3d_free (f);
	return val;
}

double
x_component_prism (const vector3d *obs, const source *src, void *data)
{
	return component_prism (obs, src, MGCAL_X_COMPONENT);
}

double
y_component_prism (const vector3d *obs, const source *src, void *data)
{
	return component_prism (obs, src, MGCAL_Y_COMPONENT);
}

double
z_component_prism (const vector3d *obs, const source *src, void *data)
{
	return component_prism (obs, src, MGCAL_Z_COMPONENT);
}

double
total_force_prism (const vector3d *obs, const source *src, void *data)
{
	return component_prism (obs, src, MGCAL_TOTAL_FORCE);
}

/*** dipole yz ***/
static double
component_dipole_xz (const vector3d *obs, const source *src, MgcalComponent comp)
{
	vector3d	*f = dipole_xz (obs, src);
	double	val = calc_component (f, src, comp);
	vector3d_free (f);
	return val;
}

double
x_component_dipole_xz (const vector3d *obs, const source *src, void *data)
{
	return component_dipole_xz (obs, src, MGCAL_X_COMPONENT);
}

double
z_component_dipole_xz (const vector3d *obs, const source *src, void *data)
{
	return component_dipole_xz (obs, src, MGCAL_Z_COMPONENT);
}

double
total_force_dipole_xz (const vector3d *obs, const source *src, void *data)
{
	return component_dipole_xz (obs, src, MGCAL_TOTAL_FORCE);
}

/*** prism yz ***/
static double
component_prism_xz (const vector3d *obs, const source *src, MgcalComponent comp)
{
	vector3d	*f = prism_xz (obs, src);
	double	val = calc_component (f, src, comp);
	vector3d_free (f);
	return val;
}

double
x_component_prism_xz (const vector3d *obs, const source *src, void *data)
{
	return component_prism_xz (obs, src, MGCAL_X_COMPONENT);
}

double
z_component_prism_xz (const vector3d *obs, const source *src, void *data)
{
	return component_prism_xz (obs, src, MGCAL_Z_COMPONENT);
}

double
total_force_prism_xz (const vector3d *obs, const source *src, void *data)
{
	return component_prism_xz (obs, src, MGCAL_TOTAL_FORCE);
}
