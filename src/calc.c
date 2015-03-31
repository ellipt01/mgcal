/*
 * calc.c
 *
 *  Created on: 2015/03/14
 *      Author: utsugi
 */

#include <stdio.h>
#include <math.h>
#include <float.h>

#include "cvector.h"
#include "source.h"
#include "calc.h"

static void
append_ (const int flag, cvector *f0, const cvector *f1)
{
	if (flag == 1) cvector_add (f0, f1);
	else if (flag == -1) cvector_sub (f0, f1);
	return;
}

static cvector *
dipole_kernel (const double x, const double y, const double z, const cvector *mag)
{
	double	fx, fy, fz;

	double	r  = sqrt (x * x + y * y + z * z);
	double	r3 = pow (r, 3.0);
	double	r5 = pow (r, 5.0);

	double	jx = mag->x;
	double	jy = mag->y;
	double	jz = mag->z;

	fx =
		- jx * (1.0 / r3 - 3.0 * pow (x, 2.0) / r5)
		+ jy * (3.0 * x * y / r5)
		+ jz * (3.0 * x * z / r5);

	fy =
		+ jx * (3.0 * y * x / r5)
		- jy * (1.0 / r3 - 3.0 * pow (y, 2.0) / r5)
		+ jz * (3.0 * y * z / r5);

	fz =
		+ jx * (3.0 * z * x / r5)
		+ jy * (3.0 * z * y / r5)
		- jz * (1.0 / r3 - 3.0 * pow (z, 2.0) / r5);

	return cvector_new (fx, fy, fz);
}

static cvector *
prism_kernel (const double x, const double y, const double z, const cvector *mag)
{
	double	fx, fy, fz;
	double	lnx, lny, lnz;

	double	r = sqrt (x * x + y * y + z * z);

	if (fabs (y) <= DBL_EPSILON && fabs (z) <= DBL_EPSILON && x < 0.0)
		lnx = - log (r - x);
	else
		lnx = log (r + x);

	if (fabs (x) <= DBL_EPSILON && fabs (z) <= DBL_EPSILON && y < 0.0)
		lny = - log (r - y);
	else
		lny = log (r + y);

	if (fabs (x) <= DBL_EPSILON && fabs (y) <= DBL_EPSILON && z < 0.0)
		lnz = - log (r - z);
	else
		lnz = log (r + z);

	{
		double	jx = mag->x;
		double	jy = mag->y;
		double	jz = mag->z;

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

	return cvector_new (fx, fy, fz);
}

cvector *
dipole (const cvector *obs, const source *s)
{
	double		x, y, z;
	double		x0, y0, z0;
	cvector		*f;
	source_item	*cur;


	if (!obs) error_and_exit ("dipole", "cvector *obs is empty.", __FILE__, __LINE__);
	if (!s) error_and_exit ("dipole", "source *s is empty.", __FILE__, __LINE__);

	x0 = obs->x;
	y0 = obs->y;
	z0 = obs->z;

	f = cvector_new (0., 0., 0.);
	cur = s->begin;
	while (cur) {
		cvector	*tmp;
		x = cur->pos->x;
		y = cur->pos->y;
		z = cur->pos->z;
		tmp = dipole_kernel (x - x0, y - y0, z - z0, cur->mgz);
		append_ (1, f, tmp);
		cvector_free (tmp);
		cur = cur->next;
	}
	return f;
}

/*********************************************
    cvector	*pos  : 観測点位置
    source	*s    : プリズムのプロパティ
*********************************************/
cvector *
prism (const cvector *obs, const source *s)
{
	int			i, j, k;
	double		a[2], b[2], c[2];
	double		x, y, z;
	double		x0, y0, z0;
	double		dx, dy, dz;
	cvector		*f;
	source_item	*cur;

	if (!obs) error_and_exit ("prism", "cvector *obs is empty.", __FILE__, __LINE__);
	if (!s) error_and_exit ("prism", "source *s is empty.", __FILE__, __LINE__);

	// 観測点座標
	x0 = obs->x;
	y0 = obs->y;
	z0 = obs->z;

	f = cvector_new (0., 0., 0.);

	cur = s->begin;
	while (cur) {
		// プリズムのディメンジョン
		if (!cur->dim) error_and_exit ("prism", "dimension of source is empty.", __FILE__, __LINE__);
		dx = cur->dim->x;
		dy = cur->dim->y;
		dz = cur->dim->z;

		// プリズムの座標
		x = cur->pos->x;	// プリズム東西中心
		y = cur->pos->y;	// 　　　　南北中心
		z = cur->pos->z;	// 　　　　上端深さ

		a[0] = x - 0.5 * dx - x0;
		b[0] = y - 0.5 * dy - y0;
		c[0] = z - 0.5 * dz - z0;

		a[1] = x + 0.5 * dx - x0;
		b[1] = y + 0.5 * dy - y0;
		c[1] = z + 0.5 * dz - z0;

		for (i = 0; i <= 1; i++) {
			for (j = 0; j <= 1; j++) {
				for (k = 0; k <= 1; k++) {
					int		flag = 1;
					cvector	*tmp;
					if (i == 0) flag *= -1;
					if (j == 0) flag *= -1;
					if (k == 0) flag *= -1;

					tmp = prism_kernel (a[i], b[j], c[k], cur->mgz);
					append_ (flag, f, tmp);
					cvector_free (tmp);
				}
			}
		}
		cur = cur->next;
	}
	return f;
}

static cvector *
dipole_yz_kernel (const double y, const double z, const cvector *mag)
{
	double	fy, fz;

	double	jy = mag->y;
	double	jz = mag->z;

	double	r2 = y * y + z * z;
	double	r4 = pow(r2, 2.);

	fy = jy * (1. / r2 - 2. * pow (y, 2.) / r4) - jz * 2. * y * z / r4;
	fz = - jy * 2. * y * z / r4 + jz * (1. / r2 - 2. * pow (z, 2.) / r4);

	return cvector_new (0., - 2. * fy, - 2. * fz);
}

static cvector *
prism_yz_kernel (const double y, const double z, const cvector *mag)
{
	double	fy, fz;

	double	jy = mag->y;
	double	jz = mag->z;

	double	r = sqrt (y * y + z * z);

	fy = jy * atan(z / y) + jz * log(r);
	fz = jy * log(r) + jz * atan(y / z);

	return cvector_new (0., - 2. * fy, - 2. * fz);
}

cvector *
dipole_yz (const cvector *obs, const source *s)
{
	double		y, z;
	double		y0, z0;
	cvector		*f;
	source_item	*cur;

	// 観測点座標
	y0 = obs->y;
	z0 = obs->z;

	f = cvector_new (0., 0., 0.);
	cur = s->begin;
	while (cur) {
		cvector	*tmp;
		// dipoleの座標
		y = cur->pos->y - y0;
		z = cur->pos->z - z0;

		tmp = dipole_yz_kernel (y, z, cur->mgz);
		append_ (1, f, tmp);
		cvector_free (tmp);
		cur = cur->next;
	}
	return f;
}

cvector *
prism_yz (const cvector *obs, const source *s)
{
	int			j, k;
	double		b[2], c[2];
	double		y, z;
	double		y0, z0;
	double		dy, dz;
	cvector		*f;
	source_item	*cur;

	// 観測点座標
	y0 = obs->y;
	z0 = obs->z;

	f = cvector_new (0., 0., 0.);

	cur = s->begin;
	while (cur) {
		// プリズムのディメンジョン
		if (!cur->dim) error_and_exit ("prism_yz", "dimension of source is empty.", __FILE__, __LINE__);
		dy = cur->dim->y;
		dz = cur->dim->z;

		// プリズムの座標
		y = cur->pos->y;	// プリズム南北中心
		z = cur->pos->z;	// 　　　　上端深さ

		b[0] = y - 0.5 * dy - y0;
		c[0] = z - 0.5 * dz - z0;

		b[1] = y + 0.5 * dy - y0;
		c[1] = z + 0.5 * dz - z0;

		for (j = 0; j <= 1; j++) {
			for (k = 0; k <= 1; k++) {
				int		flag = 1;
				cvector	*tmp;
				if (j == 0) flag *= -1;
				if (k == 0) flag *= -1;

				tmp = prism_yz_kernel (b[j], c[k], cur->mgz);
				append_ (flag, f, tmp);
				cvector_free (tmp);
			}
		}
		cur = cur->next;
	}
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
	double	val;
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
