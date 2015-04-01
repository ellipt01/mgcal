/*
 * kernel.c
 *
 *  Created on: 2015/03/15
 *      Author: utsugi
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "cvector.h"
#include "source.h"
#include "data_array.h"
#include "grid.h"
#include "kernel.h"
#include "util.h"

static mgcal_func *
mgcal_func_alloc (void)
{
	mgcal_func	*f = (mgcal_func *) malloc (sizeof (mgcal_func));
	f->function = NULL;
	f->parameter = NULL;
	return f;
}

mgcal_func *
mgcal_func_new (const theoretical func, void *data)
{
	mgcal_func	*f = mgcal_func_alloc ();
	f->function = func;
	f->parameter = data;
	return f;
}

void
mgcal_func_free (mgcal_func *f)
{
	if (f) free (f);
	return;
}

void
kernel_matrix_set (double *a, const data_array *array, const grid *g, const cvector *mgz, const cvector *exf, const mgcal_func *f)
{
	int		m;
	int		nx;
	int		ny;
	int		nz;
	int		nh;

	if (!a) error_and_exit ("kernel_matrix_set", "double *a is empty.", __FILE__, __LINE__);

	m = array->n;
	nx = g->nx;
	ny = g->ny;
	nz = g->nz;
	nh = g->nh;

#pragma omp parallel
	{
		int		i, j, k, l;
		double	*z1 = NULL;
		cvector	*obs = cvector_new (0., 0., 0.);
		source	*src = source_new (0., 0.);
		if (exf) src->exf = cvector_copy (exf);
		source_append_item (src);
		src->begin->pos = cvector_new (0., 0., 0.);
		src->begin->dim = cvector_new (0., 0., 0.);
		if (mgz) src->begin->mgz = cvector_new (0., 0., 0.);

#pragma omp for
		for (k = 0; k < nz; k++) {
			double	*zk = g->z + k;	// for parallel calculation
			double	*yj = g->y;
			for (j = 0; j < ny; j++) {
				double	*xi = g->x;
				if (g->z1) z1 = g->z1 + j * nx;
				for (i = 0; i < nx; i++) {
					double	*xl = array->x;
					double	*yl = array->y;
					double	*zl = array->z;
					double	z1k = *zk;
					if (z1) z1k += z1[i];
					cvector_set (src->begin->pos, *xi, *yj, z1k);
					cvector_set (src->begin->dim, g->dx[i], g->dy[j], g->dz[k]);
					if (mgz) cvector_set (src->begin->mgz, mgz->x, mgz->y, mgz->z);
					for (l = 0; l < array->n; l++) {
						int	index = (k * nh + j * nx + i) * m + l;	// for parallel calculation
						cvector_set (obs, *xl, *yl, *zl);
						a[index] = f->function (obs, src, f->parameter);
						xl++;
						yl++;
						zl++;
					}
					xi++;
				}
				yj++;
			}
		}
		cvector_free (obs);
		free (src);
	}
	return;
}

double *
kernel_matrix (const data_array *array, const grid *g, const cvector *mgz, const cvector *exf, const mgcal_func *f)
{
	int		m, n;
	double	*a;

	m = array->n;
	n = g->n;
	a = (double *) malloc (m * n * sizeof (double));
	if (!a) error_and_exit ("kernel_matrix", "failed to allocate memory of *a.", __FILE__, __LINE__);
	kernel_matrix_set (a, array, g, mgz, exf, f);
	return a;
}

void
kernel_set_grid_volume (double *a, const data_array *array, const grid *g)
{
	int		m, n;
	int		nx, ny, nz, nh;
	cvector	*dim;

	dim = cvector_new (0., 0., 0.);

	m = array->n;
	nx = g->nx;
	ny = g->ny;
	nz = g->nz;
	nh = g->nh;

#pragma omp parallel for
	for (n = 0; n < g->n; n++) {
		grid_get_nth (g, n, NULL, dim);
		double	dx = (g->nx > 1) ? dim->x : 1.;
		double	dy = (g->ny > 1) ? dim->y : 1.;
		double	dz = (g->nz > 1) ? dim->z : 1.;
		double	dv = fabs (dx * dy * dz);
		array_scal (m, dv, a + m * n);
	}

	return;
}
