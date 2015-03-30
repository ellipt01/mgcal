/*
 * kernel.c
 *
 *  Created on: 2015/03/15
 *      Author: utsugi
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

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
mgcal_func_new (theoretical func, void *data)
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
kernel_matrix_set (double *a, data_array *array, grid *g, cvector *mgz, cvector *exf, mgcal_func *f)
{
	int		m;
	int		nx;
	int		ny;
	int		nz;


	if (!a) error_and_exit ("kernel_matrix_set", "double *a is empty.", __FILE__, __LINE__);

	m = array->n;
	nx = g->nx;
	ny = g->ny;
	nz = g->nz;

#pragma omp parallel
	{
		int		i, j, k, l;
		double	*z1 = NULL;
		cvector	*obs = cvector_new (0., 0., 0.);
		source	*src = source_new ();
		src->pos = cvector_new (0., 0., 0.);
		src->dim = cvector_new (0., 0., 0.);
		if (exf) src->exf = exf;
		if (mgz) src->mgz = mgz;

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
					cvector_set (src->pos, *xi, *yj, z1k);
					cvector_set (src->dim, g->dx[i], g->dy[j], g->dz[k]);
					for (l = 0; l < array->n; l++) {
						int	index = (k * nx * ny + j * nx + i) * m + l;	// for parallel calculation
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
kernel_matrix (data_array *array, grid *g, cvector *mgz, cvector *exf, mgcal_func *f)
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
