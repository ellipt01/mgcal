/*
 * example_dipole.c
 *
 *  Created on: 2015/03/15
 *      Author: utsugi
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "mgcal.h"

extern const double	mgz_int;
extern const double	dec;
extern const double	inc;

void
example_dipole (FILE *stream, const int nx, const int ny, const int nz, const double x[], const double y[], const double z[])
{
	double		*a;
	grid		*g;
	source		*s;
	mgcal_func	*f;

	double		t;

	g = grid_new (nx, ny, nz, x, y, z);
	s = source_new (dec, inc);
	source_append_item (s);
	source_set_position (s, 0., 0., -2.);
	source_set_magnetization (s, 5. * mgz_int, dec, inc);

	f = mgcal_func_new (total_force_dipole, NULL);
	a = (double *) malloc (g->n * sizeof (double));

	t = omp_get_wtime ();
#pragma omp parallel
	{
		int		n;
		cvector	*obs = cvector_new (0., 0., 0.);

#pragma omp for
		for (n = 0; n < g->n; n++) {
			grid_get_nth (g, n, obs, NULL);
			a[n] = f->function (obs, s, f->parameter);
		}
		cvector_free (obs);
	}
	fprintf (stderr, "time(1) = %.4e\n", omp_get_wtime () - t);

	fwrite_grid_with_data (stream, g, a, NULL);

	grid_free (g);
	source_free (s);
	mgcal_func_free (f);
	free (a);

	return;
}

static double *
irregular_surface (const int nx, const int ny, const double x[], const double y[])
{
	int		i, j, k;
	int		n = nx * ny;
	double	dx = (x[1] - x[0]) / (double) nx;
	double	dy = (y[1] - y[0]) / (double) ny;
	double	*z1 = (double *) malloc (n * sizeof (double));
	k = 0;
	for (j = 0; j < ny; j++) {
		double	yj = y[0] + dy * (double) j;
		for (i = 0; i < nx; i++) {
			double	xi = x[0] + dx * (double) i;
			z1[k++] = 1. / sqrt (pow (2. * xi, 2.) + pow (yj, 2.) + 4.);
		}
	}
	return z1;
}

void
example_dipole_irregular_surface (FILE *stream, const int nx, const int ny, const int nz, const double x[], const double y[], const double z[])
{
	double		*z1;
	double		*a;
	grid		*g;
	source		*s;
	mgcal_func	*f;

	double		t;

	z1 = irregular_surface (nx, ny, x, y);
	g = grid_new_full (nx, ny, nz, x, y, z, NULL, NULL, NULL, z1);
	free (z1);
	s = source_new (dec, inc);
	source_append_item (s);
	source_set_position (s, 0., 0., -2.);
	source_set_magnetization (s, 5. * mgz_int, dec, inc);

	f = mgcal_func_new (total_force_dipole, NULL);
	a = (double *) malloc (g->n * sizeof (double));

	t = omp_get_wtime ();
#pragma omp parallel
	{
		int		n;
		cvector	*obs = cvector_new (0., 0., 0.);

#pragma omp for
		for (n = 0; n < g->n; n++) {
			grid_get_nth (g, n, obs, NULL);
			a[n] = f->function (obs, s, f->parameter);
		}
		cvector_free (obs);
	}
	fprintf (stderr, "time(2) = %.4e\n", omp_get_wtime () - t);

	fwrite_grid_with_data (stream, g, a, NULL);

	grid_free (g);
	source_free (s);
	mgcal_func_free (f);
	free (a);

	return;
}

void
example_dipole_multi_sources (FILE *stream, const int nx, const int ny, const int nz, const double x[], const double y[], const double z[])
{
	double		*a;
	grid		*g;
	source		*s;
	source		*cur;
	mgcal_func	*f;

	double		t;

	g = grid_new (nx, ny, nz, x, y, z);
	s = source_new (dec, inc);
	source_append_item (s);
	source_set_position (s, 0., 0., -1.);
	source_set_magnetization (s, 10., dec, inc);

	source_append_item (s);
	source_set_position (s, -2., -2., -1.5);
	source_set_magnetization (s, 15., dec, inc);

	source_append_item (s);
	source_set_position (s, 2., 4., -2.);
	source_set_magnetization (s, 25., dec, inc);

	f = mgcal_func_new (total_force_dipole, NULL);
	a = (double *) malloc (g->n * sizeof (double));

	t = omp_get_wtime ();
#pragma omp parallel
	{
		int		n;
		cvector	*obs = cvector_new (0., 0., 0.);

#pragma omp for
		for (n = 0; n < g->n; n++) {
			grid_get_nth (g, n, obs, NULL);
			a[n] = f->function (obs, s, f->parameter);
		}
		cvector_free (obs);
	}
	fprintf (stderr, "time(3) = %.4e\n", omp_get_wtime () - t);

	fwrite_grid_with_data (stream, g, a, NULL);

	grid_free (g);
	source_free (s);
	mgcal_func_free (f);
	free (a);

	return;
}
