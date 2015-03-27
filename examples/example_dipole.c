/*
 * example_dipole.c
 *
 *  Created on: 2015/03/15
 *      Author: utsugi
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "mgcal.h"

void
example_dipole (FILE *stream, const int nx, const int ny, const int nz, const double x[], const double y[], const double z[])
{
	double		*a;
	grid		*g;
	source		*s;
	mgcal_func	*f;

	g = grid_new (nx, ny, nz, x, y, z);
	s = source_new ();
	source_set_position (s, 0., 0., -2.);
	source_set_magnetization (s, 10., 45., -7.);
	source_set_external_field (s, 45., -7.);

	f = mgcal_func_new (total_force_dipole, NULL);


	a = (double *) malloc (g->n * sizeof (double));

	{
		int			i, j;
		cvector		*obs = cvector_new (0., 0., 0);

		for (j = 0; j < g->ny; j++) {
			for (i = 0; i < g->nx; i++) {
				cvector_set (obs, g->x[i], g->y[j], g->z[0]);
				a[i + j * g->nx] = f->function (obs, s, f->parameter);
			}
		}
		cvector_free (obs);
	}

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

	z1 = irregular_surface (nx, ny, x, y);
	g = grid_new_full (nx, ny, nz, x, y, z, NULL, NULL, NULL, z1);
	free (z1);
	s = source_new ();
	source_set_position (s, 0., 0., -2.);
	source_set_magnetization (s, 10., 45., -7.);
	source_set_external_field (s, 45., -7.);

	f = mgcal_func_new (total_force_dipole, NULL);
	a = (double *) malloc (g->n * sizeof (double));

	{
		int		i, j;
		cvector	*obs = cvector_new (0., 0., 0);

		for (j = 0; j < g->ny; j++) {
			for (i = 0; i < g->nx; i++) {
				double	zk = g->z[0];
				if (g->z1) zk += g->z1[i];
				cvector_set (obs, g->x[i], g->y[j], zk);
				a[i + j * g->nx] = f->function (obs, s, f->parameter);
			}
		}
		cvector_free (obs);
	}

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

	g = grid_new (nx, ny, nz, x, y, z);
	s = source_new ();

	cur = s;
	source_set_position (cur, 0., 0., -2.);
	source_set_magnetization (cur, 10., 45., -7.);
	source_set_external_field (cur, 45., -7.);

	cur->next = source_new ();
	cur = cur->next;
	source_set_position (cur, -2., -2., -2.);
	source_set_magnetization (cur, 5., 45., -7.);
	source_set_external_field (cur, 45., -7.);

	cur->next = source_new ();
	cur = cur->next;
	source_set_position (cur, 2., -4., -2.);
	source_set_magnetization (cur, 15., 45., -7.);
	source_set_external_field (cur, 45., -7.);

	f = mgcal_func_new (total_force_dipole, NULL);
	a = (double *) malloc (g->n * sizeof (double));

	{
		int		i, j;
		cvector	*obs = cvector_new (0., 0., 0);

		for (j = 0; j < g->ny; j++) {
			for (i = 0; i < g->nx; i++) {
				double	zk = g->z[0];
				if (g->z1) zk += g->z1[i];
				cvector_set (obs, g->x[i], g->y[j], zk);
				a[i + j * g->nx] = f->function (obs, s, f->parameter);
			}
		}
		cvector_free (obs);
	}

	fwrite_grid_with_data (stream, g, a, NULL);

	grid_free (g);
	source_free (s);
	mgcal_func_free (f);
	free (a);

	return;
}
