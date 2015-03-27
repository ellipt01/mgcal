/*
 * prism.c
 *
 *  Created on: 2015/03/19
 *      Author: utsugi
 */

#include <stdio.h>
#include <stdlib.h>

#include "mgcal.h"

void
example_prism (FILE *stream, const int nx, const int ny, const int nz, const double x[], const double y[], const double z[])
{
	int			i, j;
	double		*a;
	grid		*g;
	source		*s;
	mgcal_func	*f;

	cvector		*obs;

	g = grid_new (nx, ny, nz, x, y, z);
	s = source_new ();
	source_set_position (s, 0., 0., -2.);
	source_set_dimension (s, 2., 2., 2.);
	source_set_magnetization (s, 2.5, 45., -7.);
	source_set_external_field (s, 45., -7.);

	f = mgcal_func_new (total_force_prism, NULL);

	obs = cvector_new (0., 0., 0);

	a = (double *) malloc (g->n * sizeof (double));
	for (j = 0; j < g->ny; j++) {
		for (i = 0; i < g->nx; i++) {
			cvector_set (obs, g->x[i], g->y[j], g->z[0]);
			a[i + j * g->nx] = f->function (obs, s, f->parameter);
		}
	}

	fwrite_grid_with_data (stream, g, a, NULL);

	grid_free (g);
	source_free (s);
	mgcal_func_free (f);
	free (a);

	return;
}
