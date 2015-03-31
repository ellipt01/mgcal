/*
 * prism.c
 *
 *  Created on: 2015/03/19
 *      Author: utsugi
 */

#include <stdio.h>
#include <stdlib.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "mgcal.h"

void
example_prism (FILE *stream, const int nx, const int ny, const int nz, const double x[], const double y[], const double z[])
{
	int			i, j;
	double		*a;
	grid		*g;
	source		*s;
	mgcal_func	*f;

	double		t;

	g = grid_new (nx, ny, nz, x, y, z);
	s = source_new (45., -7.);
	source_append_item (s);
	source_set_position (s, 0., 0., -2.);
	source_set_dimension (s, 2., 2., 2.);
	source_set_magnetization (s, 2.5, 45., -7.);

	f = mgcal_func_new (total_force_prism, NULL);
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
	fprintf (stderr, "time(4) = %.4e\n", omp_get_wtime () - t);

	fwrite_grid_with_data (stream, g, a, NULL);

	grid_free (g);
	source_free (s);
	mgcal_func_free (f);
	free (a);

	return;
}
