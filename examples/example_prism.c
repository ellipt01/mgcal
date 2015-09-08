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

extern const double	mgz_int;
extern const double	dec;
extern const double	inc;

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
	s = source_new (dec, inc);
	source_append_item (s);
	source_set_position (s, 0., 0., -2.);
	source_set_dimension (s, 5., 5., 2.);
	source_set_magnetization (s, 0.25 * mgz_int, dec, inc);

	f = mgcal_func_new (total_force_prism, NULL);
	a = (double *) malloc (g->n * sizeof (double));

	t = omp_get_wtime ();
#pragma omp parallel
	{
		int			n;
		vector3d	*obs = vector3d_new (0., 0., 0.);

#pragma omp for
		for (n = 0; n < g->n; n++) {
			grid_get_nth (g, n, obs, NULL);
			a[n] = f->function (obs, s, f->parameter);
		}
		vector3d_free (obs);
	}
	fprintf (stderr, "time(4) = %.4e\n", omp_get_wtime () - t);

	fwrite_grid_with_data (stream, g, a, NULL);

	grid_free (g);
	source_free (s);
	mgcal_func_free (f);
	free (a);

	return;
}
