/*
 * io.c
 *
 *  Created on: 2015/03/14
 *      Author: utsugi
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>

#include "cvector.h"
#include "data_array.h"
#include "grid.h"

typedef struct s_datalist	datalist;

struct s_datalist {
	double		x;
	double		y;
	double		z;
	double		data;
	datalist	*next;
};

static datalist *
datalist_alloc (void)
{
	datalist	*list = (datalist *) malloc (sizeof (datalist));
	list->next = NULL;
	return list;
}

static datalist *
datalist_push_back (datalist *list, const double x, const double y, const double z, const double data)
{
	datalist	*p = list;
	while (p->next) p = p->next;
	p->next = datalist_alloc ();
	p = p->next;
	p->x = x;
	p->y = y;
	p->z = z;
	p->data = data;
	return p;
}

static datalist *
fread_datalist (FILE *stream, int *n)
{
	char		buf[BUFSIZ];
	datalist	*list = datalist_alloc ();
	datalist	*pp = list;
	int			count = 0;

	while (fgets (buf, BUFSIZ, stream) != NULL) {
		double	x, y, z, data;
		char	*p = buf;
		if (p[0] == '#' || p[0] == '\n') continue;
		while (p[0] == ' ' || p[0] == '\t') p++;
		sscanf (p, "%lf\t%lf\t%lf\t%lf", &x, &y, &z, &data);
		pp = datalist_push_back (pp, x, y, z, data);
		count++;
	}
	*n = count;
	return list;
}

static void
datalist_free (datalist *list)
{
	datalist	*cur = list;
	while (cur) {
		datalist	*p = cur;
		cur = cur->next;
		free (p);
	}
	return;
}

data_array *
fread_data_array (FILE *stream)
{
	data_array	*array;
	int			n, k;
	datalist	*list, *cur, *prev;

	list = fread_datalist (stream, &n);
	array = data_array_new (n);

	k = 0;
	prev = list;
	cur = list->next;
	while (cur) {
		array->x[k] = cur->x;
		array->y[k] = cur->y;
		array->z[k] = cur->z;
		array->data[k] = cur->data;
		free (prev);
		prev = cur;
		cur = cur->next;
		if (++k >= n) break;
	}
	if (cur) datalist_free (cur);
	return array;
}

void
fwrite_data_array_with_data (FILE *stream, const data_array *array, const double *data, const char *format)
{
	int		i;
	char	fm[BUFSIZ];
	if (!format) strcpy (fm, "%f %f %f %f\n");
	else sprintf (fm, "%s\n", format);
	for (i = 0; i < array->n; i++) fprintf (stream, fm, array->x[i], array->y[i], array->z[i], data[i]);
	return;
}

void
fwrite_data_array (FILE *stream, const data_array *array, const char *format)
{
	fwrite_data_array_with_data (stream, array, array->data, format);
	return;
}

static void
fprintf_array (FILE *stream, const int noneline, const int n, const double *array, const char *format)
{
	int		i;
	for (i = 0; i < n; i++) {
		fprintf (stream, format, array[i]);
		if (i < n - 1) {
			if ((i + 1) % noneline == 0) fprintf (stream, "\n");
			else fprintf (stream, " ");
		}
	}
	fprintf (stream, "\n");
	return;
}

static char *
read_valid_line (FILE *stream)
{
	char	buf[BUFSIZ];
	char	*p = NULL;
	while (1) {
		if (fgets (buf, BUFSIZ, stream) == NULL) return NULL;
		p = buf;
		while (p[0] == ' ' || p[0] == '\t') p++;
		if (p[0] == '#' || p[0] == '\n') continue;
		break;
	}
	return p;
}

static int
read_one_line (char *buf, double *x)
{
	int		i;
	char	*p;
	if (!buf) return 0;
	for (i = 0, p = strtok (buf, " "); p; p = strtok (NULL, " ")) x[i++] = (double) atof (p);
	return i;
}

grid *
fread_grid (FILE *stream)
{
	int		i;
	char	*p;
	grid	*g;

	g = (grid *) malloc (sizeof (grid));

	// read dimensions
	p = read_valid_line (stream);
	sscanf (p, "%d %d %d", &g->nx, &g->ny, &g->nz);
	g->nh = g->nx * g->ny;
	g->n = g->nh * g->nz;

	// read positions
	p = read_valid_line (stream);
	g->pos0 = cvector_new (0., 0., 0.);
	sscanf (p, "%lf %lf %lf", &g->pos0->x, &g->pos0->y, &g->pos0->z);
	p = read_valid_line (stream);
	g->pos1 = cvector_new (0., 0., 0.);
	sscanf (p, "%lf %lf %lf", &g->pos1->x, &g->pos1->y, &g->pos1->z);

	// read x
	g->x = (double *) malloc (g->nx * sizeof (double));
	i = 0;
	while (i < g->nx) {
		int	k = read_one_line (read_valid_line (stream), g->x + i);
		i += k;
	}
	// read dx
	g->dx = (double *) malloc (g->nx * sizeof (double));
	i = 0;
	while (i < g->nx) {
		int	k = read_one_line (read_valid_line (stream), g->dx + i);
		i += k;
	}
	// read y
	g->y = (double *) malloc (g->ny * sizeof (double));
	i = 0;
	while (i < g->ny) {
		int	k = read_one_line (read_valid_line (stream), g->y + i);
		i += k;
	}
	// read dy
	g->dy = (double *) malloc (g->ny * sizeof (double));
	i = 0;
	while (i < g->ny) {
		int	k = read_one_line (read_valid_line (stream), g->dy + i);
		i += k;
	}
	// read z
	g->z = (double *) malloc (g->nz * sizeof (double));
	i = 0;
	while (i < g->nz) {
		int	k = read_one_line (read_valid_line (stream), g->z + i);
		i += k;
	}
	// read dz
	g->dz = (double *) malloc (g->nz * sizeof (double));
	i = 0;
	while (i < g->nz) {
		int	k = read_one_line (read_valid_line (stream), g->dz + i);
		i += k;
	}
	// read z1
	g->z1 = (double *) malloc (g->nh * sizeof (double));
	i = 0;
	while (i < g->nh) {
		int	k = read_one_line (read_valid_line (stream), g->z1 + i);
		i += k;
	}
	if (i == 0) {
		free (g->z1);
		g->z1 = NULL;
	}

	return g;
}

const int	n_oneline = 10;

void
fwrite_grid (FILE *stream, const grid *g)
{
	fprintf (stream, "# dimension [nx, ny, nz]\n");
	fprintf (stream, "%d %d %d\n\n", g->nx, g->ny, g->nz);

	fprintf (stream, "# South-West (left-bottom) position [x, y, z]\n");
	fprintf (stream, "%f %f %f\n\n", g->pos0->x, g->pos0->y, g->pos0->z);

	fprintf (stream, "# North-East (right-top) position [x, y, z]\n");
	fprintf (stream, "%f %f %f\n\n", g->pos1->x, g->pos1->y, g->pos1->z);

	fprintf (stream, "# grid x\n");
	fprintf_array (stream, n_oneline, g->nx, g->x, "%f");
	fprintf (stream, "# grid inc x: dx\n");
	fprintf_array (stream, n_oneline, g->nx, g->dx, "%f");
	fprintf (stream, "\n");

	fprintf (stream, "# grid y\n");
	fprintf_array (stream, n_oneline, g->ny, g->y, "%f");
	fprintf (stream, "# grid inc y: dy\n");
	fprintf_array (stream, n_oneline, g->ny, g->dy, "%f");
	fprintf (stream, "\n");

	fprintf (stream, "# grid z\n");
	fprintf_array (stream, n_oneline, g->nz, g->z, "%f");
	fprintf (stream, "# grid inc z: dz\n");
	fprintf_array (stream, n_oneline, g->nz, g->dz, "%f");
	fprintf (stream, "\n");

	fprintf (stream, "# grid z1\n");
	fprintf_array (stream, n_oneline, g->nh, g->z1, "%f");
	fprintf (stream, "\n");

	return;
}

void
fwrite_grid_with_data (FILE *stream, const grid *g, const double *data, const char *format)
{
	int		n;
	cvector	*pos;
	char	fm[BUFSIZ];
	if (!format) strcpy (fm, "%f %f %f %f\n");
	else sprintf (fm, "%s\n", format);

	pos = cvector_new (0., 0., 0.);
	for (n = 0; n < g->n; n++) {
		double	val = (data) ? data[n] : 0.;
		grid_get_nth (g, n, pos, NULL);
		fprintf (stream, fm, pos->x, pos->y, pos->z, val);
	}
	cvector_free (pos);
	return;
}
