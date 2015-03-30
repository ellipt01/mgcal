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
	char		buf[BUFSIZ];
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
fwrite_data_array (FILE *stream, const data_array *array, const char *format)
{
	int		i;
	char	fm[BUFSIZ];
	if (!format) strcpy (fm, "%f %f %f %f\n");
	else sprintf (fm, "%s\n", format);
	for (i = 0; i < array->n; i++) fprintf (stream, fm, array->x[i], array->y[i], array->z[i], array->data[i]);
	return;
}

void
fwrite_grid (FILE *stream, const grid *g, const char *format)
{
	int		n;
	cvector	*pos;
	char	fm[BUFSIZ];
	if (!format) strcpy (fm, "%f %f %f\n");
	else sprintf (fm, "%s\n", format);

	pos = cvector_new (0., 0., 0.);
	for (n = 0; n < g->n; n++) {
		grid_get_nth (g, n, pos, NULL);
		fprintf (stream, fm, pos->x, pos->y, pos->z);
	}
	cvector_free (pos);
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
		grid_get_nth (g, n, pos, NULL);
		fprintf (stream, fm, pos->x, pos->y, pos->z, data[n]);
	}
	cvector_free (pos);
	return;
}
