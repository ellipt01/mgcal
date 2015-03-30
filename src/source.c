/*
 * source.c
 *
 *  Created on: 2015/03/14
 *      Author: utsugi
 */

#include <stdlib.h>

#include "cvector.h"
#include "source.h"
#include "util.h"

static source *
source_alloc (void)
{
	source	*s = (source *) malloc (sizeof (source));
	s->exf = NULL;
	s->mgz = NULL;
	s->pos = NULL;
	s->dim = NULL;
	s->next = NULL;
	return s;
}

source *
source_new (void)
{
	source	*s = source_alloc ();
	return s;
}

static void
s_source_free (source *s)
{
	if (s) {
		if (s->exf) free (s->exf);
		if (s->mgz) free (s->mgz);
		if (s->pos) free (s->pos);
		if (s->dim) free (s->dim);
		free (s);
	}
	return;
}

void
source_free (source *s)
{
	source	*cur = s;
	while (cur) {
		source	*next = cur->next;
		s_source_free (cur);
		cur = next;
	}
	return;
}

void
source_set_position (source *src, const double x, const double y, const double z)
{
	src->pos = cvector_new (x, y, z);
	return;
}

void
source_set_dimension (source *src, const double dx, const double dy, const double dz)
{
	src->dim = cvector_new (dx, dy, dz);
	return;
}

void
source_set_external_field (source *src, const double inc, const double dec)
{
	src->exf = cvector_new_with_geodesic_poler (1., inc, dec);
	return;
}

void
source_set_magnetization (source *src, const double mgz, const double inc, const double dec)
{
	src->mgz = cvector_new_with_geodesic_poler (mgz, inc, dec);
	return;
}
