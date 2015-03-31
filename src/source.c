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

static source_item *
source_item_alloc (void)
{
	source_item	*item = (source_item *) malloc (sizeof (source_item));
	item->mgz = NULL;
	item->pos = NULL;
	item->dim = NULL;
	item->next = NULL;
	return item;
}

static void
source_item_free (source_item *si)
{
	if (si) {
		if (si->mgz) free (si->mgz);
		if (si->pos) free (si->pos);
		if (si->dim) free (si->dim);
		free (si);
	}
	return;
}

int
source_append_item (source *src)
{
	int			i = 0;
	source_item	*cur = src->item;
	while (cur->next) {
		cur = cur->next;
		i++;
	}
	cur->next = source_item_alloc ();
	src->end = cur->next;
	return i;
}

static source *
source_alloc (void)
{
	source	*src = (source *) malloc (sizeof (source));
	src->exf = NULL;
	src->item = source_item_alloc ();
	src->begin = NULL;
	src->end = NULL;
	return src;
}

source *
source_new (const double inc, const double dec)
{
	source	*src = source_alloc ();
	src->exf = cvector_new_with_geodesic_poler (1., inc, dec);
	source_append_item (src);
	src->begin = src->item->next;
	return src;
}

void
source_free (source *src)
{
	if (src) {
		source_item	*cur;
		cvector_free (src->exf);
		cur = src->item;
		while (cur) {
			source_item	*next = cur->next;
			source_item_free (cur);
			cur = next;
		}
		free (src);
	}
	return;
}

void
source_set_position (source *src, const double x, const double y, const double z)
{
	src->end->pos = cvector_new (x, y, z);
	return;
}

void
source_set_dimension (source *src, const double dx, const double dy, const double dz)
{
	src->end->dim = cvector_new (dx, dy, dz);
	return;
}

void
source_set_magnetization (source *src, const double mgz, const double inc, const double dec)
{
	src->end->mgz = cvector_new_with_geodesic_poler (mgz, inc, dec);
	return;
}
