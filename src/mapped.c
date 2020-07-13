/*
 * mapped.c
 *
 *  Created on: 2020/07/04
 *      Author: utsugi
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#include <unistd.h>
#include <fcntl.h>
#include <sys/mman.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "../include/vector3d.h"
#include "source.h"
#include "data_array.h"
#include "grid.h"
#include "scattered.h"
#include "kernel.h"
#include "private/util.h"
#include "mapped.h"

static mapped *
mapped_alloc ()
{
	mapped	*m = (mapped *) malloc (sizeof (mapped));
	m->fd = -1;
	m->m = -1;
	m->n = -1;
	m->nnz = -1;
	m->size = -1;
	m->data = NULL;
	m->fn[0] = '\0';
	return m;
}

mapped *
mapped_new (int id, int nobs, int nsrc, const char *dir)
{
	int		stat;
	long	psize;

	mapped	*map = mapped_alloc ();

	map->m = nobs;
	map->n = nsrc;
	map->nnz = ((long) nobs) * ((long) nsrc);

#ifdef BSD
	psize = getpagesize ();
#else
	psize = sysconf (_SC_PAGE_SIZE);
#endif

	map->size = (map->nnz * sizeof (double) / psize + 1) * psize;

	/* open mmap file */
	if (dir) sprintf (map->fn, "%s/data%05d.map", dir, id);
	else sprintf (map->fn, "./data%05d.map", id);

	map->fd = open (map->fn, O_RDWR | O_CREAT | O_TRUNC, 0666);
	if (map->fd == -1) {
		char	msg[BUFSIZ];
		sprintf (msg, "cannot open file %s.", map->fn);
		error_and_exit ("mapped_new", msg, __FILE__, __LINE__);
	}
	stat = lseek (map->fd, map->size, SEEK_SET);
	if (stat == -1) error_and_exit ("mapped_new", "lseek failed.", __FILE__, __LINE__);

	{
		char	c;
		if (read (map->fd, &c, sizeof (char)) == -1) c ='\0';
		if (write (map->fd, &c, sizeof (double)) == -1) error_and_exit ("mapped_new", "write failed.", __FILE__, __LINE__);
	}

	map->data = (double *) mmap (NULL, map->size, PROT_READ | PROT_WRITE, MAP_SHARED, map->fd, 0);
	return map;
}

void
mapped_free (mapped *map)
{
	if (map) {
		if (map->fd != -1) {
			int	stat = munmap (map->data, map->size);
			if (stat == -1) error_and_exit ("mapped_free", "msync failed.", __FILE__, __LINE__);
			close (map->fd);
		}
		free (map);
	}
	return;
}

mapped *
kernel_matrix_mapped (int id, const data_array *array, const grid *g, const vector3d *mgz, const vector3d *exf, const mgcal_func *f,
	const char *dir)
{
	int		m;
	int		n;
	int		nx;
	int		ny;
	int		nz;
	int		nh;

	mapped	*map;
	double	*a;

	m = array->n;
	nx = g->nx;
	ny = g->ny;
	nz = g->nz;
	nh = g->nh;
	n = nh * nz;

	map = mapped_new (id, m, n, dir);
	a = map->data;

#pragma omp parallel
	{
		int			stat;

		size_t		i, j, k, l;
		double		*z1 = NULL;
		vector3d	*obs = vector3d_new (0., 0., 0.);
		source		*src = source_new (0., 0.);
		if (exf) src->exf = vector3d_copy (exf);
		source_append_item (src);
		src->begin->pos = vector3d_new (0., 0., 0.);
		src->begin->dim = vector3d_new (0., 0., 0.);
		if (mgz) src->begin->mgz = vector3d_copy (mgz);

#pragma omp for
		for (k = 0; k < nz; k++) {
			double			*zk = g->z + k;	// for parallel calculation
			double			*dzk = g->dz + k;	// for parallel calculation
			double			*yj = g->y;
			double			*dyj = g->dy;
			unsigned long	offsetk = ((unsigned long) k) * ((unsigned long) nh) * ((unsigned long) m);
			double			*ak = a + offsetk;
			for (j = 0; j < ny; j++) {
				double			*xi = g->x;
				double			*dxi = g->dx;
				unsigned long	offsetj = ((unsigned long) j) * ((unsigned long) nx) * ((unsigned long) m);
				double			*aj = ak + offsetj;
				if (g->z1) z1 = g->z1 + j * nx;
				for (i = 0; i < nx; i++) {
					double	*al = aj + i * m;
					double	*xl = array->x;
					double	*yl = array->y;
					double	*zl = array->z;
					double	z1k = *zk;
					if (z1) z1k += z1[i];
					vector3d_set (src->begin->pos, *xi, *yj, z1k);
					vector3d_set (src->begin->dim, *dxi, *dyj, *dzk);
					for (l = 0; l < m; l++) {
						vector3d_set (obs, *xl, *yl, *zl);
						*al = f->function (obs, src, f->parameter);
						al++;
						xl++;
						yl++;
						zl++;
					}
					xi++;
					dxi++;
				}
				yj++;
				dyj++;
			}

			/* synclonize file */
			stat = msync (map->data, map->size, 0);
			if (stat == -1) error_and_exit ("kernel_matrix_map", "msync failed.", __FILE__, __LINE__);
		}
		vector3d_free (obs);
		source_free (src);
	}
	return map;
}

