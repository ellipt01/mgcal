/*
 * kernel.h
 *
 *  Created on: 2015/03/14
 *      Author: utsugi
 */

#ifndef KERNEL_H_
#define KERNEL_H_

#ifdef __cplusplus
extern "C" {
#endif

typedef double	(*theoretical) (const cvector *pos, const source *src, void *data);

typedef struct s_mgcal_func	mgcal_func;

struct s_mgcal_func
{
	theoretical	function;
	void		*parameter;
};

mgcal_func	*mgcal_func_new (const theoretical func, void *data);
void		mgcal_func_free (mgcal_func *f);
void		kernel_matrix_set (double *a, const data_array *array, const grid *g, const cvector *mgz, const cvector *exf, const mgcal_func *f);
double		*kernel_matrix (const data_array *array, const grid *g, const cvector *mgz, const cvector *exf, const mgcal_func *f);
void		kernel_set_grid_volume (double *a, const data_array *array, const grid *g);

#ifdef __cplusplus
}
#endif

#endif /* KERNEL_H_ */
