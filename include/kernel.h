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

typedef double	(*theoretical) (const cvector *pos, source *src, void *data);

typedef struct s_mgcal_func	mgcal_func;

struct s_mgcal_func
{
	theoretical	function;
	void		*parameter;
};

mgcal_func	*mgcal_func_new (theoretical func, void *data);
void		mgcal_func_free (mgcal_func *f);
void		kernel_matrix_set (double *a, data_array *array, grid *g, cvector *mgz, cvector *exf, mgcal_func *f);
double		*kernel_matrix (data_array *array, grid *g, cvector *mgz, cvector *exf, mgcal_func *f);

#ifdef __cplusplus
}
#endif

#endif /* KERNEL_H_ */
