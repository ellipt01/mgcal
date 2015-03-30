/*
 * calc.h
 *
 *  Created on: 2015/03/14
 *      Author: utsugi
 */

#ifndef CALC_H_
#define CALC_H_

#ifdef __cplusplus
extern "C" {
#endif

cvector	*dipole (const cvector *obs, const source *s);
cvector	*prism (const cvector *obs, const source *s);

cvector	*dipole_yz (const cvector *obs, const source *s);
cvector	*prism_yz (const cvector *obs, const source *s);

double	total_force (const cvector *exf, const cvector *f);

double	x_component_dipole (const cvector *obs, const source *src, void *data);
double	y_component_dipole (const cvector *obs, const source *src, void *data);
double	z_component_dipole (const cvector *obs, const source *src, void *data);
double	total_force_dipole (const cvector *obs, const source *src, void *data);

double	x_component_prism (const cvector *obs, const source *src, void *data);
double	y_component_prism (const cvector *obs, const source *src, void *data);
double	z_component_prism (const cvector *obs, const source *src, void *data);
double	total_force_prism (const cvector *obs, const source *src, void *data);

double	y_component_prism_yz (const cvector *obs, const source *src, void *data);
double	z_component_prism_yz (const cvector *obs, const source *src, void *data);
double	total_force_dipole_yz (const cvector *obs, const source *src, void *data);

double	y_component_prism_yz (const cvector *obs, const source *src, void *data);
double	z_component_prism_yz (const cvector *obs, const source *src, void *data);
double	total_force_prism_yz (const cvector *obs, const source *src, void *data);

#ifdef __cplusplus
}
#endif

#endif /* CALC_H_ */
