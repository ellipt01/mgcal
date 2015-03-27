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

cvector	*dipole (const cvector *obs, source *s);
cvector	*prism (const cvector *obs, source *s);

cvector	*dipole_yz (const cvector *obs, source *s);
cvector	*prism_yz (const cvector *obs, source *s);

double	total_force (const cvector *exf, cvector *f);
double	total_force_dipole (const cvector *obs, source *src, void *data);
double	total_force_prism (const cvector *obs, source *src, void *data);
double	total_force_dipole_yz (const cvector *obs, source *src, void *data);
double	total_force_prism_yz (const cvector *obs, source *src, void *data);

#ifdef __cplusplus
}
#endif

#endif /* CALC_H_ */
