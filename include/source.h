/*
 * source.h
 *
 *  Created on: 2015/03/14
 *      Author: utsugi
 */

#ifndef SOURCE_H_
#define SOURCE_H_

#ifdef __cplusplus
extern "C" {
#endif

typedef struct s_source	source;

struct s_source {
	cvector		*exf;	// external field
	cvector		*mgz;	// magnetization
	cvector		*pos;	// center of the magnetized body
	cvector		*dim;	// dimension of source
	source		*next;
};

source	*source_new (void);
void	source_free (source *src);
void	source_set_position (source *src, const double x, const double y, const double z);
void	source_set_dimension (source *src, const double dx, const double dy, const double dz);
void	source_set_external_field (source *src, const double inc, const double dec);
void	source_set_magnetization (source *src, const double mgz, const double inc, const double dec);

#ifdef __cplusplus
}
#endif

#endif /* SOURCE_H_ */
