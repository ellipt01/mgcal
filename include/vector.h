/*
 * vector.h
 *
 *  Created on: 2015/03/14
 *      Author: utsugi
 */

#ifndef VECTOR_H_
#define VECTOR_H_

#ifdef __cplusplus
extern "C" {
#endif

#define	deg2rad(d)	((d) * M_PI / 180.)
#define	rad2deg(r)	((r) * 180. / M_PI)

/* cartesian vector */
typedef struct s_cvector	cvector;
/* poler vector */
typedef struct s_pvector	pvector;

struct s_cvector {
	double	x;
	double	y;
	double	z;
};

typedef enum {
	XCOMP,
	YCOMP,
	ZCOMP
} CCOMPONENT;

cvector	*cvector_new (const double x, const double y, const double z);
cvector	*cvector_new_with_geodesic_poler (const double r, const double inc, const double dec);
void	cvector_set (cvector *cv, const double x, const double y, const double z);
void	cvector_free (cvector *cv);

cvector	*cvector_copy (const cvector *src);
void	cvector_add (cvector *v1, const cvector *v2);
void	cvector_sub (cvector *v1, const cvector *v2);

#ifdef __cplusplus
}
#endif

#endif /* VECTOR_H_ */
