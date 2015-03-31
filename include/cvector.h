/*
 * vector.h
 *
 *  Created on: 2015/03/14
 *      Author: utsugi
 */

#ifndef CVECTOR_H_
#define CVECTOR_H_

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
void	cvector_axpy (const double alpha, const cvector *x, cvector *y);

#ifdef __cplusplus
}
#endif

#endif /* CVECTOR_H_ */
