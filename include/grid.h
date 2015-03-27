#ifndef grid_H
#define grid_H

#ifdef __cplusplus
extern "C" {
#endif

typedef struct s_grid	grid;

struct s_grid {
	int			n;
	int			nh;
	int			nx;
	int			ny;
	int			nz;

	cvector		*pos0;
	cvector		*pos1;

	double		*x;
	double		*y;
	double		*z;
	double		*z1;	// irregular surface

	double		*dx;
	double		*dy;
	double		*dz;

};

grid	*grid_new (const int nx, const int ny, const int nz, const double x[], const double y[], const double z[]);
grid	*grid_new_full (const int nx, const int ny, const int nz, const double x[], const double y[], const double z[], const double *dx, const double *dy, const double *dz, const double *z1);
bool	grid_set_surface (grid *g, const double *z1);
void	grid_free (grid *g);

void	grid_get_nth (const grid *g, const int n, cvector *center, cvector *dim);

#ifdef __cplusplus
}
#endif

#endif
