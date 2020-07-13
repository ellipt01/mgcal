/*
 * examples.h
 *
 *  Created on: 2015/03/15
 *      Author: utsugi
 */

#ifndef EXAMPLES_H_
#define EXAMPLES_H_

void	example_dipole (FILE *stream, const int nx, const int ny, const int nz, const double x[], const double y[], const double z[]);
void	example_dipole_irregular_surface (FILE *stream, const int nx, const int ny, const int nz, const double x[], const double y[], const double z[]);
void	example_dipole_multi_sources (FILE *stream, const int nx, const int ny, const int nz, const double x[], const double y[], const double z[]);
void	example_prism (FILE *stream, const int nx, const int ny, const int nz, const double x[], const double y[], const double z[]);
double	*example_kernel (const int nx, const int ny, const int nz, const double x[], const double y[], const double z[], const double zobs);
double	*example_kernel_nth_grid (const int nx, const int ny, const int nz, const double x[], const double y[], const double z[], const double zobs);
double	*example_kernel_mth_site (const int nx, const int ny, const int nz, const double x[], const double y[], const double z[], const double zobs);

#endif /* EXAMPLES_H_ */
