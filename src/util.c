/*
 * util.c
 *
 *  Created on: 2015/03/14
 *      Author: utsugi
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

void
error_and_exit (const char *funcname, const char *msg, const char *file, const int line)
{
	fprintf (stderr, "ERROR: %s : %s : %s %d\n", funcname, msg, file, line);
	exit (1);
}

bool
array_set_all (const int n, double *x, const double val)
{
	int		i;
	int		mod;
	double	*xi;

	if (!x) return false;

	mod = n % 4;
	xi = x;
	for (i = 0; i < mod; i++) {
		xi[0] = val;
		xi++;
	}

	for (; i < n; i += 4) {
		xi[0] = val;
		xi[1] = val;
		xi[2] = val;
		xi[3] = val;
		xi += 4;
	}
	return true;
}

bool
array_copy (const int n, double *dist, const double *src)
{
	int		i;
	int		mod;
	double	*di;
	double	*si;

	if (!dist || !src) return false;

	mod = n % 4;
	di = dist;
	si = (double *) src;
	for (i = 0; i < mod; i++) {
		di[0] = si[0];
		di++;
		si++;
	}

	for (; i < n; i += 4) {
		di[0] = si[0];
		di[1] = si[1];
		di[2] = si[2];
		di[3] = si[3];
		di += 4;
		si += 4;
	}
	return true;
}

bool
array_scal (const int n, double alpha, const double *src)
{
	int		i;
	int		mod;
	double	*si;

	if (!src) return false;

	mod = n % 4;
	si = (double *) src;
	for (i = 0; i < mod; i++) {
		si[0] *= alpha;
		si++;
	}

	for (; i < n; i += 4) {
		si[0] *= alpha;
		si[1] *= alpha;
		si[2] *= alpha;
		si[3] *= alpha;
		si += 4;
	}
	return true;
}
