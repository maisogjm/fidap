#include <math.h>
double dcos(double theta);
double dsin(double theta);

/* ************************************************************************** */
double inv3x3(double mat[3][3], double inverse[3][3])
{
	double det;

	det = mat[0][0] * (mat[1][1]*mat[2][2] - mat[1][2]*mat[2][1]) +
	      mat[0][1] * (mat[1][2]*mat[2][0] - mat[1][0]*mat[2][2]) +
	      mat[0][2] * (mat[1][0]*mat[2][1] - mat[1][1]*mat[2][0]);

	if  (det == 0.) return det;	/* Singular matrix, no inverse */

	inverse[0][0] = (mat[1][1]*mat[2][2] - mat[2][1]*mat[1][2])/det;
	inverse[1][0] = (mat[2][0]*mat[1][2] - mat[1][0]*mat[2][2])/det;
	inverse[2][0] = (mat[1][0]*mat[2][1] - mat[2][0]*mat[1][1])/det;
	inverse[0][1] = (mat[2][1]*mat[0][2] - mat[0][1]*mat[2][2])/det;
	inverse[1][1] = (mat[0][0]*mat[2][2] - mat[2][0]*mat[0][2])/det;
	inverse[2][1] = (mat[2][0]*mat[0][1] - mat[0][0]*mat[2][1])/det;
	inverse[0][2] = (mat[0][1]*mat[1][2] - mat[1][1]*mat[0][2])/det;
	inverse[1][2] = (mat[1][0]*mat[0][2] - mat[0][0]*mat[1][2])/det;
	inverse[2][2] = (mat[0][0]*mat[1][1] - mat[1][0]*mat[0][1])/det;

	return det;
}

/* ************************************************************************** */
void mkrotmat(double x_theta, double y_theta, double z_theta,
							float rotmat[3][3])
{
	rotmat[0][0] = (float) dcos(y_theta)*dcos(z_theta);
	rotmat[0][1] = (float) -dsin(z_theta);
	rotmat[0][2] = (float) dsin(y_theta)*dcos(z_theta);
	rotmat[1][0] = (float) dcos(x_theta)*dcos(y_theta)*dsin(z_theta)
						+ dsin(x_theta)*dsin(y_theta);
	rotmat[1][1] = (float) dcos(x_theta)*dcos(z_theta);
	rotmat[1][2] = (float) dcos(x_theta)*dsin(y_theta)*dsin(z_theta)
						- dsin(x_theta)*dcos(y_theta);
	rotmat[2][0] = (float) dsin(x_theta)*dcos(y_theta)*dsin(z_theta)
						- dcos(x_theta)*dsin(y_theta);
	rotmat[2][1] = (float) dsin(x_theta)*dcos(z_theta);
	rotmat[2][2] = (float) dsin(x_theta)*dsin(y_theta)*dsin(z_theta)
						+ dcos(x_theta)*dcos(y_theta);
	return;
}

/* ************************************************************************** */
double dcos(double theta)
{
	return cos(theta*asin(1.)/90.);
}

double dsin(double theta)
{
	return sin(theta*asin(1.)/90.);
}

/* ************************************************************************** */
void multmatvect(double mat[3][3], double in_vect[3], double outvect[3])
{
	int i,j;

	for (i = 0; i <= 2; i++) {
	    outvect[i] = 0.;
	    for (j = 0; j <= 2; j++)
		outvect[i] += mat[i][j] * in_vect[j];
	}
}
