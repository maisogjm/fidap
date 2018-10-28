#include <stdio.h>

/* ************************************************************************** */
void centergrav16bit(unsigned short int *mat, int xdim, int ydim, int zdim, double *center)
/* For use with fread method of reading in image data. */
/* 07.20.93 : JMM : Modified from centergrav2 to handle 16-bit data. */
{
	int x,y,z;
	double xsum,ysum,zsum,num;

	xsum = ysum = zsum = num = 0.;
/*	printf("centergrav: calculating center of gravity...\n"); */
	for (z = 0; z <= zdim-1; z++)
	    for (y = 0; y <= ydim-1; y++)
		for (x = 0; x <= xdim-1; x++)
		    if (*(mat + z*xdim*ydim + y*xdim + x) != 0) {
			xsum += (float) x;
			ysum += (float) y;
			zsum += (float) z;
			++num;
		    }
	
	*(center+0) = xsum/(double) num;
	*(center+1) = ysum/(double) num;
	*(center+2) = zsum/(double) num;
}
