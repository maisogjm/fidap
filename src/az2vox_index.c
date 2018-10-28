#include <stdio.h>
#include <stdlib.h>

/* ************************************************************************** */
main(int argc, char *argv[])
/*
*/
{
	int xdim,ydim,zdim;	/* Image dimensions.				*/
	int vox_index;		/* Voxel index number.				*/
	int x,y,z;		/* XYZ coordinates.				*/

	/* ****************************************************************** */
	/* CHECK COMMAND LINE ARGUMENTS*/
        if (argc != 7) {
            printf("az2vox_index: Need XYZ coordinates, xdim, ydim, and zdim.\n");
            exit(1);
        }
	x = atoi(argv[1]);
	y = atoi(argv[2]);
	z = atoi(argv[3]);
	xdim = atoi(argv[4]);
	ydim = atoi(argv[5]);
	zdim = atoi(argv[6]);
	if (x<1) {
	    printf("x is too small.\n");
	    exit(2);
	}
	if (x>xdim) {
	    printf("x is too big.\n");
	    exit(3);
	}
	if (y<1) {
	    printf("y is too small.\n");
	    exit(4);
	}
	if (y>ydim) {
	    printf("y is too big.\n");
	    exit(5);
	}
	if (z<1) {
	    printf("z is too small.\n");
	    exit(6);
	}
	if (z>zdim) {
	    printf("z is too big.\n");
	    exit(7);
	}

	/* ******************************************************************** */
	/* DETERMINE vox_index.							*/

	vox_index=((z-1)*xdim*ydim+(y-1)*xdim)+x-1;
	printf("%d\n",vox_index);
}
