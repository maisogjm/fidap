#include <stdio.h>
#include <stdlib.h>
#include "analyze6.h"
#include <float.h>

struct analyze_struct readhdr(char *argv[]);

/* ************************************************************************** */
main(int argc, char *argv[])
{
    /* 09.08.93 : JM Maisog
	Skips a specified number of bytes from an input file, and
	then passes a specified number of bytes to an output file,
	APPENDING to that file if it already exist.                */

    FILE *fIN;
    struct analyze_struct img_info;
    float *mat;
    int read_status;
    int xdim,ydim,zdim;
    int x,y,z;

    if (argc < 6) {
	printf("dopixsampfp: need input headerfile, input imagefile,\n");
	printf("             AND X-, Y-, & Z-coordinate.\n");
	exit(1);
    }

    img_info = readhdr(argv+1);
    xdim = img_info.dime.dim[1];
    ydim = img_info.dime.dim[2];
    zdim = img_info.dime.dim[3];

    mat = (float *) malloc(xdim*ydim*zdim*sizeof(float
));
    if (mat==NULL) {
	printf("dopixsampfp16bit: malloc failed for mat.\n");
	exit(2);
    }

    /* OPEN FILES FOR INPUT AND OUTPUT.                                 */
    if ((fIN=fopen(argv[2],"r")) == NULL) {
        printf("dopixsampfp: can't open %s\n",argv[2]);
        exit(3);
    }
    read_status=fread(mat,sizeof(float),xdim*ydim*zdim,fIN);
    fclose(fIN);
    if (read_status!=xdim*ydim*zdim) {
	printf("dopixsampfp: only %d voxels read in, should have been %d.\n",read_status);
	exit(4);
    }

    x = atoi(argv[3]);
    y = atoi(argv[4]);
    z = atoi(argv[5]);
    if ((x<1) || (x>xdim) || (y<1) || (y>ydim) || (z<1) || (z>zdim)) {
	printf("dopixsampfp: bad coordinates for given dimensions.\n");
	printf("             Cannot do dopixsampfp.\n");
	exit(5);
    }

    /* Correct coordinates for use with pointer */
    x--;
    y--;
    z--;

    printf("%g\n",*(mat+(z*ydim+y)*xdim+x));
}
