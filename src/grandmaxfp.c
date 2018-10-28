#include <stdio.h>
#include <stdlib.h>
#include "analyze6.h"
#include <float.h>

struct analyze_struct readhdr(char *argv[]);

/* ************************************************************************** */
main(int argc, char *argv[])
{
    FILE *fp;
    struct analyze_struct img_info;
    int image_index,xdim, ydim, zdim;
    float *image, max;
    int num_vox, vox_index, x, y, z;
    int read_status;

/* -------------------------------------------------------------------- */
/*	Check to see if correct number of arguments were passed.	*/

    if (argc < 2) {
	printf("grandmaxfp: need headerfile, and names of images file(s).\n");
	exit(1);
    }

/* -------------------------------------------------------------------- */
/*	Read header file.						*/

    img_info = readhdr(argv+1);
    xdim = img_info.dime.dim[1];
    ydim = img_info.dime.dim[2];
    zdim = img_info.dime.dim[3];

/*    printf("xdim = %d\n",xdim);
    printf("ydim = %d\n",ydim);
    printf("zdim = %d\n",zdim);
*/
    num_vox=xdim*ydim*zdim;
    image = (float *) malloc(num_vox*sizeof(float));

/* -------------------------------------------------------------------- */
/*	Read in images, find global max.				*/

    max=-FLT_MAX;

    for (image_index=2; image_index<argc; image_index++) {

	fp = fopen(*(argv+image_index),"r");
	    read_status=fread(image, sizeof(float), num_vox, fp);
	    if (read_status==num_vox)
		printf("Read in %d voxels from %s.\n",read_status,
							*(argv+image_index));
	    else
		printf("ERROR: read in only %d voxels from %s!\n",read_status,
							*(argv+image_index));
	fclose(fp);

	for (z=0; z<zdim; z++)
	   for (y=0; y<ydim; y++)
		for (x=0; x<xdim; x++) {
		    vox_index = (z*ydim + y)*xdim + x;

		    if (max<*(image+vox_index)) {
			printf("New max = %g at x=%d, y=%d, z=%d.\n",
			    *(image+vox_index),x+1,y+1,z+1);
			max = *(image+vox_index);
		    }
		}

    }

    free(image);
    printf("\nThe grand global maximum floating point number was %g.\n",max);

    exit(0);
}

