#include <stdio.h>
#include <stdlib.h>
#include "analyze6.h"

struct analyze_struct readhdr(char *argv[]);
/* ************************************************************************** */
main(int argc, char *argv[])
{
	struct analyze_struct img_info;
	int i;

/* -------------------------------------------------------------------- */
/*	Check to see if correct number of arguments were passed.	*/

	if (argc < 2) {
	    printf("dordhdr4: need headerfile(s).\n");
	    exit(1);
	}

/* -------------------------------------------------------------------- */
/*	Read header file.						*/

	for (i=1; i<argc; i++) {
	    img_info = readhdr(argv+i);

	    if (img_info.hk.session_error != 0)
		printf("dordhdr2: ?error reading header\n");

	    printf("\nContents of %s:",*(argv+i));
	    printf("\nX dim = %d\n",img_info.dime.dim[1]);
	    printf("Y dim = %d\n",img_info.dime.dim[2]);
	    printf("Z dim = %d\n",img_info.dime.dim[3]);
	    printf("X pixdim = %f\n",img_info.dime.pixdim[1]);
	    printf("Y pixdim = %f\n",img_info.dime.pixdim[2]);
	    printf("Z pixdim = %f\n",img_info.dime.pixdim[3]);
	    printf("Volumes = %d\n",img_info.dime.dim[4]);
	    printf("Bits per pixel = %d\n",img_info.dime.bitpix);
	    printf("Pixel maximum = %d\n",img_info.dime.glmax);
	    printf("Pixel minimum = %d\n",img_info.dime.glmin);
	    printf("Calibration max = %f\n",img_info.dime.cal_max);
	    printf("Calibration min = %f\n",img_info.dime.cal_min);
	    printf("Data type = %d\n",img_info.dime.datatype);
	    printf("Scale factor = %f\n",img_info.dime.funused1);
	}

	exit(0);
}
