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
	    printf("get_dims_from: need headerfile(s).\n");
	    exit(1);
	}

/* -------------------------------------------------------------------- */
/*	Read header file.						*/

	for (i=1; i<argc; i++) {
	    img_info = readhdr(argv+i);

	    if (img_info.hk.session_error != 0)
		printf("get_dims_from: error reading header.\n");

	    printf("%d %d %d %g %g %g\n",img_info.dime.dim[1],
	    img_info.dime.dim[2],
	    img_info.dime.dim[3],
            img_info.dime.pixdim[1],
            img_info.dime.pixdim[2],
            img_info.dime.pixdim[3]);
	}

	exit(0);
}
