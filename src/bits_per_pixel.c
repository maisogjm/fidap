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
	    printf("bits_per_pixel: need headerfile(s).\n");
	    exit(1);
	}

/* -------------------------------------------------------------------- */
/*	Read header file.						*/

	for (i=1; i<argc; i++) {
	    img_info = readhdr(argv+i);

	    if (img_info.hk.session_error != 0)
		printf("bits_per_pixel: error reading header.\n");

	    printf("%d\n",img_info.dime.bitpix);
	}

	exit(0);
}
