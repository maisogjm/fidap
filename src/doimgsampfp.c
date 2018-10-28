#include <stdio.h>
#include <stdlib.h>
#include "analyze6.h"


struct analyze_struct readhdr(char *argv[]);

float imgsampfp(char *imgfile, long offset);
main(int argc, char *argv[])
{
    struct analyze_struct img_info;
    int fset;                   /* Value returned from fseek function.          */
    int xdim,ydim,zdim;
    long offset;
    int x,y,z;
    float voxval;


    /* ************************************************************************ */
    /* CHECK COMMAND LINE ARGUMENTS.                                            */
    if (argc != 6) {
        printf("doimgsampfp: need input headerfile, input imagefile,\n");
        printf("             AND X-, Y-, & Z-coordinate.\n");
        exit(1);
    }

    /* ************************************************************************ */
    /* GET HEADER INFO.                                                         */
    img_info = readhdr(argv+1);
    xdim = img_info.dime.dim[1];
    ydim = img_info.dime.dim[2];
    zdim = img_info.dime.dim[3];

    /* ************************************************************************ */
    /* READ X, Y, AND Z COORDINATES FROM COMMAND LINE.                          */
    x = atoi(argv[3]);
    y = atoi(argv[4]);
    z = atoi(argv[5]);
    if ((x<1) || (x>xdim) || (y<1) || (y>ydim) || (z<1) || (z>zdim)) {
        printf("dopixsamp: bad coordinates for given dimensions.\n");
        exit(2);
    }

    /* ************************************************************************ */
    /* CORRECT COORDINATES FOR POINTER ARITHMETIC.                              */
    /* THEN DETERMINE OFFSET.                                                   */
    x--;
    y--;
    z--;
    offset = (z*ydim+y)*xdim+x;

    /* ************************************************************************ */
    /* CALL imgsampfp.                                                          */
    voxval = imgsampfp(argv[2],4*offset);
    printf("%g\n",voxval);
    }
