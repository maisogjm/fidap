#include <stdio.h>
#include <stdlib.h>

/* Function to extract a voxel value from an unsigned 8-bit file, given
   an offset.									*/

unsigned char imgsamp(char *imgfile, long offset)

{
    FILE *fIMGFILE;		/* Pointer to image file.			*/
    char name[80];              /* Name of an input data file.                  */
    int fset;			/* Value returned from fseek function.		*/
    unsigned char voxval;	/* Voxel value from image file.			*/

    /* ************************************************************************ */
    /* OPEN IMAGE FILE.								*/
    if ((fIMGFILE=fopen(imgfile,"r"))==NULL) {
	printf("imgsamp: can't open %s.\n",imgfile);
	exit(1);
    }

    /* ************************************************************************ */
    /* SET FILE POINTER.							*/
    fset = fseek(fIMGFILE,offset,0);
    if (fset!=0) {
	printf("imgsamp: error setting file pointer for %s.\n",imgfile);
	exit(2);
    }

    /* ************************************************************************ */
    /* EXTRACT VOXEL VALUE.							*/
    fread(&voxval,sizeof(unsigned char),1,fIMGFILE);
    fclose(fIMGFILE);
    if (ferror(fIMGFILE)) {
	printf("imgsamp: error sampling %s.\n",imgfile);
	exit(3);
    }
    return voxval;
}
