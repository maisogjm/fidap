#include <stdio.h>
#include <stdlib.h>
#include "analyze6.h"
main(int argc, char *argv[])
{
    FILE *fhdr;
    FILE *fimgin;
    FILE *fimgout;
    char headerfile[30];
    struct analyze_struct img_info;
    char imgfile[30];
    int x,y,z;
    int undersampling_factor;
    int num_vox;
    int vox_index;
    int xdim, ydim, zdim;
    signed short int *imgin;
    signed short int *imgout;
    signed short int min, max;
    int read_status;
    int write_status;
    
    /* ------------------------------------------------------------------------ */
    /* CHECK TO SEE IF CORRECT NUMBER OF ARGUMENTS WERE PASSED.			*/

    if (argc < 4) {
	printf("undersample: need rootname of input image file, output\n");
	printf("             image file, AND undersampling factor.\n");
	exit(1);
    }

    /* ------------------------------------------------------------------------ */
    /* DETERMINE FILE NAMES AND UNDERSAMPLING FACTOR.				*/

    sprintf(headerfile,"%s.hdr",argv[1]);
    sprintf(imgfile,"%s.img",argv[1]);
    undersampling_factor = atoi(argv[3]);

    printf("Header file = %s\n",headerfile);
    printf(".img file = %s\n",imgfile);
    printf("Undersampling Factor = %d\n",undersampling_factor);

    /* ------------------------------------------------------------------------ */
    /* READ IN xdim, ydim & zdim FROM ANALYZE HEADER.                           */

    if ((fhdr = fopen(headerfile,"r")) == NULL) {
        printf("globminmaxfp: can't open %s\n",headerfile);
        exit(1);
    }
    if(fread(&img_info, sizeof(struct analyze_struct), 1, fhdr) != 1){
        printf("globminmaxfp: error reading %s\n",headerfile);
        img_info.hk.session_error = 2;
        exit(2);
    }
    fclose(fhdr);

    xdim = img_info.dime.dim[1];
    ydim = img_info.dime.dim[2];
    zdim = img_info.dime.dim[3];
    printf("xdim = %d\n",xdim);
    printf("ydim = %d\n",ydim);
    printf("zdim = %d\n",zdim);
    num_vox = xdim*ydim*zdim;

    /* ------------------------------------------------------------------------ */
    /* ALLOCATE MEMORY FOR *scalefactor, *offset, *imgin, AND *imgout		*/

    imgin = (signed short int *) malloc(num_vox*sizeof(signed short int));
    if (imgin == NULL) {
	printf("mri2anlz: failed to malloc for imgin.\n");
	exit(4);
    }
    imgout = (signed short int *) malloc((num_vox/(undersampling_factor*undersampling_factor))*sizeof(signed short int));
    if (imgout == NULL) {
	printf("mri2anlz: failed to malloc for imgout.\n");
	exit(5);
    }
     
    /* ------------------------------------------------------------------------ */
    /* READ IN 16-BIT IMAGE.							*/

    if ((fimgin=fopen(imgfile,"r"))==NULL) {
	printf("undersample: unable to open %s.\n",imgfile);
	exit(7);
    }
    printf("undersample: reading in %s...\n",imgfile);
    read_status=fread(imgin, sizeof(signed short int), num_vox, fimgin);
    fclose(fimgin);
    if (read_status != num_vox) {
	printf("%d voxels read in.  Should have been %d voxels.\n",read_status,
								num_vox);
	exit(8);
    }

    /* ------------------------------------------------------------------------ */
    /* UNDERSAMPLE.								*

    for (x=0; x<xdim; x+=undersampling_factor)
	for (y=0; y<ydim; y+=undersampling_factor)
	    for (z=0; z<zdim; z++) {
		imgout[(z*(ydim/undersampling_factor)+(y/undersampling_factor))
			*(xdim/undersampling_factor)+(x+undersampling_factor)] =
						imgin[(z*ydim+ydim-1-y)*xdim+x];
		min = (min<imgout[(z*ydim+y)*xdim+x]) ? min : imgout[(z*ydim+y)*xdim+x];
		max = (max>imgout[(z*ydim+y)*xdim+x]) ? max : imgout[(z*ydim+y)*xdim+x];
	    }

    /* ------------------------------------------------------------------------ */
    /* WRITE OUT FLOATING POINT IMAGE.						*/

    if ((fimgout=fopen(argv[2],"w"))==NULL) {
	printf("undersample: unable to open %s.\n",argv[2]);
	exit(9);
    }
    printf("undersample: writing to %s...\n",argv[2]);
    write_status=fwrite(imgout, sizeof(signed short int), num_vox, fimgout);
    fclose(fimgout);
    if (write_status != num_vox) {
	printf("%d voxels written.  Should have been %d voxels.\n",write_status,
								num_vox);
	exit(10);
    }

    /* ------------------------------------------------------------------------ */
    /* DEALLOCATE MEMORY.							*/

    free(imgin);
    free(imgout);
}
