#include <stdio.h>
#include <stdlib.h>
#include "analyze6.h"
main(int argc, char *argv[])
{
    FILE *fhdr;
    FILE *fimgin;
    FILE *fimgout;
    struct analyze_struct img_info;
    char headerfile[30];
    char in_name[30];
    char out_name[30];
    int x,y,z;
    int xdim, ydim, zdim;
    int num_vox;
    signed short int *imgin;
    signed short int *imgout;
    int read_status;
    int write_status;
    
    /* ------------------------------------------------------------------------ */
    /* CHECK TO SEE IF CORRECT NUMBER OF ARGUMENTS WERE PASSED.			*/

    if (argc < 2) {
	printf("flipY16bit: need rootname of input ANALYZE .img file\n");
	printf("            and rootname of output ANALYZE .img file.\n");
	exit(1);
    }

    /* ------------------------------------------------------------------------ */
    /* DETERMINE FILE NAMES.							*/

    sprintf(headerfile,"%s.hdr",argv[1]);
    sprintf(in_name,"%s.img",argv[1]);
    sprintf(out_name,"%s.img",argv[2]);

    printf(".hdr file = %s\n",headerfile);
    printf("input .img file = %s\n",in_name);
    printf("output .img file = %s\n",out_name);

    /* ------------------------------------------------------------------------ */
    /* READ IN xdim, ydim & zdim FROM ANALYZE HEADER.				*/

    if ((fhdr = fopen(headerfile,"r")) == NULL) {
	printf("flipY16bit: can't open %s\n",headerfile);
	exit(1);
    }
    if(fread(&img_info, sizeof(struct analyze_struct), 1, fhdr) != 1){
	printf("flipY16bit: error reading %s\n",headerfile);
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
    /* ALLOCATE MEMORY FOR *imgin and *imgout.					*/

    imgin = (signed short int *) malloc(num_vox*sizeof(signed short int));
    if (imgin == NULL) {
	printf("flipY16bit: failed to malloc for imgin.\n");
	exit(4);
    }
     
    imgout = (signed short int *) malloc(num_vox*sizeof(signed short int));
    if (imgout == NULL) {
	printf("flipY16bit: failed to malloc for imgout.\n");
	exit(4);
    }
     
    /* ------------------------------------------------------------------------ */
    /* READ IN IMAGE.								*/

    if ((fimgin=fopen(in_name,"r"))==NULL) {
	printf("flipY16bit: unable to open %s.\n",in_name);
	exit(9);
    }
    printf("flipY16bit: reading in %s...\n",in_name);
    read_status=fread(imgin, sizeof(signed short int), num_vox, fimgin);
    fclose(fimgin);
    if (read_status != num_vox) {
        printf("%d voxels read in.  Should have been %d voxels.\n",read_status,
                                                                num_vox);
        exit(8);
    }


    /* ------------------------------------------------------------------------ */
    /* FLIP IN Y.								*/

    for (z=0; z<zdim; z++)
	for (y=0; y<ydim; y++)
	    for (x=0; x<xdim; x++)
		imgout[(z*ydim+y)*xdim+x]=imgin[(z*ydim+ydim-1-y)*xdim+x];


    /* ------------------------------------------------------------------------ */
    /* WRITE OUT TO DISK.							*/

    if ((fimgout=fopen(out_name,"w")) == NULL) {
        printf("flipY16bit: can't open %s\n",out_name);
        exit(3);
    }
    printf("Writing flipped %s to %s...\n",in_name,out_name);
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
