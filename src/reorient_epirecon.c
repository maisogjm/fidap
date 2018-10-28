#include <stdio.h>
#include <stdlib.h>
#include "analyze6.h"

main (int argc, char *argv[])
{
    /* 02.01.93 : JM Maisog
	Reorients output images of Peter Jezzard's offline
	reconstruction program.					*/

    FILE *fHDR, *fIN, *fOUT;
    struct analyze_struct img_info;
    unsigned short int *imgin, *imgout;
    char in_header[60];
    char in_name[60];
    char out_name[60];
    int num_voxX2;
    int num_vox;
    int x,y,z,vox_index;
    int xdim,ydim,zdim;
    int read_status, write_status;

    /* **************************************************************** */
    /* VERIFY NUMBER OF COMMAND LINE ARGUMENTS.				*/
    if (argc < 3) {
	printf("reorient_epirecon: need root name input file and output file.\n");
	exit(1);
    }

    /* **************************************************************** */
    /* DETERMINE FILE NAMES.						*/

    sprintf(in_header,"%s.hdr",argv[1]);
    sprintf(in_name,"%s.img",argv[1]);
    sprintf(out_name,"%s.img",argv[2]);
/*    printf("  Reorienting %s, output to %s.\n",in_name,out_name); */

    /* **************************************************************** */
    /*  READ IN xdim, ydim & zdim FROM HEADER.				*/
    if ((fHDR=fopen(in_header,"r"))==NULL) {
	printf("reorient_epirecon: unable to open %s.\n",in_header);
	exit(2);
    }
    if(fread(&img_info, sizeof(struct analyze_struct), 1, fHDR) != 1){
	printf("reorient_epirecon: error reading %s\n",in_header);
	exit(3);
    }
    fclose(fHDR);
    xdim = img_info.dime.dim[1];
    ydim = img_info.dime.dim[2];
    zdim = img_info.dime.dim[3];
/*    printf("INPUT DIMENSIONS ARE:\n");
    printf("  xdim = %d\n",xdim);
    printf("  ydim = %d\n",ydim);
    printf("  zdim = %d\n",zdim);
    printf("OUTPUT DIMENSIONS WILL BE:\n");
    printf("  xdim = %d\n",ydim);
    printf("  ydim = %d\n",xdim);
    printf("  zdim = %d\n",zdim);
*/
    num_vox=xdim*ydim*zdim;

    /* **************************************************************** */
    /* ASSIGN MEMORY.							*/
    imgin = (unsigned short int *) malloc(num_vox*sizeof(unsigned short int));
    if (imgin == NULL) {
	printf("reorient_epirecon: failed to malloc for imgin.\n");
	exit(4);
    }
    imgout = (unsigned short int *) malloc(num_vox*sizeof(unsigned short int));
    if (imgout == NULL) {
	printf("reorient_epirecon: failed to malloc for imgout.\n");
	exit(5);
    }


    /* **************************************************************** */
    /* OPEN INPUT FILE, AND  READ IN WORD DATA.				*/
    if ((fIN=fopen(in_name,"r")) == NULL) {
	printf("reorient_epirecon: can't open %s\n",in_name);
	exit(6);
    }
/*    printf("  reorient_epirecon: reading in voxels from %s...\n",in_name); */
    read_status=fread(imgin, sizeof(unsigned short int), num_vox, fIN);
    fclose(fIN);
    if (read_status != num_vox) {
        printf("%d voxels read in.  Should have been %d voxels.\n",
						read_status, num_vox);
        exit(7);
    }

    /* **************************************************************** */
    /* REORIENT INPUT IMAGE INTO OUTPUT IMAGE ARRAY.			*/
/*    printf("  reorient_epirecon: reorienting %s...\n",in_name); */
    for (x=0; x<xdim; x++)
	for (y=0; y<ydim; y++)
	    for (z=0; z<zdim; z++)
		imgout[(z*xdim+(xdim-x))*ydim+(ydim-y)]=imgin[(z*ydim+y)*xdim+x];

    /* **************************************************************** */
    /* WRITE OUT REORIENTED IMAGE TO OUTPUT FILE.			*/
/*    printf("  reorient_epirecon: writing to %s...\n",out_name); */
    if ((fOUT=fopen(out_name,"w")) == NULL) {
	printf("reorient_epirecon: can't open %s\n",out_name);
	exit(8);
    }
    write_status=fwrite(imgout, sizeof(unsigned short), num_vox, fOUT);
    fclose(fOUT);
    if (write_status != num_vox) {
        printf("%d voxels written.  Should have been %d voxels.\n",
					write_status,num_vox);
        exit(9);
    }


    /* **************************************************************** */
    /* DEALLOCATE MEMORY.						*/
    free(imgin);
    free(imgout);
}
