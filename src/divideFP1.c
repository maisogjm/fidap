#include <stdio.h>
#include <stdlib.h>
#include <math.h>

main (int argc, char *argv[])
{
    /* 12.03.93 : JM Maisog : Divides floating point images.	*/

    FILE *fIN, *fOUT;
    float *imgin1, *imgin2, *imgout;
    int cmd_lin_arg;
    int num_voxX4;
    int num_vox;
    int vox_index;
    int read_status;
    int write_status;


    /* **************************************************************** */
    /* VERIFY NUMBER OF COMMAND LINE ARGUMENTS.				*/
    if (argc < 3) {
	printf("divideFP1: need two input files and one output file.\n");
	exit(1);
    }

    /* **************************************************************** */
    /* COUNT NUMBER OF BYTES IN FIRST INPUT FILE			*/
    printf("divideFP1: opening file %s...\n",argv[1]);
    if ((fIN=fopen(argv[1],"r")) == NULL) {
	printf("divideFP1: can't open %s\n",argv[1]);
	exit(2);
    }
    num_voxX4 = 0;
    while (fgetc(fIN)!=EOF)
	num_voxX4++;

    fclose(fIN);
    num_vox = num_voxX4/4;
    printf("divideFP1: Number of voxels = %d\n",num_vox);


    /* **************************************************************** */
    /* ASSIGN MEMORY.							*/
    imgin1 = (float *) malloc(num_vox*sizeof(float));
    if (imgin1 == NULL) {
	printf("divideFP1: failed to malloc for imgin1.\n");
	exit(4);
    }
    imgin2 = (float *) malloc(num_vox*sizeof(float));
    if (imgin2 == NULL) {
	printf("divideFP1: failed to malloc for imgin2.\n");
	exit(4);
    }
    imgout = (float *) malloc(num_vox*sizeof(float));
    if (imgout == NULL) {
	printf("divideFP1: failed to malloc for imgout.\n");
	exit(5);
    }

    /* **************************************************************** */
    /* ZERO OUT *imgout.						*/
    for (vox_index=0; vox_index<num_vox; vox_index++)
	imgout[vox_index] = 0.;

    /* **************************************************************** */
    /* READ IN VOXELS FROM FIRST INPUT FILE.				*/
    if ((fIN=fopen(argv[1],"r")) == NULL) {
	printf("divideFP1: unable to open %s...\n",argv[1]);
	exit(6);
    }
    printf("divideFP1: reading in voxels from %s...\n",argv[1]);
    read_status=fread(imgin1, sizeof(float), num_vox, fIN);
    fclose(fIN);
    if (read_status != num_vox) {
	printf("%d voxels read in.  Should have been %d voxels.\n",
						read_status, num_vox);
	exit(8);
    }

    /* **************************************************************** */
    /* READ IN VOXELS FROM SECOND INPUT FILE.				*/
    if ((fIN=fopen(argv[2],"r")) == NULL) {
	printf("divideFP1: unable to open %s...\n",argv[2]);
	exit(6);
    }
    printf("divideFP1: reading in voxels from %s...\n",argv[2]);
    read_status=fread(imgin2, sizeof(float), num_vox, fIN);
    fclose(fIN);
    if (read_status != num_vox) {
	printf("%d voxels read in.  Should have been %d voxels.\n",
						read_status, num_vox);
	exit(8);
    }

    /* **************************************************************** */
    /* PUT DIVIDEND IN *imgout.						*/
    printf("Dividing %s by %s...\n",argv[1], argv[2]);
    printf("If any voxel in %s is zero, the dividend will be set to zero.\n",
								argv[2]);
    for (vox_index=0; vox_index<num_vox; vox_index++)
	if (imgin2[vox_index]==0)
	    imgout[vox_index] = 0.;
	else
	    imgout[vox_index] = imgin1[vox_index] / imgin2[vox_index];


    /* **************************************************************** */
    /* WRITE TO OUTPUT FILE, AND CLOSE FILES.				*/
    printf("divideFP1: opening file %s...\n",argv[3]);
    if ((fOUT=fopen(argv[3],"w")) == NULL) {
	printf("divideFP1: can't open %s\n",argv[3]);
	exit(3);
    }
    printf("divideFP1: writing to %s...\n",argv[3]);
    write_status=fwrite(imgout, sizeof(float), num_vox, fOUT);
    fclose(fOUT);
    if (write_status != num_vox) {
        printf("%d voxels written.  Should have been %d voxels.\n",
						write_status,num_vox);
        exit(10);
    }
    exit(0);
}
