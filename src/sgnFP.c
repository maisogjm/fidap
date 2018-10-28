#include <stdio.h>
#include <stdlib.h>

main (int argc, char *argv[])
{
    /* 12.03.93 : JM Maisog : Returns pixel-by-pixel signs of floating point images.	*/

    FILE *fIN, *fOUT;
    float *imgin, *imgout;
    int cmd_lin_arg;
    int num_voxX4;
    int num_vox;
    int vox_index;
    int read_status;
    int write_status;


    /* **************************************************************** */
    /* VERIFY NUMBER OF COMMAND LINE ARGUMENTS.				*/
    if (argc < 3) {
	printf("sgnFP: need input file and output file.\n");
	exit(1);
    }

    /* **************************************************************** */
    /* COUNT NUMBER OF BYTES IN FIRST INPUT FILE			*/
    printf("sgnFP: opening file %s...\n",argv[1]);
    if ((fIN=fopen(argv[1],"r")) == NULL) {
	printf("sgnFP: can't open %s\n",argv[1]);
	exit(2);
    }
    num_voxX4 = 0;
    while (fgetc(fIN)!=EOF)
	num_voxX4++;

    fclose(fIN);
    num_vox = num_voxX4/4;
    printf("sgnFP: Number of voxels = %d\n",num_vox);


    /* **************************************************************** */
    /* ASSIGN MEMORY.							*/
    imgin = (float *) malloc(num_vox*sizeof(float));
    if (imgin == NULL) {
	printf("sgnFP: failed to malloc for imgin.\n");
	exit(4);
    }
    imgout = (float *) malloc(num_vox*sizeof(float));
    if (imgout == NULL) {
	printf("sgnFP: failed to malloc for imgout.\n");
	exit(5);
    }


    /* **************************************************************** */
    /* READ IN VOXELS FROM INPUT FILE.					*/
    if ((fIN=fopen(argv[1],"r")) == NULL) {
	printf("sgnFP: unable to open %s...\n",argv[1]);
	exit(6);
    }
    printf("sgnFP: reading in voxels from %s...\n",argv[1]);
    read_status=fread(imgin, sizeof(float), num_vox, fIN);
    fclose(fIN);
    if (read_status != num_vox) {
	printf("%d voxels read in.  Should have been %d voxels.\n",
						read_status, num_vox);
	exit(8);
    }

    /* **************************************************************** */
    /* PUT SIGN IN *imgout.					        */
    printf("Taking signs of %s...\n",argv[1]);
    for (vox_index=0; vox_index<num_vox; vox_index++)
	if (imgin[vox_index]==0.)
	    imgout[vox_index] = 0.;
	else if (imgin[vox_index]>0.)
            imgout[vox_index] = 1.;
	else
	    imgout[vox_index] = -1.;


    /* **************************************************************** */
    /* WRITE TO OUTPUT FILE, AND CLOSE FILES.				*/
    printf("sgnFP: opening file %s...\n",argv[2]);
    if ((fOUT=fopen(argv[2],"w")) == NULL) {
	printf("sgnFP: can't open %s\n",argv[2]);
	exit(3);
    }
    printf("sgnFP: writing to %s...\n",argv[2]);
    write_status=fwrite(imgout, sizeof(float), num_vox, fOUT);
    fclose(fOUT);
    if (write_status != num_vox) {
        printf("%d voxels written.  Should have been %d voxels.\n",
						write_status,num_vox);
        exit(10);
    }
    fclose(fOUT);
}
