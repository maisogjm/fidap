#include <stdio.h>
#include <stdlib.h>
#include <math.h>

main (int argc, char *argv[])
{
    /* 12.03.93 : JM Maisog : Multiplies floating point images.	*/

    FILE *fIN, *fOUT;
    float *imgin, *product;
    int cmd_lin_arg;
    int num_voxX4;
    int num_vox;
    int vox_index;
    int read_status;
    int write_status;


    /* **************************************************************** */
    /* VERIFY NUMBER OF COMMAND LINE ARGUMENTS.				*/
    if (argc < 4) {
	printf("multFP1: need input files (at least 2) and output file.\n");
	exit(1);
    }

    /* **************************************************************** */
    /* OPEN FILES FOR INPUT AND OUTPUT.					*/
    printf("multFP1: opening file %s...\n",argv[1]);
    if ((fIN=fopen(argv[1],"r")) == NULL) {
	printf("multFP1: can't open %s\n",argv[1]);
	exit(2);
    }

    /* **************************************************************** */
    /* COUNT NUMBER OF BYTES IN FIRST INPUT FILE			*/
    num_voxX4 = 0;
    while (fgetc(fIN)!=EOF)
	num_voxX4++;

    fclose(fIN);
    num_vox = num_voxX4/4;
    printf("multFP1: Number of voxels = %d\n",num_vox);


    /* **************************************************************** */
    /* ASSIGN MEMORY.							*/
    imgin = (float *) malloc(num_vox*sizeof(float));
    if (imgin == NULL) {
	printf("multFP1: failed to malloc for imgin.\n");
	exit(4);
    }

    product = (float *) malloc(num_vox*sizeof(float));
    if (product == NULL) {
	printf("multFP1: failed to malloc for product.\n");
	exit(5);
    }

    /* **************************************************************** */
    /* SET *product TO 1.						*/
    for (vox_index=0; vox_index<num_vox; vox_index++)
	product[vox_index] = 1.;

    /* **************************************************************** */
    /* LOOP OVER ALL INPUT FILES....					*/

    for (cmd_lin_arg = 1; cmd_lin_arg<=argc-2 ;cmd_lin_arg++) {

	/* ************************************************************ */
	/* READ IN VOXELS.						*/
	fIN=fopen(argv[cmd_lin_arg],"r");
	read_status=fread(imgin, sizeof(float), num_vox, fIN);
	fclose(fIN);
	if (read_status != num_vox) {
            printf("%d voxels read in.  Should have been %d voxels.\n",
						read_status, num_vox);
            exit(8);
	}

	/* ************************************************************ */
	/* MAKE CUMULATIVE PRODUCT.					*/
	printf("multFP1: multiplying data from %s...\n",argv[cmd_lin_arg]);
	for (vox_index=0; vox_index<num_vox; vox_index++)
    	    product[vox_index] *= imgin[vox_index];

    }

    /* **************************************************************** */
    /* WRITE OUT FLOATING POINT DATA TO OUTPUT FILE, AND CLOSE FILES.	*/
    printf("multFP1: opening file %s...\n",argv[argc-1]);
    if ((fOUT=fopen(argv[argc-1],"w")) == NULL) {
	printf("multFP1: can't open %s\n",argv[argc-1]);
	exit(3);
    }
    printf("multFP1: writing to %s...\n",argv[argc-1]);
    write_status=fwrite(product, sizeof(float), num_vox, fOUT);
    fclose(fOUT);
    if (write_status != num_vox) {
        printf("%d voxels written.  Should have been %d voxels.\n",
						write_status,num_vox);
        exit(10);
    }
}
