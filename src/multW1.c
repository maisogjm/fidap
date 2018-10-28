#include <stdio.h>
#include <stdlib.h>

main (int argc, char *argv[])
{
    /* 12.03.93 : JM Maisog : Multiplies short int images.	*/

    FILE *fIN, *fOUT;
    short int *imgin, *product;
    int cmd_lin_arg;
    int num_voxX2;
    int num_vox;
    int vox_index;
    int read_status;
    int write_status;


    /* **************************************************************** */
    /* VERIFY NUMBER OF COMMAND LINE ARGUMENTS.				*/
    if (argc < 4) {
	printf("multW1: need input files (at least 2) and output file.\n");
	exit(1);
    }

    /* **************************************************************** */
    /* OPEN FILES FOR INPUT AND OUTPUT.					*/
    printf("multW1: opening file %s...\n",argv[1]);
    if ((fIN=fopen(argv[1],"r")) == NULL) {
	printf("multW1: can't open %s\n",argv[1]);
	exit(2);
    }

    /* **************************************************************** */
    /* COUNT NUMBER OF BYTES IN FIRST INPUT FILE			*/
    num_voxX2 = 0;
    while (fgetc(fIN)!=EOF)
	num_voxX2++;

    fclose(fIN);
    num_vox = num_voxX2/2;
    printf("multW1: Number of voxels = %d\n",num_vox);


    /* **************************************************************** */
    /* ASSIGN MEMORY.							*/
    imgin = (short int *) malloc(num_vox*sizeof(short int));
    if (imgin == NULL) {
	printf("multW1: failed to malloc for imgin.\n");
	exit(4);
    }

    product = (short int *) malloc(num_vox*sizeof(short int));
    if (product == NULL) {
	printf("multW1: failed to malloc for product.\n");
	exit(5);
    }

    /* **************************************************************** */
    /* SET *product TO 1.						*/
    for (vox_index=0; vox_index<num_vox; vox_index++)
	product[vox_index] = 1;

    /* **************************************************************** */
    /* LOOP OVER ALL INPUT FILES....					*/

    for (cmd_lin_arg = 1; cmd_lin_arg<=argc-2 ;cmd_lin_arg++) {

	/* ************************************************************ */
	/* READ IN VOXELS.						*/
	fIN=fopen(argv[cmd_lin_arg],"r");
	read_status=fread(imgin, sizeof(short int), num_vox, fIN);
	fclose(fIN);
	if (read_status != num_vox) {
            printf("%d voxels read in.  Should have been %d voxels.\n",
						read_status, num_vox);
            exit(8);
	}

	/* ************************************************************ */
	/* MAKE CUMULATIVE PRODUCT.					*/
	printf("Multiplying data from %s...\n",argv[cmd_lin_arg]);
	for (vox_index=0; vox_index<num_vox; vox_index++)
    	    product[vox_index] *= imgin[vox_index];

    }

    /* **************************************************************** */
    /* WRITE OUT UNSIGNED 16 BIT TO OUTPUT FILE, AND CLOSE FILES.	*/
    printf("multW1: opening file %s...\n",argv[argc-1]);
    if ((fOUT=fopen(argv[argc-1],"w")) == NULL) {
	printf("multW1: can't open %s\n",argv[argc-1]);
	exit(3);
    }
    printf("multW1: writing to %s...\n",argv[argc-1]);
    write_status=fwrite(product, sizeof(short int), num_vox, fOUT);
    fclose(fOUT);
    if (write_status != num_vox) {
        printf("%d voxels written.  Should have been %d voxels.\n",
						write_status,num_vox);
        exit(10);
    }
    /* fclose(fOUT); Joe - 25 October 2002 - commenting out extra fclose() */
}
