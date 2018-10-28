#include <stdio.h>
#include <stdlib.h>

unsigned short *unsigned_short_malloc(int num_items);
int read_unsigned_short_int_data(char *filename, unsigned short int *array, int num_items);
int unsigned_short_fwrite(char *filename, unsigned short *array, int num_items);

main (int argc, char *argv[])
{
    /* 12.03.93 : JM Maisog : Adds floating point images.	*/

    FILE *fIN, *fOUT;
    unsigned short *imgin, *sum;
    int cmd_lin_arg;
    int num_voxX2;
    int num_vox;
    int vox_index;
    int read_status;
    int write_status;


    /* **************************************************************** */
    /* VERIFY NUMBER OF COMMAND LINE ARGUMENTS.				*/
    if (argc < 3) {
	printf("add16bit1: need input file(s) and output file.\n");
	exit(1);
    }


    /* **************************************************************** */
    /* COUNT NUMBER OF BYTES IN FIRST INPUT FILE			*/
    if ((fIN=fopen(argv[1],"r")) == NULL) {
	printf("add16bit1: can't open %s\n",argv[1]);
	exit(2);
    }
    num_voxX2 = 0;
    while (fgetc(fIN)!=EOF)
	num_voxX2++;

    fclose(fIN);
    num_vox = num_voxX2/2;
    printf("add16bit1: Number of voxels = %d\n",num_vox);


    /* **************************************************************** */
    /* ASSIGN MEMORY.							*/
    imgin = unsigned_short_malloc(num_vox);
    sum = unsigned_short_malloc(num_vox);


    /* **************************************************************** */
    /* ZERO OUT *sum.							*/
    for (vox_index=0; vox_index<num_vox; vox_index++)
	sum[vox_index] = 0;


    /* **************************************************************** */
    /* LOOP OVER ALL INPUT FILES....					*/

    for (cmd_lin_arg = 1; cmd_lin_arg<=argc-2 ;cmd_lin_arg++) {

	/* ************************************************************ */
	/* READ IN VOXELS.						*/

	read_unsigned_short_int_data(argv[cmd_lin_arg], imgin, num_vox);


	/* ************************************************************ */
	/* CUMULATIVE SUM INPUT IMAGE INTO *sum.			*/
	printf("Summing data from %s...\n",argv[cmd_lin_arg]);
	for (vox_index=0; vox_index<num_vox; vox_index++)
    	    sum[vox_index] += imgin[vox_index];

    }

    /* **************************************************************** */
    /* WRITE TO OUTPUT FILE, AND CLOSE FILES.				*/

    unsigned_short_fwrite(argv[argc-1], sum, num_vox);
}
