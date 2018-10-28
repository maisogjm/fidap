#include <stdio.h>
#include <stdlib.h>

main (int argc, char *argv[])
{
    /* 12.23.93 : JM Maisog : Finds mean of floating point images.	*/

    FILE *fIN, *fOUT;
    float *imgin, *sum;
    int cmd_lin_arg;
    int num_voxX4;
    int num_vox;
    int vox_index;
    int read_status;
    int write_status;


    /* **************************************************************** */
    /* VERIFY NUMBER OF COMMAND LINE ARGUMENTS.				*/
    if (argc < 3) {
	printf("meanFP: need input file(s) and output file.\n");
	exit(1);
    }

    /* **************************************************************** */
    /* COUNT NUMBER OF BYTES IN FIRST INPUT FILE			*/
    if ((fIN=fopen(argv[1],"r")) == NULL) {
	printf("meanFP: can't open %s\n",argv[1]);
	exit(2);
    }

    num_voxX4 = 0;
    while (fgetc(fIN)!=EOF)
	num_voxX4++;

    fclose(fIN);
    num_vox = num_voxX4/4;
    printf("meanFP: Number of voxels = %d\n",num_vox);


    /* **************************************************************** */
    /* ASSIGN MEMORY.							*/
    imgin = (float *) malloc(num_vox*sizeof(float));
    if (imgin == NULL) {
	printf("meanFP: failed to malloc for imgin.\n");
	exit(4);
    }

    sum = (float *) malloc(num_vox*sizeof(float));
    if (sum == NULL) {
	printf("meanFP: failed to malloc for sum.\n");
	exit(5);
    }

    /* **************************************************************** */
    /* ZERO OUT *sum.							*/
    for (vox_index=0; vox_index<num_vox; vox_index++)
	sum[vox_index] = 0.;

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
	/* CUMULATIVE SUM INPUT IMAGE INTO *sum.			*/
	printf("Summing data from %s...\n",argv[cmd_lin_arg]);
	for (vox_index=0; vox_index<num_vox; vox_index++)
    	    sum[vox_index] += imgin[vox_index];

    }

    /* **************************************************************** */
    /* FIND MEAN BY DIVIDING SUM BY NUMBER OF INPUT IMAGES.		*/
    printf("Making mean image...\n");
    for (vox_index=0; vox_index<num_vox; vox_index++)
	sum[vox_index] /= (argc-2);

    /* **************************************************************** */
    /* WRITE TO OUTPUT FILE, AND CLOSE FILES.				*/
    printf("meanFP: opening file %s...\n",argv[argc-1]);
    if ((fOUT=fopen(argv[argc-1],"w")) == NULL) {
	printf("meanFP: can't open %s\n",argv[argc-1]);
	exit(3);
    }
    printf("meanFP: writing to %s...\n",argv[argc-1]);
    write_status=fwrite(sum, sizeof(float), num_vox, fOUT);
    fclose(fOUT);
    if (write_status != num_vox) {
        printf("%d voxels written.  Should have been %d voxels.\n",
						write_status,num_vox);
        exit(10);
    }
    free(imgin);
    free(sum);
}
