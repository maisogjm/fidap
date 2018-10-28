#include <stdio.h>
#include <stdlib.h>

main (int argc, char *argv[])
{
    /* 01.03.93 : JM Maisog : Adds a constant to floating point images.	*/

    FILE *fIN, *fOUT;
    float *imgin, *sum;
    float offset;
    int num_voxX4;
    int num_vox;
    int vox_index;
    int read_status;
    int write_status;


    /* **************************************************************** */
    /* VERIFY NUMBER OF COMMAND LINE ARGUMENTS.				*/
    if (argc < 4) {
	printf("addFP2: need input file, number to add to input, and output file.\n");
	exit(1);
    }

    /* **************************************************************** */
    /* OPEN FILES FOR INPUT AND OUTPUT.					*/
    printf("addFP2: opening file %s...\n",argv[1]);
    if ((fIN=fopen(argv[1],"r")) == NULL) {
	printf("addFP2: can't open %s\n",argv[1]);
	exit(2);
    }

    /* **************************************************************** */
    /* COUNT NUMBER OF BYTES IN FIRST INPUT FILE			*/
    num_voxX4 = 0;
    while (fgetc(fIN)!=EOF)
	num_voxX4++;

    fclose(fIN);
    num_vox = num_voxX4/4;
    printf("addFP2: Number of voxels = %d\n",num_vox);


    /* **************************************************************** */
    /* ASSIGN MEMORY.							*/
    imgin = (float *) malloc(num_vox*sizeof(float));
    if (imgin == NULL) {
	printf("addFP2: failed to malloc for imgin.\n");
	exit(4);
    }

    sum = (float *) malloc(num_vox*sizeof(float));
    if (sum == NULL) {
	printf("addFP2: failed to malloc for sum.\n");
	exit(5);
    }

    offset = atof(argv[2]);

    /* **************************************************************** */
    /* READ IN DATA.							*/
    fIN=fopen(argv[1],"r");
    read_status=fread(imgin, sizeof(float), num_vox, fIN);
    fclose(fIN);
    if (read_status != num_vox) {
	printf("%d voxels read in.  Should have been %d voxels.\n",
							read_status, num_vox);
	exit(8);
    }


    /* **************************************************************** */
    /* MAKE SUM.							*/
    printf("Adding %f to data from %s...\n",offset,argv[1]);
    for (vox_index=0; vox_index<num_vox; vox_index++)
	sum[vox_index] = imgin[vox_index] + offset;


    /* **************************************************************** */
    /* WRITE OUT UNSIGNED 8 BIT TO OUTPUT FILE, AND CLOSE FILES.	*/
    printf("addFP2: opening file %s...\n",argv[argc-1]);
    if ((fOUT=fopen(argv[argc-1],"w")) == NULL) {
	printf("addFP2: can't open %s\n",argv[argc-1]);
	exit(3);
    }
    printf("addFP2: writing to %s...\n",argv[argc-1]);
    write_status=fwrite(sum, sizeof(float), num_vox, fOUT);
    fclose(fOUT);
    if (write_status != num_vox) {
        printf("%d voxels written.  Should have been %d voxels.\n",
						write_status,num_vox);
        exit(10);
    }
}
