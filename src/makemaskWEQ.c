#include <stdio.h>
#include <stdlib.h>
#include <math.h>

main (int argc, char *argv[])
{
    /* 03.09.94 : JM Maisog
	Makes a mask out of a floating point image.

       10.11.94 : JMM : Modified to work on word images. */

    FILE *fIN, *fOUT;
    signed short int *imgin, blobnumber;
    unsigned char *imgout;
    int num_voxX2;
    int num_vox;
    int vox_index;
    int read_status;
    int write_status;

    /* **************************************************************** */
    /* VERIFY NUMBER OF COMMAND LINE ARGUMENTS.				*/
    if (argc < 4) {
	printf("makemaskWEQ: need 16-bit input file, output file, and\n");
	printf("              a 16-bit value.\n");
	printf("An 8-bit mask will be made from the 16-bit input.\n");
	printf("All word values equal to the input 16-bit value will be\n");
	printf("set to 1 in the mask.\n");
	exit(1);
    }

    /* **************************************************************** */
    /* OPEN FILES FOR INPUT AND OUTPUT.					*/
    if ((fIN=fopen(argv[1],"r")) == NULL) {
	printf("makemaskWEQ: can't open %s\n",argv[1]);
	exit(2);
    }
    if ((fOUT=fopen(argv[2],"w")) == NULL) {
	printf("makemaskWEQ: can't open %s\n",argv[2]);
	exit(3);
    }

    /* **************************************************************** */
    /* ASSIGN VALUE TO blobnumber.					*/
    blobnumber = atoi(argv[3]);
    printf("blobnumber = %d\n",blobnumber);

    /* **************************************************************** */
    /* COUNT NUMBER OF BYTES IN INPUT FILE				*/
    num_voxX2 = 0;
    while (fgetc(fIN)!=EOF)
	num_voxX2++;

    num_vox = num_voxX2/2;
    printf("makemaskWEQ: Number of voxels = %d\n",num_vox);

    /* **************************************************************** */
    /* ASSIGN MEMORY.							*/
    imgin = (signed short int *) malloc(num_vox*sizeof(signed short int));
    if (imgin == NULL) {
	printf("makemaskWEQ: failed to malloc for imgin.\n");
	exit(4);
    }

    imgout = (unsigned char *) malloc(num_vox*sizeof(unsigned char));
    if (imgout == NULL) {
	printf("makemaskWEQ: failed to malloc for imgout.\n");
	exit(5);
    }

    /* **************************************************************** */
    /* READ IN FP DATA.							*/
    rewind(fIN);
    read_status=fread(imgin, sizeof(signed short int), num_vox, fIN);
    fclose(fIN);
    if (read_status != num_vox) {
        printf("%d voxels read in.  Should have been %d voxels.\n",read_status,
                                                                num_vox);
        exit(8);
    }

    /* **************************************************************** */
    /* MASK FP DATA TO UNSIGNED 8 BIT.					*/

    for (vox_index=0; vox_index<num_vox; vox_index++)
	imgout[vox_index] = (imgin[vox_index]==blobnumber) ? 1 : 0;

    /* **************************************************************** */
    /* WRITE OUT UNSIGNED 8 BIT TO OUTPUT FILE				*/
    write_status=fwrite(imgout, sizeof(unsigned char), num_vox, fOUT);
    fclose(fOUT);
    if (write_status != num_vox) {
        printf("%d voxels written.  Should have been %d voxels.\n",write_status,
                                                                num_vox);
        exit(10);
    }
    exit(0);
}
