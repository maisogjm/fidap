#include <stdio.h>
#include <stdlib.h>
#include <limits.h>

main (int argc, char *argv[])
{
    /* 03.07.94 : JM Maisog
	Returns the minimum of a signed 16-bit binary file.	*/

    FILE *fIN;
    signed short int *imgin, min, max;
    int num_voxX2;
    int num_vox;
    int vox_index;
    int read_status;

    /* **************************************************************** */
    /* VERIFY NUMBER OF COMMAND LINE ARGUMENTS.				*/
    if (argc < 2) {
	printf("min16bit2: need name of 16-bit input file.\n");
	exit(1);
    }

    /* **************************************************************** */
    /* OPEN FILES FOR INPUT AND OUTPUT.					*/
    if ((fIN=fopen(argv[1],"r")) == NULL) {
	printf("min16bit2: can't open %s\n",argv[1]);
	exit(2);
    }

    /* **************************************************************** */
    /* COUNT NUMBER OF BYTES IN INPUT FILE				*/
    num_voxX2 = 0;
    while (fgetc(fIN)!=EOF)
	num_voxX2++;

    num_vox = num_voxX2/2;

    /* **************************************************************** */
    /* ASSIGN MEMORY.							*/
    imgin = (signed short int *) malloc(num_vox*sizeof(signed short int));
    if (imgin == NULL) {
	printf("min16bit2: failed to malloc for imgin.\n");
	exit(4);
    }

    /* **************************************************************** */
    /* SET FILE POINTER BACK TO BEGINNING AND READ IN DATA.		*/

    rewind(fIN);
    read_status=fread(imgin, sizeof(signed short int), num_vox, fIN);
    fclose(fIN);
    if (read_status != num_vox) {
        printf("%d voxels read in.  Should have been %d voxels.\n",
						read_status, num_vox);
        exit(8);
    }

    /* **************************************************************** */
    /* DETERMINE MAX OF FLOAT DATA.					*/

    min = SHRT_MAX;
    for (vox_index=0; vox_index<num_vox; vox_index++)
	min = (min<imgin[vox_index]) ? min : imgin[vox_index];
    printf("%d\n",min);

    /* **************************************************************** */
    /* CLOSE FILES							*/

    free(imgin);
    exit(0);
}
