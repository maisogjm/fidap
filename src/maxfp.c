#include <stdio.h>
#include <stdlib.h>

main (int argc, char *argv[])
{
    /* 03.07.94 : JM Maisog
	Returns the maximum of a floating point binary file.	*/

    FILE *fIN;
    float *imgin, min, max;
    int num_voxX4;
    int num_vox;
    int vox_index;
    int read_status;

    /* **************************************************************** */
    /* VERIFY NUMBER OF COMMAND LINE ARGUMENTS.				*/
    if (argc < 2) {
	printf("maxfp: need name of floating point input file.\n");
	exit(1);
    }

    /* **************************************************************** */
    /* OPEN FILES FOR INPUT AND OUTPUT.					*/
    if ((fIN=fopen(argv[1],"r")) == NULL) {
	printf("maxfp: can't open %s\n",argv[1]);
	exit(2);
    }

    /* **************************************************************** */
    /* COUNT NUMBER OF BYTES IN INPUT FILE				*/
    num_voxX4 = 0;
    while (fgetc(fIN)!=EOF)
	num_voxX4++;

    num_vox = num_voxX4/4;

    /* **************************************************************** */
    /* ASSIGN MEMORY.							*/
    imgin = (float *) malloc(num_vox*sizeof(float));
    if (imgin == NULL) {
	printf("maxfp: failed to malloc for imgin.\n");
	exit(4);
    }

    /* **************************************************************** */
    /* SET FILE POINTER BACK TO BEGINNING AND READ IN DATA.		*/

    rewind(fIN);
    read_status=fread(imgin, sizeof(float), num_vox, fIN);
    fclose(fIN);
    if (read_status != num_vox) {
        printf("%d voxels read in.  Should have been %d voxels.\n",
						read_status, num_vox);
        exit(8);
    }

    /* **************************************************************** */
    /* DETERMINE MAX OF FLOAT DATA.					*/

    max = imgin[0];
    for (vox_index=0; vox_index<num_vox; vox_index++)
	max = (max>imgin[vox_index]) ? max : imgin[vox_index];
    printf("%g\n",max);

    /* **************************************************************** */
    /* CLOSE FILES							*/
/*    fclose(fIN); Joe - 25 October 2002 - commenting out extra fclose() */
    free(imgin);
}
