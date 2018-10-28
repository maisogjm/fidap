#include <stdio.h>
#include <stdlib.h>

main (int argc, char *argv[])
{
    /* 03.07.94 : JM Maisog
	Returns the minimum of a floating point binary file.	*/

    FILE *fIN;
    float *imgin, min, max;
    int num_voxX4;
    int num_vox;
    int vox_index;
    int read_status;

    /* **************************************************************** */
    /* VERIFY NUMBER OF COMMAND LINE ARGUMENTS.				*/
    if (argc < 2) {
	printf("minfp: need name of floating point input file.\n");
	exit(1);
    }

    /* **************************************************************** */
    /* OPEN FILES FOR INPUT AND OUTPUT.					*/
    if ((fIN=fopen(argv[1],"r")) == NULL) {
	printf("minfp: can't open %s\n",argv[1]);
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
	printf("minfp: failed to malloc for imgin.\n");
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
    /* DETERMINE MIN OF FLOAT DATA.					*/

    min = imgin[0];
    for (vox_index=0; vox_index<num_vox; vox_index++)
	min = (min<imgin[vox_index]) ? min : imgin[vox_index];
    printf("%f\n",min);

    /* **************************************************************** */
    /* CLOSE FILES							*/
/*    fclose(fIN); Joe - 25 October 2002 - commenting out extra fclose() */
    free(imgin);
}
