#include <stdio.h>
#include <stdlib.h>

main (int argc, char *argv[])
{
    /* 04.25.94 : JM Maisog
	Returns the maximum of a signed 8-bit binary file.	*/

    FILE *fIN;
    unsigned char *imgin, min, max;
    int num_vox;
    int vox_index;
    int read_status;

    /* **************************************************************** */
    /* VERIFY NUMBER OF COMMAND LINE ARGUMENTS.				*/
    if (argc < 2) {
	printf("max8bit2: need name of 8-bit input file.\n");
	exit(1);
    }

    /* **************************************************************** */
    /* OPEN FILES FOR INPUT AND OUTPUT.					*/
    if ((fIN=fopen(argv[1],"r")) == NULL) {
	printf("max8bit2: can't open %s\n",argv[1]);
	exit(2);
    }

    /* **************************************************************** */
    /* COUNT NUMBER OF BYTES IN INPUT FILE				*/
    num_vox = 0;
    while (fgetc(fIN)!=EOF)
	num_vox++;

    /* **************************************************************** */
    /* ASSIGN MEMORY.							*/
    imgin = (unsigned char *) malloc(num_vox*sizeof(unsigned char));
    if (imgin == NULL) {
	printf("max8bit2: failed to malloc for imgin.\n");
	exit(4);
    }

    /* **************************************************************** */
    /* SET FILE POINTER BACK TO BEGINNING AND READ IN DATA.		*/

    rewind(fIN);
    read_status=fread(imgin, sizeof(unsigned char), num_vox, fIN);
    fclose(fIN);
    if (read_status != num_vox) {
        printf("%d voxels read in.  Should have been %d voxels.\n",
						read_status, num_vox);
        exit(8);
    }

    /* **************************************************************** */
    /* DETERMINE MAX OF FLOAT DATA.					*/

    max = 0;
    for (vox_index=0; vox_index<num_vox; vox_index++)
	max = (max>imgin[vox_index]) ? max : imgin[vox_index];
    printf("%d\n",max);

    /* **************************************************************** */
    /* CLOSE FILES							*/

    free(imgin);
}
