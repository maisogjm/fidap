#include <stdio.h>
#include <stdlib.h>

main (int argc, char *argv[])
{
    /* 10.13.93 : JM Maisog
	Returns the maximum of a signed 16-bit binary file.	*/

    FILE *fIN;
    signed short int *imgin, min, max;
    char in_name[30];
    int num_voxX2;
    int num_vox;
    int vox_index;
    int read_status;

    /* **************************************************************** */
    /* VERIFY NUMBER OF COMMAND LINE ARGUMENTS.				*/
    if (argc < 2) {
	printf("max16bit: need root name input file.\n");
	exit(1);
    }

    /* **************************************************************** */
    /* DETERMINE FILE NAME.						*/

    sprintf(in_name,"%s.img",argv[1]);

    /* **************************************************************** */
    /* OPEN FILES FOR INPUT AND OUTPUT.					*/
    if ((fIN=fopen(in_name,"r")) == NULL) {
	printf("max16bit: can't open %s\n",in_name);
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
	printf("max16bit: failed to malloc for imgin.\n");
	exit(4);
    }

    /* **************************************************************** */
    /* CLOSE, THEN RE-OPEN INPUT FILE, AND  READ IN WORD DATA		*/

    fclose(fIN);
    fIN=fopen(in_name,"r");
    read_status=fread(imgin, sizeof(signed short int), num_vox, fIN);
    fclose(fIN);
    if (read_status != num_vox) {
        printf("%d voxels read in.  Should have been %d voxels.\n",
						read_status, num_vox);
        exit(8);
    }

    /* **************************************************************** */
    /* DETERMINE MAX OF WORD DATA.					*/

    max = -32767;
    for (vox_index=0; vox_index<num_vox; vox_index++)
	max = (max>imgin[vox_index]) ? max : imgin[vox_index];
    printf("%d\n",max);

    /* **************************************************************** */
    /* CLOSE FILES							*/
    fclose(fIN);
    free(imgin);
}
