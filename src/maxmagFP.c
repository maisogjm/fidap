#include <stdio.h>
#include <stdlib.h>
#include <math.h>

float *float_malloc(int num_items);
int read_float_data(char *filename, float *array, int num_items);
int float_fwrite(char *filename, float *array, int num_items);

main (int argc, char *argv[])
{
    /* 03.12.96 : JM Maisog : Finds value of maximum absolute value of
       each voxel across several floating point images.                */

    FILE *fIN;
    float *imgin, *vox_max;
    int cmd_lin_arg;
    int num_voxX4;
    int num_vox;
    int vox_index;
    int read_status;
    int write_status;


    /* **************************************************************** */
    /* VERIFY NUMBER OF COMMAND LINE ARGUMENTS.                         */
    if (argc < 3) {
        printf("maxmagFP: need input file(s) and output file.\n");
	printf("Returns the value having the max absolute value,\n");
	printf("voxel-by-voxel.\n");
        exit(1);
    }

    /* **************************************************************** */
    /* COUNT NUMBER OF BYTES IN FIRST INPUT FILE                        */
    if ((fIN=fopen(argv[1],"r")) == NULL) {
        printf("maxmagFP: can't open %s\n",argv[1]);
        exit(2);
    }

    num_voxX4 = 0;
    while (fgetc(fIN)!=EOF)
        num_voxX4++;

    fclose(fIN);
    num_vox = num_voxX4/4;
    printf("maxmagFP: Number of voxels = %d\n",num_vox);


    /* **************************************************************** */
    /* ASSIGN MEMORY.                                                   */
    imgin   = float_malloc(num_vox);
    vox_max = float_malloc(num_vox);

    /* **************************************************************** */
    /* INITIALIZE *vox_max TO FLOATING POINT MINIMUM VALUE.             */
    for (vox_index=0; vox_index<num_vox; vox_index++)
        vox_max[vox_index] = 0;

    /* **************************************************************** */
    /* LOOP OVER ALL INPUT FILES....                                    */

    for (cmd_lin_arg = 1; cmd_lin_arg<=argc-2 ;cmd_lin_arg++) {

        /* ************************************************************ */
        /* READ IN VOXELS.                                              */
        printf("Reading in data from %s...\n",argv[cmd_lin_arg]);
        read_float_data(argv[cmd_lin_arg], imgin, num_vox);

        /* ************************************************************ */
        /* SET vox_max EQUAL TO THE MAX BETWEEN vox_max AND imgin.      */
        for (vox_index=0; vox_index<num_vox; vox_index++)
            vox_max[vox_index]=((float)fabs((double)imgin[vox_index])>vox_max[vox_index]) 
                                ? imgin[vox_index] : vox_max[vox_index];

    }


    /* **************************************************************** */
    /* WRITE TO OUTPUT FILE, AND CLOSE FILES.                           */
    printf("maxmagFP: opening file %s...\n",argv[argc-1]);
    float_fwrite(argv[argc-1], vox_max, num_vox);

    free(imgin);
    free(vox_max);
}
