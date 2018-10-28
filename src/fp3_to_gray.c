#include <stdio.h>
#include <stdlib.h>

float *float_malloc(int num_items);
int read_float_data(char *filename, float *array, int num_items);
int read_unsigned_char_data(char *filename, unsigned char *array, int num_items);
unsigned char *unsigned_char_malloc(int num_items);
int unsigned_char_fwrite(char *filename, unsigned char *array, int num_items);

main (int argc, char *argv[])
{

    FILE *fIN, *fOUT;
    float *imgin1, *imgin2, *imgin3;
    unsigned char *mask, *output;
    int num_voxX4;
    int num_vox;
    int vox_index;


    /* **************************************************************** */
    /* VERIFY NUMBER OF COMMAND LINE ARGUMENTS.                         */
    if (argc != 9) {
        printf("fp3_to_gray: need 3 input files, 3 numbers to represent each, 8-bit mask, and output file.\n");
        printf("Will output an 8-bit file whose voxels will assume one of the three colors,\n");
        printf("depending on which file had the greatest voxel value at that location.\n");
        exit(1);
    }

    /* **************************************************************** */
    /* OPEN FILES FOR INPUT AND OUTPUT.                                 */
    if ((fIN=fopen(argv[1],"r")) == NULL) {
        printf("fp3_to_gray: can't open %s\n",argv[1]);
        exit(2);
    }

    /* **************************************************************** */
    /* COUNT NUMBER OF BYTES IN FIRST INPUT FILE                        */
    /* THEN READ THEM IN.                                               */

    num_voxX4 = 0;
    while (fgetc(fIN)!=EOF)
        num_voxX4++;

    num_vox = num_voxX4/4;
    printf("fp3_to_gray: Number of voxels = %d\n",num_vox);


    /* **************************************************************** */
    /* ASSIGN MEMORY, AND READ IN DATA.                                 */
    imgin1 = float_malloc(num_vox);
    imgin2 = float_malloc(num_vox);
    imgin3 = float_malloc(num_vox);
    mask = unsigned_char_malloc(num_vox);
    output = unsigned_char_malloc(num_vox);

    printf("fp3_to_gray: opening file %s...\n",argv[1]);
    read_float_data(argv[1],imgin1,num_vox);
    printf("fp3_to_gray: opening file %s...\n",argv[2]);
    read_float_data(argv[2],imgin2,num_vox);
    printf("fp3_to_gray: opening file %s...\n",argv[3]);
    read_float_data(argv[3],imgin3,num_vox);
    printf("fp3_to_gray: opening file %s...\n",argv[7]);
    read_unsigned_char_data(argv[7],mask,num_vox);


    /* **************************************************************** */
    /* MAKE GRAY-LEVEL MAP.                                             */

    for (vox_index=0; vox_index<num_vox; vox_index++) {
        if (mask[vox_index]!=0) {
            if ((imgin1[vox_index]>imgin2[vox_index]) && (imgin1[vox_index]>imgin3[vox_index]))
                output[vox_index] = atoi(argv[4]);
            else if ((imgin2[vox_index]>imgin1[vox_index]) && (imgin2[vox_index]>imgin3[vox_index]))
                output[vox_index] = atoi(argv[5]);
            else if ((imgin3[vox_index]>imgin2[vox_index]) && (imgin3[vox_index]>imgin1[vox_index]))
                output[vox_index] = atoi(argv[6]);
            else
                output[vox_index] = 0;
        }
        else
            output[vox_index] = 0;

    }


    /* **************************************************************** */
    /* WRITE OUT UNSIGNED 8 BIT TO OUTPUT FILE, AND CLOSE FILES.        */

    unsigned_char_fwrite(argv[8], output, num_vox);
    fclose(fOUT);
}
