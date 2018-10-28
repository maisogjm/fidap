#include "analyze6.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

unsigned char *unsigned_char_malloc(int num_items);
int read_unsigned_char_data(char *filename, unsigned char *array, int num_items);
int unsigned_char_fwrite(char *filename, unsigned char *array, int num_items);

main (int argc, char *argv[])
{
    /* 06.07.96 : JM Maisog : Performs XOR operation on 8-bit maps.     */

    FILE *fIN;
    unsigned char *imgin1, *imgin2, *xor;
    int cmd_lin_arg;
    int num_vox;
    int vox_index;

    /* **************************************************************** */
    /* VERIFY NUMBER OF COMMAND LINE ARGUMENTS.                         */
    if (argc != 4) {
        printf("XOR8bit: need two input files and name for output file.\n");
        printf("         Program will perform EXCLUSIVE OR.\n");
        exit(1);
    }

    /* **************************************************************** */
    /* COUNT NUMBER OF BYTES IN FIRST INPUT FILE                        */
    if ((fIN=fopen(argv[1],"r")) == NULL) {
        printf("XOR8bit: can't open %s\n",argv[1]);
        exit(2);
    }
    num_vox= 0;
    while (fgetc(fIN)!=EOF)
        num_vox++;
    fclose(fIN);

    /* **************************************************************** */
    /* ASSIGN MEMORY, INITIALIZE xor[*]                                 */
    imgin1 = unsigned_char_malloc(num_vox);
    imgin2 = unsigned_char_malloc(num_vox);
    xor    = unsigned_char_malloc(num_vox);

    /* **************************************************************** */
    /* READ IN INPUT DATA FILES.                                        */
    read_unsigned_char_data(argv[1], imgin1, num_vox);
    read_unsigned_char_data(argv[2], imgin2, num_vox);
    for (vox_index=0; vox_index<num_vox; vox_index++)
        xor[vox_index] = (((imgin1[vox_index]!=0) || (imgin2[vox_index]==0))
			    || (imgin1[vox_index]==0) || (imgin2[vox_index]!=0))
								? 1 : 0;

    /* **************************************************************** */
    /* WRITE OUT UNSIGNED 8 BIT TO OUTPUT FILE, AND CLOSE FILES.        */
    unsigned_char_fwrite(argv[3], xor, num_vox);
    free(imgin1);
    free(imgin2);
    free(xor);
}
