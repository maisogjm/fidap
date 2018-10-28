#include "analyze6.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

unsigned char *unsigned_char_malloc(int num_items);
int read_unsigned_char_data(char *filename, unsigned char *array, int num_items);
int unsigned_char_fwrite(char *filename, unsigned char *array, int num_items);

main (int argc, char *argv[])
{
    /* 06.07.96 : JM Maisog : Performs AND operation on 8-bit maps.     */

    FILE *fIN;
    unsigned char *imgin, *product;
    int cmd_lin_arg;
    int num_vox;
    int vox_index;

    /* **************************************************************** */
    /* VERIFY NUMBER OF COMMAND LINE ARGUMENTS.                         */
    if (argc < 4) {
        printf("AND8bit: need input files (at least 2) and output file.\n");
        exit(1);
    }

    /* **************************************************************** */
    /* COUNT NUMBER OF BYTES IN FIRST INPUT FILE                        */
    if ((fIN=fopen(argv[1],"r")) == NULL) {
        printf("AND8bit: can't open %s\n",argv[1]);
        exit(2);
    }
    num_vox= 0;
    while (fgetc(fIN)!=EOF)
        num_vox++;
    fclose(fIN);

    /* **************************************************************** */
    /* ASSIGN MEMORY, INITIALIZE product[*]                             */
    imgin = unsigned_char_malloc(num_vox);
    product = unsigned_char_malloc(num_vox);
    for (vox_index=0; vox_index<num_vox; vox_index++)
        product[vox_index] = 1.;

    /* **************************************************************** */
    /* LOOP OVER ALL INPUT FILES....                                    */
    for (cmd_lin_arg = 1; cmd_lin_arg<=argc-2 ;cmd_lin_arg++) {
        read_unsigned_char_data(argv[cmd_lin_arg], imgin, num_vox);
        for (vox_index=0; vox_index<num_vox; vox_index++)
                product[vox_index] = (product[vox_index]&&imgin[vox_index]) ? 1 : 0;

    }

    /* **************************************************************** */
    /* WRITE OUT UNSIGNED 8 BIT TO OUTPUT FILE, AND CLOSE FILES.        */
    unsigned_char_fwrite(argv[argc-1], product, num_vox);
}
