#include <stdio.h>
#include <stdlib.h>

float *float_malloc(int num_items);
int read_unsigned_char_data(char *filename, unsigned char *array, int num_items);
unsigned char *unsigned_char_malloc(int num_items);
int unsigned_char_fwrite(char *filename, unsigned char *array, int num_items);

main (int argc, char *argv[])
{
    /* 12.23.93 : JM Maisog : Finds mean of floating point images.        */

    FILE *fIN, *fOUT;
    unsigned char *data8bit;
    float *sum, *denominator;
    int cmd_lin_arg;
    int num_vox;
    int vox_index;
    int read_status;
    int write_status;


    /* **************************************************************** */
    /* VERIFY NUMBER OF COMMAND LINE ARGUMENTS.                         */
    if (argc < 3) {
        printf("softmeanB: need input file(s) and output file.\n");
        exit(1);
    }

    /* **************************************************************** */
    /* COUNT NUMBER OF BYTES IN FIRST INPUT FILE                        */
    if ((fIN=fopen(argv[1],"r")) == NULL) {
        printf("softmeanB: can't open %s\n",argv[1]);
        exit(2);
    }

    num_vox = 0;
    while (fgetc(fIN)!=EOF)
        num_vox++;

    fclose(fIN);
    printf("softmeanB: Number of voxels = %d\n",num_vox);


    /* **************************************************************** */
    /* ASSIGN MEMORY.                                                   */
    data8bit = unsigned_char_malloc(num_vox);
    sum = float_malloc(num_vox);
    denominator = float_malloc(num_vox);

    for (vox_index=0; vox_index<num_vox; vox_index++) {
	denominator[vox_index]=0.;
        sum[vox_index] = 0.;
    }

    /* **************************************************************** */
    /* LOOP OVER ALL INPUT FILES....                                    */

    for (cmd_lin_arg = 1; cmd_lin_arg<=argc-2 ;cmd_lin_arg++) {

        /* ************************************************************ */
        /* READ IN VOXELS, AND CUMULATIVE SUM INPUT IMAGE INTO *sum.    */
        read_unsigned_char_data(argv[cmd_lin_arg], data8bit, num_vox);

        printf("Summing data from %s...\n",argv[cmd_lin_arg]);
        for (vox_index=0; vox_index<num_vox; vox_index++) {
            sum[vox_index] += (float)data8bit[vox_index];
            if (data8bit[vox_index]!=0) denominator[vox_index]+=1.;
        }

    }

    /* **************************************************************** */
    /* FIND MEAN BY DIVIDING SUM BY NUMBER OF INPUT IMAGES, THEN WRITE  */
    /* TO OUTPUT FILE, AND CLOSE FILES.                                 */
    for (vox_index=0; vox_index<num_vox; vox_index++) {
	if (denominator[vox_index]!=0) {
            sum[vox_index] /= denominator[vox_index];
	}
    }

    unsigned_char_fwrite(argv[argc-1], data8bit, num_vox);
    free(data8bit);
    free(denominator);
    free(sum);
}
