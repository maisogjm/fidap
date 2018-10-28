#include <stdio.h>
#include <stdlib.h>
#include <math.h>

float *float_malloc(int num_items);
int read_unsigned_short_int_data(char *filename, unsigned short int *array, int num_items);
unsigned short *unsigned_short_malloc(int num_items);
int unsigned_short_fwrite(char *filename, unsigned short *array, int num_items);


main (int argc, char *argv[])
{
    /* 09.16.96 : JM Maisog : Finds mean of 16-bit images.        */

    FILE *fIN, *fOUT;
    unsigned short *data16bit;
    float *sum, *denominator;
    int cmd_lin_arg;
    int num_vox;
    int vox_index;
    int read_status;
    int write_status;


    /* **************************************************************** */
    /* VERIFY NUMBER OF COMMAND LINE ARGUMENTS.                         */
    if (argc < 3) {
        printf("softmean16bitB: need input file(s) and output file.\n");
        exit(1);
    }

    /* **************************************************************** */
    /* COUNT NUMBER OF BYTES IN FIRST INPUT FILE                        */
    if ((fIN=fopen(argv[1],"r")) == NULL) {
        printf("softmean16bitB: can't open %s\n",argv[1]);
        exit(2);
    }

    num_vox = 0;
    while (fgetc(fIN)!=EOF)
        num_vox++;

    fclose(fIN);
    num_vox/=2;
    printf("softmean16bitB: Number of voxels = %d\n",num_vox);


    /* **************************************************************** */
    /* ASSIGN MEMORY.                                                   */
    data16bit = unsigned_short_malloc(num_vox);
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
        read_unsigned_short_int_data(argv[cmd_lin_arg], data16bit, num_vox);

        printf("Summing data from %s...\n",argv[cmd_lin_arg]);
        for (vox_index=0; vox_index<num_vox; vox_index++) {
            sum[vox_index] += (float)data16bit[vox_index];
            if (data16bit[vox_index]!=0) denominator[vox_index]+=1.;
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

    for (vox_index=0; vox_index<num_vox; vox_index++)
        data16bit[vox_index] = (unsigned short)floor(sum[vox_index]+0.5);

    unsigned_short_fwrite(argv[argc-1], data16bit, num_vox);
    free(data16bit);
    free(denominator);
    free(sum);
}
