#include <stdio.h>
#include <stdlib.h>

main (int argc, char *argv[])
{
    /* 06.07.95 : JM Maisog : Sums voxels in a floating point image.      */

    FILE *fIN;
    float *imgin,sum;
    unsigned char *mask;
    int num_voxX4;
    int num_vox;
    int vox_index;
    int read_status;


    /* **************************************************************** */
    /* VERIFY NUMBER OF COMMAND LINE ARGUMENTS.                         */
    if ((argc!=2) && (argc!=3)) {
        printf("sum_sqrFP: need input floating point image file.\n");
        printf("       Will return the sum of the squares all voxels in this file.\n");
        printf("       If you include a second argument, name of an 8-bit\n");
        printf("       mask, only voxels in the mask will be used in the sum.\n");
        exit(1);
    }

    /* **************************************************************** */
    /* COUNT NUMBER OF BYTES IN FIRST INPUT FILE                        */
    if ((fIN=fopen(argv[1],"r")) == NULL) {
        printf("sum_sqrFP: can't open %s\n",argv[1]);
        exit(2);
    }

    num_voxX4 = 0;
    while (fgetc(fIN)!=EOF)
        num_voxX4++;

    fclose(fIN);
    num_vox = num_voxX4/4;
/*    printf("sum_sqrFP: Number of voxels = %d\n",num_vox); */


    /* **************************************************************** */
    /* ASSIGN MEMORY.                                                   */
    imgin = (float *) malloc(num_vox*sizeof(float));
    if (imgin == NULL) {
        printf("sum_sqrFP: failed to malloc for imgin.\n");
        exit(4);
    }

    /* **************************************************************** */
    /* ASSIGN MEMORY FOR MASK, AND READ IT IN.                          */
    if (argc==3) {
        mask = (unsigned char *) malloc(num_vox*sizeof(unsigned char));
        if (mask == NULL) {
            printf("mask: failed to malloc for sum.\n");
            exit(5);
        }
        if ((fIN=fopen(argv[2],"r"))==NULL) {
            printf("sum_sqrFP: unable to open %s\n",argv[2]);
            exit(5);
        }
        else {
            read_status=fread(mask, sizeof(unsigned char), num_vox, fIN);
            if (read_status<num_vox) {
                printf("sum_sqrFP: only %d voxels read in from %s\n",read_status,argv[2]);
                printf("       Should have been %d.\n",num_vox);
                exit(6);
            }
        }
    }


    /* ************************************************************ */
    /* READ IN INPUT IMAGE.                                         */
    fIN=fopen(argv[1],"r");
    read_status=fread(imgin, sizeof(float), num_vox, fIN);
    fclose(fIN);
    if (read_status != num_vox) {
        printf("%d voxels read in.  Should have been %d voxels.\n",
                                            read_status, num_vox);
        exit(8);
    }

    /* ************************************************************ */
    /* CUMULATIVE SUM VOXELS IN INPUT IMAGE                         */
/*    printf("Summing squares of voxels in %s...\n",argv[1]); */
    sum=0.;
    for (vox_index=0; vox_index<num_vox; vox_index++)
        if(((argc==3)&&mask[vox_index]!=0) || (argc==2))
            sum += (imgin[vox_index]*imgin[vox_index]);

    printf("%g\n",sum);
    free(imgin);
    if (argc==3) free(mask);
}
