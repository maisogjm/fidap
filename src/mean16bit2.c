#include <stdio.h>
#include <stdlib.h>

main (int argc, char *argv[])
{
    /* 06.07.95 : JM Maisog : Sums voxels in a floating point image.
    
       10.16.95 : JMM : Modified sumFP to generate mean16bit2.  Program
       outputs the mean voxel intensity of a floating point image.      */

    FILE *fIN;
    signed short *imgin;
    float sum;
    unsigned char *mask;
    int num_vox_in_mask;
    int num_voxX2;
    int num_vox;
    int vox_index;
    int read_status;


    /* **************************************************************** */
    /* VERIFY NUMBER OF COMMAND LINE ARGUMENTS.                         */
    if ((argc!=2) && (argc!=3)) {
        printf("mean16bit2: need input 16-bit image file.\n");
        printf("         Will return the mean of all voxels in this file.\n");
        printf("         If you include a second argument, name of an 8-bit\n");
        printf("         mask, only voxels in the mask will be used in the mean.\n");
        exit(1);
    }

    /* **************************************************************** */
    /* COUNT NUMBER OF BYTES IN FIRST INPUT FILE                        */
    if ((fIN=fopen(argv[1],"r")) == NULL) {
        printf("mean16bit2: can't open %s\n",argv[1]);
        exit(2);
    }

    num_voxX2 = 0;
    while (fgetc(fIN)!=EOF)
        num_voxX2++;

    fclose(fIN);
    num_vox = num_voxX2/2;

    /* **************************************************************** */
    /* ASSIGN MEMORY.                                                   */
    imgin = (signed short *) malloc(num_vox*sizeof(signed short));
    if (imgin == NULL) {
        printf("mean16bit2: failed to malloc for imgin.\n");
        exit(4);
    }

    /* **************************************************************** */
    /* ASSIGN MEMORY FOR MASK, AND READ IT IN.                          */
    /* Count number of nonzero voxels in mask.                          */
    if (argc==3) {
        mask = (unsigned char *) malloc(num_vox*sizeof(unsigned char));
        if (mask == NULL) {
            printf("mask: failed to malloc for sum.\n");
            exit(5);
        }
        if ((fIN=fopen(argv[2],"r"))==NULL) {
            printf("mean16bit2: unable to open %s\n",argv[2]);
            exit(5);
        }
        else {
            read_status=fread(mask, sizeof(unsigned char), num_vox, fIN);
            fclose(fIN);
            if (read_status<num_vox) {
                printf("mean16bit2: only %d voxels read in from %s\n",read_status,argv[2]);
                printf("         Should have been %d.\n",num_vox);
                exit(6);
            }
            num_vox_in_mask=0;
            for (vox_index=0; vox_index<num_vox; vox_index++)
                if (mask[vox_index]!=0)
                    num_vox_in_mask++;
        }
    }
    else num_vox_in_mask=num_vox;

    /* ************************************************************ */
    /* READ IN INPUT IMAGE.                                         */
    fIN=fopen(argv[1],"r");
    read_status=fread(imgin, sizeof(signed short), num_vox, fIN);
    fclose(fIN);
    if (read_status != num_vox) {
        printf("%d voxels read in.  Should have been %d voxels.\n",
                                            read_status, num_vox);
        exit(8);
    }

    /* ************************************************************ */
    /* CUMULATIVE SUM VOXELS IN INPUT IMAGE                         */
    sum=0.;
    for (vox_index=0; vox_index<num_vox; vox_index++)
        if(((argc==3)&&mask[vox_index]!=0) || (argc==2))
            sum += (float) imgin[vox_index];

    printf("%g\n",sum/(float)num_vox_in_mask);
    free(imgin);
    if (argc==3) free(mask);
}
