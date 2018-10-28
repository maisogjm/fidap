#include <stdio.h>
#include <stdlib.h>

main (int argc, char *argv[])
{
    /* 06.07.95 : JM Maisog : Sums voxels in a floating point image.
    
       10.16.95 : JMM : Modified sumFP to generate meanFP2.  Program
       outputs the mean voxel intensity of a floating point image.
       
       10.20.95: JMM : Modified meanFP2 to generate varFP2.  Program
       outputs the variance of a floating point image.                  */

    FILE *fIN;
    float *imgin,sum,mean;
    unsigned char *mask;
    int num_vox_in_mask;
    int num_voxX4;
    int num_vox;
    int vox_index;
    int read_status;
    float ss;


    /* **************************************************************** */
    /* VERIFY NUMBER OF COMMAND LINE ARGUMENTS.                         */
    if ((argc!=2) && (argc!=3)) {
        printf("varFP2: need input floating point image file.\n");
        printf("        Will return the variance of all voxels in this file.\n");
        printf("        If you include a second argument, name of an 8-bit\n");
        printf("        mask, only voxels in the mask will be used.\n");
        exit(1);
    }

    /* **************************************************************** */
    /* COUNT NUMBER OF BYTES IN FIRST INPUT FILE                        */
    if ((fIN=fopen(argv[1],"r")) == NULL) {
        printf("varFP2: can't open %s\n",argv[1]);
        exit(2);
    }

    num_voxX4 = 0;
    while (fgetc(fIN)!=EOF)
        num_voxX4++;

    fclose(fIN);
    num_vox = num_voxX4/4;


    /* **************************************************************** */
    /* ASSIGN MEMORY.                                                   */
    imgin = (float *) malloc(num_vox*sizeof(float));
    if (imgin == NULL) {
        printf("varFP2: failed to malloc for imgin.\n");
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
            printf("varFP2: unable to open %s\n",argv[2]);
            exit(5);
        }
        else {
            read_status=fread(mask, sizeof(unsigned char), num_vox, fIN);
            fclose(fIN);
            if (read_status<num_vox) {
                printf("varFP2: only %d voxels read in from %s\n",read_status,argv[2]);
                printf("        Should have been %d.\n",num_vox);
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
    read_status=fread(imgin, sizeof(float), num_vox, fIN);
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
            sum += imgin[vox_index];

    mean=sum/(float)num_vox_in_mask;

    /* ************************************************************ */
    /* CALCULATE VARIANCE.                                          */

    ss=0.;
    for (vox_index=0; vox_index<num_vox; vox_index++)
        if(((argc==3)&&mask[vox_index]!=0) || (argc==2))
            ss += ((imgin[vox_index]-mean)*(imgin[vox_index]-mean));

    printf("%g\n",ss/(float)(num_vox_in_mask-1));
    free(imgin);
    free(mask);
}
