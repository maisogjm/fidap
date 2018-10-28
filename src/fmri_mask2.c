#include <stdio.h>
#include <stdlib.h>
#include "analyze6.h"
#include <math.h>

                       /* FUNCTION DECLARATIONS.                                */
struct  analyze_struct readhdr(char *argv[]);
float betai(float a, float b, float x);

/* **************************************************************************** */
/* PROGRAM fmri_mask.c                                                          */
/*                                                                              */
/* ANSI C code completed 08/24/94 by Jose' Ma. Maisog, M.D.                     */
/*                                                                              */
/*                                                                              */
/* Building 10, Rm 4C110                                                        */
/* NIH/NIMH/LPP/SFBI                                                            */
/* 9000 Rockville Pike                                                          */
/* Bethesda, Maryland  20892                                                    */
/*                                                                              */
/* Tel: (301) 402-0416                                                          */
/* Fax: (301) 402-0921                                                          */
/* email: joem@alw.nih.gov                                                      */
/*                                                                              */
/* Program reads in names of image volumes from a text file, then reads in      */
/* those image volumes, and generates the mean image of        those volumes.   */
/* Mean image is thresholded with a cutoff specified by the user, generating a  */
/* preliminary mask.  All voxels in this preliminary mask for which one or more */
/* original images have a zero value are eliminated from the mask, to avoid     */
/* missing cells problem.  This final mask is then written to disk.             */
/*                                                                              */
/* HISTORY: 11.24.93 : JMM : started modifying fmri_pairedt_ratio.c to generate */
/*                       ultra_anova.c.  ultra_anova reads in 16-bit images and */
/*                       performs an N-WAY ANOVA.                               */
/*            01.19.94 : JMM : added code to ultra_anova to output in addition  */
/*                       the error SS image.                                    */
/*            03.04.94 : JMM : final debugging of ultra_anova.                  */
/*            03.07.94 : JMM : final touches on output.                         */
/*            08.24.94 : JMM : Modified ultra_anova to just generate a mask,    */
/*                           new program is fmri_mask.                          */

main(int argc, char *argv[])
{
                            /* READING IN RAW DATA.                             */
    struct analyze_struct   /* Header information variable.                     */
                img_info;   /*                                                  */
    FILE *fimgin;           /* Pointer for first input image file.              */
    char choice;            /* For keyboard input to yes/no question.           */
    char dummy[1024];         /* Dummy pointer for counting number of entries     */
                            /* in a text file.                                  */
    int total_n;            /* Total number of images.                          */
    int n_index;            /* Index indicating which image in group.           */
                            /* Range: [0,total_n-1]                             */
    char **img_name;        /* Pointer to pointers to input file names.         */
    int num_vox;            /* Number of voxels in an image = xdim*ydim*zdim.   */
    int vox_index;          /* Index into imgin[*][*][], extracting a specific  */
                            /* voxel value.  Range: [0,vox_index-1]             */
    int xdim,ydim,zdim;     /* Image dimensions.                                */
    signed short *temp;     /* Temporary buffer for reading in images.          */
    signed short **imgin;   /* Pointer to pointers to input image data.         */

                            /* INTENSITY SCALING & GLOBAL MEAN NORMALIZATION    */
    FILE *fmean;            /* Pointer for output mean file.                    */
    float *mean;            /* Pointer to mean images.                          */
    unsigned char *outmask; /* Unsigned char copy of mask, for writing to       */
                            /* disk.                                            */
    float max_vox_val;      /* Maximum voxel value in the mean image.           */
    float min_vox_val;      /* Minimum voxel value in the mean image.           */
    float threshold;        /* Threshold used to generate a mask from the       */
                            /* average image.                                   */
    float prob_thresh;      /* Probability threshold for significance.          */
    float *mask;            /* Pointer to mask image.                           */
    char outmask_name[80];  /* Name of output mask image.                       */
    int num_mask_vox;       /* Number of non-zero voxels in *mask.              */

    int x,y,z;              /* Indices into ANALYZE format image arrays.        */

/* ---------------------------------------------------------------------------- */
/*  CHECK TO SEE IF CORRECT NUMBER OF ARGUMENTS WERE PASSED.                    */

    if (argc != 3) {
        printf("fmri_mask: need header file and input textfile\n\n");

        printf("fmri_mask will generate an 8-bit binary mask, binary in the\n");
        printf("sense that voxels will be set to either 1 or 0.\n\n");

        printf("This is done in the following procedure: a mean image is made\n");
        printf("of all the input 16-bit images.  This mean image is then\n");
        printf("thresholded using a cutoff specified by the user, generating\n");
        printf("a preliminary mask.  Then, all voxels in this preliminary mask\n");
        printf("for which even one of the image volumes has a zero value are\n");
        printf("eliminated from the mask, generating the final mask.\n");

        printf("\nThe input textfile should contain the pathnames\n");
        printf("(either absolute or relative) of the input images.\n") ;

        exit(1);
    }

/* ---------------------------------------------------------------------------- */
/*  READ HEADER FILE; ASSUME IMAGES ARE SAME SIZE, SO USE SAME HEADER FILE      */
/*  FOR ALL IMAGES.                                                             */

    img_info = readhdr(argv+1);

    xdim = img_info.dime.dim[1];
    ydim = img_info.dime.dim[2];
    zdim = img_info.dime.dim[3];
    printf("xdim = %d, ydim = %d, zdim = %d\n",xdim,ydim,zdim);fflush(NULL);

    num_vox = xdim*ydim*zdim;

/* ---------------------------------------------------------------------------- */
/*  COUNT NUMBER OF ENTRIES IN INPUT TEXT FILE.                                 */

    if ((fimgin = fopen(argv[2],"r")) == NULL) {
        printf ("fmri_mask: can't open %s\n",argv[2]);
        exit(4);
    }
    printf("Just opened %s.\n",argv[2]);fflush(NULL);
    total_n = 1;
    fscanf(fimgin,"%s",dummy);
    while(fscanf(fimgin,"%s",dummy)!=EOF)
        total_n++;
    printf("There are %d images.\n",total_n);fflush(NULL);

/* ---------------------------------------------------------------------------- */
/*  READ IN FILENAMES AND GLOBAL MEAN WEIGHTS IN INPUT TEXT FILE.               */

    img_name = (char **) malloc(total_n*sizeof(char *));
    if (img_name==NULL) {
        printf("fmri_mask: error malloc'ing for img_name.\n");
        exit(6);
    }

    rewind(fimgin);
    for (n_index=0; n_index<=total_n-1; n_index++) {
        img_name[n_index] = (char *) malloc(1024*sizeof(char));
        if (img_name[n_index]==NULL) {
            printf("fmri_mask: error malloc'ing for img_name.\n");
            exit(7);
        }
        fscanf(fimgin,"%s", img_name[n_index]);
    }
    fclose(fimgin);

/* ---------------------------------------------------------------------------- */
/*  SHOW USER INFORMATION READ IN SO FAR.                                       */
/*  Allow user to quit program if incorrect information has been input.         */

    printf("\nThese are the %d image volumes:\n", total_n);
    for (n_index=0; n_index<=total_n-1; n_index++)
        printf("  %20s\n", img_name[n_index]);

    printf("\nImages have dimensions:\n");
    printf("  xdim = %d\n  ydim = %d\n  zdim = %d\n",xdim,ydim,zdim);

    printf("\nIt is assumed that the input images are registered.  A mean\n");
    printf("image will be calculated.  A mask will be made from this mean\n");
    printf("image by thresholding it.  Below, you will need to specify the\n");
    printf("number T.  All voxels in the mean image which have pixel\n");
    printf("intensity T %% the maximum voxel intensity in the mean image or\n");
    printf("greater will be set to 1 in the mask.\n");
    printf("E.g., if you want a 10 %% threshold, type in <10>.\n");
    printf("This mask will be used to calculate the global mean\n");
    printf("intensities of the individual images.\n");
    printf("\nENTER the threshold T for the mask: ");
    scanf("%f",&threshold);
    printf("\nThreshold = %3.2f %%.\n",threshold);
    printf("\n");

/* ---------------------------------------------------------------------------- */
/*  READ IN IMAGES DATA.                                                        */

    temp = (signed short *) malloc(num_vox*sizeof(signed short));
    imgin = (signed short **) malloc(total_n*sizeof(signed short *));

    if (imgin==NULL) {
        printf("fmri_mask: error malloc'ing for imgin.\n");
        exit(21);
    }

    for (n_index=0; n_index<total_n; n_index++) {
        imgin[n_index] = (signed short *) malloc(num_vox*sizeof(signed short));
        if (imgin[n_index]==NULL) {
            printf("fmri_mask: error malloc'ing for imgin.\n");
            exit(22);
        }
        printf("fmri_mask: reading in data from %s....\n",img_name[n_index]);
        if ((fimgin=fopen(img_name[n_index],"r")) == NULL) {
            printf("fmri_mask: can't open %s\n",img_name[n_index]);
            exit(23);
        }
        fread(temp, sizeof(signed short), num_vox, fimgin);
        for (vox_index=0; vox_index<num_vox; vox_index++)
            imgin[n_index][vox_index]=temp[vox_index];
    }
    free(temp);

/* ---------------------------------------------------------------------------- */
/*  ALLOCATE MEMORY FOR OUTPUT IMAGES.                                          */

    /* ------------------------------------------------------------------------ */
    /* ALLOCATE MEMORY FOR mean AND outmask.                                    */

    mean = (float *) malloc(num_vox*sizeof(float));
    if (mean==NULL) {
        printf("fmri_mask: error malloc'ing for mean.\n");
        exit(24);
    }
    outmask = (unsigned char *) malloc(num_vox*sizeof(unsigned char));
    if (outmask==NULL) {
        printf("fmri_mask: error malloc'ing for outmask.\n");
        exit(26);
    }

    /* ------------------------------------------------------------------------ */
    /* CALCULATE MEAN IMAGE.                                                    */

        /* -------------------------------------------------------------------- */
        /* INITIALIZE MEAN IMAGE.                                               */

        for (vox_index=0; vox_index<num_vox; vox_index++)
            mean[vox_index] = 0.;

        /* -------------------------------------------------------------------- */
        /* CALCULATE MEAN IMAGE.                                                */
        /* This mean image is used to generate a mask, which in turn is         */
        /* used to calculate global means.  Later, the mean is recalculated     */
        /* using the normalized values, and then used to calculate the t-test.  */

        printf("fmri_mask: calculating mean image...\n");
        for (n_index=0; n_index<total_n; n_index++)
            for (vox_index=0; vox_index<num_vox; vox_index++)
                mean[vox_index] += (float)fabs((double)imgin[n_index][vox_index]);

        for (vox_index=0; vox_index<num_vox; vox_index++)
            mean[vox_index] /= ((float) total_n);

    /* ------------------------------------------------------------------------ */
    /* MAKE MASK FROM mean.                                                     */

    printf("fmri_mask: making mask image.\n");

        /* -------------------------------------------------------------------- */
        /* DETERMINE MAXIMUM AND MINIMUM PIXEL VALUES OF MEAN IMAGE.            */

        max_vox_val = 0.;
        min_vox_val = 100000.;
        for (vox_index=0; vox_index<num_vox; vox_index++) {
            max_vox_val = (mean[vox_index] > max_vox_val)
                                        ? mean[vox_index] : max_vox_val;
            min_vox_val = (mean[vox_index] < min_vox_val)
                                        ? mean[vox_index] : min_vox_val;
        }
	printf("Threshold will be %f.\n\n",(max_vox_val*threshold/100.));

        /* -------------------------------------------------------------------- */
        /* THRESHOLD MEAN IMAGE USING threshold, MAKING *mask.                  */

        mask = (float *) malloc(num_vox*sizeof(float));
        if (mask==NULL) {
            printf("fmri_mask: error malloc'ing for mask.\n");
            exit(27);
        }
        for (vox_index=0; vox_index<num_vox; vox_index++)
            mask[vox_index] =
                (mean[vox_index] >= (max_vox_val*threshold/100.)) ? 1. : 0.;

        /* -------------------------------------------------------------------- */
        /* FURTHER CLEANING UP OF MASK.                                         */
        /* If any of the original data is equal to 0, set corresponding entry   */
        /* in mask to zero.                                                     */

        for (n_index=0; n_index<total_n; n_index++)
            for (vox_index=0; vox_index<num_vox; vox_index++)
                if (imgin[n_index][vox_index]==0) mask[vox_index]=0.;

	/*
        for (n_index=0; n_index<total_n; n_index++)
            for (vox_index=0; vox_index<num_vox; vox_index++)
                if (imgin[n_index][vox_index]<=0) mask[vox_index]=0.;
	*/

        /* -------------------------------------------------------------------- */
        /* MASK OUT mean.                                                       */

        for (vox_index=0; vox_index<num_vox; vox_index++)
            mean[vox_index] *= mask[vox_index];

        /* -------------------------------------------------------------------- */
        /* SCALE float mean TO unsigned char outmask.                           */

        for (vox_index=0; vox_index<num_vox; vox_index++)
            outmask[vox_index] = (unsigned char) mask[vox_index];

        /* -------------------------------------------------------------------- */
        /*  WRITE outmask TO DISK.                                              */

        sprintf(outmask_name,"mask_%d.img\0",(int) threshold);
        if ((fmean = fopen(outmask_name,"w")) == NULL) {
            printf("fmri_mask: can't open %s.img.\n",outmask_name);
            exit(28);
        }
        printf("fmri_mask: writing mask image to %s.\n",outmask_name);
        fwrite(outmask,sizeof(unsigned char),num_vox,fmean);
        fclose(fmean);

        /* -------------------------------------------------------------------- */
        /* COUNT NUMBER OF NON-ZERO VOXELS IN *mask.                            */

        num_mask_vox = 0.;
        for (vox_index=0; vox_index<num_vox; vox_index++)
            if (mask[vox_index] != 0.)
                num_mask_vox++;

/* ---------------------------------------------------------------------------- */
/*  FREE MEMORY.                                                                */

    printf("\nfmri_mask: releasing memory...\n");
    free(img_name);
    free(imgin);
    free(mean);
    free(outmask);
    free(mask);
}
