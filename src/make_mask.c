#include <stdio.h>
#include <stdlib.h>
#include "analyze6.h"


                       /* FUNCTION DECLARATIONS.                                */
struct  analyze_struct readhdr(char *argv[]);
char **char_malloc2(int num_items1, int num_items2);
unsigned char *unsigned_char_malloc(int num_items);
float *float_malloc(int num_items);
signed short *signed_short_malloc(int num_items);
int read_unsigned_char_data(char *filename, unsigned char *array, int num_items);
int read_signed_short_int_data(char *filename, signed short int *array, int num_items);
int read_float_data(char *filename, float *array, int num_items);
int unsigned_char_fwrite(char *filename, unsigned char *array, int num_items);

/* **************************************************************************** */
/* PROGRAM make_mask.c                                                          */
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
/*            11.15.96 : JMM : Modified fmri_mask to generate make_mask         */
/*                       Program now runs on 8-bit, 16-bit, or floating point   */
/*                       data and doesn't load all data into memory at once.    */

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
    int vox_index;          /* Voxel index.  Range: [0,vox_index-1]             */
    int xdim,ydim,zdim;     /* Image dimensions.                                */
    int num_bits_per_pixel; /* Bits per pixel in image.                         */
    unsigned char *Bbuffer; /* Buffer to hold 8-bit image as read in            */
                            /* from disk.                                       */
    signed short int *Wbuffer;/* Buffer to hold 16-bit image as read in         */
                            /* from disk.                                       */
    float *Fbuffer;         /* Buffer to hold floating point image as           */
                            /* read in from disk.                               */
    unsigned char *zero_vox;/* Set to 1 if any voxel in any of the scans was 0. */

                            /* INTENSITY SCALING & GLOBAL MEAN NORMALIZATION    */
    FILE *fmean;            /* Pointer for output mean file.                    */
    float *mean;            /* Pointer to mean images.                          */
    unsigned char *mask;    /* Unsigned char copy of mask, for writing to       */
                            /* disk.                                            */
    float max_vox_val;      /* Maximum voxel value in the mean image.           */
    float threshold;        /* Threshold used to generate a mask from the       */
                            /* average image.                                   */
    float prob_thresh;      /* Probability threshold for significance.          */
    char mask_name[80];     /* Name of output mask.                             */
    int num_mask_vox;       /* Number of non-zero voxels in mask.               */

    int x,y,z;              /* Indices into ANALYZE format image arrays.        */

/* ---------------------------------------------------------------------------- */
/*  CHECK TO SEE IF CORRECT NUMBER OF ARGUMENTS WERE PASSED.                    */

    if (argc != 4) {
        printf("make_mask: need header file input textfile, and percentage threshold.\n\n");

        printf("make_mask will generate an 8-bit binary mask, binary in the\n");
        printf("sense that voxels will be set to either 1 or 0.\n\n");

        printf("Input scans can be 8-bit or 16-bit integer, or single-precision\n");
        printf("floating point.\n\n");

        printf("The mask is created in the following procedure: a mean scan is\n");
        printf("made using a cutoff specified by the user, generating a preliminary\n");
        printf("mask.  Then, all voxels in this preliminary mask for which even\n");
        printf("one of the image volumes has a zero value are eliminated from\n");
        printf("the mask, generating the final mask.\n");

        printf("\nThe input textfile should contain the pathnames\n");
        printf("(either absolute or relative) of the input scans.\n") ;

        exit(1);
    }

/* ---------------------------------------------------------------------------- */
/*  READ HEADER FILE; ASSUME IMAGES ARE SAME SIZE, SO USE SAME HEADER FILE      */
/*  FOR ALL IMAGES.                                                             */

    printf("Reading AVW header...\n");fflush(NULL);
    img_info = readhdr(argv+1);

    xdim = img_info.dime.dim[1];
    ydim = img_info.dime.dim[2];
    zdim = img_info.dime.dim[3];
    num_bits_per_pixel=img_info.dime.bitpix;
    num_vox = xdim*ydim*zdim;
    printf("xdim = %d\n",xdim);fflush(NULL);
    printf("ydim = %d\n",ydim);fflush(NULL);
    printf("zdim = %d\n",zdim);fflush(NULL);

/* ---------------------------------------------------------------------------- */
/*  COUNT NUMBER OF ENTRIES IN INPUT TEXT FILE.                                 */

    printf("Reading input text file...\n");fflush(NULL);
    if ((fimgin = fopen(argv[2],"r")) == NULL) {
        printf ("make_mask: can't open %s\n",argv[2]);
        exit(4);
    }
    total_n = 1;
    fscanf(fimgin,"%s",dummy);
    while(fscanf(fimgin,"%s",dummy)!=EOF)
        total_n++;
    printf("There are %d image volumes.\n",total_n);fflush(NULL);

/* ---------------------------------------------------------------------------- */
/*  READ IN FILENAMES AND GLOBAL MEAN WEIGHTS IN INPUT TEXT FILE.               */

    img_name = char_malloc2(total_n, 1024);

    rewind(fimgin);
    for (n_index=0; n_index<=total_n-1; n_index++)
        fscanf(fimgin,"%s", img_name[n_index]);
    fclose(fimgin);

/* ---------------------------------------------------------------------------- */
/*  SHOW USER INFORMATION READ IN SO FAR.                                       */

/*    printf("\nThese are the %d image volumes:\n", total_n);
    for (n_index=0; n_index<=total_n-1; n_index++)
        printf("  %20s\n", img_name[n_index]);fflush(NULL);
*/

    printf("There are %d bits per pixel.\n",num_bits_per_pixel);
    printf("Image volumes have dimensions:\n");
    printf("  xdim = %d\n  ydim = %d\n  zdim = %d\n",xdim,ydim,zdim);

    threshold=atof(argv[3]);
    printf("\nThreshold = %3.2f %%.\n\n",threshold);fflush(NULL);

/* ---------------------------------------------------------------------------- */
/* INITIALIZE mean[] AND zero_vox[].                                            */

    zero_vox=unsigned_char_malloc(num_vox);
    mean=float_malloc(num_vox);
    for (vox_index=0; vox_index<num_vox; vox_index++) {
        mean[vox_index] = 0.;
        zero_vox[vox_index] = 0.;
    }
    mask=unsigned_char_malloc(num_vox);

/* ---------------------------------------------------------------------------- */
/*  READ IN IMAGE DATA, CUMULATE INTO mean[*].                                  */
/*  SET zero_vox[vox_index] TO ZERO IF ANY VOXEL AT vox_index IS ZERO.          */

    if (num_bits_per_pixel==8)
        Bbuffer = unsigned_char_malloc(num_vox);
    else if (num_bits_per_pixel==16)
        Wbuffer = signed_short_malloc(num_vox);
    else if (num_bits_per_pixel==32)
        Fbuffer = float_malloc(num_vox);


    printf("make_mask: reading in data....\n");fflush(NULL);
    for (n_index=0; n_index<total_n; n_index++) {
/*        printf("make_mask: reading in data from %s....\n",img_name[n_index]); */
        if (num_bits_per_pixel==8) {
            read_unsigned_char_data(img_name[n_index], Bbuffer, num_vox);
            for (vox_index=0; vox_index<num_vox; vox_index++) {
                mean[vox_index]+=(float)Bbuffer[vox_index];
                zero_vox[vox_index]=(Bbuffer[vox_index]==0)?1:zero_vox[vox_index];
            }
        }
        else if (num_bits_per_pixel==16) {
            read_signed_short_int_data(img_name[n_index], Wbuffer, num_vox);
            for (vox_index=0; vox_index<num_vox; vox_index++) {
                mean[vox_index]+=(float)Wbuffer[vox_index];
                zero_vox[vox_index]=(Wbuffer[vox_index]==0)?1:zero_vox[vox_index];
            }
        }
        else if (num_bits_per_pixel==32) {
            read_float_data(img_name[n_index], Fbuffer, num_vox);
            for (vox_index=0; vox_index<num_vox; vox_index++) {
                mean[vox_index]+=(float)Fbuffer[vox_index];
                zero_vox[vox_index]=(Fbuffer[vox_index]==0.)?1:zero_vox[vox_index];
            }
        }
    }
    for (vox_index=0; vox_index<num_vox; vox_index++)
        mean[vox_index]/=(float)total_n;

    if (num_bits_per_pixel==8)
        free(Bbuffer);
    else if (num_bits_per_pixel==16)
        free(Wbuffer);
    else if (num_bits_per_pixel==32)
        free(Fbuffer);

/* ---------------------------------------------------------------------------- */
/* DETERMINE MAXIMUM AND MINIMUM PIXEL VALUES OF MEAN IMAGE.                    */

    max_vox_val = 0.;
    for (vox_index=0; vox_index<num_vox; vox_index++)
        max_vox_val = (mean[vox_index] > max_vox_val) ? mean[vox_index] : max_vox_val;
    printf("Threshold will be %f.\n\n",(max_vox_val*threshold/100.));fflush(NULL);

/* ---------------------------------------------------------------------------- */
/* THRESHOLD MEAN SCAN USING threshold, MAKING MASK.                           */

    for (vox_index=0; vox_index<num_vox; vox_index++)
        mask[vox_index] = (mean[vox_index]>=(max_vox_val*threshold/100.))?1:0;

/* ---------------------------------------------------------------------------- */
/* FURTHER CLEANING UP OF MASK.                                                 */
/* If any of the original data is equal to 0, set corresponding entry           */
/* in mask to zero.                                                             */

    for (vox_index=0; vox_index<num_vox; vox_index++)
        if (zero_vox[vox_index]!=0) mask[vox_index]=0;

/* ---------------------------------------------------------------------------- */
/* COUNT NUMBER OF NON-ZERO VOXELS IN mask.                                     */

    num_mask_vox = 0.;
    for (vox_index=0; vox_index<num_vox; vox_index++)
        if (mask[vox_index] != 0.)
            num_mask_vox++;

/* ---------------------------------------------------------------------------- */
/*  WRITE mask TO DISK.                                                         */

    sprintf(mask_name,"mask_%d.img\0",(int) threshold);
    printf("There are %d voxels set in %s.\n\n",num_mask_vox,mask_name);fflush(NULL);
    unsigned_char_fwrite(mask_name, mask, num_vox);

/* ---------------------------------------------------------------------------- */
/*  FREE MEMORY.                                                                */

    printf("\nmake_mask: releasing memory...\n");
    free(img_name);
    free(mean);
    free(mask);
}
