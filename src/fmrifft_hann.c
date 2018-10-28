#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "analyze6.h"

                       /* FUNCTION DECLARATIONS.                                */

void four1(float data[], unsigned long nn, int isign);
struct  analyze_struct readhdr(char *argv[]);

/* **************************************************************************** */
/* Program fmrifft_hann.c                                                       */
/* ANSI C Code by Jose' Ma. Maisog.                                             */
/* Theory developed by Jack Van Horn.                                           */
/*                                                                              */
/*                                                                              */
/*                                                                              */
/* For correspondence, please refer to Dr. Van Horn.                            */
/* Program reads in sequence of time series images and performs FFT analysis.   */
/*                                                                              */

main(int argc, char *argv[])
{
    FILE *ftxt;             /* Pointer for input text file containing names of  */
                            /* input images.                                    */
    struct analyze_struct   /* Header information variable.                     */
                        img_info;
    FILE *fimgin;           /* Pointer for input image file.                    */
    char choice;            /* For keyboard input to yes/no question.           */
    char window;            /* For keyboard input to yes/no question.           */
    char normalize;         /* For keyboard input to yes/no question.           */

    char dummy[80];         /* Dummy pointer for counting number of entries     */
                            /* in a text file.                                  */
    float dummy_float;      /* Dummy floating point number for counting number  */
                            /* of entries in a text file.                       */
    char **header_name;     /* Pointer to pointer to name of header file.       */
    char **img_name;        /* Pointer to pointers to input image file names.   */
    float *input_waveform;  /* Pointer to data for input waveform.              */
    int num_vox;            /* Number of voxels in an image.                    */
    char *outmat_filename;  /* Pointer to output .mat file name.                */
    int num_read_in;        /* Number of objects read in in fread operation.    */
    unsigned
        short int *temp;    /* Pointer to temporary short int array, for        */
                            /* reading in 16-bit fMRI data.                     */
                            /* byte data.                                       */
    float float_thresh;     /* Temporary variable for reading in threshold.     */
    float psig;             /* Significance threshold.                          */
    float pi;               /* Pi.                                              */

    int x,y,z;              /* Indices into image arrays.                       */
    int dat_index2;         /* Index into images in time series.                */
    FILE *fintrmdt;         /* Pointer for output intermediate values text file.*/
    FILE *fOUTPUT_IMG;      /* Pointer for output image files.                  */
    float series_mean;      /* Mean value of a time series.                     */
    float series_mag;       /* Magnitued of a time series expressed as a vector */
                            /* after the mean has been set to zero.             */
    float *hann;            /* Values of the Hann window.                       */
    int x_int,y_int,z_int;  /* Coordinates of voxel of interest.                */
    int vox_index_int;      /* Voxel index of voxel of interest.                */
    float threshold;        /* Threshold for masking output.                    */
    FILE *fMASK;            /* File pointer to mask file.                       */
    unsigned char *mask;    /* Mask to separate brain vs. non-brain.            */

int     xdim,ydim,zdim;/* Image dimensions.                                     */
int     vox_index1;    /* Index into imgin[*], extracting a specific voxel.     */
int     vox_index2;    /* Another voxel index.                                  */
float  *imgin;         /* Pointer to pointers to input image data.              */
float  *Mimgout;       /* Pointer to Mag output data, ANALYZE form.             */
float  *Pimgout;       /* Pointer to Phase output data, ANALYZE form.           */
float **Rimgout;       /* Pointer to pointers to Real output data, ANALYZE form.*/
float **Iimgout;       /* Pointer to pointers to Imag output data, ANALYZE form.*/
int     num_dat_pts;   /* The number of images in time series.                  */
int     num_dat_pts2;  /* Next integer power of 2 >= num_dat_pts.               */
float   log2;          /* Integer power of 2.                                   */
int     dat_index;     /* Index into images in time series.                     */
int     num_harm;      /* The number of harmonics to be output.                 */
int     harm_index;    /* Index for parameters.                                 */

/* ---------------------------------------------------------------------------- */
/*  CHECK TO SEE IF CORRECT NUMBER OF ARGUMENTS WERE PASSED.                    */

    if (argc < 7) {
        printf("You typed in %d arguments.\n",argc);
        printf("fmrifft_hann: need input textfile, number of harmonics desired,\n");
        printf("              name of 8-bit .img mask image, and  X, Y, and Z\n");
        printf("              coordinates of voxel of interest.\n\n");

        printf("The first line in the input text file should be the name of\n");
        printf("the header file to be used with the input images.  The rest of\n") ;
        printf("the input text file should contain the names of the input .img\n") ;
        printf("files, with one filename per line.\n");

        printf("Output will be written to files named mN.img and pN.img,\n");
        printf("where N is the index of the harmonic.  The m files are\n");
        printf("the magnitudes, and the p files are the phases.\n");
        exit(1);
    }

/* ---------------------------------------------------------------------------- */
/*  COUNT NUMBER OF DATA POINTS IN INPUT TEXT FILE.                             */

    if ((ftxt = fopen(argv[1],"r")) == NULL) {
        printf ("fmrifft_hann: can't open %s\n",*(argv+1));
        exit(2);
    }
    printf("Counting number of data points in %s...\n",*(argv+1));
    num_dat_pts = 0;
    fscanf(ftxt,"%s",dummy);
    while(fscanf(ftxt,"%s %f",dummy,&dummy_float)!=EOF)
        num_dat_pts++;
    log2 = (int)(log((double)num_dat_pts)/log(2.));
    if ((int)pow(2.,log2)==num_dat_pts)
        num_dat_pts2=num_dat_pts;
    else
        num_dat_pts2=(int)pow(2,log2+1);


/* ---------------------------------------------------------------------------- */
/*  READ IN NAMES OF INPUT IMAGE FILES.                                         */

    header_name = (char **) malloc(sizeof(char *));
    *header_name = (char *) malloc(80*sizeof(char));
    img_name = (char **) malloc(num_dat_pts*sizeof(char *));
    input_waveform = (float *) malloc(2*num_dat_pts*sizeof(float));

    rewind(ftxt);
    fscanf(ftxt,"%s",*header_name);
    printf("Filename            Input waveform\n");
    for (dat_index=0; dat_index<num_dat_pts; dat_index++) {
        img_name[dat_index] = (char *) malloc(80*sizeof(char));
        fscanf(ftxt,"%s %f",img_name[dat_index],input_waveform+2*dat_index);
        input_waveform[2*dat_index+1] = 0.;
        printf("%20s %f\n",img_name[dat_index],input_waveform[2*dat_index+0]);
    }
    fclose(ftxt);

/* ---------------------------------------------------------------------------- */
/*  READ HEADER FILE; ASSUME IMAGES ARE SAME SIZE, SO USE SAME HEADER FILE      */
/*  FOR ALL IMAGES.                                                             */

    img_info = readhdr(header_name);

    xdim = img_info.dime.dim[1];
    ydim = img_info.dime.dim[2];
    zdim = img_info.dime.dim[3];

    num_vox = xdim*ydim*zdim;

/* ---------------------------------------------------------------------------- */
/*  SHOW USER INFORMATION READ IN SO FAR.                                       */

    printf("There are %d input images.\n", num_dat_pts);
    printf("Next integer power of 2 >= %d is %d\n",num_dat_pts,num_dat_pts2);
    printf("Will pad end of time series with %d zeros.\n",num_dat_pts2-num_dat_pts);
    printf("\nImages have dimensions:\n");
    printf("  xdim = %d\n  ydim = %d\n  zdim = %d\n",xdim,ydim,zdim);

    num_harm = atoi(argv[2]);
    if (num_harm>num_dat_pts) {
        printf("Number of harmonics requested, %d, exceeds\n",num_harm);
        printf("number of data points, %d.  No can do.\n",num_dat_pts);
        exit(3);
    }
    printf("\n%d harmonics will be output.\n", num_harm);

    x_int=atoi(argv[4]);
    y_int=atoi(argv[5]);
    z_int=atoi(argv[6]);
    printf("The voxel of interest is [%d, %d, %d].\n",x_int,y_int,z_int);
    x_int--;
    y_int--;
    z_int--;
    vox_index_int=(z_int*ydim+y_int)*xdim+x_int;
    printf("vox_index_int = %d\n",vox_index_int);
    printf("\nThere will be two output files per harmonic.\n");
    printf("Harmonic output will go to floating point ANALYZE .img files:\n");
    for (harm_index=0; harm_index<num_harm; harm_index++)
        printf("  mag%d.img and phs%d.img\n",harm_index,harm_index);

    printf("\nDo you wish to apply a Hann window to the time series <y/n>? ");
    scanf("%s",&window);
    while((window!='y') && (window!='n')) {
        printf("Type 'y' or 'n': ");
        scanf("%s",&window);
    }
    if (window=='y') printf("Will apply a Hann window to time series.\n");
    else printf("Will NOT apply a Hann window to time series.\n\n");

    printf("The mean of a time series may vary from voxel to voxel.\n");
    printf("Likewise, after setting the mean to zero, the total power\n");
    printf("of a time series may vary from voxel to voxel.  It may be\n");
    printf("desirable to set the means of all time series to zero, and\n");
    printf("power of all time series to 1, to be able to compare the\n");
    printf("spectral *distribution* from one voxel to that of another.\n\n");
    printf("Do you wish to normalize each time series <y/n>? ");
    scanf("%s",&normalize);
    while((normalize!='y') && (normalize!='n')) {
        printf("Type 'y' or 'n': ");
        scanf("%s",&normalize);
    }
    if (normalize=='y') printf("Will normalize time series.\n");
    else printf("Will NOT normalize time series.\n");

    printf("\nDoes the above information look okay <y/n>? ");
    scanf("%s",&choice);
    while((choice!='y') && (choice!='n')) {
        printf("Type 'y' or 'n': ");
        scanf("%s",&choice);
    }
    if (choice=='n') exit(4);

/* ---------------------------------------------------------------------------- */
/*  READ IN IMAGES.  ALLOCATE MEMORY FOR INPUT AND OUTPUT IMAGES.               */

    temp = (unsigned short int *) malloc(num_vox*sizeof(unsigned short int));
    if (temp==NULL) {
        printf("fmrifft_hann: malloc failed for temp.\n");
        exit(5);
    }
    imgin  = (float *) malloc(2*num_vox*num_dat_pts2*sizeof(float));
    if (imgin==NULL) {
        printf("fmrifft_hann: malloc failed for imgin.\n");
        exit(6);
    }
    for (vox_index1=0; vox_index1<2*num_vox*num_dat_pts2; vox_index1++)
        imgin[vox_index1]=0.;

    for (dat_index=0; dat_index<num_dat_pts; dat_index++) {
        if ((fimgin = fopen(img_name[dat_index],"r")) == NULL) {
            printf("fmrifft_hann: can't open %s\n",img_name[dat_index]);
            exit(8);
        }
        printf("Reading in data from %s...\n",img_name[dat_index]);
        num_read_in=fread(temp,sizeof(short int),num_vox,fimgin);
        fclose(fimgin);
        if (num_read_in<num_vox) {
            printf("Read in only %d objects, not %d.\n",num_read_in,num_vox);
            exit(9);
        }
        printf("Transfering data from %s to memory array...\n",img_name[dat_index]);
        for (vox_index1=0; vox_index1<num_vox; vox_index1++) {
            vox_index2=vox_index1*num_dat_pts2+dat_index;
            imgin[2*vox_index2+0]=(float)(*(temp+vox_index1));
            imgin[2*vox_index2+1]=0.;
        }
    }

    printf("Allocating memory for Mimgout, Pimgout, Rimgout, Iimgout,...\n");
    Mimgout = (float *) malloc(num_vox*sizeof(float));
    Pimgout = (float *) malloc(num_vox*sizeof(float));
    Rimgout = (float **) malloc(num_harm*sizeof(float *));
    Iimgout = (float **) malloc(num_harm*sizeof(float *));
    for (harm_index=0; harm_index<num_harm; harm_index++) {
        Rimgout[harm_index]=(float *) malloc(num_vox*sizeof(float));
        if (Rimgout[harm_index]==NULL) {
            printf("fmrifft_hann: malloc failed for Rimgout at harm_index = %d.\n",
                                                                        harm_index);
            exit(10);
        }
        Iimgout[harm_index]=(float *) malloc(num_vox*sizeof(float));
        if (Iimgout[harm_index]==NULL) {
            printf("fmrifft_hann: malloc failed for Iimgout at harm_index = %d.\n",
                                                                        harm_index);
            exit(11);
        }
    }

/* ---------------------------------------------------------------------------- */
/* READ IN 8-BIT MASK.                                                          */

    mask = (unsigned char *) malloc(num_vox*sizeof(unsigned char));
    if (mask == NULL) {
        printf("fmri_lagNdisp4: malloc failed for mask.\n");
        exit(28);
    }
    if ((fMASK = fopen(argv[3],"r")) == NULL) {
        printf("fmri_lagNdisp4: can't open %s\n",argv[3]);
        exit(29);
    }
    printf("Reading in mask from %s...\n",argv[3]);
    num_read_in=fread(mask,sizeof(unsigned char),num_vox,fMASK);
    fclose(fMASK);
    if(num_read_in!=num_vox) {
        printf("fmri_lagNdisp4: %d voxels read in; should have been %d.\n",
                                                num_read_in,num_vox);
        exit(30);
    }


/* ---------------------------------------------------------------------------- */
/*  CALCULATE FFT OF INPUT WAVE FORM.                                           */
/*                                                                              */
    four1(input_waveform-1,(unsigned long) num_dat_pts,1);

/* ---------------------------------------------------------------------------- */
/* ALLOCATE MEMORY FOR HANNING WINDOW, AND CALCULATE ITS VALUES.                */

    hann = (float *) malloc(num_dat_pts*sizeof(float));
    if (hann==NULL) {
        printf("fmrifft_hann: error malloc'ing for hann.\n");
        exit(14);
    }
    pi = 2.*acos(0.);
/*    printf("Hann window:\n"); */
    for (dat_index=0; dat_index<num_dat_pts; dat_index++) {
        hann[dat_index]=0.5*(1.-cos(2*pi*dat_index/(num_dat_pts-1)));
/*          printf("  hann[%3d] = %f\n",dat_index+1,hann[dat_index]); */
    }


/* ---------------------------------------------------------------------------- */
/*  NOW CALCULATE FFT, VOXEL-BY-VOXEL.                                          */
/*                                                                              */

    printf("\nNow generating harmonic images...\n\n");

    for (vox_index1=0; vox_index1<num_vox; vox_index1++) {

        /* -------------------------------------------------------------------- */
        /* IF normalize IS 'y', NORMALIZE TIME SERIES MEAN AND POWER.           */

        if (normalize=='y') {

            /* -------------------------------------------------------------------- */
            /* CALCULATE MEAN VALUE OF THIS VOXEL'S TIME SERIES.                    */

            series_mean = 0.;
            for (dat_index=0; dat_index<num_dat_pts2; dat_index++) {
                vox_index2=vox_index1*num_dat_pts2+dat_index;
                series_mean += imgin[2*vox_index2];
            }
            series_mean /= num_dat_pts2;


            /* -------------------------------------------------------------------- */
            /* SUBTRACT MEAN VALUE FROM TIME SERIES, GIVING IT A ZERO MEAN.         */

            for (dat_index=0; dat_index<num_dat_pts2; dat_index++) {
                vox_index2=vox_index1*num_dat_pts2+dat_index;
                imgin[2*vox_index2] -= series_mean;
            }


            /* -------------------------------------------------------------------- */
            /* NOW CALCULATE MAGNITUDE OF THIS VOXEL'S TIME SERIES.                 */

            series_mag = 0.;
            for (dat_index=0; dat_index<num_dat_pts2; dat_index++) {
                vox_index2=vox_index1*num_dat_pts2+dat_index;
                series_mag += (imgin[2*vox_index2] * imgin[2*vox_index2]);
            }
            series_mag = (float)sqrt((double)series_mag);


            /* -------------------------------------------------------------------- */
            /* DIVIDE TIME SERIES BY MAGNITUDE, GIVING IT A UNIT MAGNITUDE.         */
            /* This is equivalent to normalizing the total power in the series.     */
            /* Consider adding mean back, but this is not done here, commented out. */
            /* It is now decommented.                                               */

            for (dat_index=0; dat_index<num_dat_pts2; dat_index++) {
                vox_index2=vox_index1*num_dat_pts2+dat_index;
                if (series_mag!=0.)
                    imgin[2*vox_index2] /= series_mag;
                imgin[2*vox_index2]+=series_mean;
            }
        }

        /* -------------------------------------------------------------------- */
        /* MULTIPLY TIME SERIES ELEMENT-BY-ELEMENT BY HANNING WINDOW IF FLAG    */
        /* window IS SET TO 'y'.                                                */

        if (window=='y')
            for (dat_index=0; dat_index<num_dat_pts2; dat_index++) {
                vox_index2=vox_index1*num_dat_pts2+dat_index;
                imgin[2*vox_index2] *= hann[dat_index];
            }

        /* -------------------------------------------------------------------- */
        /* DO FFT....                                                           */
        /* Note that Numerical Recipes in C follows the physicists' convention, */
        /* where the sign of the exponent in the Fourier transform is positive, */
        /* rather than negative, as in the engineering convention.  This means  */
        /* that the sign of the imaginary component of the FFT will be opposite */
        /* that which is expected under the engineering convention.             */

        if (vox_index1==vox_index_int) {
            printf("Voxel of interest, [%d, %d, %d], time series:\n",
                                                        x_int+1,y_int+1,z_int+1);
            for (dat_index=0; dat_index<num_dat_pts2; dat_index++) {
                vox_index2=vox_index_int*num_dat_pts2+dat_index;
                printf("%f\n",imgin[2*vox_index2]);
            }
        }
        vox_index2=2.*vox_index1*num_dat_pts2;

        four1(imgin+vox_index2-1,num_dat_pts2,1);

        if (vox_index1==vox_index_int) {
            printf("Voxel of interest, [%d, %d, %d], FFT, real & imaginary components:\n",
                                                        x_int+1,y_int+1,z_int+1);
            for (dat_index=0; dat_index<num_dat_pts2; dat_index++) {
                vox_index2=vox_index_int*num_dat_pts2+dat_index;
                printf("%f %f\n",imgin[2*vox_index2],imgin[2*vox_index2+1]);
            }
        }

        /* -------------------------------------------------------------------- */
        /* TRANSFER RESULTS TO OUTPUT ARRAY....                                 */
        /* Extract Real and Imaginary components.                               */

        for (dat_index=0; dat_index<num_harm; dat_index++) {
            vox_index2=vox_index1*num_dat_pts2+dat_index;
            Rimgout[dat_index][vox_index1] = imgin[2*vox_index2+0];
            Iimgout[dat_index][vox_index1] = imgin[2*vox_index2+1];
        }
    }


    /* ------------------------------------------------------------------------ */
    /* LOOP OVER HARMONICS, AND OUTPUT REAL, IMAG, MAG, AND PHASE IMAGES.       */

    outmat_filename = (char *) malloc(80*sizeof(char));

    for (harm_index=0; harm_index<num_harm; harm_index++) {

        /* -------------------------------------------------------------------- */
        /* CONVERT REAL/IMAG DATA TO MAG/PHASE.                                 */
        /* Mag = sqrt(R*R + I*I)                                                */
        /* Phase = atan(I/R)                                                    */

        for (vox_index1=0; vox_index1<num_vox; vox_index1++) {

            vox_index2=vox_index1*num_dat_pts2+harm_index;

            Mimgout[vox_index1] = (float)
                sqrt((double) ((imgin[2*vox_index2+0]*imgin[2*vox_index2+0])
                                                   +
                               (imgin[2*vox_index2+1]*imgin[2*vox_index2+1])))*(float)mask[vox_index1];

            if (imgin[2*vox_index2+0]!=0.)
                Pimgout[vox_index1]=(float)
                        atan((double)(imgin[2*vox_index2+1]/imgin[2*vox_index2+0]))*(float)mask[vox_index1];
            else
                Pimgout[vox_index1]=0.;
        }


    /* ------------------------------------------------------------------------ */
    /*  WRITE OUT PARAMETERS TO OUTPUT FILES.                                   */
    /*  FOR EACH OUTPUT ARRAY...                                                */

        /* -------------------------------------------------------------------- */
        /* Mag image.                                                           */
        sprintf(outmat_filename,"mag%d.img\0",harm_index);
        if ((fOUTPUT_IMG = fopen(outmat_filename,"w")) == NULL) {
            printf ("fmrifft_hann: can't open %s\n",outmat_filename);
            exit(15);
        }
        printf("  Writing output to %s...\n",outmat_filename);
        fwrite(Mimgout,sizeof(float),num_vox,fOUTPUT_IMG);
        fclose(fOUTPUT_IMG);

        /* -------------------------------------------------------------------- */
        /* Phase image.                                                         */
        sprintf(outmat_filename,"phs%d.img\0",harm_index);
        if ((fOUTPUT_IMG = fopen(outmat_filename,"w")) == NULL) {
            printf ("fmrifft_hann: can't open %s\n",outmat_filename);
            exit(16);
        }
        printf("  Writing output to %s...\n",outmat_filename);
        fwrite(Pimgout,sizeof(float),num_vox,fOUTPUT_IMG);
        fclose(fOUTPUT_IMG);

        /* -------------------------------------------------------------------- */
        /* Real image.                                                          */
/*        sprintf(outmat_filename,"r%d.img\0",harm_index);
        if ((fOUTPUT_IMG = fopen(outmat_filename,"w")) == NULL) {
            printf ("fmrifft_hann: can't open %s\n",outmat_filename);
            exit(17);
        }
        printf("  Writing output to %s...\n",outmat_filename);
        fwrite(Rimgout[harm_index],sizeof(float),num_vox,fOUTPUT_IMG);
        fclose(fOUTPUT_IMG);
*/
        /* -------------------------------------------------------------------- */
        /* Imag image.                                                          */
/*      sprintf(outmat_filename,"i%d.img\0",harm_index);
        if ((fOUTPUT_IMG = fopen(outmat_filename,"w")) == NULL) {
            printf ("fmrifft_hann: can't open %s\n",outmat_filename);
            exit(18);
        }
        printf("  Writing output to %s...\n",outmat_filename);
        fwrite(Iimgout[harm_index],sizeof(float),num_vox,fOUTPUT_IMG);
        fclose(fOUTPUT_IMG);
*/
    } /* CLOSE harm_index LOOP */


/* ---------------------------------------------------------------------------- */
/*  FREE MEMORY.                                                                */

    printf("Releasing memory...\n");
    free(header_name);
    free(img_name);
    free(outmat_filename);
    free(imgin);
    free(Mimgout);
    free(Pimgout);
    free(Rimgout);
    free(Iimgout);
    free(temp);
    free(mask);
}

