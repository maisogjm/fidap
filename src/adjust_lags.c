#include <math.h>
#include "nrutil.h"
#include <stdio.h>
#include <stdlib.h>
#include "analyze6.h"

                            /* ************************************************ */
                            /* FUNCTIONS CALLED.                                */
struct analyze_struct readhdr(char *argv[]);

unsigned char *unsigned_char_malloc(int num_items);
int *int_malloc(int num_items);
unsigned short *unsigned_short_malloc(int num_items);
char *char_malloc(int num_items);
char ***char_malloc3(int num_items1, int num_items2, int num_items3);
float *float_malloc(int num_items);
float **float_malloc2(int num_items1, int num_items2);
unsigned short ***unsigned_short_malloc3(int num_items1, int num_items2, int num_items3);

int read_unsigned_char_data(char *filename, unsigned char *array, int num_items);
int read_unsigned_short_int_data(char *filename, unsigned short int *array, int num_items);
int float_fwrite(char *filename, float *array, int num_items);


/* **************************************************************************** */
/* Program adjust_lags.c                                                        */
/* ANSI C Code by Jose' Ma. Maisog, M.D.                                        */
/* Section on Functional Brain Imaging,                                         */
/* Laboratory of Psychology and Psychopathology                                 */
/* National Institute of Mental Health                                          */
/*                                                                              */
/* Program reads in RELATIVE lag maps and adjusts the lags for the slice        */
/* acquisition sequences, outputing an ABSOLUTE lag map.                        */
/*                                                                              */
/* OCTOBER 1994: Started changes to allow cross-run simultaneous optimization   */
/* of lag and dispersion.  Previous version was only for within-run.            */
/*                                                                              */
/* NOVEMBER 1994: Started changes to perform a blob-wise modeling of the        */
/* hemodynamic response function.  This will effectively be a selective spatial */
/* averaging, which should improve signal-to-noise.  Additional inputs are a    */
/* 16-bit cluster/blobs image which defines clusters of spatially contiguous    */
/* voxels, a floating point lag image, and a floating point dispersion image.   */
/* Outputs will be a new 16-bit cluster image, and new floating point lag and   */
/* dispersion images.  Now named fmri_lagNdisp3.                                */
/*                                                                              */
/* JANUARY 1995: Started changes to have program reassess the spatial extent of */
/* clusters, after having estimated the hemodynamic response of a cluster.      */
/* This requires checking whether a cluster's spatial extent has changed, and   */
/* whether it has split into two clusters or become contiguous with another     */
/* cluster.                                                                     */
/*                                                                              */
/* FEBRUARY 1995: Changing from correlation to ANCOVA.                          */
/*                                                                              */
/* APRIL 1995: Correcting misunderstanding regarding calculation of slope.      */
/* Before, slope was calculated by pooling data across voxels in a blob.        */
/* Now slope is calculated on a per voxel basis.                                */
/*                                                                              */
/* 04.19.95: Adding code to take into account the difference in time acquisi-   */
/* tions for different slices.                                                  */
/*                                                                              */
/* 04.24.95: Changed output waveforms to simpler blob mean time series, not     */
/* taking into account difference in time lag between slices.  The mathematics, */
/* however, has not been changed and still takes into account inter-slice time  */
/* lags.                                                                        */
/*                                                                              */
/* 04.25.95: Cleaned up code, removing debugging code & made final adjustments  */
/* to optimization routines.                                                    */
/*                                                                              */
/* 05.02.95: Modified indexing into imgin.  Before, indexing was                */
/* imgin[run][voxel][scan]; now it will be imgin[voxel][run][scan]              */
/* It is hoped that this will improve memory performance, and make the program  */
/* run faster.                                                                  */
/*                                                                              */
/* 07.11.95: Modified code so that instead of an iterative procedure for esti-  */
/* mating the hemodynamic response, analytic formulas are used to estimate the  */
/* temporal dispersion and lag, generating fmri_lagNdisp8.  A Gaussian curve is */
/* still used to model the hemodynamic response.                                */
/*                                                                              */
/* 12.22.95 : JMM : Modified fmri_lagNdisp8 to generate extract_time_series.    */
/* Program just outputs fMRI time series for specified voxel indices of         */
/* interest                                                                     */
/*                                                                              */
/* 12.26.95 : JMM : Modified extract_time_series to generate                    */
/* multivar_multiple_regression.  Program reads in multiple-run fMRI data and   */
/* regressors, and performs a voxel-by-voxel multivariate multiple regression.  */
/*                                                                              */
/* 01.11.96 : JMM : Added code to multivar_multiple_regression to do within-    */
/* voxel estimate of lag and dispersion.                                        */
/*                                                                              */
/* 03.04.96 : JMM : Added code to make an estimate of effective degrees of      */
/* freedom, using method of Worsley KJ, Friston KJ, "Analysis of fMRI           */
/* Time-Series Revisited -- Again," Neuroimage, submitted.  Also, added         */
/* subroutines to do memory allocation, file I/O, and matrix multiplication, in */
/* an attempt to modularize code and clean it up a little.                      */
/*                                                                              */
/* 03.16.96 : JMM : Modified Dmultivar_multiple_regression to create            */
/* voxelwise_hemodynamics.  Basically chopped out code which estimates          */
/* hemodynamic response for each voxel and each independent varaible, making it */
/* a stand-alone program.                                                       */
/*                                                                              */
/* 04.02.96 : JMM : Added option to voxelwise_hemodynamics, allowing user to    */
/* have lag maps with absolute or relative times.                               */
/*                                                                              */
/* 09.16.96 : JMM : Modified voxelwise_hemodynamics to generate adjust_lags.    */
/* Program converts RELATIVE lag maps into ABSOLUTE lag maps.                   */

main(int argc, char *argv[])
{
                            /* INDICES AND CONSTANTS.                           */
    int num_vox;            /* Number of voxels in an image.                    */
    int vox_index;          /* Voxel index.                                     */
    int xdim,ydim,zdim;     /* Image dimensions.                                */
    int x,y,z;              /* XYZ coordinates.                                 */
    int *slice_seq1;        /* Slice acquisition sequence.                      */
    int *slice_seq2;        /* Slice acquisition sequence.                      */
    float TR;               /* TR in seconds.                                   */
    float pi;               /* Pi.                                              */
    struct analyze_struct
                img_info;   /* Header information variable.                     */
    FILE *ftxt;             /* File pointer to textfile.                        */
    int check_fscanf;       /* Check on result of fscanf operation.             */

                            /* LAG MAPS.                                        */
    float *rel_lag_map;     /* Relative lag map.                                */
    float *abs_lag_map;     /* Relative lag map.                                */

                            /* MASK                                             */
    unsigned char *mask;    /* Pointer to mask image.                           */



/* ---------------------------------------------------------------------------- */
/*  CHECK TO SEE IF CORRECT NUMBER OF ARGUMENTS WERE PASSED.                    */

    if (argc<7) {
        printf("adjust_lags: need input header file, floating point map of RELATIVE\n");
        printf("             lags, name of 8-bit mask, text file containing slice\n");
        printf("             acquisition sequence, TR in seconds, and output file name.\n\n");

        printf("All arguments are of course to be separated by spaces.\n\n");
        printf("This program will subtract an adjustment factor from the lags\n");
        printf("to convert them from relative lags to absolute lags.  Lags less\n");
        printf("than zero will be truncated to zero.\n");
        exit(1);
    }


/* ---------------------------------------------------------------------------- */
/*  READ HEADER FILE; ASSUME IMAGES ARE SAME SIZE, SO USE SAME HEADER FILE      */
/*  FOR ALL IMAGES.  READ IN TR.                                                */

    printf("Reading header info from %s...\n",argv[1]);
    img_info = readhdr(argv+1);

    xdim = img_info.dime.dim[1];
    ydim = img_info.dime.dim[2];
    zdim = img_info.dime.dim[3];
    num_vox = xdim*ydim*zdim;
    printf("Images have dimensions: [%d %d %d]\n\n",xdim,ydim,zdim);

    TR=atof(argv[5]);
    printf("TR = %g seconds.\n\n",TR);


/* ---------------------------------------------------------------------------- */
/* READ IN SLICE ACQUISITION SEQUENCE.                                          */
/* Convert to numbers from 0 to zdim-1.  If slice k was the nth slice acquired, */
/* set slice_seq2[k-1] equal to (n-1).                                          */

    printf("\nReading in slice acquisition sequence from %s...\n",argv[4]);
    if ((ftxt = fopen(argv[4],"r")) == NULL) {
        printf("voxelwise_hemodynamics: can't open %s\n",argv[4]);
        exit(2);
    }
    slice_seq1=int_malloc(zdim);
    slice_seq2=int_malloc(zdim);

    printf("Slices were acquired in this sequence: ");
    for (z=0; z<zdim; z++) {
        check_fscanf=fscanf(ftxt,"%d",slice_seq1+z);
        if (check_fscanf==EOF) {
            printf("adjust_lags: error reading slice sequence from %s.\n",
                                                                        argv[4]);
            exit(3);
        }
        printf("%d ",slice_seq1[z]);
    }
    fclose(ftxt);
    printf("\n");
/*        for (z=0; z<zdim; z++)
        printf("slice_seq1[%d] = %d\n",z,slice_seq1[z]);
*/
    for (z=0; z<zdim; z++)
        if ((slice_seq1[z]<1) || (slice_seq1[z]>zdim)) {
            printf("slice_seq1[%d] = %d\n",z,slice_seq1[z]);
            printf("adjust_lags: weird values for slice sequence.\n");
            printf("Have you forgotten to put this information in the input text file?\n");
            exit(4);
         }
    printf("In other words,\n");
    for (z=0; z<zdim; z++) {
        slice_seq2[slice_seq1[z]-1]=z;
        printf("  slice # %d was acquired # %d in sequence\n",
                slice_seq1[z],slice_seq2[slice_seq1[z]-1]+1);
    }
    printf("\n\n");


/* ---------------------------------------------------------------------------- */
/* READ IN 8-BIT MASK.                                                          */

    mask = unsigned_char_malloc(num_vox);
    read_unsigned_char_data(argv[3],mask,num_vox);


/* ---------------------------------------------------------------------------- */
/*  READ IN RELATIVE LAG MAP.                                                   */

    rel_lag_map = float_malloc(num_vox);
    read_float_data(argv[2],rel_lag_map,num_vox);


/* ---------------------------------------------------------------------------- */
/* LOOP OVER VOXELS, CONVERT RELATIVE LAGS TO ABSOLUTE LAGS.                    */

    abs_lag_map = float_malloc(num_vox);
    pi = 2.*acos(0.);

    printf("\n\nNow converting relative lags into absolute lags...\n");
    for (z=0; z<zdim; z++) {
        printf("Subtracting %6g seconds to lags in slice %d...\n",
                                TR*(float)slice_seq2[z]/(float)zdim,z+1);
        for (y=0; y<ydim; y++)
            for (x=0; x<xdim; x++) {

                vox_index=((z*ydim)+y)*xdim+x;

                if (mask[vox_index]!=0) {
                    abs_lag_map[vox_index]=rel_lag_map[vox_index]
                                -(TR*(float)slice_seq2[z]/(float)zdim);

                }
                else
                    abs_lag_map[vox_index]=rel_lag_map[vox_index];

                if (abs_lag_map[vox_index]<0.)
                    abs_lag_map[vox_index]=0.;
            } /* END x LOOP. */
    } /* END z LOOP. */

    printf("\n\n");

/* ---------------------------------------------------------------------------- */
/*  WRITE ESTIMATED REGRESSION COEFFICIENT MAPS AND WILKS LAMBDA MAPS TO        */
/*  FLOATING POINT FILES.                                                       */

    float_fwrite(argv[6],abs_lag_map,num_vox);


/* ---------------------------------------------------------------------------- */
/*  FREE MEMORY.                                                                */

    free(slice_seq1);
    free(slice_seq2);
    free(rel_lag_map);
    free(abs_lag_map);
    free(mask);
}

