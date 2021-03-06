#include <string.h>
#include <math.h>
#include "nrutil.h"
#include <stdio.h>
#include <stdlib.h>
#include "analyze6.h"


/* **************************************************************************** */
/* FUNCTIONS FOR MATRIX INVERSION.                                              */

void Dludcmp(double **a, int n, int *indx, double *d);
void Dlubksb(double **a, int n, int *indx, double b[]);


/* **************************************************************************** */
/* FUNCTIONS FOR MEMORY ALLOCATION                                              */

unsigned char *unsigned_char_malloc(int num_items);
int *int_malloc(int num_items);
signed short *signed_short_malloc(int num_items);
signed short ***signed_short_malloc3(int num_items1, int num_items2,
                                                                int num_items3);
char *char_malloc(int num_items);
char **char_malloc2(int num_items1, int num_items2);
char ***char_malloc3(int num_items1, int num_items2, int num_items3);
float *float_malloc(int num_items);
float **float_malloc2(int num_items1, int num_items2);
float ***float_malloc3(int num_items1, int num_items2, int num_items3);


/* **************************************************************************** */
/* FUNCTIONS FOR FILE I/O.                                                      */

struct analyze_struct readhdr(char *argv[]);
int read_unsigned_char_data(char *filename, unsigned char *array, int num_items);
int read_signed_short_int_data(char *filename, signed short int *array,
                                                                int num_items);
int read_float_data(char *filename, float *array, int num_items);
int float_fwrite(char *filename, float *array, int num_items);
int signed_short_fwrite(char *filename, signed short *array, int num_items);
int write_double_matrix2textfile(char *filename, double **array, int num_rows,
                                                                int num_cols);
int write_NRCdouble_matrix2textfile(char *filename, double **array, int num_rows,
                                                                int num_cols);
int write_float_matrix2textfile(char *filename, float **array, int num_rows,
                                                                int num_cols);

int mult_NRCdouble_matrices(double **input1, double **input2, int num_rows1,
        int num_cols1, int num_cols2, double **output);
int mult_NRCdouble_matrices_tr_both(double **input1, double **input2, int num_rows1,
        int num_cols1, int num_rows2, double **output);
int mult_NRCdouble_matrices_tr1(double **input1, double **input2, int num_rows1,
        int num_cols1, int num_rows2, double **output);
int mult_NRCdouble_matrices_tr2(double **input1, double **input2, int num_rows1,
        int num_cols1, int num_rows2, double **output);


/* **************************************************************************** */
/* Program multivar_multiple_regression.c                                       */
/* ANSI C Code by Jose' Ma. Maisog, M.D.                                        */
/* Section on Functional Brain Imaging,                                         */
/* Laboratory of Psychology and Psychopathology                                 */
/* National Institute of Mental Health                                          */
/*                                                                              */
/* Program reads in fMRI data for several runs and independent variables, and   */
/* performs a multivariate multiple regression.  Matrix inversion is done using */
/* LU decomposition method.  Many thanks to Dr. John D. Van Horn for his        */
/* assistance with matrix inversion, and for many discussions about             */
/* multivariate methods.                                                        */
/*                                                                              */
/* Reference: Rencher AC, Methods of Multivariate Analysis, New York: John      */
/* Wiley & Sons, 1995, pp. 366-380.  Note that this implementation is for       */
/* "fixed", not "random", X's.  An attempt will be made in this code to follow  */
/* the naming convention of the regression variables in this reference.         */
/*                                                                              */
/* See also :                                                                   */
/* Worsley KJ, Friston KJ, "Analysis of fMRI Time-Series Revisited -- Again,"   */
/* Neuroimage, 2:173-181, for the method of estimating the effective degrees of */
/* freedom;                                                                     */
/*                                                                              */
/* and:                                                                         */
/* Press WH, Teukolsky SA, Vetterling WT, Flannery BP, Numerical Recipes in C,  */
/* 2nd edition, Cambridge:Cambridge University Press (1992), pp. 43-49, for     */
/* code used to perform matrix inversion.                                       */
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
/* 03.19.96 : JMM : Stripped out code to estimate hemodynamic response on a     */
/* voxel-wise basis, and made it into a stand-alone program named               */
/* voxelwise_hemodynamics.  Dmultivar_multiple_regression will now read in      */
/* lag and dispersion maps, possibly generated by voxelwise_hemodynamics.       */
/*                                                                              */
/* 03.26.96 : JMM : Changed datatype of imgin from float *** to                 */
/* unsigned short int ***, to avoid problem of slow memory access (probably     */
/* swapping scans in and out of memory).                                        */
/*                                                                              */
/* 04.04.96 : JMM : Modified array for hemodynamic response to sample time more */
/* finely, thus ameliorating aliasing problem with very small estimates of      */
/* temporal dispersion.  Time is sampled every 1/(TR/zdim) seconds, rather than */
/* every 1/TR, and is still a discrete representation of time.  The aliasing    */
/* will still be problematic if dispersions equal to or less than 1/(TR/zdim)   */
/* are used.                                                                    */
/*                                                                              */
/* 05.02.96 : JMM : Changed from the Singular Value Decomposition to the LU     */
/* decomposition, in an effort to improve speed.  Also, if no within-voxel      */
/* estimate of effective df is desired, then I use the method at the top of     */
/* page 49 of Numerical Recipes in C (Press WH, Teukolsky SA, Vetterling WT,    */
/* Flannery BP, Numerical Recipes in C, 2nd edition, Cambridge:Cambridge        */
/* University Press (1992)) to efficiently calculate inv(X'*X)*X', without      */
/* explicitly calculating inv(X'*X); this saves time and is more accurate.      */
/* However, in the case the an estimate of effective df is desired, then        */
/* inv(X'*X) will need to be explicitly calculated but might need inv(A) itself */
/* if calculating effective df, which is an option.                             */
/*                                                                              */
/* 05.18.96 : JMM : Implemented Jim's idea for faster convolution -- calculate  */
/* terms only for values of hemodynamic response greater than some threshold.   */
/*                                                                              */
/* 09.23.96 : JMM : Added code so that a textfile containing estimated effect-  */
/* ive df is output in a voxel of interest's directory.                         */
/*                                                                              */
/* 12.03.96 : JMM : Modified code to be able to handle input stimulus waveforms */
/* which are not time-locked with the TR.  For example, the input stimulus      */
/* may be 2 seconds on, 3 seconds off, repeating, with a TR of 3 seconds.       */
/* Input textfiles for effects of interest must now be two-columns, the first   */
/* to indicate duration in seconds, the second to indicate effect size          */
/* (e.g., second column can be 1's and 0's to indicate on-off squarewave).      */
/*                                                                              */
/* 08.22.05 : JMM : Added lines to check for "bad" matrices.                    */
/*                                                                              */
/* 09.20.2005 : JMM : Changed unsigned chars to signed chars to allow signed    */
/* voxels values.                                                               */


main(int argc, char *argv[])
{
                            /* INDICES AND CONSTANTS.                           */
    int num_vox;            /* Number of voxels in an image.                    */
    int vox_index;          /* Index into imgin[][][*], extracting a specific   */
                            /* voxel.                                           */
    int vox_index2;         /* Secondary index into voxels.                     */
    int *vox_int;           /* List of voxel indices of interest.               */
    int yes_voxel_of_int;   /* Flag indicating whether current voxel is a voxel */
                            /* of interest.                                     */
    int num_vox_int;        /* Number of voxel indices of interest.             */
    int xdim,ydim,zdim;     /* Image dimensions.                                */
    int x,y,z;              /* XYZ coordinates.                                 */
    int num_dep_var;        /* Number of dependent variables.  This is the      */
                            /* number of columns in the matrix of observations, */
                            /* Y.  If you're concatenating fMRI runs, this will */
                            /* be equal to 1.                                   */
    int col_index;          /* Index into columns in the matrix of observations.*/
    int num_rows;           /* The number of rows in the matrix of observations.*/
    int zdimXnum_rows;      /* zdim*num_rows, number of elements in over-       */
                            /* sampled time series.                             */
    int num_dep_varXnum_rows;
                            /* num_dep_var times num_rows                       */
    register int row_index; /* Index into time points.                          */
    int tot_num_scans;      /* Total number of scans.  Should be equal to       */
                            /* num_dep_var * num_rows                           */
    int *slice_seq1;        /* Slice acquisition sequence.                      */
    int *slice_seq2;        /* Slice acquisition sequence.                      */
    float TR;               /* TR in seconds.                                   */
    int tick;               /* Fraction of a second represented in oversampled  */
                            /* arrays X_os and SX_os.  E.g., with tick=4, time  */
                            /* is represented in units of 1/4 seconds.          */
    char line[10];

                            /* INPUT IMAGES.                                    */
    FILE *ftxt;             /* Pointer for input text file containing filenames */
    char ***img_name;       /* Pointer to pointers to input image file names.   */
                            /* Indexing:                                        */
                            /*      img_name[column][timepoint][ASCII character]*/
    char dummy[80];         /* Dummy pointer for counting number of entries     */
                            /* in a text file.                                  */
    int check_fscanf;       /* Flag to check error status of fscanf operation.  */
    struct analyze_struct
                img_info;   /* Header information variable.                     */
    unsigned                /* Pointer to temporary short int array, for        */
        short int *temp;    /* reading in 16-bit fMRI data.                     */
                            /* byte data.                                       */
    signed short ***imgin;/* Pointer to input image data.                     */
                            /* Indexing: imgin[voxel][run][scan]                */
    char normalize[2];      /* Flag indicating whether normalization is desired */
                            /* not.                                             */
    float **global_mean;    /* Global mean intensity of imgin[*][run][scan].    */
                            /* Indexing: global_mean[run][scan].                */
    float grand_global_mean;/* Grand global mean intensity across scans.        */
    int num_files_read_in;  /* Number of files read in so far.                  */

                            /* MASK                                             */
    unsigned char *mask;    /* Pointer to mask image.                           */
    int num_mask_vox;       /* Number of non-zero voxels in *mask.              */
    int count;              /* Number of voxels processed thus far.             */
    int *mask_voxel;        /* Voxel indices of nonzero mask voxels.            */
                            /* mask_voxel[N] is the index of the (N+1)th        */
                            /* nonzero voxel in the mask.                       */

                            /* MULTIVARIATE MULTIPLE REGRESSION.                */
    FILE *f_indvar;         /* Pointer to file containing independent variable. */
    int num_ind_var;        /* Number of independent variables.                 */
    int num_ind_var_int;    /* Number of independent variables of interest.     */
    int var_index;          /* Index over independent variables.                */
    int argv_index;         /* Index over command line arguments.               */
    int *normalize_flag;    /* Flag to indicate whether or not an independent   */
                            /* variable is to be Z-score normalized (mean set   */
                            /* to zero, variance set to one.)                   */
    double max,min;         /* Max and min values in an independent variable.   */
    double peak2peak;       /* max-min, peak-to-peak amplitude.                 */
    int *lag_adjust_flag;   /* Flag to indicate whether or not the lag of an    */
                            /* independent variable is to be adjusted for the   */
                            /* slight lag due to the slice-by-slice acquisition.*/
                            /* for constant and ramp functions.                 */
    float var_mean,         /* Mean of an independent variable, to be set to 0. */
          var_var;          /* Variance of an ind. variable, to be set to 1.    */
    float float_buffer;     /* Buffer for reading in floating point numbers.    */
    char only_vox_of_int[2];/* Flag indicating whether only voxels of interst   */
                            /* are to be analyzed.                              */
    float tot_time;         /* Total time for one fMRI run = TR * num_scans     */
    float duration;         /* Duration of an epoch in the input stimulus       */
                            /* waveform.                                        */
    int t_pointer;          /* Index recording the start of a new epoch in the  */
                            /* input stimulus waveform.                         */
    int t_index;            /* Index stepping through time within an epoch in   */
                            /* the input stimulus waveform.                     */
    double **X_os,          /* Matrix of independent variables, oversampled in  */
                            /* time (zdim*num_rows elements, not zdim).         */
          **SX_os,          /* Matrix of smoothed independent variables,        */
                            /* oversampled in time.                             */
          **SX,             /* Matrix of smoothed independent variables,        */
                            /* zdim elements.                                   */
          **X_r1,           /* Reduced X, for testing dependence on each        */
                            /* independent variable.                            */
          **X_r2,           /* Reduced X, for testing dependence on independent */
                            /* variables of interest.                           */
          **Y,              /* Matrix of dependent variables, the observed fMRI */
                            /* signal, one column per run, one row per scan.    */
          **XtrX,           /* Matrix to be inverted.                           */
          **inv_Xtr_X_Xtr,  /* inv(X'*X)*X, When multiplied by Y, this matrix   */
                            /* gives the least squares estimate for regression  */
                            /* coefficients.                                    */
            d,              /* +/-1 depending on number of row interchanges     */
                            /* in matrix inversion.                             */
           *col,            /* Column vector, used in inverting a matrix.       */
          **B_hat2,         /* Local copy of B_hat, to accomodate possibility   */
                            /* of averaging results across runs.                */
           *B_hat_tr_bar_tr,/* Average column of B_hat, (mean(B_hat'))'         */
          **B_hat_r1,       /* Reduced B_hat matrix, for testing dependence on  */
                            /* each independent variable.                       */
          **B_hat_r2,       /* Reduced B_hat matrix, for testing dependence on  */
                            /* independent variables of interest.               */
          **YtrY,           /* Y'*Y                                             */
          **B_hattr_Xtr,    /* B_hat'*X'                                        */
          **B_hattr_Xtr_Y,  /* B_hat'*X'*Y                                      */
          **B_hat_r1tr_X_r1tr,
                            /* B_hat_r1'*X_r1'                                  */
          **B_hat_r1tr_X_r1tr_Y,
                            /* B_hat_r1'*X_r1'*Y                                */
          **B_hat_r2tr_X_r2tr,
                            /* B_hat_r2'*X_r2'                                  */
          **B_hat_r2tr_X_r2tr_Y,
                            /* B_hat_r2'*X_r2'*Y                                */
           *y_bar,          /* Vector of the column means of Y.                 */
          **y_bar_y_bartr,  /* y_bar*y_bar'                                     */
          **E,              /* Error matrix.                                    */
            bigE,           /* Check for zero column in error matrix.           */
            bigEplusH,      /* Check for zero column in E+H matrix.             */
            bigEplusH_r1,   /* Same, for testing individual variables.          */
            bigEplusH_r2,   /* Same, for testing variables of interest.         */
          **EplusH,         /* E + H                                            */
          **EplusH_r1,      /* E + H for the reduced model.                     */
          **EplusH_r2,      /* E + H for the reduced model.                     */
            detE,           /* Determinant of E.                                */
            detEplusH,      /* Determinant of [E + H].                          */
            detEplusH_r1,   /* Same, for testing individual variables.          */
            detEplusH_r2;   /* Same, for testing variables of interest.         */
    float ***B_hat,         /* Least squares estimate of regression             */
                            /* coefficients B.                                  */
           *Wilks_Lambda,   /* Wilks' lambda map.                               */
          **Wilks_Lambda_r1,/* Same, for testing individual variables.          */
           *Wilks_Lambda_r2;/* Same, for testing variables of interest.         */
    float *SE;              /* Sum of errors map.                               */
    float *SSE;             /* SSE map.                                         */
    int E_EplusH_okay,      /* Check for singular E or [E+H] matrices.          */
        EplusH_r1_okay,     /* Same, for checking individual variables.         */
        EplusH_r2_okay;     /* Same, for variables of interest.                 */
    int *index;             /* Index used by Dludcmp routine recording the row  */
                            /* permutation -- whatever that means               */
    int i, j, k;            /* Indices into matrices.                           */
    int first_analysis;     /* Flag to indicate whether or not this is the first*/
                            /* voxel to have multivariate multiple regression   */
                            /* analysis run on it.                              */

                            /* ESTIMATE OF EFFECTIVE DEGREES OF FREEDOM.        */
    char effective_df[2];   /* Flag indicating whether estimating effective     */
                            /* degrees of freedom is desired.                   */
    char output_no_int[2];  /* Flag indicating whether it is desired to output  */
                            /* maps for variables of no interest, in addition   */
                            /* to the maps for variables of interest.           */
    char adjusted_fmri_data[2];/* Flag indicating whether it is desired to      */
                            /* output fMRI data with effects of no interest     */
                            /* removed.                                         */
    float sum_dispersions,  /* Weighted sum of estimates of temporal smoothness.*/
          sum_weights,      /* Sum of weights in weighted average.              */
          max_SD,           /* Greatest temporal dispersion across all          */
                            /* independent variables.                           */
          *overall_HDR;     /* Effective overall hemodynamic response.          */
    double **K,             /* Smoothing matrix.                                */
           **V,             /* K*K'                                             */
           **G_inv_G_tr_G_G_tr,
                            /* G*inv(G'*G)*G'                                   */
           **R,             /* Residual-forming matrix, I-G*inv(G'*G)*G'        */
           **RV,            /* R*V                                              */
           **RVRV;          /* R*V*R*V                                          */
    float trace_RV,         /* trace(R*V)                                       */
          trace_RVRV,       /* trace(R*V*R*V)                                   */
          *nu;              /* Map of effective degrees of freedom.             */

                            /* EXTERNAL VARIABLES, SMOOTHED INPUT SQUAREWAVE    */
                            /* AND HEMODYNAMIC RESPONSE FUNCTION.               */
    FILE *fimgin;           /* File pointer to lag/dispersion map.              */
    int num_read_in;        /* Number of voxels read in.                        */
    float pi;               /* Pi.                                              */
    float sqrt_2pi;         /* sqrt(2.*pi)                                      */
    float SDxSqrt_2pi;      /* SD*sqrt(2.*pi)                                   */
    register int row_index2;/* Index used in convolving of hemodynamic          */
                            /* response with input square wave.                 */
    register int kernel_index; 
                            /* Index used to smooth time series with kernel.    */
    double kernel_index2;   /* Index used to calculate hemodynamic response.    */
    float **HDresponse;     /* Gaussian kernel modeling hemodynamic response    */
                            /* with zdim*num_rows elements.  This number of     */
                            /* elements is needed rather than zdim because we   */
                            /* need to take into account the fact that different*/
                            /* slices within a volume are acquired at different */
                            /* time offsets.                                    */
    float thresh;           /* Threshold value.  Convolution calculations for   */
                            /* values of the hemodynamic response below this    */
                            /* value will not be done, resulting in a signifi-  */
                            /* savings in time, at little expense of floating   */
                            /* point accuracy.                                  */
    int low1,high1,low2;    /* Indices indicating which parts of the hemodyna-  */
                            /* mic response are to be included in calculations. */
    float **lag_map,        /* Lag and Dispersion maps.                         */
        **dispersion_map;
    float *last_lag,        /* Lag of last voxel analyzed.                      */
          *last_dispersion; /* Dispersion of last voxel analyzed.               */
    int hdr_changed;        /* Flag indicating whether the hemodynamic response */
                            /* for at least one independent variable differs    */
                            /* from the last voxel done.  If not, don't need to */
                            /* invert matrix yet again -- it will be exactly    */
                            /* the same.                                        */
    float mean;             /* Mean of Gaussian lagging/smoothing filter.       */
                            /* Will function as the temporal lag.               */
    float SD;               /* SD of Gaussian lagging/smoothing filter.         */
    float SDxSDx2;          /* SD * SD * 2.                                     */

                            /* OUTPUT FILES AND DIRECTORIES.                    */
    float *mean_map;        /* Mean image.                                      */
    FILE *fOUTFILE;         /* File pointer to output text file.                */
    signed short int *Wbuffer;
                            /* fMRI data with effects of no interest removed.   */
    int length;             /* Length of directory portion of pathname.         */
    int str_index;          /* Index for counting characters in a string array. */
    int index2;             /* Indexinto character strings.                     */
    char path[80];          /* Directory portion of pathname.                   */
    char file[80];          /* Filename portion of pathname.                    */

    char outfilename[80];   /* Array to contain name of output text file.       */
    char command[256];      /* UNIX command to be sent to system call.          */

/* ---------------------------------------------------------------------------- */
/*  CHECK TO SEE IF CORRECT NUMBER OF ARGUMENTS WERE PASSED.                    */

    if (argc<15) {
        printf("Dmultivar_multiple_regressionC: need input header file, input textfile\n");
        printf("of file names, number of runs, number of scans per run, name of 8-bit\n");
        printf("mask, total number of independent variables including any constant or\n");
        printf("ramp function), number of independent variables of interest not including\n");
        printf("\"nuisance variables\"), text file containing slice acquisition sequence,\n");
        printf("TR in seconds, and a list of input text files containing independent\n");
        printf("variables, each input filename being followed by a normalization flag\n");
        printf("and information for lag and dispersion.\n");

        printf("All arguments are of course to be separated by spaces.\n\n");

        printf("The input textfile of file names should contain the names of the\n");
        printf("fMRI image volume files, one name per line.  If you have more\n");
        printf("than one fMRI run, just list all the files for the first run,\n");
        printf("then list all the files for the second run, etc., all in the\n");
        printf("same text file, one name per line.\n\n");

        printf("The input text files containing independent variables should\n");
        printf("have two numbers per line.  There should be as many lines in each\n");
        printf("of these text files as there are scans per run.  There can be any\n");
        printf("number of independent variables, but there should be at least one.\n");
        printf("Input textfiles for effects of interest must be two-columns, the first\n");
        printf("to indicate duration in seconds, the second to indicate effect size\n");
        printf("(e.g., second column can be 1's and 0's to indicate on-off squarewave).\n\n");

        printf("If you wish to estimate an intercept (constant or DC offset), include\n");
        printf("a series of 1's as an independent variable.  Consider including a ramp\n");
        printf("function to account for linear drift in the fMRI signal.\n\n");

        printf("If you type in, e.g., 3 as the number of independent variables of,\n");
        printf("interest, the LAST 3 independent variables listed on the command line\n");
        printf("will be tested for overall contribution to the fMRI signal variance.\n\n");

        printf("After the name of each independent variable's input text file,\n");
        printf("type in a \"normalization flag.\"  Set this to 0 if no normalization\n");
        printf("is desired; to 1 if normalizing the mean to zero is desired but no\n");
        printf("normalization of the peak-to-peak amplitude is desired; or to 2 if normalizing both\n");
        printf("the mean to zero and the peak-to-peak amplitude to 1 is desired.\n\n");

        printf("After the normalization flag, type in a lag adjustment flag.\n");
        printf("Set this to 0 if the lag values are to be used as is (\"relative time\");\n");
        printf("to 1 if the lag values input still need to be adjusted for the slight\n");
        printf("time difference due to the slice-by-slice acquisition (\"absolute time\").\n");
        printf("For example, suppose you have two brain areas which are in perfect synchrony,\n");
        printf("but happen to appear in different slices.  If you have typed in a number\n");
        printf("which is to be applied to all brain areas, you will probably want this to\n");
        printf("be adjusted for the slight time difference due to the slice-by-slice\n");
        printf("acquisition, so you'd set this flag to \"1\".  On the other hand, if\n");
        printf("you have typed in the name of a floating point lag map which already\n");
        printf("has been adjusted for the slice-time differences, you'd set this flag\n");
        printf("to \"0\".\n\n");

        printf("After the lag adjustment flag, type in the name of a floating point\n");
        printf("lag map stored as a .img file, and the name of a floating point\n");
        printf("dispersion map.  Alternatively, you can type in a lag (in seconds)\n");
        printf("and a dispersion (also in seconds) to be used for all brain voxels;\n");
        printf("in this case, the time lag due to the slice acquisition sequence\n");
        printf("will be taken into account.\n\n");

        printf("The normalization flag and lag and dispersion should all be 0 for a\n");
        printf("constant or ramp function.  Otherwise, wierd things may happen.\n\n");

        printf("After the last independent variable input text file and normalization\n");
        printf("and lag and dispersion maps, you can optionally type in a list of voxel\n");
        printf("indices of interest.  Intermediate matrices for the calculation of the\n");
        printf("multivariate multiple regression will be output to individual subdirectories\n");
        printf("for each of these voxel indices.  Voxel indices must be zero offset.\n\n");

        exit(1);
    }


/* ---------------------------------------------------------------------------- */
/*  READ IN num_dep_var AND num_rows OFF OF COMMAND LINE.                       */

    num_dep_var          = atoi(argv[3]);
    num_rows             = atoi(argv[4]);
    num_dep_varXnum_rows = num_dep_var*num_rows;
    num_ind_var          = atoi(argv[6]);
    num_ind_var_int      = atoi(argv[7]);
    TR                   = atof(argv[9]);
    pi                   = 2.*acos(0.);
    sqrt_2pi             = (float)sqrt((double)(2.*pi));

    printf("*******************************************************************\n");
    printf("TR = %g seconds.\n",TR);
    printf("There are %d dependent variables.\n", num_dep_var);
    printf("There are %d rows per dependent variable.\n", num_rows);
    printf("This makes a total of %d input scans.\n",num_dep_varXnum_rows);
    printf("Will read in input scan filenames from %s.\n",argv[2]);
    printf("\nThere are to be %d independent variables in the regression\n",num_ind_var);
    printf("(including any constant or ramp function).\n\n");


/* ---------------------------------------------------------------------------- */
/*  READ IN VOXEL INDICES OF INTEREST FROM COMMAND LINE.                        */

    num_vox_int= argc-10-5*num_ind_var;
    vox_int=int_malloc(num_vox_int);

    if (num_vox_int>0) {
        printf("These are the %d voxel indices of interest:\n",num_vox_int);
        for (argv_index=0; argv_index<num_vox_int; argv_index++) {
            vox_int[argv_index]=atoi(argv[10+5*num_ind_var+argv_index]);
            printf("    %d\n",vox_int[argv_index]);
        }
        printf("\n");
    }


/* ---------------------------------------------------------------------------- */
/*  INPUT TIME SAMPLING FACTOR, tick.                                           */

    printf("*******************************************************************\n");
    printf("\nNow you need to type in an integer indicating the time interval\n");
    printf("you'd like to use.  E.g., if you want to represent time in intervals\n");
    printf("in quarter seconds, type in a '4'.\n\n");
    printf("    Temporal up-sampling factor? ");
    fgets(line,10,stdin);
    sscanf(line,"%d",&tick);
/*     scanf("%d\n",&tick); */
    printf("\n\n");
    printf("Will up-sample by a factor of %d\n\n",tick);

/* ---------------------------------------------------------------------------- */
/*  READ HEADER FILE; ASSUME IMAGES ARE SAME SIZE, SO USE SAME HEADER FILE      */
/*  FOR ALL IMAGES.                                                             */

    printf("Reading header info from %s...\n",argv[1]);
    img_info      = readhdr(argv+1);

    xdim          = img_info.dime.dim[1];
    ydim          = img_info.dime.dim[2];
    zdim          = img_info.dime.dim[3];
    zdimXnum_rows = TR*tick*num_rows;
    num_vox       = xdim*ydim*zdim;
    printf("Images have dimensions: [%d %d %d]\n\n",xdim,ydim,zdim);


/* ---------------------------------------------------------------------------- */
/* READ IN SLICE ACQUISITION SEQUENCE.                                          */
/* Convert to numbers from 0 to zdim-1.  If slice k was the nth slice acquired, */
/* set slice_seq2[k-1] equal to (n-1).                                          */

    if ((ftxt = fopen(argv[8],"r")) == NULL) {
        printf("Dmultivar_multiple_regressionC: can't open %s\n",argv[8]);
        exit(2);
    }
    slice_seq1=int_malloc(zdim);
    slice_seq2=int_malloc(zdim);

    printf("*******************************************************************\n");
    printf("Slices were acquired in this sequence: ");
    for (z=0; z<zdim; z++) {
        check_fscanf=fscanf(ftxt,"%d",slice_seq1+z);
        if (check_fscanf==EOF) {
            printf("Dmultivar_multiple_regressionC: error reading slice sequence from %s.\n",
                                                                        argv[8]);
            exit(3);
        }
        printf("%d ",slice_seq1[z]);
    }
    printf("\n");
    for (z=0; z<zdim; z++)
        if ((slice_seq1[z]<1) || (slice_seq1[z]>zdim)) {
            printf("Dmultivar_multiple_regressionC: weird values for slice sequence.\n");
            printf("Have you forgotten to put this information in the input text file?\n");
            exit(4);
         }
    printf("In other words,\n");
    for (z=0; z<zdim; z++) {
        slice_seq2[slice_seq1[z]-1]=z;
        printf("  slice # %d was acquired # %d in sequence\n",
                slice_seq1[z],slice_seq2[slice_seq1[z]-1]+1);
    }
    fclose(ftxt);
    printf("\n\n");


/* ---------------------------------------------------------------------------- */
/*  COUNT NUMBER OF INPUT IMAGE FILENAMES.                                      */

    if ((ftxt = fopen(argv[2],"r")) == NULL) {
        printf("Dmultivar_multiple_regressionC: can't open %s\n",argv[2]);
        exit(5);
    }
    tot_num_scans = 0;
    while(fscanf(ftxt,"%s",dummy)!=EOF)
        tot_num_scans++;

    if (tot_num_scans != num_dep_varXnum_rows) {
        printf("Dmultivar_multiple_regressionC: number of entries in %s not equal to number\n",argv[2]);
        printf("of runs times number of scans per run (%d versus %d).\n",
                                        num_dep_varXnum_rows,tot_num_scans);
        exit(6);
    }
    rewind(ftxt);

/* ---------------------------------------------------------------------------- */
/*  ALLOCATE MEMORY FOR img_name; THEN READ IN NAMES OF INPUT IMAGE FILES.      */

    printf("*******************************************************************\n");
    printf("Total Number of Scans = %d\n",tot_num_scans);

    printf("\nAllocating memory for filenames...\n");
    img_name = char_malloc3(num_dep_var,num_rows,80);
    for (col_index=0; col_index<num_dep_var; col_index++)
        for (row_index=0; row_index<num_rows; row_index++) {
            check_fscanf=fscanf(ftxt,"%s\0",img_name[col_index][row_index]);
            if (check_fscanf==EOF) {
                printf("Dmultivar_multiple_regressionC: error reading run %d timepoint %d filename from %s.\n",
                                                col_index+1,row_index+1,argv[2]);
                exit(7);
            }
        }
    fclose(ftxt);


/* ---------------------------------------------------------------------------- */
/* QUERY USER ABOUT ESTIMATION OF EFFECTIVE DF, AND ABOUT POOLING REGRESSION    */
/* COEFFICIENTS ACROSS RUNS.                                                    */

    printf("\nYou have the option of generating a map of estimated effective\n");
    printf("df, using the method of Worsley and Friston.  This may be very\n");
    printf("computationally intensive!)\n\n");
    printf("Do you want to generate a map of estimated effective df <y/n>? ");
    scanf("%2c",effective_df);
    while((effective_df[0]!='y') && (effective_df[0]!='n') 
                        && (effective_df[0]!='Y') && (effective_df[0]!='N')) {
        printf("Type 'y' or 'n': ");
        scanf("%2c",effective_df);
    }
    printf("\n\n");

    printf("This program will output regression coefficient maps and Wilks'\n");
    printf("Lambda maps for the independent variables of interest.  If you want,\n");
    printf("it will also output regression coefficient maps and Wilks' Lambda maps\n");
    printf("for all other independent variables; this will of course require more\n");
    printf("disk space.\n\n");
    printf("IN ADDITION TO the maps for the variables of interest, do you want to\n");
    printf("output the maps for the variables of NO interest <y/n>? ");
    scanf("%2c",output_no_int);
    while((output_no_int[0]!='y') && (output_no_int[0]!='n') 
                        && (output_no_int[0]!='Y') && (output_no_int[0]!='N')) {
        printf("Type 'y' or 'n': ");
        scanf("%2c",output_no_int);
    }
    printf("\n\n");

    printf("If you want, after all voxels are regressed, this program will output\n");
    printf("the fMRI data with all effects of no interest removed, rounded to the\n");
    printf("nearest 16-bit integer.  The grand global mean will be added back to\n");
    printf("the voxel values to avoid negative numbers.\n\n");
    printf("Do you want to output fMRI data with all effects of no interest removed <y/n>? ");
    scanf("%2c",adjusted_fmri_data);
    while((adjusted_fmri_data[0]!='y') && (adjusted_fmri_data[0]!='n') 
                        && (adjusted_fmri_data[0]!='Y') && (adjusted_fmri_data[0]!='N')) {
        printf("Type 'y' or 'n': ");
        scanf("%2c",adjusted_fmri_data);
    }
    printf("\n\n");

    printf("Do you want scans to be ratio normalized <y/n>? ");
    scanf("%2c",normalize);
    while((normalize[0]!='y') && (normalize[0]!='n')
                        && (normalize[0]!='Y') && (normalize[0]!='N')) {
        printf("Type 'y' or 'n': ");
        scanf("%2c",normalize);
    }
    printf("\n\n");

    printf("Sometimes it is desirable to rerun an analysis, adding more voxels of\n");
    printf("interest, but without reanalyzing all other voxels.\n");
    printf("\nDo you want only voxels of interest to be analyzed <y/n>? ");
    scanf("%2c",only_vox_of_int);
    while((only_vox_of_int[0]!='y') && (only_vox_of_int[0]!='n')
                        && (only_vox_of_int[0]!='Y') && (only_vox_of_int[0]!='N')) {
        printf("Type 'y' or 'n': ");
        scanf("%2c",only_vox_of_int);
    }
    printf("\n\n");



/* ---------------------------------------------------------------------------- */
/* ALLOCATE MEMORY FOR MATRICES.                                                */

    printf("*** zdimXnum_rows = %d ***\n",zdimXnum_rows);
    X_os                       = dmatrix(1,zdimXnum_rows,1,num_ind_var);
    SX                         = dmatrix(1,num_rows,1,num_ind_var);
    SX_os                      = dmatrix(1,zdimXnum_rows,1,num_ind_var);
    Y                          = dmatrix(1,num_rows,1,num_dep_var);
    XtrX                       = dmatrix(1,num_ind_var,1,num_ind_var);
    inv_Xtr_X_Xtr              = dmatrix(1,num_ind_var,1,num_rows);
    YtrY                       = dmatrix(1,num_dep_var,1,num_dep_var);
    B_hat2                     = dmatrix(1,num_ind_var,1,num_dep_var);
    B_hat_tr_bar_tr            = dvector(1,num_ind_var);
    B_hattr_Xtr                = dmatrix(1,num_dep_var,1,num_rows);
    B_hattr_Xtr_Y              = dmatrix(1,num_dep_var,1,num_dep_var);
    y_bar                      = dvector(1,num_dep_var);
    y_bar_y_bartr              = dmatrix(1,num_dep_var,1,num_dep_var);
    E                          = dmatrix(1,num_dep_var,1,num_dep_var);
    EplusH                     = dmatrix(1,num_dep_var,1,num_dep_var);
    index                      = ivector(1,num_ind_var);
    col                        = dvector(1,num_ind_var);

    X_r1                       = dmatrix(1,num_rows,1,num_ind_var-1);
    B_hat_r1                   = dmatrix(1,num_ind_var-1,1,num_dep_var);
    B_hat_r1tr_X_r1tr          = dmatrix(1,num_dep_var,1,num_rows);
    B_hat_r1tr_X_r1tr_Y        = dmatrix(1,num_dep_var,1,num_dep_var);
    EplusH_r1                  = dmatrix(1,num_dep_var,1,num_dep_var);

    X_r2                       = dmatrix(1,num_rows,1,num_ind_var-num_ind_var_int);
    B_hat_r2                   = dmatrix(1,num_ind_var-num_ind_var_int,1,num_dep_var);
    B_hat_r2tr_X_r2tr          = dmatrix(1,num_dep_var,1,num_rows);
    B_hat_r2tr_X_r2tr_Y        = dmatrix(1,num_dep_var,1,num_dep_var);
    EplusH_r2                  = dmatrix(1,num_dep_var,1,num_dep_var);

    printf("\nAllocating memory for regression coefficients maps...\n");
    B_hat                      = float_malloc3(num_ind_var,num_dep_var,num_vox);


    /* ------------------------------------------------------------------------ */
    /* MATRICES FOR ESTIMATE OF EFFECTIVE DEGREES OF FREEDOM.                   */

    K                          = dmatrix(1,num_rows,1,num_rows);
    V                          = dmatrix(1,num_rows,1,num_rows);
    G_inv_G_tr_G_G_tr          = dmatrix(1,num_rows,1,num_rows);
    R                          = dmatrix(1,num_rows,1,num_rows);
    RV                         = dmatrix(1,num_rows,1,num_rows);
    RVRV                       = dmatrix(1,num_rows,1,num_rows);
    overall_HDR                = float_malloc(num_rows);
    nu                         = float_malloc(num_vox);

    for (vox_index=0; vox_index<num_vox; vox_index++)
        nu[vox_index]=0.;


/* ---------------------------------------------------------------------------- */
/* READ IN 8-BIT MASK.                                                          */

    mask = unsigned_char_malloc(num_vox);
    read_unsigned_char_data(argv[5],mask,num_vox);

    num_mask_vox=0;
    for (vox_index=0; vox_index<num_vox; vox_index++)
        if (mask[vox_index]!=0)
            num_mask_vox++;
    printf("\nNumber of nonzero voxels in mask %s = %d.\n\n",argv[5],num_mask_vox);

    mask_voxel=int_malloc(num_mask_vox);
    vox_index2=0;
    for (vox_index=0; vox_index<num_vox; vox_index++)
        if (mask[vox_index]!=0) {
            mask_voxel[vox_index2]=vox_index;
            vox_index2++;
        }

/* ---------------------------------------------------------------------------- */
/* READ IN INDEPENDENT VARIABLES FROM TEXT FILES.                               */

    printf("*******************************************************************\n");
    printf("\nWill read in independent variables from the files:\n");
    for (var_index=0; var_index<num_ind_var; var_index++) 
        printf("  #%d : %d : %d : %s : lag: %s, dispersion: %s\n",
                var_index+1,atoi(argv[11+5*var_index]),atoi(argv[12+5*var_index]),
                argv[10+5*var_index],argv[13+5*var_index],argv[14+5*var_index]);
    printf("\nwhere the number after the first colon is the normalization flag,\n");
    printf("and the number after the second colon is the lag adjustment flag.\n");

    sprintf(outfilename,"wilks_lambda_int_p%dvH%dvE%d.img",num_dep_var,num_ind_var_int,
                                                        num_rows-num_ind_var);
    printf("The last %d of these will be the variables of interest,\n",num_ind_var_int);
    printf("and %s will be the map testing them simultaneously.\n\n",outfilename);

    lag_map                    = float_malloc2(num_ind_var,num_vox);
    dispersion_map             = float_malloc2(num_ind_var,num_vox);
    last_lag                   = float_malloc(num_ind_var);
    last_dispersion            = float_malloc(num_ind_var);
    normalize_flag             = int_malloc(num_ind_var);
    lag_adjust_flag            = int_malloc(num_ind_var);

    for (var_index=0; var_index<num_ind_var; var_index++) {
        printf("  Determining hemodynamic response for independent variable #%d...\n",var_index+1);
        normalize_flag[var_index]  = atoi(argv[11+5*var_index]);
        lag_adjust_flag[var_index] = atoi(argv[12+5*var_index]);

        /* -------------------------------------------------------------------- */
        /* IF PROGRAM IS ABLE TO OPEN LAG MAP FILE, READ IT IN.  IF NOT, ASSUME */
        /* USER HAS TYPED IN A NUMBER INSTEAD, WHICH IS TO BE USED FOR LAG.     */
        /* SAME GOES FOR DISPERSION.                                            */

        if ((fimgin = fopen(argv[13+5*var_index],"r")) == NULL) {
            mean=atof(argv[13+5*var_index]);
            printf("    lag = %g milliseconds\n",mean);
            for (vox_index=0; vox_index<num_vox; vox_index++)
                if (mask[vox_index]!=0) {
                    z=vox_index/(xdim*ydim);
                    lag_map[var_index][vox_index]=mean;
                }
                else
                    lag_map[var_index][vox_index]=0.;
        }
        else {
            printf("    Reading in input floating point lag map from %s...\n",
                                                        argv[13+5*var_index]);
            num_read_in=fread(lag_map[var_index],sizeof(float),num_vox,fimgin);
            fclose(fimgin);
            if (num_read_in!=num_vox) {
                printf("    %d voxels read in from %s, should have been %d.\n",
                                        num_read_in,argv[13+5*var_index],num_vox);
                exit(8);
            }
        }
        if (lag_adjust_flag[var_index]==0) {
/*
            printf("    Converting lag %s from seconds to TR's...\n",argv[13+5*var_index]);
            for (vox_index=0; vox_index<num_vox; vox_index++)
                lag_map[var_index][vox_index]/=TR;
/*
	    ;
        }
        else {
/*
            printf("    Converting lag %s from seconds to TR's,\n",argv[13+5*var_index]);
            printf("    AND adjusting for slice-time acquisition effect.\n");
            for (vox_index=0; vox_index<num_vox; vox_index++)
                lag_map[var_index][vox_index]=(lag_map[var_index][vox_index]/TR)+((float)slice_seq2[z]/(float)zdim);
*/
            printf("    Adjusting lag %s for slice-time acquisition effect.\n",argv[13+5*var_index]);
            for (vox_index=0; vox_index<num_vox; vox_index++)
                lag_map[var_index][vox_index]=(lag_map[var_index][vox_index])+((float)slice_seq2[z]*TR/(float)zdim);
        }

        if ((fimgin = fopen(argv[14+5*var_index],"r")) == NULL) {
            SD=atof(argv[14+5*var_index]);
            printf("    dispersion = %g milliseconds\n",SD);
            for (vox_index=0; vox_index<num_vox; vox_index++)
                if (mask[vox_index]!=0)
/*
                    dispersion_map[var_index][vox_index]=SD/TR;
*/
		    dispersion_map[var_index][vox_index]=SD;
                else
                    dispersion_map[var_index][vox_index]=0.;
        }
        else {
            printf("    Reading in input floating point dispersion map from %s...\n",
                                                        argv[14+5*var_index]);
            num_read_in=fread(dispersion_map[var_index],sizeof(float),num_vox,fimgin);
            fclose(fimgin);
            if (num_read_in!=num_vox) {
                printf("    %d voxels read in from %s, should have been %d.\n",
                                        num_read_in,argv[14+5*var_index],num_vox);
                exit(9);
            }
/*
            for (vox_index=0; vox_index<num_vox; vox_index++)
                dispersion_map[var_index][vox_index]/=TR;
*/
        }
        printf("\n");

        /* ------------------------------------------------------------------------ */
        /* Scale lags and dispersion up by tick; in anticipation of sampling time   */
        /* more finely every TR/tick seconds rather than every TR.                  */
/*
        for (vox_index=0; vox_index<num_vox; vox_index++) {
            lag_map[var_index][vox_index]*=(float)tick;
            dispersion_map[var_index][vox_index]*=(float)tick;
        }
*/
    } /* END LOOP OVER INDEPENDENT VARIABLES. */

    /* ------------------------------------------------------------------------ */
    /* Fill in X_os[][].  X_os is the input stimulus waveform, oversampled in   */
    /* time.  The smoothed version of the oversampled independent variables     */
    /* will be placed in SX_os.  This will in turn be downsampled back to time  */
    /* series of length num_rows, and stored in SX.                             */

    tot_time      = TR*num_rows;
    printf("\n");
    printf("Constructing matrix X of independent variables...\n");
    for (var_index=0; var_index<num_ind_var; var_index++) {
        if ((f_indvar = fopen(argv[10+5*var_index],"r")) == NULL) {
            printf("Dmultivar_multiple_regressionC: can't open %s\n",argv[10+5*var_index]);
            exit(10);
        }
        printf("  Reading in data from %s (independent variable #%d)...\n",
                                                argv[10+5*var_index],var_index+1);
        t_pointer = 0;
        while (((check_fscanf=fscanf(f_indvar,"%g %f",&duration, &float_buffer))!=EOF)
            && (t_pointer+1<zdimXnum_rows)){
            for (t_index=t_pointer;(t_index<zdimXnum_rows) &&
                ((float)t_index<(float)t_pointer+((float)1*duration)); t_index++) {
                X_os[t_index+1][var_index+1]=(double)float_buffer;
            }
            t_pointer=t_index;
        }
        fclose(f_indvar);
        if (t_index!=TR*num_rows*tick) {
            printf("t_index = %d\n",t_index);
            printf("TR*num_rows*tick = %d\n",(int)TR*num_rows*tick);
            printf("Dmultivar_multiple_regressionC: not enough data in %s to cover %g seconds.\n",
                                                                argv[10+5*var_index],tot_time);
            exit(11);
        }
    } /* END LOOP OVER INDEPENDENT VARIABLES. */
    printf("\n");

/*
printf("Writing output to X_os.dat...\n");
sprintf(outfilename,"X_os.dat",x,y,z);
write_NRCdouble_matrix2textfile(outfilename,X_os,zdimXnum_rows,num_ind_var);
exit(222);
*/

/* ---------------------------------------------------------------------------- */
/* ALLOCATE MEMORY FOR HEMODYNAMIC ESTIMATION AND INITIALIZE TO KRONECKER DELTA */

    HDresponse = float_malloc2(num_ind_var,zdimXnum_rows);
    for (var_index=0; var_index<num_ind_var; var_index++) {
        HDresponse[var_index][0]=1.;
        for (row_index=1; row_index<zdimXnum_rows; row_index++)
            HDresponse[var_index][row_index]=0.;
    }


/* ---------------------------------------------------------------------------- */
/*  ALLOCATE MEMORY FOR WILKS' LAMBDA AND SSE MAPS.                             */

    Wilks_Lambda=float_malloc(num_vox);
    for (vox_index=0; vox_index<num_vox; vox_index++)
        Wilks_Lambda[vox_index]=0.;

    Wilks_Lambda_r2=float_malloc(num_vox);
    for (vox_index=0; vox_index<num_vox; vox_index++)
        Wilks_Lambda_r2[vox_index]=0.;

    Wilks_Lambda_r1=float_malloc2(num_ind_var,num_vox);
    for (var_index=0; var_index<num_ind_var; var_index++)
        for (vox_index=1; vox_index<num_vox; vox_index++)
            Wilks_Lambda_r1[var_index][vox_index]=0.;

    SE=float_malloc(num_vox);
    for (vox_index=0; vox_index<num_vox; vox_index++)
        SE[vox_index]=0.;

    SSE=float_malloc(num_vox);
    for (vox_index=0; vox_index<num_vox; vox_index++)
        SSE[vox_index]=0.;


/* ---------------------------------------------------------------------------- */
/*  READ IN IMAGES.                                                             */

    temp = signed_short_malloc(num_vox);
    printf("*******************************************************************\n");
    printf("Allocating memory for fMRI data...\n");
    imgin = signed_short_malloc3(num_mask_vox,num_dep_var,num_rows);

/*     printf("\nReading in fMRI scans...\n       "); */
    printf("\nReading in fMRI scans...\n");
    num_files_read_in=0;
    for (col_index=0; col_index<num_dep_var; col_index++)
        for (row_index=0; row_index<num_rows; row_index++) {
            num_files_read_in++;
/*            printf("\b\b\b\b\b\b\b%5.1f %%",(float)num_files_read_in*100./(float)tot_num_scans); */
            read_signed_short_int_data(img_name[col_index][row_index],temp,num_vox);

            for (vox_index=0; vox_index<num_mask_vox; vox_index++)
                imgin[vox_index][col_index][row_index]=temp[mask_voxel[vox_index]];
        }
    printf("\n\n");
    free(temp);


/* ---------------------------------------------------------------------------- */
/* CALCULATE GLOBAL INTENSITIES AND GRAND GLOBAL MEAN.                          */

    global_mean=float_malloc2(num_dep_var,num_rows);
    printf("*******************************************************************\n");

    if ((normalize[0]=='y') || (normalize[0]=='Y')) {
        printf("Calculating scan global means and grand global mean...\n");
        grand_global_mean=0.;
        count=0;
        for (col_index=0; col_index<num_dep_var; col_index++)
            for (row_index=0; row_index<num_rows; row_index++) {
                count++;
/*            printf("\b\b\b\b\b\b\b%5.1f %%",(float)count*100./(float)(num_dep_var*num_rows)); */
                global_mean[col_index][row_index]=0.;
                for (vox_index=0; vox_index<num_mask_vox; vox_index++)
                    global_mean[col_index][row_index]
                                    +=(float)imgin[vox_index][col_index][row_index];
                global_mean[col_index][row_index]/=(float)num_mask_vox;
                grand_global_mean+=global_mean[col_index][row_index];
            }
        grand_global_mean/=(float)tot_num_scans;

        printf("\nGrand global mean = %g...\n\n",grand_global_mean);
    }
    else {
        for (col_index=0; col_index<num_dep_var; col_index++)
            for (row_index=0; row_index<num_rows; row_index++)
                global_mean[col_index][row_index]=1.;
        grand_global_mean=1.;
    }


/* ---------------------------------------------------------------------------- */
/* CALCULATE MEAN ACROSS ALL SCANS.                                             */

    mean_map=float_malloc(num_mask_vox);
    for (vox_index=0; vox_index<num_mask_vox; vox_index++)
        mean_map[vox_index]=0.;

    for (col_index=0; col_index<num_dep_var; col_index++)
        for (row_index=0; row_index<num_rows; row_index++)
            for (vox_index=0; vox_index<num_mask_vox; vox_index++)
                mean_map[vox_index]+=(float)imgin[vox_index][col_index][row_index];

    for (vox_index=0; vox_index<num_mask_vox; vox_index++)
        mean_map[vox_index]/=(float)(num_dep_var*num_rows);

/* ---------------------------------------------------------------------------- */
/* INITIALIZE B_hat AND Wilks_Lambda TO ZERO.                                   */

    for (vox_index=0; vox_index<num_vox; vox_index++) {
        for (var_index=0; var_index<num_ind_var; var_index++) {
            Wilks_Lambda_r1[var_index][vox_index]=0.;
            for (col_index=0; col_index<num_dep_var; col_index++)
                B_hat[var_index][col_index][vox_index]=0.;
        }
        Wilks_Lambda[vox_index]=0.;
        Wilks_Lambda_r2[vox_index]=0.;
    } /* LOOP OVER VOXELS. */


/* ---------------------------------------------------------------------------- */
/* LOOP OVER VOXELS, PERFORM VOXEL-WISE MULTIVARIATE MULTIPLE REGRESSION.       */

    printf("*******************************************************************\n");
/*    printf("Now doing multivariate multiple regression on each voxel...\n       "); */
    printf("Now doing multivariate multiple regression on each voxel...\n");
    first_analysis=1;
    for (vox_index2=0; vox_index2<num_mask_vox; vox_index2++) {

        vox_index=mask_voxel[vox_index2];

        /* ---------------------------------------------------------------- */
        /* IF IT IS DESIRED TO ANALYZE ONLY VOXELS OF INTEREST, SKIP THIS   */
        /* VOXEL IF IT IS NOT A VOXEL OF INTEREST.                          */

        yes_voxel_of_int=0;
        for (argv_index=0; argv_index<num_vox_int; argv_index++)
            if (vox_int[argv_index]==vox_index)
                yes_voxel_of_int=1;

        if ((((only_vox_of_int[0]=='y') || (only_vox_of_int[0]=='Y'))
            && (yes_voxel_of_int==1))
            || ((only_vox_of_int[0]=='n') || (only_vox_of_int[0]=='N'))) {


/*        printf("\b\b\b\b\b\b\b%5.1f %%",(float)vox_index*100./(float)num_vox); */
/*            printf("\n\nvox_index = %d\n\n       ",vox_index); */

            /* ---------------------------------------------------------------- */
            /* DETERMINE SLICE NUMBER OF CURRENT VOXEL.  THEN LOOP OVER         */
            /* INDEPENDENT VARIABLES, SMOOTHING AND LAGGING THEM AS SPECIFIED.  */

            z=vox_index/(xdim*ydim);

            hdr_changed=0;
            for (var_index=0; var_index<num_ind_var; var_index++) {

                /* ------------------------------------------------------------ */
                /* IF LAG AND DISPERSION ARE EXACTLY THE SAME AS PREVIOUS       */
                /* VOXEL'S, DON'T NEED TO RECALCULATE IDEALIZED RESPONSE.       */
                /* IF, HOWEVER, THIS IS THE FIRST VOXEL TO BE ANALYZED WITH THE */
                /* MULTIVARIATE MULTIPLE REGRESSION, WE'LL NEED TO GENERATE THE */
                /* SMOOTHED/LAGGED VERSION OF THE INDEPENDENT VARIABLES.        */

                if ((first_analysis==1) || 
                                        ((lag_map[var_index][vox_index]
                                                    !=last_lag[var_index])
                                        || (dispersion_map[var_index][vox_index]
                                                    !=last_dispersion[var_index]))) {

                    hdr_changed=1;
                    last_lag[var_index]        = lag_map[var_index][vox_index];
                    last_dispersion[var_index] = dispersion_map[var_index][vox_index];

                    /* -------------------------------------------------------- */
                    /* IF LAG AND DISPERSION ARE NOT BOTH ZERO, NEED TO         */
                    /* CALCULATE HEMODYNAMIC RESPONSE.                          */

                    if ((lag_map[var_index][vox_index]!=0.)
                                || (dispersion_map[var_index][vox_index]!=0.)) {

                        /* ---------------------------------------------------- */
                        /* CALCULATE GAUSSIAN LAGGING/SMOOTHING KERNEL.         */
                        /* Gaussian smoothing kernel used to used to lag and    */
                        /* disperse  the independent variable waveforms.  The   */
                        /* calculation of the smoothing kernel and the          */
                        /* convolution will be done in       oversampled time.  */
                        /* The resultant lagged and dispersed input waveforms   */
                        /* will then be downsampled back to the original time   */
                        /* series length, and then passed to the multivariate   */
                        /* multiple regression.                                 */

                        mean=lag_map[var_index][vox_index];
                        SD=dispersion_map[var_index][vox_index];
printf("%d lag = %g, disp = %g\n",var_index,mean,SD);
                        SDxSDx2=SD*SD*2.;
                        SDxSqrt_2pi=SD*sqrt_2pi;

                        for (row_index=0; row_index<zdimXnum_rows; row_index++) {
                            kernel_index2=(float)row_index;
                            if ((float)row_index>((float)(zdimXnum_rows)-1.)/2.)
                                kernel_index2=(float)(zdimXnum_rows-row_index);
                            if (SD>0.1)
                                HDresponse[var_index][row_index]= (1/(SDxSqrt_2pi))
                                                    *
                                    exp((double)-((kernel_index2*kernel_index2))/SDxSDx2);
                                                else {
                                if (row_index==0)
                                    HDresponse[var_index][row_index]=1.;
                                else
                                    HDresponse[var_index][row_index]=0.;
                            }
                        }

                        /* -------------------------------------------------------- */
                        /* CONVOLVE GAUSSIAN LAGGING/SMOOTHING KERNEL WITH INPUT    */
                        /* WAVEFORM, GENERATING IDEALIZED RESPONSE FUNCTION, PUT    */
                        /* INTO SX[][].  Done in time domain.  Could have done in   */
                        /* spectral domain, but that would require doing FT's.      */
                        /* Circular convolution is being done, although this could  */
                        /* be changed.  Using Jim's ideas of incrementing row_index */
                        /* by zdim, for a factor of zdim times increase in speed,   */
                        /* and of "windowing" the convolution to do calculations    */
                        /* only for values of the hemodynamic response significantly*/
                        /* greater than some threshold.  Here, that threshold is    */
                        /* to be 0.0001 times the max value in the hemodynamic      */
                        /* response.                                                */

                        thresh=HDresponse[var_index][0];
                        for (row_index=0; row_index<zdimXnum_rows; row_index++)
                            thresh=(HDresponse[var_index][row_index]>thresh)
                                        ? HDresponse[var_index][row_index] : thresh;
                        thresh*=0.0001;

                        low1      = 0;
                        row_index = 0;
                        while ((row_index<zdimXnum_rows)
                                        && (HDresponse[var_index][row_index]>=thresh))
                            row_index=row_index+1;
                        high1=row_index-1;
                    
                        if (high1<zdimXnum_rows) {
                            row_index=high1+1;
                            while((row_index<zdimXnum_rows)
                                        && (HDresponse[var_index][row_index]<thresh))
                                row_index = row_index+1;
                            low2=row_index;
                        }
                        else {
                            low2 = zdimXnum_rows+1;
                        }

                        for (row_index=0; row_index<zdimXnum_rows; row_index+=TR*tick) {
                            SX_os[row_index+1][var_index+1]=0.;
                            for (row_index2=low1; row_index2<=high1; row_index2++) {
                                kernel_index = row_index - row_index2 -mean ;
                                while (kernel_index<0)
                                    kernel_index += zdimXnum_rows;
                                while (kernel_index>=zdimXnum_rows)
                                    kernel_index -= zdimXnum_rows;
                                SX_os[row_index+1][var_index+1] +=
                                    X_os[kernel_index+1][var_index+1]
                                                *HDresponse[var_index][row_index2];
                            }
                            for (row_index2=low2; row_index2<zdimXnum_rows; row_index2++) {
                                kernel_index = row_index - row_index2 -mean ;
                                while (kernel_index<0)
                                    kernel_index += zdimXnum_rows;
                                while (kernel_index>=zdimXnum_rows)
                                    kernel_index -= zdimXnum_rows;
                                SX_os[row_index+1][var_index+1] +=
                                    X_os[kernel_index+1][var_index+1]
                                                *HDresponse[var_index][row_index2];
                            }

                        }
/*
                        for (row_index=0; row_index<zdimXnum_rows; row_index+=TR*tick) {
                            SX_os[row_index+1][var_index+1]=0.;
                            for (row_index2=0; row_index2<zdimXnum_rows; row_index2++) {
                                kernel_index = row_index - row_index2;
                                if (kernel_index<0) kernel_index += zdimXnum_rows;
                                SX_os[row_index+1][var_index+1] +=
                                    X_os[row_index2+1][var_index+1]
                                            *HDresponse[var_index][kernel_index];
                            }
                        }
*/
                    } /* END IF BLOCK TESTING WHETHER HDR ESTIMATE IS NOT A         */
                      /* KRONECKER DELTA FUNCTION.                                  */
                    else {
                        /* -------------------------------------------------------- */
                        /* HEMODYNAMIC RESPONSE IS NOT TO BE ESTIMATED FOR CURRENT  */
                        /* INDEPENDENT VARIABLE, LOAD DIRECTLY INTO SX_os[][] WITH  */
                        /* NO SMOOTHING.                                            */

                        for (row_index=1; row_index<=zdimXnum_rows; row_index++)
                            SX_os[row_index][var_index+1]=X_os[row_index][var_index+1];

                    } /* END ELSE BLOCK TESTING FOR WHETHER HDR IS DELTA FUNCTION. */

                    /* ------------------------------------------------------------ */
                    /* DOWN SAMPLE FROM SX_os[][] INTO SX[][].                      */

                    for (row_index=0; row_index<num_rows; row_index++)
                        SX[row_index+1][var_index+1]=SX_os[(int)TR*tick*row_index+1][var_index+1];

                } /* END IF BLOCK TESTING WHETHER HEMODYNAMIC RESPONSE DIFFERS  */
                  /* FROM LAST VOXEL.                                           */

            } /* END LOOP OVER INDEPENDENT VARIABLES. */

/*
sprintf(outfilename,"SX_os.dat",x,y,z);
write_NRCdouble_matrix2textfile(outfilename,SX_os,zdimXnum_rows,num_ind_var);

sprintf(outfilename,"SX.dat",x,y,z);
write_NRCdouble_matrix2textfile(outfilename,SX,num_rows,num_ind_var);

sprintf(outfilename,"hdr.dat");
if ((fOUTFILE = fopen(outfilename,"w")) == NULL) {
    printf("Dmultivar_multiple_regressionC: can't open %s.\n",outfilename);
    exit(15);
}
for(row_index=0;row_index<zdimXnum_rows;row_index++)
    fprintf(fOUTFILE,"%10f\n",HDresponse[num_ind_var-1][row_index]);
fclose(fOUTFILE);
exit(47);
*/
            /* ---------------------------------------------------------------- */
            /* IF HEMODYNAMIC RESPONSE ESTIMATE HAS NOT CHANGED FOR ANY         */
            /* INDEPENDENT VARIABLE SINCE THE LAST VOXEL, DON'T NEED TO         */
            /* RECALCULATE INVERSE OF (X'*X).  IN THAT CASE, CAN SKIP THIS      */
            /* SECTION, SAVING LOTS OF COMPUTATIONS.                            */

            if ((hdr_changed==1) || (first_analysis==1)) {

                /* -------------------------------------------------------- */
                /* NORMALIZE MEAN AND VARIANCE OF INDEPENDENT VARIABLES.    */

                for (var_index=0; var_index<num_ind_var; var_index++) {
                    if ((normalize_flag[var_index]==1) || (normalize_flag[var_index]==2)) {
                        var_mean=0.;
                        for (row_index=0; row_index<num_rows; row_index++)
                            var_mean+=SX[row_index+1][var_index+1];
                        var_mean/=(float)num_rows;
                        for (row_index=0; row_index<num_rows; row_index++)
                            SX[row_index+1][var_index+1]-=var_mean;
                    }

                    if (normalize_flag[var_index]==2) {
                        max=SX[1][var_index+1];
                        min=SX[1][var_index+1];
                        for (row_index=0; row_index<num_rows; row_index++) {
                            max=(SX[row_index+1][var_index+1]>max)
                                        ? SX[row_index+1][var_index+1] : max;
                            min=(SX[row_index+1][var_index+1]<min)
                                        ? SX[row_index+1][var_index+1] : min;
                        }
                        peak2peak=max-min;
                        for (row_index=0; row_index<num_rows; row_index++)
                            SX[row_index+1][var_index+1]/=peak2peak;
                    }
                }

                /* -------------------------------------------------------- */
                /* CALCULATE inv_Xtr_X_Xtr.  THIS WILL INVOLVE INVERTING    */
                /* XtrX.  WHEN MULTIPLIED BY Y, inv_Xtr_X_Xtr WILL GIVE     */
                /* B_HAT, THE LEAST SQUARES ESTIMATE OF THE REGRESSION      */
                /* COEFFICIENTS, B.                                         */

                /* ******************************************************** */
                /* CALCULATE X'*X (XtrX).                                   */

                for(i=1;i<=num_ind_var;i++)
                    for(j=1;j<=num_ind_var;j++) {
                        XtrX[i][j]=0.;
                        for(k=1;k<=num_rows;k++)
                            XtrX[i][j]+=SX[k][i]*SX[k][j];
                    }


                /* ******************************************************** */
                /* PERFORM LU DECOMPOSITION ON (X'*X).                      */

                Dludcmp(XtrX,num_ind_var,index,&d);

                /* ******************************************************** */
                /* CALCULATE inv(X'*X)*X', inv_Xtr_X_Xtr.  WE CAN CALCULATE */
                /* inv(X'*X)*X' DIRECTLY BY BACKSUBSTITUTING WITH THE       */
                /* COLUMNS OF X', RATHER THAN WITH UNIT VECTORS.  THIS      */
                /* SAVES COMPUTATIONS INVOLVING MATRIX MULTIPLICATION, AND  */
                /* IS MORE ACCURATE.  inv(X'*X) IS NOT EXPLICITLY           */
                /* CALCULATED, BUT IS NOT NEEDED.  Reference: Press WH,     */
                /* Teukolsky SA, Vetterling WT, Flannery BP, Numerical      */
                /* Recipes in C, 2nd edition, Cambridge:Cambridge           */
                /* University Press (1992), top of page 49.                 */

                for(j=1;j<=num_rows;j++) {
                    for(i=1;i<=num_ind_var;i++) col[i]=SX[j][i];
                    Dlubksb(XtrX,num_ind_var,index,col);
                    for(i=1;i<=num_ind_var; i++) inv_Xtr_X_Xtr[i][j]=col[i];
                }
            }
            first_analysis=0;


            /* ---------------------------------------------------------------- */
            /* LOAD CURRENT VOXEL'S fMRI SIGNAL INTO MATRIX Y.                  */

            for (col_index=0; col_index<num_dep_var; col_index++)
                for (row_index=0; row_index<num_rows; row_index++)
                    Y[row_index+1][col_index+1]
                        =(float)imgin[vox_index2][col_index][row_index]
                                *grand_global_mean/global_mean[col_index][row_index];


            /* ---------------------------------------------------------------- */
            /* MULTIPLY inv_Xtr_X_Xtr BY Y TO DO REGRESSION.                    */

            for (var_index=0; var_index<num_ind_var; var_index++) {
                for (col_index=0; col_index<num_dep_var; col_index++) {
                    B_hat[var_index][col_index][vox_index]=0.;
                    for (row_index=0; row_index<num_rows; row_index++)
                        B_hat[var_index][col_index][vox_index]+=
                                inv_Xtr_X_Xtr[var_index+1][row_index+1]
                                                    *Y[row_index+1][col_index+1];
/*                  if (isinf(B_hat[var_index][col_index][vox_index]))
                        printf("B_hat[%d][%d][%d] is Inf!\n    ",
                                                var_index,col_index,vox_index);
                    if (isnan(B_hat[var_index][col_index][vox_index]))
                        printf("B_hat[%d][%d][%d] is NaN!\n    ",
                                                var_index,col_index,vox_index);
*/
                } /* END col_index LOOP */
            } /* END var_index LOOP */


            /* ---------------------------------------------------------------- */
            /* IF AVERAGING REGRESSION COEFFICIENTS ACROSS RUN IS DESIRED, PUT  */
            /* AVERAGE COEFFICIENTS IN B_hat2; OTHERWISE PUT "RAW" COEFFICIENTS */
            /* IN B_hat2.

            COMMENTED THIS OUT.  This is unnecessary.  To perform this function,
            just concatenate fMRI runs into one long fMRI run.

            printf("Check for pooling across runs...\n");

            if ((pool_across_runs=='y') || (pool_across_runs=='Y')) {
                for (var_index=0; var_index<num_ind_var; var_index++) {
                    B_hat_tr_bar_tr[var_index+1]=0.;
                    for (col_index=0; col_index<num_dep_var; col_index++)
                        B_hat_tr_bar_tr[var_index+1]+=
                                        B_hat[var_index][col_index][vox_index];
                    B_hat_tr_bar_tr[var_index+1]/=(float)num_dep_var;
                }
                for (var_index=1; var_index<=num_ind_var; var_index++)
                    for (col_index=1; col_index<=num_dep_var; col_index++)
                        B_hat2[var_index][col_index]=B_hat_tr_bar_tr[var_index];
            } 
            ******************************************************************* */

                for (var_index=1; var_index<=num_ind_var; var_index++)
                    for (col_index=1; col_index<=num_dep_var; col_index++)
                        B_hat2[var_index][col_index]=
                                        B_hat[var_index-1][col_index-1][vox_index];


            /* ---------------------------------------------------------------- */
            /* CALCULATE E, and EplusH.                                         */
            /* Commented out calls to routines which perform matrix             */
            /* multiplication and doing the multiplication explicitly, for      */
            /* speed.

            mult_NRCdouble_matrices_tr1(Y,Y,num_rows,num_dep_var,num_dep_var,YtrY);
            mult_NRCdouble_matrices_tr_both(B_hat2,SX,num_ind_var,num_dep_var,
                                                        num_rows,B_hattr_Xtr);
            mult_NRCdouble_matrices(B_hattr_Xtr,Y,num_dep_var,num_rows,
                                                        num_dep_var,B_hattr_Xtr_Y);
              ---------------------------------------------------------------- */
            for(i=1;i<=num_dep_var;i++)
                for(j=1;j<=num_dep_var;j++) {
                    YtrY[i][j]=0.;
                    for (row_index=1; row_index<=num_rows; row_index++)
                        YtrY[i][j]+=Y[row_index][i]*Y[row_index][j];
                }

            for(i=1;i<=num_dep_var;i++)
                for(j=1;j<=num_rows;j++) {
                    B_hattr_Xtr[i][j]=0.;
                    for (var_index=1; var_index<=num_ind_var; var_index++)
                        B_hattr_Xtr[i][j]+=B_hat2[var_index][i]*SX[j][var_index];
                }

            for(i=1;i<=num_dep_var;i++)
                for(j=1;j<=num_dep_var;j++) {
                    B_hattr_Xtr_Y[i][j]=0.;
                    for (row_index=1; row_index<=num_rows; row_index++)
                        B_hattr_Xtr_Y[i][j]+=B_hattr_Xtr[i][row_index]*Y[row_index][j];
                }

            for(i=1;i<=num_dep_var;i++) {
                y_bar[i]=0.;
                for (row_index=0; row_index<num_rows; row_index++)
                    y_bar[i]+=Y[row_index+1][i];
                y_bar[i]/=(float)num_rows;
            }

            for(i=1;i<=num_dep_var;i++)
                for(j=1;j<=num_dep_var;j++)
                    y_bar_y_bartr[i][j]=y_bar[i]*y_bar[j];

            for(i=1;i<=num_dep_var;i++)
                for(j=1;j<=num_dep_var;j++)
                    E[i][j]=YtrY[i][j]-B_hattr_Xtr_Y[i][j];

            for(i=1;i<=num_dep_var;i++)
                for(j=1;j<=num_dep_var;j++)
                    EplusH[i][j]=YtrY[i][j]-((float)num_rows*y_bar_y_bartr[i][j]);




            /* ---------------------------------------------------------------- */
            /* CALCULATE det(E) AND det(EplusH).                                */

            E_EplusH_okay=1;
            for(i=1;i<=num_dep_var;i++) {
                bigE=0.;
                for(j=1;j<=num_dep_var;j++)
                    bigE=(fabs(E[i][j])>bigE) ? fabs((double)E[i][j]) : bigE;
                if (bigE==0.) {
                    E_EplusH_okay=0;
                    printf("\n  Singular error matrix at vox_index=%d\n",vox_index);
                    printf("  Will set corresponding Wilks' lambda to zero.\n\n       ");
                }
            }

            for(i=1;i<=num_dep_var;i++) {
                bigEplusH=0.;
                for(j=1;j<=num_dep_var;j++)
                    bigEplusH=(fabs(EplusH[i][j])>bigEplusH) ?
                                        fabs((double)EplusH[i][j]) : bigEplusH;
                if (bigEplusH==0.) {
                    E_EplusH_okay=0;
                    printf("\n  Singular [E+H] matrix at vox_index=%d\n",vox_index);
                    printf("  Will set corresponding Wilks' lambda to zero.\n\n       ");
                }
            }
            if (E_EplusH_okay) {
                Dludcmp(E,num_dep_var,index,&detE);
                for (col_index=1; col_index<=num_dep_var; col_index++)
                    detE*=E[col_index][col_index];
                if (detE<0.) {
                    printf("\n  Negative det(E) generated at vox_index = %d\n",vox_index);
                    printf("  Will use absolute value.\n\n       ");
                    detE=fabs(detE);
                }
                Dludcmp(EplusH,num_dep_var,index,&detEplusH);
                for (col_index=1; col_index<=num_dep_var; col_index++)
                    detEplusH*=EplusH[col_index][col_index];
                if (detEplusH<0.) {
                    printf("\n\n  Negative det(E+H) generated at vox_index = %d\n",vox_index);
                    printf("  Will use absolute value.\n\n       ");
                    detEplusH=fabs(detEplusH);
                }
            }

            SSE[vox_index]=detE;


            /* ---------------------------------------------------------------- */
            /* CALCULATE WILKS' LAMBDA (Rencher 1995, section 10.5.1)           */

            if ((E_EplusH_okay!=0.) && (detEplusH!=0.))
                Wilks_Lambda[vox_index]=detE/detEplusH;
            else
                Wilks_Lambda[vox_index]=0.;


/*            if (isinf(Wilks_Lambda[vox_index])) {
                printf("\n\n  Infinite Wilks' Lambda generated at vox_index = %d...\n",
                                                                        vox_index);
                printf("  Will set it to zero.\n");
                Wilks_Lambda[vox_index]=0.;

                printf("  Writing Y'*Y to YtrY.dat...\n");
                write_NRCdouble_matrix2textfile("YtrY.dat",YtrY,num_dep_var,num_dep_var);

                printf("  Writing mean(Y)'*mean(Y) to y_bar_y_bartr.dat...\n");
                write_NRCdouble_matrix2textfile("y_bar_y_bartr.dat",y_bar_y_bartr,
                                                                num_dep_var,num_dep_var);

                printf("  Writing E to E.dat...\n");
                write_NRCdouble_matrix2textfile("E.dat",E,num_dep_var,num_dep_var);

                printf("  Writing ^B to B_hat.dat...\n");
                if ((fOUTFILE = fopen("B_hat.dat","w")) == NULL) {
                    printf("Dmultivar_multiple_regressionC: can't open B_hat.dat.\n");
                    exit(12);
                }
                for(i=0;i<num_ind_var;i++) {
                    for(j=0;j<num_dep_var;j++)
                        fprintf(fOUTFILE,"%10f ",B_hat[i][j][vox_index]);
                    fprintf(fOUTFILE,"\n");
                }
                fclose(fOUTFILE);

                printf("  Writing ^B mean to B_hat2.dat...\n");
                write_NRCdouble_matrix2textfile("B_hat2.dat",
                                                        B_hat2,num_ind_var,num_dep_var);

                printf("  Writing ^B'*X'*Y to B_hattr_Xtr_Y.dat...\n");
                write_NRCdouble_matrix2textfile("B_hattr_Xtr_Y.dat",B_hattr_Xtr_Y,
                                                            num_dep_var,num_dep_var);

                printf("  Writing ^B'*X' to B_hattr_Xtr.dat...\n\n       ");
                write_NRCdouble_matrix2textfile("B_hattr_Xtr.dat",B_hattr_Xtr,
                                                                num_dep_var,num_rows);
            }
*/



            /* ---------------------------------------------------------------- */
            /* TEST DEPENDENCE ON EACH INDEPENDENT VARIABLE                     */
            /* (Rencher 1995, section 10.5.2)                                   */
            /* Outer loop over var_index, independent variable being tested.    */

            for (var_index=1; var_index<=num_ind_var; var_index++) {

                /* ------------------------------------------------------------ */
                /* GENERATE B_hat_r1 AND X_r1, MATRICES FOR THE REDUCED MODEL.  */

/*                printf("\nGenerating B_hat_r1 and X_r1...\n"); */
                for (i=1; i<=num_ind_var-1; i++) {
                    if (i<var_index) {
                        for (col_index=1; col_index<=num_dep_var; col_index++)
                            B_hat_r1[i][col_index]=B_hat2[i][col_index];
                        for (row_index=1; row_index<=num_rows; row_index++)
                            X_r1[row_index][i]=SX[row_index][i];
                    }
                    else {
                        for (col_index=1; col_index<=num_dep_var; col_index++)
                            B_hat_r1[i][col_index]=B_hat2[i+1][col_index];
                        for (row_index=1; row_index<=num_rows; row_index++)
                            X_r1[row_index][i]=SX[row_index][i+1];
                    }
                } /* END i LOOP, LOOP OVER ROWS IN REDUCED B MATRIX. */

                /* ------------------------------------------------------------ */
                /* CALCULATE EplusH_r1.                                         */
                /* Commented out calls to routines which perform matrix         */
                /* multiplication and doing the multiplication explicitly, for  */
                /* speed.
                mult_NRCdouble_matrices_tr_both(B_hat_r1,X_r1,num_ind_var-1,num_dep_var,
                                                        num_rows,B_hat_r1tr_X_r1tr);
                mult_NRCdouble_matrices(B_hat_r1tr_X_r1tr,Y,num_dep_var,num_rows,
                                                        num_dep_var,B_hat_r1tr_X_r1tr_Y);
                   ------------------------------------------------------------ */
/*                printf("Calculating B_hat_r1tr_X_r1tr...\n"); */
                for(i=1;i<=num_dep_var;i++)
                    for(j=1;j<=num_rows;j++) {
                        B_hat_r1tr_X_r1tr[i][j]=0.;
                        for (k=1; k<=num_ind_var-1; k++)
                            B_hat_r1tr_X_r1tr[i][j]+=B_hat_r1[k][i]*X_r1[j][k];
                    }

/*               printf("Calculating B_hat_r1tr_X_r1tr_Y...\n"); */
                for(i=1;i<=num_dep_var;i++)
                    for(j=1;j<=num_dep_var;j++) {
                        B_hat_r1tr_X_r1tr_Y[i][j]=0.;
                        for (row_index=1; row_index<=num_rows; row_index++)
                            B_hat_r1tr_X_r1tr_Y[i][j]+=
                                        B_hat_r1tr_X_r1tr[i][row_index]*Y[row_index][j];
                    }

/*                printf("Calculating EplusH_r1...\n"); */
                for(i=1;i<=num_dep_var;i++)
                    for(j=1;j<=num_dep_var;j++)
                        EplusH_r1[i][j]=YtrY[i][j]-B_hat_r1tr_X_r1tr_Y[i][j];


                /* ------------------------------------------------------------- */
                /* CALCULATE det(EplusH_r1)                                      */

                EplusH_r1_okay=1;
                for(i=1;i<=num_dep_var;i++) {
                    bigEplusH_r1=0.;
                    for(j=1;j<=num_dep_var;j++)
                        bigEplusH_r1=(fabs(EplusH_r1[i][j])>bigEplusH_r1) ?
                                            fabs((double)EplusH_r1[i][j]) : bigEplusH_r1;
                    if (bigEplusH_r1==0.) {
                        E_EplusH_okay=0;
			/* JOE MAISOG : 08-22-2005 : Need to unset EplusH_r1_okay also! */
                        EplusH_r1_okay=0;
                        printf("\n  Singular reduced [E+H] matrix in independent variable #%d, at vox_index=%d\n",
                                                                var_index+1,vox_index);
                        printf("  Will set corresponding Wilks' lambda to zero.\n\n       ");
                    }
                }
                if (EplusH_r1_okay) {
                    Dludcmp(EplusH_r1,num_dep_var,index,&detEplusH_r1);
/*                    printf("Calculating detEplusH_r1...\n"); */
                    for (col_index=1; col_index<=num_dep_var; col_index++)
                        detEplusH_r1*=EplusH_r1[col_index][col_index];
                    if (detEplusH_r1<0.) {
/*                         printf("\n\n  Negative reduced det(E+H) generated at vox_index = %d\n",vox_index);
                        printf("  Will use absolute value.\n\n       ");
*/
                        detEplusH_r1=fabs(detEplusH_r1);
                    }
                }

                /* ------------------------------------------------------------ */
                /* CALCULATE WILKS' LAMBDA FOR REDUCED MODEL.                   */

/*                printf("Calculating Wilks_Lambda_r1[%d]...\n",var_index+1); */
                if ((EplusH_r1_okay!=0) && (detEplusH_r1!=0.))
/*                        && (!isinf(detE/detEplusH_r1)) && (!isnan(detE/detEplusH_r1))) */
                    Wilks_Lambda_r1[var_index-1][vox_index]=detE/detEplusH_r1;
                else
                    Wilks_Lambda_r1[var_index-1][vox_index]=0.;


                if (Wilks_Lambda_r1[var_index-1][vox_index]>1.)
                    Wilks_Lambda_r1[var_index-1][vox_index]=1.;


            } /* END LOOP OVER INDEPENDENT VARIABLES. */

            /* ---------------------------------------------------------------- */
            /* TEST DEPENDENCE ON INDEPENDENT VARIABLES OF INTEREST.            */
            /* (Rencher 1995, section 10.5.2)                                   */

                /* ------------------------------------------------------------ */
                /* GENERATE B_hat_r2 AND X_r2, MATRICES FOR THE REDUCED MODEL.  */

                for (i=1; i<=num_ind_var-num_ind_var_int; i++) {
                    for (col_index=1; col_index<=num_dep_var; col_index++)
                        B_hat_r2[i][col_index]=B_hat2[i][col_index];
                    for (row_index=1; row_index<=num_rows; row_index++)
                        X_r2[row_index][i]=SX[row_index][i];
                } /* END i LOOP, LOOP OVER ROWS IN REDUCED B MATRIX. */

                /* ------------------------------------------------------------ */
                /* CALCULATE EplusH_r2.                                         */

                for(i=1;i<=num_dep_var;i++)
                    for(j=1;j<=num_rows;j++) {
                        B_hat_r2tr_X_r2tr[i][j]=0.;
                        for (k=1; k<=num_ind_var-num_ind_var_int; k++)
                            B_hat_r2tr_X_r2tr[i][j]+=B_hat_r2[k][i]*X_r2[j][k];
                    }

                for(i=1;i<=num_dep_var;i++)
                    for(j=1;j<=num_dep_var;j++) {
                        B_hat_r2tr_X_r2tr_Y[i][j]=0.;
                        for (row_index=1; row_index<=num_rows; row_index++)
                            B_hat_r2tr_X_r2tr_Y[i][j]+=
                                        B_hat_r2tr_X_r2tr[i][row_index]*Y[row_index][j];
                    }

                for(i=1;i<=num_dep_var;i++)
                    for(j=1;j<=num_dep_var;j++)
                        EplusH_r2[i][j]=YtrY[i][j]-B_hat_r2tr_X_r2tr_Y[i][j];


                /* ------------------------------------------------------------ */
                /* CALCULATE det(EplusH_r2)                                     */

                EplusH_r2_okay=1;
                for(i=1;i<=num_dep_var;i++) {
                    bigEplusH_r2=0.;
                    for(j=1;j<=num_dep_var;j++)
                        bigEplusH_r2=(fabs(EplusH_r2[i][j])>bigEplusH_r2) ?
                                            fabs((double)EplusH_r2[i][j]) : bigEplusH_r2;
                    if (bigEplusH_r2==0.) {
                        E_EplusH_okay=0;
			/* JOE MAISOG : 08-22-2005 : Need to unset EplusH_r2_okay also! */
                        EplusH_r2_okay=0;
                        printf("\n  Singular reduced [E+H] matrix in variables of interest, at vox_index=%d\n",
                                                                vox_index);
                        printf("  Will set corresponding Wilks' lambda to zero.\n\n       ");
                    }
                }
                if (EplusH_r2_okay) {
                    Dludcmp(EplusH_r2,num_dep_var,index,&detEplusH_r2);
                    for (col_index=1; col_index<=num_dep_var; col_index++)
                        detEplusH_r2*=EplusH_r2[col_index][col_index];
                    if (detEplusH_r2<0.) {
/*                         printf("\n\n  Negative reduced det(E+H) generated at vox_index = %d\n",vox_index);
                        printf("  Will use absolute value.\n\n       ");
*/
                        detEplusH_r2=fabs(detEplusH_r2);
                    }
                }

                /* ------------------------------------------------------------ */
                /* CALCULATE WILKS' LAMBDA FOR REDUCED MODEL.                   */

                if ((EplusH_r2_okay!=0) && (detEplusH_r2!=0.))
/*                        && (!isinf(detE/detEplusH_r2)) && (!isnan(detE/detEplusH_r2))) */
                    Wilks_Lambda_r2[vox_index]=detE/detEplusH_r2;
                else
                    Wilks_Lambda_r2[vox_index]=0.;


                if (Wilks_Lambda_r2[vox_index]>1.)
                    Wilks_Lambda_r2[vox_index]=1.;


            /* ---------------------------------------------------------------- */
            /* ESTIMATION OF EFFECTIVE DEGREES OF FREEDOM.                      */
            /* Using method of Worsley KJ, Friston KJ, Analysis of fMRI Time-   */
            /* Series Revisited - Again, Neuroimage.                            */

            if ((effective_df[0]=='y') || (effective_df[0]=='Y')
                                                        || (yes_voxel_of_int==1)) {

                /* ------------------------------------------------------------ */
                /* DETERMINE MAXIMUM ESTIMATED TEMPORAL DISPERSION.             */
                /* Base calculation of effective df on the maximum estimated    */
                /* temporal dispersion, just in case there are different        */
                /* estimates of temporal dispersion for different independent   */
                /* variables.  This seems to be the conservative thing to do.   */

                max_SD=dispersion_map[0][vox_index];
                for (var_index=1; var_index<num_ind_var; var_index++)
                    max_SD = (max_SD>dispersion_map[var_index][vox_index])
                                ? max_SD : dispersion_map[var_index][vox_index];
                max_SD=max_SD/(float)tick;

                /* ------------------------------------------------------------ */
                /* GENERATE OVERALL HEMODYNAMIC RESPONSE FROM THIS AVERAGE.     */
                /* Don't need to include lag estimate because lag doesn't       */
                /* affect the estimate of effective degrees of freedom, only    */
                /* dispersion does.  Here, we'll set lag to zero.               */

                    for (row_index=0; row_index<num_rows; row_index++) {
                        kernel_index2=(float)row_index;
                        if ((float)row_index>((float)(num_rows)-1.)/2.)
                            kernel_index2=(float)(num_rows-row_index);
                        if (max_SD>0.1)
                            overall_HDR[row_index]=(1/(max_SD*sqrt_2pi))
                                                *
                                exp((double)-((kernel_index2*kernel_index2))
                                                /(2.*max_SD*max_SD));
                        else {
                            if (row_index==0)
                                overall_HDR[row_index]=1.;
                            else
                                overall_HDR[row_index]=0.;
                        }
                    }

                /* ------------------------------------------------------------ */
                /* MAKE MATRICES K AND V.                                       */
                /* Construct K, the smoothing matrix.  The first row is the     */
                /* hemodynamic response associated with the first time point.   */
                /* Each row thereafter is a lagged version of the hemodynamic   */
                /* response, corresponding to later time points.                */
                /* V is just K*K'.                                              */

                for (row_index=1;row_index<=num_rows;row_index++) {

                    for (row_index2=num_rows-row_index+2;row_index2<=num_rows;row_index2++)
                        K[row_index][row_index2-num_rows+row_index-1]=
                                                        overall_HDR[row_index2-1];

                    for (row_index2=1;row_index2<=num_rows-row_index+1;row_index2++)
                        K[row_index][row_index2+row_index-1]=
                                                        overall_HDR[row_index2-1];

                }

                for(i=1;i<=num_rows;i++)
                    for(j=1;j<=num_rows;j++) {
                        V[i][j]=0.;
                        for(k=1;k<=num_rows;k++)
                            V[i][j]+=K[i][k]*K[j][k];
                    }

                /* ------------------------------------------------------------ */
                /* CALCULATE G*inv(G'*G)*G', G_inv_G_tr_G_G_tr.                 */
                /* Actually, matrix G corresponds to matrix SX, the smoothed    */
                /* and lagged version of X; and inv_Xtr_X_Xtr, which            */
                /* corresponds to inv(G'*G)*G', has already been calculated.    */
                /* So, to calculate G*inv(G'*G)*G', all we need to do is        */
                /* multiply SX by inv_Xtr_X_Xtr.                                */

                for(i=1;i<=num_rows;i++)
                    for(j=1;j<=num_rows;j++) {
                        G_inv_G_tr_G_G_tr[i][j]=0.;
                        for(var_index=1; var_index<=num_ind_var; var_index++)
                            G_inv_G_tr_G_G_tr[i][j]+=
                                        SX[i][var_index]*inv_Xtr_X_Xtr[var_index][j];
                    }


                /* ------------------------------------------------------------ */
                /* MAKE MATRICES R, RV, AND RVRV.                               */
                /* R = residual-forming matrix = I-G*inv(G'*G)*G'               */

                for(i=1;i<=num_rows;i++)
                    for(j=1;j<=num_rows;j++) {
                        R[i][j]=-G_inv_G_tr_G_G_tr[i][j];
                        if (i==j) R[i][j]+=1.;
                    }

                for(i=1;i<=num_rows;i++)
                    for(j=1;j<=num_rows;j++) {
                        RV[i][j]=0.;
                        for(row_index=1; row_index<=num_rows; row_index++)
                            RV[i][j]+=R[i][row_index]*V[row_index][j];
                    }

                for(i=1;i<=num_rows;i++)
                    for(j=1;j<=num_rows;j++) {
                        RVRV[i][j]=0.;
                        for(row_index=1; row_index<=num_rows; row_index++)
                            RVRV[i][j]+=RV[i][row_index]*RV[row_index][j];
                    }

                /* ------------------------------------------------------------ */
                /* CALCULATE trace(R*V) AND trace(R*V*R*V).                     */

                trace_RV=0.;
                trace_RVRV=0.;
                for(i=1;i<=num_rows;i++) {
                    trace_RV+=RV[i][i];
                    trace_RVRV+=RVRV[i][i];
                }

                /* ------------------------------------------------------------ */
                /* GENERATE ESTIMATE OF EFFECTIVE DEGREES OF FREEDOM.           */

                nu[vox_index]=trace_RV*trace_RV/trace_RVRV;

            } /* END IF CHECKING WHETHER ESTIMATE OF EFFECTIVE DF IS DESIRED. */


            /* ---------------------------------------------------------------- */
            /* WRITE OUT INTERMEDIATE MATRICES TO DISK FOR VOXELS OF INTEREST.  */

            for (argv_index=0; argv_index<num_vox_int; argv_index++) {

                if (vox_int[argv_index]==vox_index) {
                    z=vox_index/(xdim*ydim);
                    y=(vox_index-(z*xdim*ydim))/xdim;
                    x=((vox_index-(z*xdim*ydim))-y*xdim);
                    x++;
                    y++;
                    z++;

                    printf("\n\n  Writing matrices for voxel index %d to directory voxel_%d_%d_%d...\n\n",
                                                vox_index,x,y,z);
                    sprintf(command,"mkdir voxel_%d_%d_%d",x,y,z);
                    system(command);
                    
                    sprintf(outfilename,"voxel_%d_%d_%d/X.dat",x,y,z);
                    write_NRCdouble_matrix2textfile(outfilename,SX,num_rows,num_ind_var);

/*
                    sprintf(outfilename,"voxel_%d_%d_%d/K.dat",x,y,z);
                    write_NRCdouble_matrix2textfile(outfilename,K,num_rows,num_rows);
*/

/*                  COMMENTED OUT THE FOLLOWING OUTPUT MATRICES BECAUSE THEY CAN BE EXTREMELY BIG.
                    sprintf(outfilename,"voxel_%d_%d_%d/X_orig_os.dat",x,y,z);
                    write_NRCdouble_matrix2textfile(outfilename,X_os,zdimXnum_rows,num_ind_var);

                    sprintf(outfilename,"voxel_%d_%d_%d/X_os.dat",x,y,z);
                    write_NRCdouble_matrix2textfile(outfilename,SX_os,zdimXnum_rows,num_ind_var);

                    sprintf(outfilename,"voxel_%d_%d_%d/hdr.dat",x,y,z);
                    write_float_matrix2textfile(outfilename,HDresponse,num_ind_var,zdimXnum_rows);
*/

                    sprintf(outfilename,"voxel_%d_%d_%d/Xtr_X.dat",x,y,z);
                    write_NRCdouble_matrix2textfile(outfilename,XtrX,num_ind_var,num_ind_var);

                    sprintf(outfilename,"voxel_%d_%d_%d/inv_Xtr_X_Xtr.dat",x,y,z);
                    write_NRCdouble_matrix2textfile(outfilename,inv_Xtr_X_Xtr,num_ind_var,num_rows);

                    sprintf(outfilename,"voxel_%d_%d_%d/df.dat",x,y,z);
                    if ((fOUTFILE = fopen(outfilename,"w")) == NULL) {
                        printf("Dmultivar_multiple_regressionC: can't open %s.\n",outfilename);
                        exit(13);
                    }
                    fprintf(fOUTFILE,"%10f\n",nu[vox_index]);
                    fclose(fOUTFILE);

                    sprintf(outfilename,"voxel_%d_%d_%d/B_hat.dat",x,y,z);
                    if ((fOUTFILE = fopen(outfilename,"w")) == NULL) {
                        printf("Dmultivar_multiple_regressionC: can't open %s.\n",outfilename);
                        exit(14);
                    }
                    for(i=0;i<num_ind_var;i++) {
                        for(j=0;j<num_dep_var;j++)
                            fprintf(fOUTFILE,"%10f ",B_hat[i][j][vox_index]);
                        fprintf(fOUTFILE,"\n");
                    }
                    fclose(fOUTFILE);
          
                    sprintf(outfilename,"voxel_%d_%d_%d/Y.dat",x,y,z);
                    write_NRCdouble_matrix2textfile(outfilename,Y,num_rows,num_dep_var);
                    
                    sprintf(outfilename,"voxel_%d_%d_%d/YtrY.dat",x,y,z);
                    write_NRCdouble_matrix2textfile(outfilename,YtrY,num_dep_var,num_dep_var);

                    sprintf(outfilename,"voxel_%d_%d_%d/B_hattr_Xtr_Y.dat",x,y,z);
                    write_NRCdouble_matrix2textfile(outfilename,B_hattr_Xtr_Y,num_dep_var,num_dep_var);

                    sprintf(outfilename,"voxel_%d_%d_%d/BrXrY.dat",x,y,z);
                    write_NRCdouble_matrix2textfile(outfilename,B_hat_r2tr_X_r2tr_Y,num_dep_var,num_dep_var);

                    sprintf(outfilename,"voxel_%d_%d_%d/E.dat",x,y,z);
                    write_NRCdouble_matrix2textfile(outfilename,E,num_dep_var,num_dep_var);

                    sprintf(outfilename,"voxel_%d_%d_%d/EplusH.dat",x,y,z);
                    write_NRCdouble_matrix2textfile(outfilename,EplusH,num_dep_var,num_dep_var);

                    sprintf(outfilename,"voxel_%d_%d_%d/y_bar_y_bartr.dat",x,y,z);
                    write_NRCdouble_matrix2textfile(outfilename,y_bar_y_bartr,num_dep_var,num_dep_var);

                    sprintf(outfilename,"voxel_%d_%d_%d/lag.dat",x,y,z);
                    if ((fOUTFILE = fopen(outfilename,"w")) == NULL) {
                        printf("Dmultivar_multiple_regressionC: can't open %s.\n",outfilename);
                        exit(15);
                    }
                    for(i=0;i<num_ind_var;i++)
                        fprintf(fOUTFILE,"%10f\n",lag_map[i][vox_index]*TR/(float)tick);
                    fclose(fOUTFILE);

                    sprintf(outfilename,"voxel_%d_%d_%d/dispersion.dat",x,y,z);
                    if ((fOUTFILE = fopen(outfilename,"w")) == NULL) {
                        printf("Dmultivar_multiple_regressionC: can't open %s.\n",outfilename);
                        exit(16);
                    }
                    for(i=0;i<num_ind_var;i++)
                        fprintf(fOUTFILE,"%10f\n",dispersion_map[i][vox_index]*TR/(float)tick);
                    fclose(fOUTFILE);

/*                printf("Continuing multivariate multiple regression over all voxels...\n\n       "); */
                    printf("Continuing multivariate multiple regression over all voxels...\n\n");
                    /* ********************
                    printf("  det(E)   = %f\n",detE);
                    printf("  det(E+H) = %f\n",detEplusH);
                    printf("  Wilks' Lambda = %f\n\n       ",detE/detEplusH);
                    ******************** */

                } /* END IF BLOCK, CONDITIONAL TESTING VARIABLES OF INTEREST. */

            } /* END LOOP OVER INDEPENDENT VARIABLES. */

            /* ---------------------------------------------------------------- */
            /* CALCULATE SUM OF ERRORS MAP.  THIS CAN BE USED TO ESTIMATE       */
            /* SPATIAL SMOOTHNESS.                                              */

            for (col_index=0; col_index<num_dep_var; col_index++)
                for (row_index=0; row_index<num_rows; row_index++) {

                     /* ------------------------------------------------------- */
                     /* SET float_buffer EQUAL TO CURRENT VOXEL VALUE.          */
                     float_buffer=(float)imgin[vox_index2][col_index][row_index];

                     /* ------------------------------------------------------- */
                     /* SUBTRACT OUT ALL ESTIMATED EXPERIMENTAL EFFECTS.        */
                     for (var_index=0; var_index<num_ind_var; var_index++)
                         float_buffer-=B_hat[var_index][col_index][vox_index]
                                                    *SX[row_index+1][var_index+1];

                     /* ------------------------------------------------------- */
                     /* STORE ADJUSTED VALUE IN SE[vox_index].                  */
                     SE[vox_index]+=(float)floor((double)float_buffer+0.5);
                }


            /* ---------------------------------------------------------------- */
            /* REMOVE EFFECTS OF NO INTEREST FROM fMRI DATA IF DESIRED.         */
            /* DESIRED.                                                         */

            if ((adjusted_fmri_data[0]=='y') || (adjusted_fmri_data[0]=='Y')) {
                for (col_index=0; col_index<num_dep_var; col_index++)
                    for (row_index=0; row_index<num_rows; row_index++) {

                         /* ------------------------------------------------------- */
                         /* SET float_buffer EQUAL TO CURRENT VOXEL VALUE.          */
                         float_buffer=(float)imgin[vox_index2][col_index][row_index];

                         /* ------------------------------------------------------- */
                         /* SUBTRACT OUT EFFECTS OF NO INTEREST.                    */
                         for (var_index=0; var_index<num_ind_var-num_ind_var_int; var_index++)
                             float_buffer-=B_hat[var_index][col_index][vox_index]
                                                        *SX[row_index+1][var_index+1];

                         /* ------------------------------------------------------- */
                         /* ADD IN grand_global_mean SO THERE'S NO NEGATIVE NUMBERS.*/
                         float_buffer+=mean_map[vox_index2];

                         /* ------------------------------------------------------- */
                         /* STORE ADJUSTED VALUE BACK INTO imgin[][][] ARRAY.       */
                         imgin[vox_index2][col_index][row_index]=
                                    (signed short int)floor((double)float_buffer+0.5);
                    }
            }
        } /* END IF STATEMENT CHECKING FOR VOXELS OF INTEREST. */
    } /* END LOOP OVER VOXELS. */
    printf("Finished looping over voxels!\n");
    printf("\n");

    /* ------------------------------------------------------------------------ */
    /* SET INDICATOR VOXELS, AND CHECK FOR BAD NUMBERS.                         */

printf("INDICATOR VOXELS AND BAD NUMBERS...\n");
    for (var_index=0; var_index<num_ind_var; var_index++) {
/*        for (col_index=0; col_index<num_dep_var; col_index++) {
            for (vox_index=0; vox_index<num_vox; vox_index++) {
                if (isinf(B_hat[var_index][col_index][vox_index]))
                    printf("B_hat[%d][%d][%d] is Inf!\n",var_index,col_index,vox_index);
                if (isnan(B_hat[var_index][col_index][vox_index]))
                    printf("B_hat[%d][%d][%d] is NaN!\n",var_index,col_index,vox_index);
            }
        }
*/
        Wilks_Lambda_r1[var_index][0]=1.;
    }
    Wilks_Lambda[0]=1.;
    Wilks_Lambda_r2[0]=1.;


/* ---------------------------------------------------------------------------- */
/*  WRITE ESTIMATED REGRESSION COEFFICIENT MAPS AND WILKS LAMBDA MAPS TO        */
/*  FLOATING POINT FILES.                                                       */

    if ((only_vox_of_int[0]=='n') || (only_vox_of_int[0]=='N')) {

        printf("*******************************************************************\n");
        printf("Writing maps for variables of interest...\n");
        for (var_index=num_ind_var-num_ind_var_int; var_index<num_ind_var; var_index++) {
            for (col_index=0; col_index<num_dep_var; col_index++) {
                sprintf(outfilename,"B_hat_%d_col%d.img",var_index+1,col_index+1);
                float_fwrite(outfilename,B_hat[var_index][col_index],num_vox);
            }
            sprintf(outfilename,"wilks_lambda_%d_p%dvH1vE%d.img",var_index+1,
                                    num_dep_var,num_rows-num_ind_var);
            float_fwrite(outfilename,Wilks_Lambda_r1[var_index],num_vox);
        }

        printf("\n*******************************************************************\n");
        printf("Writing Wilks' Lambda maps for variables of interest, and for all variables, and SSE map...\n");
        sprintf(outfilename,"wilks_lambda_omni_p%dvH%dvE%d.img",
                                    num_dep_var,num_ind_var-1,num_rows-num_ind_var);
        float_fwrite(outfilename,Wilks_Lambda,num_vox);

        sprintf(outfilename,"wilks_lambda_int_p%dvH%dvE%d.img",
                                    num_dep_var,num_ind_var_int,num_rows-num_ind_var);
        float_fwrite(outfilename,Wilks_Lambda_r2,num_vox);

        float_fwrite("SSE.img",SSE,num_vox);

        if ((output_no_int[0]=='y') || (output_no_int[0]=='Y')) {
            printf("\n*******************************************************************\n");
            printf("Writing maps for variables of NO interest...\n");
            for (var_index=0; var_index<num_ind_var-num_ind_var_int; var_index++) {
                for (col_index=0; col_index<num_dep_var; col_index++) {
                    sprintf(outfilename,"B_hat_%d_col%d.img",var_index+1,col_index+1);
                    float_fwrite(outfilename,B_hat[var_index][col_index],num_vox);
                }
                sprintf(outfilename,"wilks_lambda_%d_p%dvH1vE%d.img",var_index+1,
                                        num_dep_var,num_rows-num_ind_var);
                float_fwrite(outfilename,Wilks_Lambda_r1[var_index],num_vox);
            }
        }

        /* -------------------------------------------------------------------- */
        /*  WRITE OUT SUM OF ERRORS MAP.                                        */

            float_fwrite("SE.img",SE,num_vox);

    }

    if ((effective_df[0]=='y') || (effective_df[0]=='Y')) {
        printf("\n*******************************************************************\n");
        printf("Writing out map of effective degrees of freedom...\n");
        float_fwrite("nu.img",nu,num_vox);
    }


/* ---------------------------------------------------------------------------- */
/*  IF DESIRED, WRITE OUT fMRI DATA WITH EFFECTS OF NO INTEREST REMOVED.        */

    if ((adjusted_fmri_data[0]=='y') || (adjusted_fmri_data[0]=='Y')) {
        printf("\n*******************************************************************\n");
        printf("Writing out fMRI data with effects of no interest removed...\n");
        Wbuffer=signed_short_malloc(num_vox);

        for (vox_index2=0; vox_index2<num_vox; vox_index2++)
            Wbuffer[vox_index2]=0;

        for (col_index=0; col_index<num_dep_var; col_index++)
            for (row_index=0; row_index<num_rows; row_index++) {

                /* ------------------------------------------------------------ */
                /* TRANSFER ADJUSTED DATA FROM imgin TO Wbuffer.                */
                for (vox_index2=0; vox_index2<num_mask_vox; vox_index2++)
                    Wbuffer[mask_voxel[vox_index2]]=imgin[vox_index2][col_index][row_index];


                /* ------------------------------------------------------------ */
                /* DETERMINE NAME OF OUTPUT FILE.  The output file will have    */
                /* the same name as the unadjusted input file, with an "A"      */
                /* prepended to the name.  E.g., if the name of the unadjusted  */
                /* input file is fmri.img, the name of the adjusted output file */
                /* will be Afmri.img.                                           */

                length = strlen(img_name[col_index][row_index]);
                str_index=length;
                while ((str_index>0) && (img_name[col_index][row_index][str_index]!='/'))
                    str_index--;
                if (str_index>0) {
                    for (index2=0;index2<=str_index; index2++)
                        path[index2]=img_name[col_index][row_index][index2];
                    path[index2]='\0';
                    for (index2=0;index2<length-str_index; index2++)
                        file[index2]=img_name[col_index][row_index][str_index+index2+1];
/*                    sprintf(outfilename,"%sA%s",path,file); */
                    sprintf(outfilename,"A%s",path,file);
                }
                else
                    sprintf(outfilename,"A%s",img_name[col_index][row_index]);


                /* ------------------------------------------------------------ */
                /* THEN WRITE TO DISK.                                          */
                signed_short_fwrite(outfilename, Wbuffer, num_vox);
            }

    }

/* ---------------------------------------------------------------------------- */
/*  FREE MEMORY.                                                                */

    free(vox_int);
    free(slice_seq1);
    free(slice_seq2);
    free(img_name);
    free(HDresponse);
    free(lag_map);
    free(dispersion_map);
    free(last_lag);
    free(last_dispersion);
/*    free(temp);
	JOE MAISOG : 08-22-2005 : Commented out "free(temp)", since temp was already freed above.*/
    free(imgin);
    free(global_mean);
    free(Wilks_Lambda);
    free(Wilks_Lambda_r1);
    free(Wilks_Lambda_r2);
    free(SE);
    free(SSE);
    free(mask);

    free_dmatrix(X_os,1,zdimXnum_rows,1,num_ind_var);
    free_dmatrix(SX,1,num_rows,1,num_ind_var);
    free_dmatrix(SX_os,1,zdimXnum_rows,1,num_ind_var);
    free_dmatrix(X_r1,1,num_rows,1,num_ind_var-1);
    free_dmatrix(X_r2,1,num_rows,1,num_ind_var-1);
    free_dmatrix(Y,1,num_rows,1,num_dep_var);
    free_dmatrix(XtrX,1,num_ind_var,1,num_ind_var);
    free_dmatrix(inv_Xtr_X_Xtr,1,num_ind_var,1,num_rows);
    free_dmatrix(B_hat_r1,1,num_ind_var-1,1,num_dep_var);
    free_dmatrix(B_hat_r2,1,num_ind_var-1,1,num_dep_var);
    free_dmatrix(B_hat2,1,num_ind_var,1,num_dep_var);
    free_dvector(B_hat_tr_bar_tr,1,num_ind_var);
    free_dmatrix(YtrY,1,num_dep_var,1,num_dep_var);
    free_dvector(y_bar,1,num_dep_var);
    free_dmatrix(y_bar_y_bartr,1,num_dep_var,1,num_dep_var);
    free_dmatrix(B_hattr_Xtr,1,num_dep_var,1,num_rows);
    free_dmatrix(B_hat_r1tr_X_r1tr,1,num_dep_var,1,num_rows);
    free_dmatrix(B_hat_r1tr_X_r1tr_Y,1,num_dep_var,1,num_dep_var);
    free_dmatrix(B_hat_r2tr_X_r2tr,1,num_dep_var,1,num_rows);
    free_dmatrix(B_hat_r2tr_X_r2tr_Y,1,num_dep_var,1,num_dep_var);
    free_dmatrix(B_hattr_Xtr_Y,1,num_dep_var,1,num_dep_var);
    free_dmatrix(E,1,num_dep_var,1,num_dep_var);
    free_dmatrix(EplusH,1,num_dep_var,1,num_dep_var);
    free_dmatrix(EplusH_r1,1,num_dep_var,1,num_dep_var);
    free_dmatrix(EplusH_r2,1,num_dep_var,1,num_dep_var);

    free_dmatrix(K,1,num_rows,1,num_rows);
    free_dmatrix(V,1,num_rows,1,num_rows);
    free_dmatrix(G_inv_G_tr_G_G_tr,1,num_rows,1,num_rows);
    free_dmatrix(R,1,num_rows,1,num_rows);
    free_dmatrix(RV,1,num_rows,1,num_rows);
    free_dmatrix(RVRV,1,num_rows,1,num_rows);

    free_dvector(col,1,num_ind_var);
    free(B_hat);
    free(normalize_flag);
    free(nu);
/*    free(img_name);
	JOE MAISOG : 08-26-2005 : Commented out "free(img_name)", since img_name was already freed above.*/
    free(mean_map);
    free(Wbuffer);
}
