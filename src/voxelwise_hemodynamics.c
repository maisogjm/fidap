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
/* Program voxelwise_hemodynamics.c                                             */
/* ANSI C Code by Jose' Ma. Maisog, M.D.                                        */
/* Section on Functional Brain Imaging,                                         */
/* Laboratory of Psychology and Psychopathology                                 */
/* National Institute of Mental Health                                          */
/*                                                                              */
/* Program reads in fMRI data for several runs and independent variables, and   */
/* estimates lag and dispersion for each voxel and each independent variable.   */
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

main(int argc, char *argv[])
{
                            /* INDICES AND CONSTANTS.                           */
    int num_vox;            /* Number of voxels in an image.                    */
    int vox_index;          /* Index into imgin[][][*], extracting a specific   */
                            /* voxel.                                           */
    int xdim,ydim,zdim;     /* Image dimensions.                                */
    int x,y,z;              /* XYZ coordinates.                                 */
    int num_runs;           /* Number of fMR runs.                              */
    int run_index;          /* Index into run number.                           */
    int num_dat_pts;        /* The number of time points in a run.              */
    int num_runsXnum_dat_pts;
                            /* num_runs times num_dat_pts                       */
    int dat_index;          /* Index into time points.                          */
    int tot_num_scans;      /* Total number of scans.  Should be equal to       */
                            /* num_runs * num_dat_pts                           */
    int *slice_seq1;        /* Slice acquisition sequence.                      */
    int *slice_seq2;        /* Slice acquisition sequence.                      */
    float TR;               /* TR in seconds.                                   */
    char absolute_time[2];  /* Flag indicating whether absolute times are       */
                            /* desired or not.                                  */
    char normalize[2];      /* Flag indicating whether normalization is desired */
                            /* not.                                             */

                            /* INPUT IMAGES.                                    */
    FILE *ftxt;             /* Pointer for input text file containing filenames */
    char ***img_name;       /* Pointer to pointers to input image file names.   */
                            /* Indexing:                                        */
                            /*      img_name[run][timepoint][ASCII character]   */
    char dummy[80];         /* Dummy pointer for counting number of entries     */
                            /* in a text file.                                  */
    int check_fscanf;       /* Flag to check error status of fscanf operation.  */
    struct analyze_struct
                img_info;   /* Header information variable.                     */
    unsigned                /* Pointer to temporary short int array, for        */
        short int *temp;    /* reading in 16-bit fMRI data.                     */
                            /* byte data.                                       */
    unsigned short ***imgin;/* Pointer to input image data.                     */
                            /* Indexing: imgin[voxel][run][scan]                */
    float *mean_fmr;        /* Mean of the fMR time series across runs.         */
    float **global_mean;    /* Global mean intensity of a scan.                 */
    int num_files_read_in;  /* Number of files read in so far.                  */

                            /* MASK                                             */
    unsigned char *mask;    /* Pointer to mask image.                           */
    int num_mask_vox;       /* Number of non-zero voxels in *mask.              */
    int count;              /* Number of voxels processed thus far.             */

                            /* INDEPENDENT VARIABLES (INPUT WAVEFORMS).         */
    FILE *f_indvar;         /* Pointer to file containing independent variable. */
    int num_ind_var;        /* Number of independent variables.                 */
    int var_index;          /* Index over independent variables.                */
    int argv_index;         /* Index over command line arguments.               */
    float var_mean,         /* Mean of an independent variable, to be set to 0. */
          var_var;          /* Variance of an ind. variable, to be set to 1.    */
    float float_buffer;     /* Buffer for reading in floating point numbers.    */
    double **X;             /* Matrix of independent variables.                 */

                            /* EXTERNAL VARIABLES, SMOOTHED INPUT SQUAREWAVE    */
                            /* AND HEMODYNAMIC RESPONSE FUNCTION.               */
    float pi;               /* Pi.                                              */
    int dat_index2;         /* Index used in calculating lookup table of        */
                            /* exponentials.                                    */
    float RFTfmr_m,IFTfmr_m;/* Real and imaginary parts of mth frequency        */
                            /* component of an fMR signal.                      */
    float RFTfmr_n,IFTfmr_n;/* Real and imaginary parts of nth frequency        */
                            /* component of an fMR signal.                      */
    float max,              /* Maximum of all magnitudes thus far considered.   */
          mag;              /* Magnitude of a frequency component.              */
    int *n,*m;              /* Indices for frequencies of interest in the FT of */
                            /* the input squarewave and the idealized response. */
                            /* The indices will be determined for each input    */
                            /* waveform.                                        */
    float **Rexp,**Iexp;    /* Real and imaginary parts of complex exponentials */
                            /* used in calculating FT's.                        */
    float *Rfft,*Ifft;      /* Real and imaginary parts of FFT of a signal      */
    float *RFTinw_m,        /* Real and imaginary parts of mth component of FT  */
                *IFTinw_m;  /* of square wave, for each run.                    */
    float *RFTinw_n,        /* Real and imaginary parts of nth component of FT  */
                *IFTinw_n;  /* of square wave, for each run.                    */
    float *phase_inw;       /* Phase of input stimulus waveform.                */
    float phase_fmr;        /* Phase of the fMR signal.                         */
    float phase_diff;       /* Difference in phase between phase_fmr and        */
                            /* phase_inw.                                       */
    float Ym,Yn;            /* Mags of mth and nth frequency components of fMR  */
                            /* signal.                                          */
    float Xm,Xn;            /* Mags of mth and nth frequency components of      */
                            /* input square wave.                               */
    float                   /* Denominator of estimator for standard deviation  */
       *denominator_factor; /* of equivalent gaussian smoothing filter.         */
    float mean;             /* Mean of Gaussian lagging/smoothing filter.       */
                            /* Will function as the temporal lag.               */
    float max_lag;          /* Maximum lag allowed, in terms of TR's.           */
    float min_lag;          /* Minimum lag allowed, in terms of TR's.           */
    float SD;               /* SD of Gaussian lagging/smoothing filter.         */
                            /* Functions as the temporal dispersion.            */
    float max_mag;          /* Max mag of mth frequency component of fMRI       */
                            /* signal, across independent variables.            */

                            /* OUTPUT FILES AND DIRECTORIES.                    */
    float **lag_map,        /* Lag and Dispersion maps.                         */
        **dispersion_map;
    float *best_lag,        /* Best estimated lag and dispersion.               */
          *best_disp;
    char outfilename[80];   /* Array to contain name of output text file.       */
    char command[256];      /* UNIX command to be sent to system call.          */
    int num_wrote_out;      /* Number of voxels written to disk.                */


/* ---------------------------------------------------------------------------- */
/*  CHECK TO SEE IF CORRECT NUMBER OF ARGUMENTS WERE PASSED.                    */

    if (argc<14) {
        printf("voxelwise_hemodynamics: need input header file, input textfile of\n");
        printf("file names, number of runs, number of scans per run, name of 8-bit\n");
        printf("mask, maximum lag to be allowed in units of seconds, minimum lag to\n");
        printf("be allowed in units of seconds, total number of independent variables\n");
        printf("including any constant or ramp function, text file containing slice\n");
        printf("acquisition sequence, TR in seconds, and a list of input text files\n");
        printf("containing independent variables, each followed by the name of the\n");
        printf("corresponding output floating point lag and dispersion map files.\n\n");

        printf("All arguments are of course to be separated by spaces.\n\n");

        printf("The input textfile of file names should contain the names of the\n");
        printf("fMRI image volume files, one name per line.  If you have more\n");
        printf("than one fMRI run, just list all the files for the first run,\n");
        printf("then list all the files for the second run, etc., all in the\n");
        printf("same text file, one name per line.\n\n");

        printf("The input text files containing independent variables should\n");
        printf("have one number per line.  There should be as many numbers in each\n");
        printf("of these text files as there are scans per run.  There can be any\n");
        printf("number of independent variables, but there should be at least one.\n\n");

        printf("In addition to the lag and dispersion maps for each independent\n");
        printf("variable, one \"best\" lag map and one \"best\" dispersion map\n");
        printf("will be generated.  These will be the best estimates, based on\n");
        printf("the power of the frequency components found in the fMRI signal.\n");
        exit(1);
    }


/* ---------------------------------------------------------------------------- */
/*  READ IN num_runs AND num_dat_pts OFF OF COMMAND LINE.                       */

    num_runs             = atoi(argv[3]);
    num_dat_pts          = atoi(argv[4]);
    num_runsXnum_dat_pts = num_runs*num_dat_pts;
    num_ind_var          = atoi(argv[8]);

    printf("There are %d runs.\n", num_runs);
    printf("There are %d input scans per run.\n", num_dat_pts);
    printf("This makes a total of %d input scans.\n",num_runsXnum_dat_pts);
    printf("Will read in input scan filenames from %s.\n",argv[2]);
    printf("\nThere are to be %d independent variables in the regression\n",num_ind_var);
    printf("(including any constant or ramp function).\n\n");


/* ---------------------------------------------------------------------------- */
/*  READ HEADER FILE; ASSUME IMAGES ARE SAME SIZE, SO USE SAME HEADER FILE      */
/*  FOR ALL IMAGES.                                                             */

    printf("Reading header info from %s...\n",argv[1]);
    img_info = readhdr(argv+1);

    xdim = img_info.dime.dim[1];
    ydim = img_info.dime.dim[2];
    zdim = img_info.dime.dim[3];
    num_vox = xdim*ydim*zdim;
    printf("Images have dimensions: [%d %d %d]\n\n",xdim,ydim,zdim);


/* ---------------------------------------------------------------------------- */
/* READ IN SLICE ACQUISITION SEQUENCE.                                          */
/* Convert to numbers from 0 to zdim-1.  If slice k was the nth slice acquired, */
/* set slice_seq2[k-1] equal to (n-1).                                          */

    printf("*********************************************************************\n");
    printf("\nNow you must indicate whether you want ABSOLUTE or RELATIVE times.\n");
    printf("In our convention, ABSOLUTE time is as follows.  Suppose two voxels\n");
    printf("are part of the same functional region and behave exactly the same,\n");
    printf("but they are in different slices.  Under \"ABSOLUTE\" time, they will\n");
    printf("have the same lag value in the lag map, but under \"RELATIVE\" time,\n");
    printf("they will have different times, because of the time difference in slice\n");
    printf("acquisition.  (The voxel in the slice acquired at a later time will\n");
    printf("have a larger lag.)\n\n");
    printf("\nDo you want ABSOLUTE times <y/n>? ");
    scanf("%2c",absolute_time);
    while((absolute_time[0]!='y') && (absolute_time[0]!='n')
                        && (absolute_time[0]!='Y') && (absolute_time[0]!='N')) {
        printf("Type 'y' or 'n': ");
        scanf("%2c",absolute_time);
    }

    if ((absolute_time[0]=='y') || (absolute_time[0]=='Y')) {
        printf("\nReading in slice acquisition sequence from %s...\n",argv[9]);
        if ((ftxt = fopen(argv[9],"r")) == NULL) {
            printf("voxelwise_hemodynamics: can't open %s\n",argv[9]);
            exit(2);
        }
        slice_seq1=int_malloc(zdim);
        slice_seq2=int_malloc(zdim);

        printf("Slices were acquired in this sequence: ");
        for (z=0; z<zdim; z++) {
            check_fscanf=fscanf(ftxt,"%d",slice_seq1+z);
            if (check_fscanf==EOF) {
                printf("voxelwise_hemodynamics: error reading slice sequence from %s.\n",
                                                                            argv[9]);
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
                printf("voxelwise_hemodynamics: weird values for slice sequence.\n");
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
    }

    printf("\nDo you want scans to be ratio normalized <y/n>? ");
    scanf("%2c",normalize);
    while((normalize[0]!='y') && (normalize[0]!='n')
                        && (normalize[0]!='Y') && (normalize[0]!='N')) {
        printf("Type 'y' or 'n': ");
        scanf("%2c",normalize);
    }

/* ---------------------------------------------------------------------------- */
/*  COUNT NUMBER OF INPUT IMAGE FILENAMES.                                      */

    if ((ftxt = fopen(argv[2],"r")) == NULL) {
        printf("voxelwise_hemodynamics: can't open %s\n",argv[2]);
        exit(5);
    }
    tot_num_scans = 0;
    while(fscanf(ftxt,"%s",dummy)!=EOF)
        tot_num_scans++;

    if (tot_num_scans != num_runsXnum_dat_pts) {
        printf("voxelwise_hemodynamics: number of entries in %s not equal to number\n",argv[2]);
        printf("of runs times number of scans per run (%d versus %d).\n",
                                        num_runsXnum_dat_pts,tot_num_scans);
        exit(6);
    }
    rewind(ftxt);

/* ---------------------------------------------------------------------------- */
/*  ALLOCATE MEMORY FOR img_name; THEN READ IN NAMES OF INPUT IMAGE FILES.      */

    printf("Total Number of Scans = %d\n",tot_num_scans);

    printf("\nAllocating memory for filenames...\n");
    img_name = char_malloc3(num_runs,num_dat_pts,80);
    for (run_index=0; run_index<num_runs; run_index++)
        for (dat_index=0; dat_index<num_dat_pts; dat_index++) {
            check_fscanf=fscanf(ftxt,"%s",img_name[run_index][dat_index]);
            if (check_fscanf==EOF) {
                printf("voxelwise_hemodynamics: error reading run %d timepoint %d filename from %s.\n",
                                                run_index+1,dat_index+1,argv[2]);
                exit(7);
            }
        }
    fclose(ftxt);


/* ---------------------------------------------------------------------------- */
/* ALLOCATE MEMORY FOR ARRAYS.                                                  */

    X              = dmatrix(1,num_dat_pts,1,num_ind_var);
    lag_map        = float_malloc2(num_ind_var,num_vox);
    dispersion_map = float_malloc2(num_ind_var,num_vox);
    best_lag       = float_malloc(num_vox);
    best_disp      = float_malloc(num_vox);

    for (var_index=0; var_index<num_ind_var; var_index++)
        for (vox_index=0; vox_index<num_vox; vox_index++) {
            lag_map[var_index][vox_index]=0.;
            dispersion_map[var_index][vox_index]=0.;
        }

    for (vox_index=0; vox_index<num_vox; vox_index++) {
        best_lag[vox_index]=0;
        best_disp[vox_index]=0;
    }

/* ---------------------------------------------------------------------------- */
/* READ IN 8-BIT MASK.                                                          */

    mask = unsigned_char_malloc(num_vox);
    read_unsigned_char_data(argv[5],mask,num_vox);

    num_mask_vox=0;
    for (vox_index=0; vox_index<num_vox; vox_index++)
        if (mask[vox_index]!=0)
            num_mask_vox++;
    printf("\nNumber of nonzero voxels in mask %s = %d.\n\n",argv[5],num_mask_vox);


/* ---------------------------------------------------------------------------- */
/* READ IN INDEPENDENT VARIABLES FROM TEXT FILES.                               */

    printf("\nWill read in independent variables from the files:\n");
    for (var_index=0; var_index<num_ind_var; var_index++) {
        printf("  #%d : %s\n",var_index+1,argv[11+3*var_index]);
        printf("       lag output will be written to %s\n",argv[12+3*var_index]);
        printf("       dispersion output will be written to %s\n\n",argv[13+3*var_index]);
    }

    printf("\n");

    printf("Constructing matrix X of regression variables...\n");
    for (var_index=0; var_index<num_ind_var; var_index++) {
        if ((f_indvar = fopen(argv[11+3*var_index],"r")) == NULL) {
            printf("voxelwise_hemodynamics: can't open %s\n",argv[11+3*var_index]);
            exit(8);
        }
        printf("  Reading in data from %s (independent variable #%d)...\n",
                                                argv[11+3*var_index],var_index+1);
        for (dat_index=0; dat_index<num_dat_pts; dat_index++) {
            if ((check_fscanf=fscanf(f_indvar,"%f",&float_buffer))==EOF) {
                printf("voxelwise_hemodynamics: error reading %s.\n",
                                                                argv[11+3*var_index]);
                exit(9);
            }
            X[dat_index+1][var_index+1]=float_buffer;
        }
        fclose(f_indvar);
    }
    printf("\n");


/* ------------------------------------------------------------------------ */
/* SET UP ARRAYS OF REAL AND IMAGINARY PARTS OF COMPLEX EXPONENTIALS.       */

    pi = 2.*acos(0.);
    Rexp = float_malloc2(num_dat_pts,num_dat_pts);
    Iexp = float_malloc2(num_dat_pts,num_dat_pts);

    for (dat_index=0; dat_index<num_dat_pts; dat_index++) {
        for (dat_index2=0; dat_index2<num_dat_pts; dat_index2++) {
            Rexp[dat_index][dat_index2]=(float)cos(
                ((2.*pi*((double)dat_index)*(double)dat_index2/(double)num_dat_pts)));
            Iexp[dat_index][dat_index2]=(float)sin(
                ((2.*pi*((double)dat_index)*(double)dat_index2/(double)num_dat_pts)));
        }
    }
/*
    printf("Complex exponentials for 12th freq component:\n");
    for (dat_index=0; dat_index<num_dat_pts; dat_index++)
        printf("%d %f %f\n",dat_index,Rexp[12][dat_index],Iexp[12][dat_index]);

    printf("Complex exponentials for 24th freq component:\n");
    for (dat_index=0; dat_index<num_dat_pts; dat_index++)
        printf("%d %f %f\n",dat_index,Rexp[24][dat_index],Iexp[24][dat_index]);
*/


/* ------------------------------------------------------------------------ */
/* Determine m and n, the indices of the frequencies of the input waveforms */
/* with the most power, not including the 0th frequency component (DC bias).*/

    m = int_malloc(num_ind_var);
    n = int_malloc(num_ind_var); 

    Rfft = float_malloc(num_dat_pts);
    Ifft = float_malloc(num_dat_pts);

    RFTinw_m = float_malloc(num_ind_var);
    IFTinw_m = float_malloc(num_ind_var);
    RFTinw_n = float_malloc(num_ind_var);
    IFTinw_n = float_malloc(num_ind_var);

    denominator_factor = float_malloc(num_ind_var);
    phase_inw = float_malloc(num_ind_var);

    TR = atof(argv[10]);
    printf("TR = %g seconds.\n",TR);
    max_lag = atof(argv[6]);
    min_lag = atof(argv[7]);
    printf("Maximum lag allowed will be %g seconds.\n\n",max_lag);
    printf("Minimum lag allowed will be %g seconds.\n\n",min_lag);
    max_lag/=TR; /* Convert from seconds to TR's */
    min_lag/=TR; /* Convert from seconds to TR's */

    printf("Calculating phases of independent variables...\n");
    for (var_index=0; var_index<num_ind_var; var_index++) {


            /* ---------------------------------------------------------------- */
            /* FFT THE INPUT WAVEFORM.                                          */
            for (dat_index=0; dat_index<num_dat_pts; dat_index++) {
                Rfft[dat_index]=0.;
                Ifft[dat_index]=0.;
                for (dat_index2=0; dat_index2<num_dat_pts; dat_index2++) {
                    Rfft[dat_index]+=X[dat_index2+1][var_index+1]*Rexp[dat_index][dat_index2];
                    Ifft[dat_index]+=X[dat_index2+1][var_index+1]*Iexp[dat_index][dat_index2];
                }
            }

            /* ---------------------------------------------------------------- */
            /* FIND mth FREQUENCY COMPONENT.  THIS IS THE FREQUENCY WITH THE    */
            /* GREATEST POWER, NOT INCLUDING 0 Hz.                              */
            max=-1.;
            m[var_index]=0;
            for (dat_index=1; dat_index<num_dat_pts/2; dat_index++) {
                mag=(float)fabs((double)(Rfft[dat_index]*Rfft[dat_index]+
                                            Ifft[dat_index]*Ifft[dat_index]));

                if (mag>max) {
                    max=mag;
                    m[var_index]=dat_index;
                }
            }

            /* ---------------------------------------------------------------- */
            /* FIND nth FREQUENCY COMPONENT.  THIS WILL BE THE FREQUENCY COMPO- */
            /* NENT WITH GREATEST POWER AFTER THE mth FREQUENCY COMPONENT, WITH */
            /* INDEX > m.                                                       */
            max=-1.;
            n[var_index]=0;
            for (dat_index=m[var_index]+1; dat_index<num_dat_pts/2; dat_index++) {
                mag=(float)fabs((double)(Rfft[dat_index]*Rfft[dat_index]+
                                                Ifft[dat_index]*Ifft[dat_index]));

                if (mag>max) {
                    max=mag;
                    n[var_index]=dat_index;
                }
            }

            printf("  The frequency components which will be used to estimate the\n");
            printf("  response for independent variable #%d are #%d and #%d.\n",
                                        var_index+1,m[var_index],n[var_index]);
            denominator_factor[var_index]=
                        pi*(n[var_index]*n[var_index]-m[var_index]*m[var_index]);


        /* ------------------------------------------------------------------------ */
        /* CALCULATE FT OF INDEPENDENT VARIABLES AT SPECIFIC FREQUENCIES m AND n.   */
        /* ALSO, CALCULATE PHASE OF mth FREQUENCY COMPONENT, FOR FIRST RUN.         */
        /* FOR THIS CALCULATION, IT IS ASSUMED THAT THE INPUT STIMULUS WAVEFORM IS  */
        /* A SQUAREWAVE AND CONSTANT ACROSS RUNS.                                   */

            RFTinw_m[var_index]=0.;
            IFTinw_m[var_index]=0.;
            RFTinw_n[var_index]=0.;
            IFTinw_n[var_index]=0.;
            for (dat_index=0; dat_index<num_dat_pts; dat_index++) {
                RFTinw_m[var_index]+=X[dat_index+1][var_index+1]
                                            *Rexp[m[var_index]][dat_index];
                IFTinw_m[var_index]+=X[dat_index+1][var_index+1]
                                            *Iexp[m[var_index]][dat_index];
                RFTinw_n[var_index]+=X[dat_index+1][var_index+1]
                                            *Rexp[n[var_index]][dat_index];
                IFTinw_n[var_index]+=X[dat_index+1][var_index+1]
                                            *Iexp[n[var_index]][dat_index];
            }
            Xm=(float)sqrt((double)(RFTinw_m[var_index]*RFTinw_m[var_index]
                                    +IFTinw_m[var_index]*IFTinw_m[var_index]));
            Xn=(float)sqrt((double)(RFTinw_n[var_index]*RFTinw_n[var_index]
                                    +IFTinw_n[var_index]*IFTinw_n[var_index]));

            /* ---------------------------------------------------------------- */
            /* CALCULATE PHASE OF INPUT STIMULUS WAVEFORM.  Using algorithm     */
	    /* kindly provided by Jagath Rajapakse, to return the correct angle */
	    /* in the complex plane.                                            */
            if (IFTinw_m[var_index] == 0.) {
                if (RFTinw_m[var_index] >= 0.)
                    phase_inw[var_index]=0.;
                else
                    phase_inw[var_index]=pi;
            }
            else if (RFTinw_m[var_index] == 0) {
                if (IFTinw_m[var_index] < 0)
                    phase_inw[var_index]=1.5*pi;
                if (IFTinw_m[var_index] > 0)
                    phase_inw[var_index]=0.5*pi;
            }
            else if (RFTinw_m[var_index] < 0 & IFTinw_m[var_index] > 0)
                phase_inw[var_index]=pi - atan(-IFTinw_m[var_index]/RFTinw_m[var_index]);
            else if (RFTinw_m[var_index] <  0 & IFTinw_m[var_index] < 0)
                phase_inw[var_index]=pi + atan(IFTinw_m[var_index]/RFTinw_m[var_index]);
            else if (RFTinw_m[var_index] > 0 & IFTinw_m[var_index] < 0)
                phase_inw[var_index]=2.0*pi - atan(-IFTinw_m[var_index]/RFTinw_m[var_index]);
            else
                phase_inw[var_index]=atan(IFTinw_m[var_index]/RFTinw_m[var_index]);

            printf("  The phase at frequency component %d is %g\n\n",        
                                    m[var_index],phase_inw[var_index]);
    }


/* ---------------------------------------------------------------------------- */
/*  READ IN IMAGES.                                                             */

    temp     = unsigned_short_malloc(num_vox);
    printf("Allocating memory for fMRI data...\n");
    imgin    = unsigned_short_malloc3(num_vox,num_runs,num_dat_pts);
    mean_fmr = float_malloc(num_dat_pts);

    printf("\nReading in fMRI scans...\n       ");
    num_files_read_in=0;
    for (run_index=0; run_index<num_runs; run_index++)
        for (dat_index=0; dat_index<num_dat_pts; dat_index++) {
            num_files_read_in++;
            printf("\b\b\b\b\b\b\b%5.1f %%",(float)num_files_read_in*100./(float)tot_num_scans);
            read_unsigned_short_int_data(img_name[run_index][dat_index],temp,num_vox);

            for (vox_index=0; vox_index<num_vox; vox_index++)
                imgin[vox_index][run_index][dat_index]=(float)temp[vox_index];
        }
    printf("\n\n");
    free(temp);
    free(img_name);


/* ---------------------------------------------------------------------------- */
/* NORMALIZE GLOBAL INTENSITY.                                                  */


    if ((normalize[0]=='y') || (normalize[0]=='Y')) {
        printf("Calculating scan global means...\n       ");
        global_mean=float_malloc2(num_runs,num_dat_pts);
        count=0;
        for (run_index=0; run_index<num_runs; run_index++)
            for (dat_index=0; dat_index<num_dat_pts; dat_index++) {
                count++;
                    printf("\b\b\b\b\b\b\b%5.1f %%",(float)count*100./(float)(num_runs*num_dat_pts));
                global_mean[run_index][dat_index]=0.;
                for (vox_index=0; vox_index<num_vox; vox_index++)
                    if (mask[vox_index]!=0)
                        global_mean[run_index][dat_index]+=imgin[vox_index][run_index][dat_index];
                global_mean[run_index][dat_index]/=(float)num_mask_vox;
            }

        /* ********************** COMMENTED THIS SECTION OUT ****************************
        This will have no effect on lag estimates, and only adds to computations.
        printf("\n\nNow normalizing global means to grand global mean (=%.3f)...\n       ",
                                                                    grand_global_mean);
        count=0;
        for (run_index=0; run_index<num_runs; run_index++)
            for (dat_index=0; dat_index<num_dat_pts; dat_index++) {
                count++;
                printf("\b\b\b\b\b\b\b%5.1f %%",(float)count*100./(float)(num_runs*num_dat_pts));
                for (vox_index=0; vox_index<num_vox; vox_index++)
                     imgin[vox_index][run_index][dat_index]*=grand_global_mean;
            }
        printf("\n\n");
        ******************************************************************************* */
    }


/* ---------------------------------------------------------------------------- */
/* LOOP OVER VOXELS, PERFORM VOXEL-WISE ESTIMATE OF LAG AND DISPERSION.         */

    printf("\n\nNow estimating lag and dispersion for each voxel and independent variable...\n");
    for (vox_index=0; vox_index<num_vox; vox_index++) {

        /* printf("\b\b\b\b\b\b\b%5.1f %%",(float)vox_index*100./(float)num_vox); */

        /* -------------------------------------------------------------------- */
        /* IF CURRENT VOXEL IS WITHIN MASK, ESTIMATE HEMODYANMIC RESPONSE.      */

        if (mask[vox_index]!=0) {

            /* ---------------------------------------------------- */
            /* CUMULUATE LAG AND DISPERSION ESTIMATES OVER ALL RUNS. */

            for (dat_index=0; dat_index<num_dat_pts; dat_index++)
                mean_fmr[dat_index]=0.;

            for (run_index=0; run_index<num_runs; run_index++)
                for (dat_index=0; dat_index<num_dat_pts; dat_index++)
                    mean_fmr[dat_index]+=(float)imgin[vox_index][run_index][dat_index]
                                        /global_mean[run_index][dat_index];

            /* -------------------------------------------------------- */
            /* LOOP OVER ALL INDEPENDENT VARIABLES.
            printf("Looping over independent variables to estimate lag and dispersion...\n"); */

            for (var_index=0; var_index<num_ind_var; var_index++) {

                /* ---------------------------------------------------- */
                /* CALCULATE mth AND nth FREQUENCY COMPONENTS.          */
                /* SUM INTO RFTfmr_m, IFTfmr_m, RFTfmr_n, IFTfmr_n.     */

                RFTfmr_m=0.;
                IFTfmr_m=0.;
                RFTfmr_n=0.;
                IFTfmr_n=0.;
                for (dat_index=0; dat_index<num_dat_pts; dat_index++) {
                    RFTfmr_m += mean_fmr[dat_index]*Rexp[m[var_index]][dat_index];
                    IFTfmr_m += mean_fmr[dat_index]*Iexp[m[var_index]][dat_index];
                    RFTfmr_n += mean_fmr[dat_index]*Rexp[n[var_index]][dat_index];
                    IFTfmr_n += mean_fmr[dat_index]*Iexp[n[var_index]][dat_index];
                }

                /* ------------------------------------------------------------ */
                /* CALCULATE LAG FROM PHASE LAG AT mth FREQUENCY COMPONENT.     */

                phase_fmr=(float)atan((double)(IFTfmr_m/RFTfmr_m));
                if ((RFTfmr_m<0.) && (IFTfmr_m<0.))
                    phase_fmr=phase_fmr-pi;
                else if ((RFTfmr_m<0.) && (IFTfmr_m>0.))
                    phase_fmr=phase_fmr+pi;

                /* ---------------------------------------------------------------- */
                /* CALCULATE PHASE OF fMRI WAVEFORM AT mth FREQUENCY COMPONENT.     */
                /*  Using algorithm kindly provided by Jagath Rajapakse, to return  */
                /* the correct angle in the complex plane.                          */
                if (IFTfmr_m == 0.) {
                    if (RFTfmr_m >= 0.)
                        phase_fmr=0.;
                    else
                        phase_fmr=pi;
                }
                else if (RFTfmr_m == 0) {
                    if (IFTfmr_m < 0)
                        phase_fmr=1.5*pi;
                    if (IFTfmr_m > 0)
                        phase_fmr=0.5*pi;
                }
                else if (RFTfmr_m < 0 & IFTfmr_m > 0)
                    phase_fmr=pi - atan(-IFTfmr_m/RFTfmr_m);
                else if (RFTfmr_m <  0 & IFTfmr_m < 0)
                    phase_fmr=pi + atan(IFTfmr_m/RFTfmr_m);
                else if (RFTfmr_m > 0 & IFTfmr_m < 0)
                    phase_fmr=2.0*pi - atan(-IFTfmr_m/RFTfmr_m);
                else
                    phase_fmr=atan(IFTfmr_m/RFTfmr_m);

                phase_diff=phase_fmr-phase_inw[var_index];

                /* Take into account slice acquisition sequence if desired.     */
                if ((absolute_time[0]=='y') || (absolute_time[0]=='Y')) {
                    z=vox_index/(xdim*ydim);
                    phase_diff-=(2.*pi*(float)slice_seq2[z]/(float)zdim);
                }

                /* Restrict phase difference to range (-pi,+pi).                */
                while (phase_diff<-pi) phase_diff+=pi;
                while (phase_diff>=pi) phase_diff-=pi;

                mean=(phase_diff*(float)num_dat_pts/(2.*pi*(float)m[var_index]));

                if (mean<min_lag)
                    lag_map[var_index][vox_index]=TR*min_lag;
                else if (mean>max_lag)
                    lag_map[var_index][vox_index]=TR*max_lag;
                else
                    lag_map[var_index][vox_index]=TR*mean;

                /* ------------------------------------------------------------ */
                /* CALCULATE DISPERSION FROM RELATIVE MAGNITUDES OF mth AND nth */
                /* FREQUENCY COMPONENTS OF POOLED fMR SIGNALS AND INPUT         */
                /* STIMULUS SQUAREWAVE.
                printf("Calculating dispersion...\n");                          */

                Ym=(float)sqrt((double)(RFTfmr_m*RFTfmr_m+IFTfmr_m*IFTfmr_m));
                Yn=(float)sqrt((double)(RFTfmr_n*RFTfmr_n+IFTfmr_n*IFTfmr_n));
                if (((Yn*Xm)==0.) || ((Ym*Xn)/(Yn*Xm)<1.))
                    SD = 0.;
                else
                    SD=(1./(float)sqrt(2.*pi))*
                        (float)num_dat_pts*
                        (float)sqrt(log((double)(Ym*Xn/(Yn*Xm)))
                                            /denominator_factor[var_index]);
                dispersion_map[var_index][vox_index]=SD*TR;

                /* ------------------------------------------------------------ */
                /* DETERMINE FREQUENCY COMPONENT WHICH HAS THE GREATEST MAG.    */
                /* THIS WILL BE USED TO DETERMINE THE BEST ESTIMATE OF LAG AND  */
                /* DISPERSION.                                                  */
                mag=(float)sqrt((double)((RFTfmr_m*RFTfmr_m)+(IFTfmr_m*IFTfmr_m)));
                if (var_index==0) {
                    max_mag=mag;
                    best_lag[vox_index]  = lag_map[var_index][vox_index];
                    best_disp[vox_index] = dispersion_map[var_index][vox_index];
                }
                else {
                    best_lag[vox_index]  = (mag>max_mag) ?
                        lag_map[var_index][vox_index] : best_lag[vox_index];
                    best_disp[vox_index] = (mag>max_mag) ?
                        dispersion_map[var_index][vox_index] : best_disp[vox_index];
                    max_mag = (mag>max_mag) ? mag : max_mag;
                }

            } /* END LOOP OVER INDEPENDENT VARIABLES. */

        } /* END IF CHECKING ON MASK. */

    } /* END LOOP OVER VOXELS. */
    printf("\n\n");

/* ---------------------------------------------------------------------------- */
/*  WRITE ESTIMATED REGRESSION COEFFICIENT MAPS AND WILKS LAMBDA MAPS TO        */
/*  FLOATING POINT FILES.                                                       */

    for (var_index=0; var_index<num_ind_var; var_index++) {
        float_fwrite(argv[12+3*var_index],lag_map[var_index],num_vox);
        float_fwrite(argv[13+3*var_index],dispersion_map[var_index],num_vox);
    }
    float_fwrite("best_lag.img",best_lag,num_vox);
    float_fwrite("best_disp.img",best_disp,num_vox);


/* ---------------------------------------------------------------------------- */
/*  FREE MEMORY.                                                                */

    free(slice_seq1);
    free(slice_seq2);
    free(img_name);
    free(lag_map);
    free(dispersion_map);
    free(temp);
    free(imgin);
    free(global_mean);
    free(mean_fmr);
    free(mask);

    free_dmatrix(X,1,num_dat_pts,1,num_ind_var);
}

