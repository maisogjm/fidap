#include <stdio.h>
#include <stdlib.h>
#include "analyze6.h"

#define MAX_STR 1024
struct analyze_struct readhdr(char *argv[]);
unsigned char *unsigned_char_malloc(int num_items);
float *float_malloc(int num_items);
float **float_malloc2(int num_items1, int num_items2);
unsigned short *unsigned_short_malloc(int num_items);
int read_unsigned_char_data(char *filename, unsigned char *array, int num_items);
int read_unsigned_short_int_data(char *filename, unsigned short int *array, int num_items);

                            /* EXTERNAL VARIABLES.                              */
int num_vox;                /* Number of voxels in an image.                    */
int vox_index;              /* Index into imgin[][][*], extracting a specific   */
                            /* voxel.                                           */
int xdim,ydim,zdim;         /* Image dimensions.                                */
unsigned short ***imgin;    /* Pointer to input image data.                     */
                            /* Indexing: imgin[voxel][run][scan]                */

                            /* EXTERNAL VARIABLES, INDICES AND CONSTANTS.       */
int num_runs;               /* Number of fMR runs.                              */
int run_index;              /* Index into run number.                           */
int num_dat_pts;            /* The number of time points in a run.              */
int num_runsXnum_dat_pts;   /* num_runs times num_dat_pts                       */
int dat_index;              /* Index into time points.                          */
int tot_num_scans;          /* Total number of scans.  Should be equal to       */
                            /* num_runs * num_dat_pts                           */


/* **************************************************************************** */
/* Program extract_time_series.c                                                */
/* ANSI C Code by Jose' Ma. Maisog, M.D.                                        */
/* Section on Functional Brain Imaging,                                         */
/* Laboratory of Psychology and Psychopathology                                 */
/* National Institute of Mental Health                                          */
/*                                                                              */
/* Program reads in fMRI data for several runs, normalizes global mean intensi- */
/* ties with a ratio normalization, and then outputs the fMRI time series for   */
/* specified voxel indices of interest.  No other processing is done.           */
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
/* 12.22.95 : JMM : Modified fmri_lagNdisp9 to generate extract_time_series.    */
/* Program just outputs fMRI time series for specified voxel indices of         */
/* interest                                                                     */

main(int argc, char *argv[])
{
                            /* INPUT IMAGES.                                    */
    FILE *ftxt;             /* Pointer for input text file containing filenames */
    char ***img_name;       /* Pointer to pointers to input image file names.   */
                            /* Indexing:                                        */
                            /*      img_name[run][timepoint][ASCII character]   */
    char choice;            /* For keyboard input to yes/no question.           */
    char dummy[MAX_STR];         /* Dummy pointer for counting number of entries     */
                            /* in a text file.                                  */
    int check_fscanf;       /* Flag to check error status of fscanf operation.  */
    struct analyze_struct
                img_info;   /* Header information variable.                     */
    FILE *fimgin;           /* Pointer for input image file.                    */
    unsigned                /* Pointer to temporary short int array, for        */
        short int *temp;    /* reading in 16-bit fMRI data.                     */
                            /* byte data.                                       */
    float **global_mean;    /* Global mean intensity of a scan.                 */
    float grand_global_mean;/* Grand global mean intensity across scans.        */
    float mag_fmr_series;   /* Magnitude of input image data.                   */
    float mag_idealized_response; 
                            /* Magnitude of input image data.                   */
    int num_objs_read_in;   /* Number of items read in by fread.                */
    int num_files_read_in;  /* Number of files read in so far.                  */

                            /* MASK                                             */
    FILE *fMASK;            /* Pointer for input mask .img file.                */
    unsigned char *mask;    /* Pointer to mask image.                           */
    int num_mask_vox;       /* Number of non-zero voxels in *mask.              */
    int count;              /* Number of voxels processed thus far.             */

                            /* OUTPUT TEXTFILE.                                 */
    int argv_index;         /* Index over command line arguments.               */
    FILE *fOUTTEXT;         /* File pointer to output text file.                */
    char outtext_name[MAX_STR];  /* Array to contain name of output text file.       */
    int x,y,z;              /* XYZ coordinates.                                 */

/* ---------------------------------------------------------------------------- */
/*  CHECK TO SEE IF CORRECT NUMBER OF ARGUMENTS WERE PASSED.                    */

    if (argc<7) {
        printf("extract_time_series: need input header file, input textfile of\n");
        printf("                     file names, number of runs, number of scans\n");
        printf("                     per run, name of 8-bit mask, and list of\n");
        printf("                     voxel indices of interest.\n\n");

        printf("The input text file should contain the names of the input .img\n") ;
        printf("files.\n");
        exit(1);
    }


/* ---------------------------------------------------------------------------- */
/*  READ IN num_runs AND num_dat_pts OFF OF COMMAND LINE.                       */
/*  COUNT NUMBER OF DATA POINTS IN INPUT TEXT FILE.                             */

    num_runs = atoi(argv[3]);
    num_dat_pts = atoi(argv[4]);
    num_runsXnum_dat_pts=num_runs*num_dat_pts;
    printf("There are %d runs.\n", num_runs);
    printf("There are %d input scans per run.\n", num_dat_pts);
    printf("This makes a total of %d input scans.\n",num_runsXnum_dat_pts);
    printf("Will read in input scan filenames from %s.\n",argv[2]);


/* ---------------------------------------------------------------------------- */
/*  READ HEADER FILE; ASSUME IMAGES ARE SAME SIZE, SO USE SAME HEADER FILE      */
/*  FOR ALL IMAGES.                                                             */

    printf("Reading header info from %s...\n",argv[1]);
    img_info = readhdr(argv+1);

    xdim = img_info.dime.dim[1];
    ydim = img_info.dime.dim[2];
    zdim = img_info.dime.dim[3];
    num_vox = xdim*ydim*zdim;
    printf("Images have dimensions:\n");
    printf("  xdim = %d\n  ydim = %d\n  zdim = %d\n\n",xdim,ydim,zdim);


/* ---------------------------------------------------------------------------- */
/*  COUNT NUMBER OF INPUT IMAGE FILENAMES.                                      */

    if ((ftxt = fopen(argv[2],"r")) == NULL) {
        printf("extract_time_series: can't open %s\n",argv[2]);
        exit(4);
    }
    tot_num_scans = 0;
    while(fscanf(ftxt,"%s",dummy)!=EOF)
        tot_num_scans++;

    if (tot_num_scans != num_runsXnum_dat_pts) {
        printf("extract_time_series: number of entries in %s not equal to number\n",argv[2]);
        printf("                of runs times number of scans per run.\n");
        exit(10);
    }


/* ---------------------------------------------------------------------------- */
/*  ALLOCATE MEMORY FOR img_name; THEN READ IN NAMES OF INPUT IMAGE FILES.      */

    img_name = (char ***) malloc(num_runs*sizeof(char **));
    if (img_name==NULL) {
        printf("extract_time_series: error malloc'ing for img_name.\n");
        exit(12);
    }

    for (run_index=0; run_index<num_runs; run_index++) {
        img_name[run_index]=(char **) malloc(num_dat_pts*sizeof(char *));
        if (img_name[run_index]==NULL) {
            printf("extract_time_series: error malloc'ing for img_name[%d].\n",
                                                                        run_index);
            exit(14);
        }
    }
    rewind(ftxt);
printf("num_dat_pts = %d\n",num_dat_pts);
    for (run_index=0; run_index<num_runs; run_index++)
        for (dat_index=0; dat_index<num_dat_pts; dat_index++) {
            img_name[run_index][dat_index] = (char *) malloc(MAX_STR*sizeof(char));
            if (img_name[run_index][dat_index]==NULL) {
                printf("extract_time_series: error malloc'ing for img_name[%d][%d].\n",
                                                        run_index,dat_index);
                exit(17);
            }
            check_fscanf=fscanf(ftxt,"%s",img_name[run_index][dat_index]);
printf("%s\n",img_name[run_index][dat_index]);
            if (check_fscanf==EOF) {
                printf("extract_time_series: error reading run %d timepoint %d filename from %s.\n",
                                                run_index+1,dat_index+1,argv[2]);
                exit(18);
            }
        }
    fclose(ftxt);
printf("CHECK: %s\n",img_name[0][0]);


/* ---------------------------------------------------------------------------- */
/*  ASK USER IF IT'S OKAY TO PROCEED.                                           */

    printf("Mask will be read in from 8-bit image %s.\n",argv[5]);

    printf("Will output fMRI signals for voxel indices:");
    for (argv_index=6; argv_index<argc; argv_index++) 
        printf(" %d",atoi(argv[argv_index]));
    printf("\n\n");
    printf("\nDoes the above information look okay <y/n>? ");
printf("CHECK: %s\n",img_name[0][0]);
    scanf("%c",&choice);
printf("CHECK: %s\n",img_name[0][0]);
    while((choice!='y') && (choice!='n') && (choice!='Y') && (choice!='N')) {
        printf("Type 'y' or 'n': ");
        scanf("%s",&choice);
    }
printf("CHECK: %s\n",img_name[0][0]);
    if ((choice=='n') || (choice=='N')) exit(19);
printf("CHECK: %s\n",img_name[0][0]);


/* ---------------------------------------------------------------------------- */
/* READ IN 8-BIT MASK.                                                          */

    mask = unsigned_char_malloc(num_vox);
    read_unsigned_char_data(argv[5],mask,num_vox);

    num_mask_vox=0;
    for (vox_index=0; vox_index<num_vox; vox_index++)
        if (mask[vox_index]!=0)
            num_mask_vox++;
    printf("Number of nonzero voxels in mask = %d.\n\n",num_mask_vox);
printf("CHECK: %s\n",img_name[0][0]);


    /* -------------------------------------------------------------------- */
    /* INITIALIZE OUTPUT FILES.                                             */
    for (argv_index=6; argv_index<argc; argv_index++) {
        vox_index=atoi(argv[argv_index]);
        z=vox_index/(xdim*ydim);
        y=(vox_index-(z*xdim*ydim))/xdim;
        x=((vox_index-(z*xdim*ydim))-y*xdim);
        x++;
        y++;
        z++;
        sprintf(outtext_name,"voxel_%d_%d_%d.dat",x,y,z);
        if ((fOUTTEXT = fopen(outtext_name,"w")) == NULL) {
                printf ("extract_time_series: can't open %s.\n",outtext_name);
                exit(2);
        }
        fclose(fOUTTEXT);
    } /* END LOOP OVER VOXEL INDICES OF INTEREST. */
printf("CHECK: %s\n",img_name[0][0]);

    /* -------------------------------------------------------------------- */
    /* LOOP OVER SCANS AND VOXEL INDICES OF INTEREST.                       */

    temp     = unsigned_short_malloc(num_vox);
    grand_global_mean=0.;
    global_mean=float_malloc2(num_runs,num_dat_pts);
printf("CHECK: %s\n",img_name[0][0]);
    printf("Now printing out time series for %d timepoints...\n",num_dat_pts);
    printf("Timepoint:\n");
    for (dat_index=0; dat_index<num_dat_pts; dat_index++) {
        printf("%d\n",dat_index+1);
        for (run_index=0; run_index<num_runs; run_index++) {
            read_unsigned_short_int_data(img_name[run_index][dat_index], temp, num_vox);
            global_mean[run_index][dat_index]=0.;
            for (vox_index=0; vox_index<num_vox; vox_index++)
                if (mask[vox_index]!=0)
                    global_mean[run_index][dat_index]+=(float)temp[vox_index];
            global_mean[run_index][dat_index]/=(float)num_mask_vox;
            grand_global_mean+=global_mean[run_index][dat_index];
        
            for (argv_index=6; argv_index<argc; argv_index++) {
                vox_index=atoi(argv[argv_index]);
                z=vox_index/(xdim*ydim);
                y=(vox_index-(z*xdim*ydim))/xdim;
                x=((vox_index-(z*xdim*ydim))-y*xdim);
                x++;
                y++;
                z++;
                sprintf(outtext_name,"voxel_%d_%d_%d.dat",x,y,z);
                if ((fOUTTEXT = fopen(outtext_name,"a")) == NULL) {
                        printf ("extract_time_series: can't open %s.\n",outtext_name);
                        exit(2);
                }
                fprintf(fOUTTEXT,"%10g ",(float)temp[atoi(argv[argv_index])]
					/global_mean[run_index][dat_index]);
                fclose(fOUTTEXT);

            } /* END LOOP OVER VOXEL INDICES OF INTEREST. */
        } /* END LOOP OVER RUNS. */

        /* -------------------------------------------------------------------- */
        /* PUT NEWLINES AFTER STRING OF NUMBERS.                                */
        for (argv_index=6; argv_index<argc; argv_index++) {
            vox_index=atoi(argv[argv_index]);
            z=vox_index/(xdim*ydim);
            y=(vox_index-(z*xdim*ydim))/xdim;
            x=((vox_index-(z*xdim*ydim))-y*xdim);
            x++;
            y++;
            z++;
            sprintf(outtext_name,"voxel_%d_%d_%d.dat",x,y,z);
            if ((fOUTTEXT = fopen(outtext_name,"a")) == NULL) {
                    printf ("extract_time_series: can't open %s.\n",outtext_name);
                    exit(2);
            }
            fprintf(fOUTTEXT,"\n");
            fclose(fOUTTEXT);
        } /* END LOOP OVER VOXEL INDICES OF INTEREST. */
    } /* END LOOP OVER TIMEPOINTS. */
    grand_global_mean/=(float)tot_num_scans;
    printf("\n\nAll finished!  Output has been normalized to a global mean intensity of 1.\n");
    printf("The grand global mean was %g, so you could multiply the output by this number\n",grand_global_mean);
    printf("to bring the output back into the usual fMRI scale in the hundreds.\n");
    if ((fOUTTEXT = fopen("grand_global_mean.dat","w")) == NULL) {
        printf ("extract_time_series: can't open grand_global_mean.dat.\n");
        exit(2);
    }
    fprintf(fOUTTEXT,"%10g ",grand_global_mean);
    fclose(fOUTTEXT);



/* ---------------------------------------------------------------------------- */
/*  FREE MEMORY.                                                                */

    free(img_name);
    free(temp);
    free(imgin);
    free(mask);
    free(global_mean);
    free(temp);
    free(img_name);
}
