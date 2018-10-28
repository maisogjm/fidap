#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "analyze6.h"


                       /* FUNCTION DECLARATIONS.                                */
struct  analyze_struct readhdr(char *argv[]);

/* **************************************************************************** */
/* Program reads in names of images from a text file, then reads in those       */
/* images.  Images are scaled, and the global mean of each image is also scaled */
/* to some weighted average of the global means of the input images.            */
/* Parametric images are calculated by voxel-by-voxel curve-fitting, and they   */
/* are output to disk as MATLAB .mat files.  These .mat files can be converted  */
/* to an ANALYZE .img file using mtz.  The MATLAB programs read_param_files,    */
/* xyz, and plot_sigA can be used to analyze the output of this program.        */
/*                                                                              */
/* HISTORY: 02.11.93 : fitcurv2voxA written by Jose' Ma. Maisog, M.D.           */
/*          02.12.93 : fitcurv2voxA modified by JMM to make fitcurv2voxB        */
/*                     fitcurv2voxB normalizes the global mean of the input     */
/*                     images to the global mean of the resting state image,    */
/*                     or to some weighted average of the global means of the   */
/*                     resting state images if there are more than one.  This   */
/*                     weighted average is specified by weights given in the    */
/*                     input text file.                                         */
/*          03.02.93 : added code to fitcurv2voxB to output in addition a       */
/*                     MATLAB .mat file containing the variables scale_factor,  */
/*                     NUMRTR and glob_mean, because these variables will be    */
/*                     needed for image analysis in MATLAB                      */
/*                     (see read_param_filesA.m).                               */
/*          03.03.93 : changed X-offset parameter into a steepnesS parameter,   */
/*                     because an X-offset parameter is redundant -- the K      */
/*                     parameter controls X-offset.  Needed a parameter to      */
/*                     describe how steep the sigmoid is, so changed X-offset   */
/*                     to a "slope" parameter.  Note that this steepness        */
/*                     parameter is NOT a slope in the sense of a change in Y   */
/*                     per change in X.                                         */
/*                     Also, note that the function being minimized is not      */
/*                     the squared error per se but the squared error plus 1.   */
/*                     This is because the Nelder-Mead Simplex search algorithm */
/*                     as currently implemented relies on termination of the    */
/*                     search by calculating an intermediate value rtol, but    */
/*                     the way this rtol is calculated doesn't work for a cost  */
/*                     function which is minimized at a value of zero; and the  */
/*                     squared error function is minimized at zero.  However,   */
/*                     minimizing the squared error plus 1 will be equivalent   */
/*                     to minimizing the squared error itself.                  */
/*                     When the squared error image is written to disk as a     */
/*                     MATLAB .mat data file, the image is decremented by 1 so  */
/*                     the what is saved to disk is the true squared error, and */
/*                     not the squared error plus 1.                            */
/*          03.15.93 : JMM: Changed ALL floats to doubles.  All calculations    */
/*                     now done in double.                                      */
/*          10.08.93 : JMM : modified fitcurv2voxB.c to generate fmri_group_t.c */
/*                     fmri_group_t reads in 16-bit images from two groups,     */
/*                     calculates a t-map, and outputs the t-map in a MATLAB    */
/*                     .mat data file.  This .mat data file can be converted    */
/*                     into an 8-bit ANALYZE image using mtz.                   */
/*        10.31.94 : JMM : modified fmri_group_t.c to generate fmri_group_tC.c  */
/*                     Program reads in mask from a file rather than creating it*/
/*                      in memory.  Also, writes output into floating point     */
/*                      ANALYZE images rather than MATLAB .mat data files.      */
/*                     Accepts coordinates of voxel of interest as command line */
/*                      arguments.                                              */
/*          04.02.95 : JMM : Modified fmri_group_tC to generate fmri_group_tD.  */
/*                     The two contrast groups are now normalized to the same   */
/*                     across-group average global mean, rather than to within- */
/*                     group means.                                             */

main(int argc, char *argv[])
{
                            /* READING IN RAW DATA.                             */
    struct analyze_struct   /* Header information variable.                     */
                img_info;   /*                                                  */
    FILE *ftext;            /* Pointer for input text file.                     */
    FILE *fimgin;           /* Pointer for input image file.                    */
    char choice;            /* For keyboard input to yes/no question.           */
    char dummy[80];         /* Dummy pointer for counting number of entries     */
                            /* in a text file.                                  */
    char **header_name;     /* Pointer to pointer to name of header file.       */
    int grp_index;          /* Index indicating group 1 vs group 2.             */
    int n[2];               /* Number of images in groups 1 and 2.              */
    int n_index;            /* Index for group 1 and group 2 data.              */
    int max_n;              /* The maximum between n0] and n[2].                */
    char ***img_name;       /* Pointer to pointers to group 1 input file names. */
    int xdim,ydim,zdim;     /* Image dimensions.                                */
    int num_vox;            /* Number of voxels in an image = xdim*ydim*zdim.   */
    signed short *temp;     /* Temporary buffer for reading in input images.    */
    float ***imgin;         /* Pointer to pointers to input image data in       */
                            /* floating point form.                             */

                            /* INTENSITY SCALING & GLOBAL MEAN NORMALIZATION    */
    int vox_index;          /* Index into imgin[0][*][] & imgin[1][*][],        */
                            /* extracting a specific voxel value.               */
    FILE *fmean;            /* Pointer for output mean file.                    */
    float **mean;           /* Pointer to pointers to mean images.              */
    float *diff;            /* Pointer to difference of group means image.      */
    float **sqrerr;         /* Pointer to pointers to squared error images.     */
    float float_thresh;     /* Temporary variable for reading in threshold.     */
    float threshold;        /* If the images are to be normalized to some       */
                            /* global mean intensity, this is a threshold used  */
                            /* to generate a mask from the average image.       */
    float float_t_thresh;   /* Temporary variable for reading in significance   */
                            /* threshold.                                       */
    float t_thresh;         /* t threshold for significance.                    */
    FILE *fmask;            /* Pointer to mask file.                            */
    int write_status;       /* Returns number of elements written to file.      */
    int read_status;        /* Returns number of elements read from file.       */
    unsigned char *mask;    /* Unsigned char copy of mask.                      */
    int num_mask_vox;       /* Number of non-zero voxels in *mask.              */
    float **glob_mean;      /* Global mean intensity for both groups, for       */
                            /* output to disk.                                  */
    float float_c;          /* Temporary variable for reading in *c.            */
    float **c;              /* Weighting factor for normalizing global means.   */
    float grand_global_mean;/* Grand global mean of all input images, used for  */
                            /* normalizing global means.                        */
    float SUM;              /* Sum used for normalizing global means.           */

                            /* VOXEL OF INTEREST.                               */
    int x_int, y_int,       /* X, Y, and Z coordinates of a voxel of interest.  */
                z_int;      /* The ANOVA for this voxel of interest will be     */
                            /* output to a text file.                           */
    int vox_index_int;      /* The voxel index for the voxel of interest,       */
                            /* equal to (z_int*ydim+y_int)*xdim+x_int, after    */
                            /* x_int, y_int, and z_int have been decremented by */
                            /* 1 (pointer arithmetic starts counting at zero).  */

                            /* t-IMAGE AND OUTPUT TO DISK.                      */
    float *timg;            /* Pointer to output image data, ANALYZE format.    */
    float t_factor;         /* Factor used to calculate paired t.               */
    char *out_name;         /* Pointer to output .mat file name.                */
    int x,y,z;              /* Indices into ANALYZE format image arrays.        */
    FILE *fOUTTEXT;         /* Pointer to output text file.                     */
    char outtext_name[80];  /* Name of output text file.                        */
    FILE *fOUT;             /* Pointer for output floating point .img file.     */

/* ---------------------------------------------------------------------------- */
/*  CHECK TO SEE IF CORRECT NUMBER OF ARGUMENTS WERE PASSED.                    */

    if (argc != 10) {
        printf("fmri_group_tD: need header file, name of 8-bit mask file, input group\n");
        printf("               1 textfile, input group 2 textfile, output .img file\n");
        printf("               name for t tests, output .img file name for difference\n");
        printf("               of group means, AND XYZ coordinates of voxel of interest.\n");

        printf("\nThe input group 1 textfile should contain the pathnames\n");
        printf("(either absolute or relative) of the group 1 input images.\n") ;
        printf("The input group 2 textfile should similarly contain the pathnames\n");
        printf("of the group 2 input images.\n") ;
        printf("These input text files should contain the names of the input .img\n") ;
        printf("files and a weighting factor c(i) used for scaling of the global means,\n");
        printf("with one filename and weighting factor pair per line.\n\n");
        printf("To correct for differences in global mean intensity, the\n");
        printf("following multiplicative scaling will be done on an image-\n");
        printf("by-image basis.  A global mean GBM(i)  will be calculated for\n");
        printf("each image IMG(i).  Then each image IMG(i) will be multiplied by\n\n");

        printf("                  grand_global_mean/GBM(i)\n\n");

        printf("where grand_global_mean is a weighted average of all the GBM's.\n");
        printf("Let c(i) be the weight assigned to GBM(i).  Let\n\n");

        printf("           SUM = abs(c(1))+abs(c(2))+...\n\n");

        printf("for all the c(i)`s.  Then\n\n");

        printf("         grand_global_mean = [c(1)*GBM(1)+c(2)*GBM(2)+...]/SUM.\n\n");

        printf("For example, if you want the numerator to be the average of\n");
        printf("GBM(1) and GBM(6) set\n\n");

        printf("          c(1) = 1.0\n");
        printf("          c(2) = 0.0\n");
        printf("          c(3) = 0.0\n");
        printf("          c(4) = 0.0\n");
        printf("          c(5) = 0.0\n");
        printf("          c(6) = 1.0\n\n");

        printf("A text file will be output which will contain the intermediate\n");
        printf("variables grand_global_mean, glob_mean, and threshold (the number used to\n");
        printf("generate the mask for calculation of the global means).\n\n");

        exit(1);
    }

/* ---------------------------------------------------------------------------- */
/*  READ HEADER FILE; ASSUME IMAGES ARE SAME SIZE, SO USE SAME HEADER FILE      */
/*  FOR ALL IMAGES.                                                             */

    img_info = readhdr(argv+1);

    xdim = img_info.dime.dim[1];
    ydim = img_info.dime.dim[2];
    zdim = img_info.dime.dim[3];

    num_vox = xdim*ydim*zdim;

/* ---------------------------------------------------------------------------- */
/*  COUNT NUMBER OF ENTRIES IN GROUP 1 INPUT TEXT FILE.                         */

    if ((ftext = fopen(argv[3],"r")) == NULL) {
        printf ("fmri_group_tD: can't open %s\n",argv[3]);
        exit(2);
    }
    n[0] = 0;
    fscanf(ftext,"%s",dummy);
    while(fscanf(ftext,"%s %f",dummy,&float_c)!=EOF)
        n[0]++;
    printf("There are %d images in group 1.\n",n[0]);

/* ---------------------------------------------------------------------------- */
/* DETERMINE VOXEL OF INTEREST COORDINATES AND NAME OF OUTPUT TEXT FILE.        */
    x_int = atoi(argv[7]);
    y_int = atoi(argv[8]);
    z_int = atoi(argv[9]);

    sprintf(outtext_name,"group_t_test_%d_%d_%d.txt\0",x_int,y_int,z_int);

    vox_index_int = ((z_int-1)*ydim+(y_int-1))*xdim+(x_int-1);

/* ---------------------------------------------------------------------------- */
/*  READ IN FILENAMES AND GLOBAL MEAN WEIGHTS IN GROUP 1 INPUT TEXT FILE.       */

    c = (float **) malloc(2*sizeof(float *));
    c[0] = (float *) malloc(n[0]*sizeof(float));
    img_name = (char ***) malloc(2*sizeof(char **));
    img_name[0] = (char **) malloc(n[0]*sizeof(char *));

    rewind(ftext);
    for (n_index=0; n_index<=n[0]-1; n_index++) {
        img_name[0][n_index] = (char *) malloc(80*sizeof(char));
        fscanf(ftext,"%s %f",img_name[0][n_index],&float_c);
        c[0][n_index] = (float) float_c;
    }
    fclose(ftext);

/* ---------------------------------------------------------------------------- */
/*  COUNT NUMBER OF ENTRIES IN GROUP 2 INPUT TEXT FILE.                         */

    if ((ftext = fopen(argv[4],"r")) == NULL) {
        printf ("fmri_group_tD: can't open %s\n",argv[4]);
        exit(3);
    }
    n[1] = 0;
    fscanf(ftext,"%s",dummy);
    while(fscanf(ftext,"%s %f",dummy,&float_c)!=EOF)
        n[1]++;
    printf("There are %d images in group 2.\n",n[1]);

    max_n = (n[0]>n[1]) ? n[0] : n[1];

/* ---------------------------------------------------------------------------- */
/*  READ IN FILENAMES AND GLOBAL MEAN WEIGHTS IN GROUP 2 INPUT TEXT FILE.       */

    c[1] = (float *) malloc(n[1]*sizeof(float));
    img_name[1] = (char **) malloc(n[1]*sizeof(char *));

    rewind(ftext);
    for (n_index=0; n_index<=n[1]-1; n_index++) {
        img_name[1][n_index] = (char *) malloc(80*sizeof(char));
        fscanf(ftext,"%s %f",img_name[1][n_index],&float_c);
        c[1][n_index] = (float) float_c;
    }
    fclose(ftext);

/* ---------------------------------------------------------------------------- */
/*  SHOW USER INFORMATION READ IN SO FAR.                                       */

    printf("\nThese are the %d group 1 input images and weighting factors:\n", n[0]);
    for (n_index=0; n_index<=n[0]-1; n_index++)
        printf("  %20s   %f\n", img_name[0][n_index],c[0][n_index]);

    printf("\nThese are the %d group 2 input images and weighting factors:\n", n[1]);
    for (n_index=0; n_index<=n[1]-1; n_index++)
        printf("  %20s   %f\n", img_name[1][n_index],c[1][n_index]);

    printf("\nImages have dimensions:\n");
    printf("  xdim = %d\n  ydim = %d\n  zdim = %d\n",xdim,ydim,zdim);

    printf("\nThe t-test file will be %s.\n", argv[4]);
    printf("It is an floating point ANALYZE image.\n");

    printf("Your t-test has %d degrees of freedom.\n",n[0]+n[1]-2);
    printf("\nENTER the threshold t for significance: ");
    scanf("%f",&float_t_thresh);
    t_thresh = (float) float_t_thresh;
    printf("\nSignificance threshold = %3.2f.\n",t_thresh);

    printf("Coordinates of voxel of interest = [%d, %d, %d]\n",x_int,y_int,z_int);
    printf("Intermediate variables will be written to the text file %s.\n",
                                                                outtext_name);

    printf("\nDoes the above information look okay <y/n>? ");
    scanf("%s",&choice);
    while((choice!='y') && (choice!='n')) {
        printf("Type 'y' or 'n': ");
        scanf("%s",&choice);
    }
    if (choice=='n') exit(4);

/* ---------------------------------------------------------------------------- */
/*  READ IN GROUP 1 AND GROUP 2 IMAGES.  ALLOCATE MEMORY FOR OUTPUT IMAGE.      */

    temp  = (signed short *) malloc(num_vox*sizeof(signed short));
    imgin  = (float ***) malloc(2*sizeof(float **));
    imgin[0]  = (float **) malloc(n[0]*sizeof(float *));
    imgin[1]  = (float **) malloc(n[1]*sizeof(float *));

    printf("Reading in group 1 images...\n");
    for (n_index=0; n_index<=n[0]-1; n_index++) {
        imgin[0][n_index] = (float *) malloc(num_vox*sizeof(float));
        printf("  Reading in image from %s....\n",img_name[0][n_index]);
        if ((fimgin=fopen(img_name[0][n_index],"r"))==NULL) {
            printf("fmri_group_tD: unable to open %s\n",img_name[0][n_index]);
            exit(6);
        }
        fread(temp, sizeof(signed short), num_vox, fimgin);
        for (vox_index=0; vox_index<num_vox; vox_index++)
            imgin[0][n_index][vox_index]=(float)temp[vox_index];
        fclose(fimgin);
    }

    printf("Reading in group 2 Images...\n");
    for (n_index=0; n_index<=n[1]-1; n_index++) {
        imgin[1][n_index] = (float *) malloc(num_vox*sizeof(float));
        printf("  Reading in image from %s....\n",img_name[1][n_index]);
        if ((fimgin=fopen(img_name[1][n_index],"r"))==NULL) {
            printf("fmri_group_tD: unable to open %s\n",img_name[0][n_index]);
            exit(6);
        }
        fread(temp, sizeof(signed short), num_vox, fimgin);
        for (vox_index=0; vox_index<num_vox; vox_index++)
            imgin[1][n_index][vox_index]=(float)temp[vox_index];
        fclose(fimgin);
    }

    timg = (float *) malloc(num_vox*sizeof(float));

/* ---------------------------------------------------------------------------- */
/*  NORMALIZE GLOBAL MEANS.                                                     */

    /* ------------------------------------------------------------------------ */
    /* ALLOCATE MEMORY FOR mean, sqrerr, AND glob_mean.                                */

    mean = (float **) malloc(3*sizeof(float *));
    mean[0] = (float *) malloc(num_vox*sizeof(float));
    mean[1] = (float *) malloc(num_vox*sizeof(float));
    mean[2] = (float *) malloc(num_vox*sizeof(float));
    if ((mean == NULL) || (mean[0] == NULL) ||
                                (mean[1] == NULL) || (mean[2] == NULL)) {
        printf("fmri_group_tD: failed to malloc for mean.\n");
        exit(5);
    }
    sqrerr = (float **) malloc(2*sizeof(float *));
    sqrerr[0] = (float *) malloc(num_vox*sizeof(float));
    sqrerr[1] = (float *) malloc(num_vox*sizeof(float));
    if ((sqrerr == NULL) || (sqrerr[0] == NULL) ||
                                (sqrerr[1] == NULL)) {
        printf("fmri_group_tD: failed to malloc for sqrerr.\n");
        exit(5);
    }
    glob_mean = (float **) malloc(2*sizeof(float *));
    glob_mean[0] = (float *) malloc(n[0]*sizeof(float));
    glob_mean[1] = (float *) malloc(n[1]*sizeof(float));
    if ((glob_mean==NULL) || (glob_mean[0]==NULL) || (glob_mean[1]==NULL) ) {
        printf("fmri_group_tD: failed to malloc for glob_mean.\n");
        exit(5);
    }
    diff = (float *) malloc(num_vox*sizeof(float));
    if (diff == NULL) {
        printf("fmri_group_tD: failed to malloc for diff.\n");
        exit(5);
    }


    /* ------------------------------------------------------------------------ */
    /* INITIALIZE OVERALL MEAN IMAGE.                                           */

    for (vox_index=0; vox_index<num_vox; vox_index++)
        mean[2][vox_index] = 0.;

    /* ------------------------------------------------------------------------ */
    /* CALCULATE OVERALL MEAN IMAGE.                                            */

    printf("Calculating overall mean image...\n");
    for (grp_index=0; grp_index<=1; grp_index++)
        for (n_index=0; n_index<n[grp_index]; n_index++)
            for (vox_index=0; vox_index<num_vox; vox_index++)
                mean[2][vox_index] += imgin[grp_index][n_index][vox_index];

    for (vox_index=0; vox_index<num_vox; vox_index++)
            mean[2][vox_index] /= ((float) (n[0] + n[1]));


    /* ------------------------------------------------------------------------ */
    /* READ IN MASK FROM DISK.                                                  */

    printf("Reading in mask from %s...\n", argv[2]);
    mask = (unsigned char *) malloc(num_vox*sizeof(unsigned char));
    if (mask == NULL) {
        printf("fmri_group_tD: Unable to open %s for writing to.\n",argv[2]);
        exit(15);
    }
    if ((fmask = fopen(argv[2],"r")) == NULL) {
        printf("fmri_group_tD: Unable to open %s for writing to.\n",argv[2]);
        exit(16);
    }
    else {
        read_status=fread(mask, sizeof(unsigned char), num_vox, fmask);
        if (read_status != num_vox) {
            printf("Error reading from %s\n",argv[2]);
            fclose(fmask);
            exit(17);
        }
    }

    /* -------------------------------------------------------------------- */
    /* MASK OUT mean[2].                                                    */

    printf("Masking out overall mean image...\n");
    for (vox_index=0; vox_index<num_vox; vox_index++)
        mean[2][vox_index] *= mask[vox_index];


    /* ------------------------------------------------------------------------ */
    /* COUNT NUMBER OF NON-ZERO VOXELS IN *mask.                                */

    num_mask_vox=0.;
    for (vox_index=0; vox_index<num_vox; vox_index++)
        if (mask[vox_index]!=0.)
            num_mask_vox++;
    printf("fmri_group_tD: number of voxels in mask = %d\n", num_mask_vox);

    /* ------------------------------------------------------------------------ */
    /* CALCULATE WITHIN-VOLUME GLOBAL MEANS OF ALL INPUT IMAGE VOLUMES.         */

    for (grp_index=0; grp_index<=1; grp_index++)
        for (n_index=0; n_index<n[grp_index]; n_index++) {
            glob_mean[grp_index][n_index] = 0.;
            for (vox_index=0; vox_index<num_vox; vox_index++)
                if (mask[vox_index]!=0)
                    glob_mean[grp_index][n_index] +=
                                        imgin[grp_index][n_index][vox_index];
            glob_mean[grp_index][n_index] /= ((float) num_mask_vox);
            printf("global mean of %s = %f\n", img_name[grp_index][n_index],
                                                        glob_mean[grp_index][n_index]);
        }

    /* -------------------------------------------------------------------- */
    /* CALCULATE grand_global_mean.                                         */

    grand_global_mean = 0.;
    SUM    =0.;
    for (grp_index=0; grp_index<=1; grp_index++)
        for (n_index=0; n_index<n[grp_index]; n_index++) {
            grand_global_mean += 
                            	(c[grp_index][n_index] * glob_mean[grp_index][n_index]);
            SUM += (fabs(c[grp_index][n_index]));
        }

    grand_global_mean /= SUM;
    printf("grand_global_mean = GRAND GLOBAL MEAN = %f\n", grand_global_mean);
            

    /* ------------------------------------------------------------------------ */
    /* LOOP OVER TWO GROUPS, NORMALIZE INDIVIDUAL VOLUMES, CALCULATE MEAN AND   */
    /* SQUARED ERRORS IMAGES.                                                   */

    for (grp_index=0; grp_index<=1; grp_index++) {

        /* -------------------------------------------------------------------- */
        /* NORMALIZE INDIVIDUAL IMAGES.  MULTIPLY EACH VOXEL BY                 */
        /* grand_global_mean/glob_mean[i].                                      */

        printf("Normalizing group %d global means...\n",grp_index+1);
        for (n_index=0; n_index<n[grp_index]; n_index++)
            for (vox_index=0; vox_index<num_vox; vox_index++)
                imgin[grp_index][n_index][vox_index] *=
                        ((float)grand_global_mean/(float)glob_mean[grp_index][n_index]);

        /* -------------------------------------------------------------------- */
        /* INITIALIZE GROUP MEAN IMAGE.                                         */

        printf("Initializing group %d mean image...\n",grp_index+1);
        for (vox_index=0; vox_index<num_vox; vox_index++)
            mean[grp_index][vox_index] = 0.;

        /* -------------------------------------------------------------------- */
        /* CALCULATE GROUP MEAN IMAGE.                                          */

        printf("Calculating group %d mean image...\n",grp_index+1);
        for (n_index=0; n_index<n[grp_index]; n_index++)
            for (vox_index=0; vox_index<num_vox; vox_index++)
                mean[grp_index][vox_index] += imgin[grp_index][n_index][vox_index];

        for (vox_index=0; vox_index<num_vox; vox_index++)
            mean[grp_index][vox_index] /= ((float) (n[grp_index]));

        /* -------------------------------------------------------------------- */
        /* INITIALIZE SQUARED ERRORS IMAGE.                                     */

        printf("Initializing group %d squared errors image...\n",grp_index+1);
        for (n_index=0; n_index<n[grp_index]; n_index++)
            for (vox_index=0; vox_index<num_vox; vox_index++)
                sqrerr[grp_index][vox_index] = 0.;

        /* -------------------------------------------------------------------- */
        /* CALCULATE SQUARED ERRORS IMAGE.                                      */

        printf("Calculating group %d squared errors image...\n",grp_index+1);
        for (n_index=0; n_index<n[grp_index]; n_index++)
            for (vox_index=0; vox_index<num_vox; vox_index++)
                sqrerr[grp_index][vox_index] +=
                    (imgin[grp_index][n_index][vox_index]-mean[grp_index][vox_index])
                                                *
                    (imgin[grp_index][n_index][vox_index]-mean[grp_index][vox_index]);

    }

    /* ------------------------------------------------------------------------ */
    /* SAVE grand_global_mean, glob_mean, & scale_factor TO OUTPUT FILE.        */

    if ((fOUTTEXT = fopen(outtext_name,"w")) == NULL) {
        printf ("fmri_group_tD: can't open %s.\n",outtext_name);
        exit(2);
    }

    printf("Writing output to %s...\n",outtext_name);
    fprintf(fOUTTEXT,"T-TEST FOR VOXEL AT COORDINATES [%d,%d,%d]\n\n",
                                                        x_int,y_int,z_int);
    fprintf(fOUTTEXT,"Coordinates of voxel of interest = [%d, %d, %d]\n",
                                                                x_int,y_int,z_int);
    fprintf(fOUTTEXT,"Voxel index of interest = %d\n",vox_index_int);
    fprintf(fOUTTEXT,"grand_global_mean = %f\n",grand_global_mean);
    fprintf(fOUTTEXT,"input text file # 1 = %s\n",argv[3]);
    fprintf(fOUTTEXT,"glob_mean1 =\n");
    for (n_index=0; n_index<n[0]; n_index++)
        fprintf(fOUTTEXT,"  %s = %f\n",
                        img_name[0][n_index],glob_mean[0][n_index]);
    fprintf(fOUTTEXT,"\n");
    fprintf(fOUTTEXT,"input text file # 2 = %s\n",argv[4]);
    fprintf(fOUTTEXT,"glob_mean2 =\n");
    for (n_index=0; n_index<n[1]; n_index++)
        fprintf(fOUTTEXT,"  %s = %f\n",
                        img_name[1][n_index],glob_mean[1][n_index]);
    fprintf(fOUTTEXT,"\n");
    fprintf(fOUTTEXT,"mask = %s\n",argv[2]);
    fprintf(fOUTTEXT,"t-score threshold = %f\n",t_thresh);


/* ---------------------------------------------------------------------------- */
/*  CALCULATE DIFFERENCE OF MEAN IMAGES.                                        */

    printf("Calculating difference of mean images...\n");
    for (vox_index=0; vox_index<num_vox; vox_index++)
        diff[vox_index]=(mask[vox_index]==1.)
                ? mean[1][vox_index]-mean[0][vox_index] : 0.;

    /* ------------------------------------------------------------------------ */
    /*  AND THEN SAVE TO DISK.                                                  */

    if ((fOUT = fopen(argv[6],"w")) == NULL) {
        printf ("fmri_group_tD: can't open %s\n",argv[6]);
        exit(18);
    }
    printf("Writing difference image to %s...\n",argv[6]);
    fwrite(diff,sizeof(float),num_vox,fOUT);
    fclose(fOUT);


/* ---------------------------------------------------------------------------- */
/*  NOW CALCULATE T-TESTS.                                                      */

    printf("Calculating t tests...\n");
    t_factor = (1./((float) (n[0]+n[1]) - 2.))
                                        * ((1./(float) n[0])+(1./(float) n[1]));
/*    printf("t_factor = %f\n",t_factor);
    fprintf(fOUTTEXT,"t_factor = %f\n",t_factor); */
    for (vox_index=0; vox_index<num_vox; vox_index++)
        if (sqrerr[0][vox_index]+sqrerr[1][vox_index] == 0.)
            timg[vox_index] = 0.;
        else
            timg[vox_index] =(float) ((float)mask[vox_index] * (mean[1][vox_index] - mean[0][vox_index]))
                                                /
        ((float) sqrt((double)(t_factor*(sqrerr[0][vox_index]+sqrerr[1][vox_index]))));

    printf("Group 1 mean at [%d, %d, %d] = %f\n",x_int,y_int,z_int,mean[0][vox_index_int]);
    printf("Group 2 mean at [%d, %d, %d] = %f\n",x_int,y_int,z_int,mean[1][vox_index_int]);
    printf("Group 1 squared error at [%d, %d, %d] = %f\n",x_int,y_int,z_int,sqrerr[0][vox_index_int]);
    printf("Group 2 squared error at [%d, %d, %d] = %f\n",x_int,y_int,z_int,sqrerr[1][vox_index_int]);
    printf("t-test at [%d, %d, %d] = %f\n",x_int,y_int,z_int,timg[vox_index_int]);
    printf("df = %d\n",n[0]+n[1]-2);
    fprintf(fOUTTEXT,"Group 1 mean at [%d, %d, %d] = %f\n",x_int,y_int,z_int,mean[0][vox_index_int]);
    fprintf(fOUTTEXT,"Group 2 mean at [%d, %d, %d] = %f\n",x_int,y_int,z_int,mean[1][vox_index_int]);
    fprintf(fOUTTEXT,"Group 1 squared error at [%d, %d, %d] = %f\n",x_int,y_int,z_int,sqrerr[0][vox_index_int]);
    fprintf(fOUTTEXT,"Group 2 squared error at [%d, %d, %d] = %f\n",x_int,y_int,z_int,sqrerr[1][vox_index_int]);
    fprintf(fOUTTEXT,"t-test at [%d, %d, %d] = %f\n",x_int,y_int,z_int,timg[vox_index_int]);
    fprintf(fOUTTEXT,"df = %d\n",n[0]+n[1]-2);

    fclose(fOUTTEXT);

/* ---------------------------------------------------------------------------- */
/*  WRITE OUT PARAMETERS TO OUTPUT FILES.                                       */

    /* ------------------------------------------------------------------------ */
    /*  PLACE INDICATOR SQUARE IN LOWER LEFT CORNER.                            */

    for (x = 0; x <= 5; x++)
        for (y = 0; y <= 5; y++)
            for (z = 0; z <= zdim-1; z++)
                timg[(z*ydim+y)*xdim + x]=(float)t_thresh;


    /* ------------------------------------------------------------------------ */
    /*  AND THEN SAVE TO DISK.                                                  */

    if ((fOUT = fopen(argv[5],"w")) == NULL) {
        printf ("fmri_group_tD: can't open %s\n",argv[5]);
        exit(18);
    }
    printf("Writing t-test image to %s...\n",argv[5]);
    fwrite(timg,sizeof(float),num_vox,fOUT);
    fclose(fOUT);

/* ---------------------------------------------------------------------------- */
/*  FREE MEMORY.                                                                */

    printf("Releasing memory...\n");
    free(header_name);
    free(c);
    free(img_name);
    free(temp);
    free(imgin);
    free(timg);
    free(mean);
    free(mask);
}
