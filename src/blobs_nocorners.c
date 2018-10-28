#include <stdio.h>
#include <stdlib.h>
#include <float.h>        /* FOR FLT_MAX                                                */
#include "analyze6.h"
#include <math.h>

                        /* GLOBAL VARIABLES.                                            */
int *clustermask;       /* Mask used to fill in clusters in recursive routine.          */
int *clustermask2;      /* Mask used to fill in outimg in recursive routine.            */
int *intensitymask;     /* Mask used to fill in clusterintensity in recursive routine.  */
int *probmask;          /* Mask used to fill inusterprob in recursive routine.          */
float *image;           /* Input image to be analyzed.                                  */
unsigned short *outimg; /* Output image of clusters.                                    */
int xdim, ydim, zdim;   /* Dimensions of image.                                         */
int num_cluster;        /* Cumulative number of clusters found.                         */
int clustersize;        /* Size of a cluster.  Initialize to zero before                */
                        /* invoking growcluster.                                        */
float threshold;        /* Upper threshold for z-scores, in single-precision            */
                        /* floating point.                                              */
double xmean, ymean,    /* Mean coordinate values and intensity value for               */
        zmean, meanval; /* cluster.                                                     */
double Probability;     /* Probability that n is >= k.                                  */
double Probability2;    /* Probability that n is >= k, based on input FWHM or W.        */
float *clusterintensity;/* Cluster intensity image.                                     */
float *clusterprob;     /* Probability cluster image.                                   */

struct analyze_struct readhdr(char *argv[]);
void growcluster(int x, int y, int z);
void filloutimg(int x, int y, int z);
void fillintensityimg(int x, int y, int z);
void fillprobimg(int x, int y, int z, float fillval);
float gammln(float xx);
double PnGEx(double beta, double x, double D);
double PnmaxGEk(double beta, double x, double D, double Em);
double PnGEx2(double beta, double x, double D);
double PnmaxGEk2(double beta, double x, double D, double Em);
float erffc(float x);

/* **************************************************************************** */
main(int argc, char *argv[])
/*
HISTORY:
DEC-91: written by Jose' Ma. Maisog, M.D.
09.30.92 : JMM : modified from original findlocalmax.c to simply find
                 local maxima in a .mat file of dimensions 65x87, so that
                 it can be used on Barry's correlational output.
10.02.92 : JMM : modified from findmaxA.c to accept threshold and below_thresh as
                 command line arguments
02.05.93 : JMM : modified from findmaxB.c to work on ANALYZE .img files
                 rather than MATLAB .mat files.
12:20.93 : JMM : Modified from findmaxE to count clusters.
03.08.94 : JMM : Started modifications from clusterGE to do cluster/particle
                 analysis as proposed by Friston/Worsley.

Given a floating pointing ANALYZE image file and a threshold, finds
clusters of voxels GREATER THAN OR EQUAL TO a threshold.
A weighted average coordinate is generated and outputed as well.

This program makes no assumption about image dimensions.
They are read in from the header file.
*/
{
    FILE *ftext;
    FILE *fIN;
    FILE *fOUT;

    struct analyze_struct img_info; /* Header information variable.        */
    int num_read_in;                /* Number of voxels read in.           */
    int num_wrote_out;              /* Number of voxels written out.       */
    int *num_size;                  /* Number of clusters of a certain     */
                                    /* size.                               */
    int vox_index,                  /* Index into voxels.                  */
        num_vox,                    /* Total number of voxels.             */
        cluster_index,              /* Index into clusters.                */
        num_suprathresh_vox;        /* Number of SUPRA-threshold voxels.   */
    unsigned char *mask;            /* Mask to be applied to input         */
                                    /* statistical image.                  */
    int x, y, z;                    /* Indices for looping over 3          */
                                    /* spatial dimensions.                 */
    int xx, yy, zz;                 /* Indices for looping over 3          */
                                    /* spatial dimensions.                 */
    int not_max;                    /* Flag indicating whether or not a    */
                                    /* given voxel is a local maximum.     */

                                    /* VARIABLES FOR DOING PARTICLE        */
                                    /* ANALYSIS.                           */
    int num_vox_in_mask;            /* Number of voxels in the mask.       */
    double pi;                      /* Constant containing value for Pi    */
    double xvoxdim, yvoxdim,        /* Voxel sizes along the 3             */
                        zvoxdim;    /*         dimensions.                 */
    double vox_vol;                 /* Volume of one voxel.                */
    double S;                       /* Total volume of voxels in mask,     */
                                    /*    = num_vox_in_mask * vox_vol      */
    int num_derivatives;            /* Number of derivatives in one row    */
                                    /* or column.                          */
    double derivative_sum;          /* Sum of derivatives.                 */
    int total_num_derivatives;      /* Number of derivatives in sum        */
                                    /* along one axis.                     */
    double D=3.;                    /* Number of dimensions.               */
    double derivative_meanX;        /* Mean of derivatives along X axis.   */
    double derivative_meanY;        /* Mean of derivatives along Y axis.   */
    double derivative_meanZ;        /* Mean of derivatives along Z axis.   */
    double residual_derivative;     /* Difference from mean derivative.    */
    double derivative_varX;         /* Variance of derivatives along X.    */
    double derivative_varY;         /* Variance of derivatives along Y.    */
    double derivative_varZ;         /* Variance of derivatives along Z.    */
    double probability_threshold;   /* Probability threshold for           */
                                    /* significance.                       */
    double W;                       /* W variable in Friston/Worsley's     */
                                    /* paper.                              */
    double Em;                      /* Expectation or mean of m.           */
    double dz;                      /* Partial of z, used in               */
                                    /* calculating EN.                     */
    double sum_term;                /* Summation term in calculation of    */
                                    /* integral in formula for EN.         */
    double max_sum_term;            /* Maximum ummation term.              */
    double EN;                      /* Expectation of N.                   */
    double u;                       /* Threshold for z-scores, in          */
                                    /* double-precision floating point.    */
    double gamma;                   /* Value of gamma function, used in    */
                                    /* evaluating beta.                    */
    double beta;                    /* Derived parameter for calculating   */
                                    /* P(n=x)                              */

    int clustersize_index;          /* Index into sizes of clusters.       */
    int xval,yval,zval;             /* Coordinates of a SUPRA-threshold    */
                                    /* voxel.                              */
    int *xc, *yc, *zc;              /* List of coordinates of SUPRA-       */
                                    /* threshold voxels.                   */
    float tx,ty,tz;                 /* Talairach coordinates.              */

                                    /* USER DEFINED FWHM.                  */
    double FWHMx, FWHMy, FWHMz;     /* FWHM along 3 axes, used if user     */
                                    /* inputs them on the command line.    */
    double W2;                      /* W calculated from user-defined      */
                                    /* FWHM's.                             */
    double Em2;                     /* E{m} calculated from user-          */
                                    /* defined FWHM's.                     */
    double beta2;                   /* beta calculated from user-          */
                                    /* defined FWHM's.                     */
    double Probability;             /* Probability that n is >= k.         */
    double Probability2;            /* Probability that n is >= k,         */
                                    /* based on input FWHM's or W.         */
    float mean, var;                /* Mean and variance of floating       */
                                    /* point image.                        */

                                    /* TALAIRACH COORDINATES.              */
    char spatially_normalized;      /* Flag set to 'y' if data is spatially*/
                                    /* normalized, to 'n' if not.          */
    float AX,AY,AZ;                 /* ANALYZE coordinates corresponding   */
                                    /* to the Talairach origin.            */




    /* ************************************************************************ */
    /* CHECK COMMAND LINE ARGUMENTS                                             */

    if (argc < 10) {
        printf("blobs_nocorners: Need name of input .hdr file, input statistical\n");
        printf("                 image file, input 8-bit mask image, output cluster\n");
        printf("                 mask file, output zscore cluster file, output\n");
        printf("                 probability cluster file, output text file, zscore\n");
        printf("                 threshold for generating clusters, and probability\n");
        printf("                 threshold for significance.\n\n");
        printf("Clusters of voxels >= threshold will be found.\n");
        printf("\nYou may in addition enter THREE numbers after the threshold\n");
        printf("for output, representing the known FWHM in X, Y, and Z in the\n");
        printf("same units as the voxel dimensions in the ANALYZE header, if\n");
        printf("this information is known.\n");
        printf("\nAlternatively, after the threshold for output you may input\n");
        printf("ONE number, representing your estimate for W for this zscore\n");
        printf("image.\n");
        exit(1);
    }


    /* ************************************************************************ */
    /* OPEN OUTPUT TEXT FILE.                                                   */

    if ((ftext = fopen(argv[7],"w")) == NULL) {
        printf ("blobs_nocorners: can't open %s\n",argv[7]);
        exit(2);
    }
    printf("Output text will be written to %s\n",argv[7]);

 
    /* ************************************************************************ */
    /*  READ HEADER FILE.                                                       */

    img_info = readhdr(argv+1);

    xdim = img_info.dime.dim[1];
    ydim = img_info.dime.dim[2];
    zdim = img_info.dime.dim[3];
    num_vox = xdim*ydim*zdim;
    xvoxdim = (double) img_info.dime.pixdim[1];
    yvoxdim = (double) img_info.dime.pixdim[2];
    zvoxdim = (double) img_info.dime.pixdim[3];
    if (argc==13) {
        FWHMx = (double) atof(argv[10]);
        FWHMy = (double) atof(argv[11]);
        FWHMz = (double) atof(argv[12]);
        printf("FWHMx = %g\n",(float)FWHMx);
        printf("FWHMy = %g\n",(float)FWHMy);
        printf("FWHMz = %g\n",(float)FWHMz);
        fprintf(ftext,"FWHMx = %g\n",(float)FWHMx);
        fprintf(ftext,"FWHMy = %g\n",(float)FWHMy);
        fprintf(ftext,"FWHMz = %g\n",(float)FWHMz);

        FWHMx/=xvoxdim;
        FWHMy/=yvoxdim;
        FWHMz/=zvoxdim;
        printf("FWHMx = %g in voxel lengths.\n",(float)FWHMx);
        printf("FWHMy = %g in voxel lengths.\n",(float)FWHMy);
        printf("FWHMz = %g in voxel lengths.\n",(float)FWHMz);
        fprintf(ftext,"FWHMx = %g in voxel lengths.\n",(float)FWHMx);
        fprintf(ftext,"FWHMy = %g in voxel lengths.\n",(float)FWHMy);
        fprintf(ftext,"FWHMz = %g in voxel lengths.\n",(float)FWHMz);
    }
    vox_vol = xvoxdim * yvoxdim * zvoxdim;
    printf("Image dimensions are %d x %d x %d.\n",xdim,ydim,zdim);
    printf("Voxel size is %g x %g x %g.\n",xvoxdim,yvoxdim,zvoxdim);
    printf("Volume of one voxel = %g\n\n",vox_vol);
    fprintf(ftext,"Image dimensions are %d x %d x %d.\n",xdim,ydim,zdim);
    fprintf(ftext,"Voxel size is %g x %g x %g.\n",xvoxdim,yvoxdim,zvoxdim);
    fprintf(ftext,"Volume of one voxel = %g\n\n",vox_vol);


    /* ************************************************************************ */
    /* READ IN TALAIRACH ORIGIN FROM STANDARD INPUT.                            */

    printf("Is this data spatially normalized <y/n>? ");
    scanf("%c",&spatially_normalized);
    while ((spatially_normalized!='y') && (spatially_normalized!='n')) {
        printf("\nType a 'y' or a 'n'.\n");
        printf("Is this data spatially normalized <y/n>? ");
        scanf("%c",&spatially_normalized);
    }
    if (spatially_normalized=='y') {
        printf("Please type in the ANALYZE coordinates of the Talairach origin.\n");
        printf("The origin for \"classic\" SPM is [33,53,8].  The slice corresponding\n");
        printf("to Talairach Z=0 may be different for SPM95 output.  The origin\n");
        printf("in functional MRI images may be different for all three axes.\n");
        printf("\n  ANALYZE X-coordinate for Talairach origin? ");
        scanf("%f",&AX);
        printf("\n  ANALYZE Y-coordinate for Talairach origin? ");
        scanf("%f",&AY);
        printf("\n  ANALYZE Z-coordinate for Talairach origin? ");
        scanf("%f",&AZ);
        printf("\nANALYZE coordinates corresponding to Talairach [0,0,0]: [%g %g %g]\n",
                                                                            AX,AY,AZ);
        fprintf(ftext,"\nANALYZE coordinates corresponding to Talairach [0,0,0]: [%g %g %g]\n",
                                                                            AX,AY,AZ);
    }

    /* ************************************************************************ */
    /* READ IN Z-SCORE IMAGE.                                                   */

    image = (float *) malloc(sizeof(float)*num_vox);
    if (image==NULL) {
        printf("blobs_nocorners: unable to malloc for image.\n");
        exit(3);
    }
    printf("\nLoading z-score image %s...\n",argv[2]);
    fprintf(ftext,"\nLoading z-score image %s...\n",argv[2]);
    if ((fIN = fopen(argv[2],"r")) == NULL) {
        printf ("blobs_nocorners: can't open %s\n",argv[2]);
        exit(4);
    }
    num_read_in=fread(image,sizeof(float),num_vox,fIN);
    fclose(fIN);
    if(num_read_in!=num_vox) {
        printf("blobs_nocorners: %d voxels read in; should have been %d.\n",
                                                num_read_in,num_vox);
        exit(5);
    }


    /* ************************************************************************ */
    /* READ IN mask, COUNT NUMBER OF VOXELS INCLUDED IN MASK.                   */

    mask      = (unsigned char *) malloc(sizeof(unsigned char)*num_vox);
    if (mask==NULL) {
        printf("blobs_nocorners: unable to malloc for mask.\n");
        exit(6);
    }
    printf("Loading mask image %s...\n",argv[3]);
    fprintf(ftext,"Loading %s...\n",argv[3]);
    if ((fIN = fopen(argv[3],"r")) == NULL) {
        printf ("blobs_nocorners: can't open %s\n",argv[3]);
        exit(7);
    }
    num_read_in=fread(mask,sizeof(unsigned char),num_vox,fIN);
    fclose(fIN);
    if(num_read_in!=num_vox) {
        printf("blobs_nocorners: %d voxels read in; should have been %d.\n",
                                                num_read_in,num_vox);
        exit(8);
    }
    num_vox_in_mask=0;
    for (vox_index=0; vox_index<num_vox; vox_index++)
        if (mask[vox_index]==1) num_vox_in_mask++;
    printf("\nNumber of voxels in mask = %d\n",num_vox_in_mask);
    fprintf(ftext,"\nNumber of voxels in mask = %d\n",num_vox_in_mask);


    /* ************************************************************************ */
    /* MASK Z-SCORE IMAGE.                                                      */;

    printf("Masking %s with %s...\n",argv[2],argv[3]);
    fprintf(ftext,"Masking %s with %s...\n",argv[2],argv[3]);
    for (vox_index=0; vox_index<num_vox; vox_index++)
        if (mask[vox_index]!=1) image[vox_index]=0.;


    /* ************************************************************************ */
    /* OPEN OUTPUT IMAGE FILE.                                                  */



    /* ************************************************************************ */
    /* ALLOCATE MEMORY FOR clustermask, clustermask2, AND outimg.               */

    clustermask = (int *) malloc(sizeof(int)*num_vox);
    if (clustermask==NULL) {
        printf("blobs_nocorners: unable to malloc for clustermask.\n");
        exit(10);
    }
    clustermask2 = (int *) malloc(sizeof(int)*num_vox);
    if (clustermask2==NULL) {
        printf("blobs_nocorners: unable to malloc for clustermask2.\n");
        exit(11);
    }
    intensitymask = (int *) malloc(sizeof(int)*num_vox);
    if (intensitymask==NULL) {
        printf("blobs_nocorners: unable to malloc for intensitymask.\n");
        exit(12);
    }
    probmask = (int *) malloc(sizeof(int)*num_vox);
    if (probmask==NULL) {
        printf("blobs_nocorners: unable to malloc for probmask.\n");
        exit(13);
    }
    outimg = (unsigned short *) malloc(sizeof(unsigned short)*num_vox);
    if (outimg==NULL) {
        printf("blobs_nocorners: unable to malloc for outimg.\n");
        exit(14);
    }
    clusterintensity = (float *) malloc(sizeof(float)*num_vox);
    if (clusterintensity==NULL) {
        printf("blobs_nocorners: unable to malloc for clusterintensity.\n");
        exit(15);
    }
    clusterprob = (float *) malloc(sizeof(float)*num_vox);
    if (clusterprob==NULL) {
        printf("blobs_nocorners: unable to malloc for clusterprob.\n");
        exit(16);
    }
    num_size = (int *) malloc(num_vox*sizeof(int));
    if (num_size==NULL) {
        printf("blobs_nocorners: unable to malloc for num_size.\n");
        exit(17);
    }


    /* ************************************************************************ */
    /* ZERO OUT outimg, clustermask2, clusterintensity, clusterprob,            */
    /* intensitymask, AND probmask.                                             */

    for (vox_index=0; vox_index<num_vox; vox_index++) {
        outimg[vox_index]=0;
        clustermask2[vox_index]=0;

        clusterintensity[vox_index]=0.;
        intensitymask[vox_index]=0;

        clusterprob[vox_index]=0.;
        probmask[vox_index]=0;
    }


    /* ************************************************************************ */
    /* ACCEPT ZSCORE THRESHOLD AND PROBABILITY THRESHOLD AS COMMAND LINE        */
    /* ARGUMENTS.                                                               */

    threshold = atof(argv[8]);
    u = (double) threshold;
    probability_threshold = (double)atof(argv[9]);


    /* ************************************************************************ */
    /* CALCULATE MEAN OF INPUT IMAGE.                                           */
    mean = 0.;
    for (vox_index=0; vox_index<num_vox; vox_index++)
        if (mask[vox_index]!=0)
            mean += image[vox_index];

    mean /= (float)num_vox_in_mask;


    /* ************************************************************************ */
    /* CALCULATE VARIANCE OF INPUT IMAGE.                                       */
    var = 0.;
    for (vox_index=0; vox_index<num_vox; vox_index++)
        if (mask[vox_index]!=0)
            var += ((image[vox_index]-mean)*(image[vox_index]-mean));

    var /= (float) (num_vox_in_mask-1);


    /* ************************************************************************ */
    /* COUNT NUMBER OF VOXELS GREATER THAN/EQUAL TO threshold.                  */

    num_suprathresh_vox = 0;
    for (vox_index=0; vox_index<num_vox; vox_index++)
        if (image[vox_index] <= threshold) num_suprathresh_vox++;

    xc = (int *) malloc(num_suprathresh_vox*sizeof(int));
    yc = (int *) malloc(num_suprathresh_vox*sizeof(int));
    zc = (int *) malloc(num_suprathresh_vox*sizeof(int));

    num_suprathresh_vox=0;
    for (zval = 0; zval <= zdim-1; zval++)
        for (yval = 0; yval <= ydim-1; yval++)
            for (xval = 0; xval <= xdim-1; xval++) {
                if (image[(zval*ydim+yval)*xdim+xval] >= threshold) {
                    xc[num_suprathresh_vox]=xval;
                    yc[num_suprathresh_vox]=yval;
                    zc[num_suprathresh_vox]=zval;
                    num_suprathresh_vox++;
                }
            }
    printf("There are %d voxels greater than/equal to Threshold = %g.\n\n",
                                                num_suprathresh_vox,threshold);
    fprintf(ftext,"There are %d voxels greater than/equal to Threshold = %g.\n\n",
                                                num_suprathresh_vox,threshold);


    /* ************************************************************************ */
    /* CALCULATE CLUSTER PROBABILITIES AS PROPOSED BY FRISTON/WORSLEY.          */

        /* ******************************************************************** */
        /* CALCULATE MEAN OF DERIVATIVES.                                       */

            /* **************************************************************** */
            /* ALONG X-AXIS, mean{SPMx}.                                        */

            total_num_derivatives = 0;
            derivative_sum = 0.;
            for (z=0; z<zdim; z++)
                for (y=0; y<ydim; y++)
                    for(x=0; x<xdim-1; x++)
                        if ((mask[(z*ydim+y)*xdim+x]==1)
                                && (mask[(z*ydim+y)*xdim+x+1]==1)) {
                            derivative_sum += (double)
                                ((image[(z*ydim+y)*xdim+x+1]
                                        -
                                  image[(z*ydim+y)*xdim+x]));
                            total_num_derivatives++;
                        }

            derivative_meanX = derivative_sum / (double) total_num_derivatives;


            /* **************************************************************** */
            /* ALONG Y-AXIS, mean{SPMy}.                                        */

            total_num_derivatives = 0;
            derivative_sum = 0.;
            for (z=0; z<zdim; z++)
                for (x=0; x<xdim; x++)
                    for(y=0; y<ydim-1; y++)
                        if ((mask[(z*ydim+y)*xdim+x]==1)
                                && (mask[(z*ydim+y+1)*xdim+x]==1)) {
                            derivative_sum += (double)
                                ((image[(z*ydim+y+1)*xdim+x]
                                        -
                                  image[(z*ydim+y)*xdim+x]));
                            total_num_derivatives++;
                        }

            derivative_meanY = derivative_sum / (double) total_num_derivatives;


            /* **************************************************************** */
            /* ALONG Z-AXIS, mean{SPMz}.                                        */

            total_num_derivatives = 0;
            derivative_sum = 0.;
            for (y=0; y<ydim; y++)
                for (x=0; x<xdim; x++)
                    for(z=0; z<zdim-1; z++)
                        if ((mask[(z*ydim+y)*xdim+x]==1)
                                && (mask[((z+1)*ydim+y)*xdim+x]==1)) {
                            derivative_sum += (double)
                                ((image[((z+1)*ydim+y)*xdim+x]
                                        -
                                  image[(z*ydim+y)*xdim+x]));
                            total_num_derivatives++;
                        }

            derivative_meanZ = derivative_sum / (double) total_num_derivatives;


        /* ******************************************************************** */
        /* CALCULATE VARIANCE OF DERIVATIVES.                                   */

            /* **************************************************************** */
            /* ALONG X-AXIS, var{SPMx}.                                         */

            derivative_varX = 0.;
            total_num_derivatives = 0;
            for (z=0; z<zdim; z++)
                for (y=0; y<ydim; y++)
                    for(x=0; x<xdim-1; x++)
                        if ((mask[(z*ydim+y)*xdim+x]==1)
                                && (mask[(z*ydim+y)*xdim+x+1]==1)) {
                            residual_derivative = (double)
                                (image[(z*ydim+y)*xdim+x+1]
                                        -
                                 image[(z*ydim+y)*xdim+x])
                                        -
                                 derivative_meanX;
                            derivative_varX +=
                                (residual_derivative
                                        *
                                 residual_derivative);
                            total_num_derivatives++;
                        }

            derivative_varX /= (double)(total_num_derivatives-1);


            /* **************************************************************** */
            /* ALONG Y-AXIS, var{SPMy}.                                         */

            derivative_varY = 0.;
            total_num_derivatives = 0;
            for (z=0; z<zdim; z++)
                for (x=0; x<xdim; x++)
                    for(y=0; y<ydim-1; y++)
                        if ((mask[(z*ydim+y)*xdim+x]==1)
                                && (mask[(z*ydim+y+1)*xdim+x]==1)) {
                            residual_derivative = (double)
                                (image[(z*ydim+y+1)*xdim+x]
                                        -
                                 image[(z*ydim+y)*xdim+x])
                                        -
                                 derivative_meanY;
                            derivative_varY +=
                                (residual_derivative
                                        *
                                 residual_derivative);
                            total_num_derivatives++;
                        }

            derivative_varY /= (double)(total_num_derivatives-1);


            /* **************************************************************** */
            /* ALONG Z-AXIS, var{SPMz}                                          */

            derivative_varZ = 0.;
            total_num_derivatives = 0;
            for (y=0; y<ydim; y++)
                for (x=0; x<xdim; x++)
                    for(z=0; z<zdim-1; z++)
                        if ((mask[(z*ydim+y)*xdim+x]==1)
                                && (mask[((z+1)*ydim+y)*xdim+x]==1)) {
                            residual_derivative = (double)
                                (image[((z+1)*ydim+y)*xdim+x]
                                        -
                                 image[(z*ydim+y)*xdim+x])
                                        -
                                 derivative_meanZ;
                            derivative_varZ +=
                                (residual_derivative
                                        *
                                 residual_derivative);
                            total_num_derivatives++;
                        }

            derivative_varZ /= (double)(total_num_derivatives-1);


        /* ******************************************************************** */
        /* CALCULATE S, TOTAL VOLUME OF D-DIMENSIONAL SPM.  HERE, D = 3.        */
        /* S = (NUMBER OF VOXELS COVERING BRAIN) * (VOLUME PER VOXEL)           */
        /* S = num_vox_in_mask * vox_vol; */

        S = (double) num_vox_in_mask;


        /* ******************************************************************** */
        /* CALCULATE E{N}.                                                      */
        /* Equation 3.                                1                         */
        /*         E{N} = S *  INTEGRAL(z=u,Inf) [--------- * exp[(-z^2)/2] dz  */ 
        /*                                        sqrt[2 Pi]                    */

        EN = 0.;
        pi = 2.*acos(0.);
        sum_term = exp(-threshold*threshold/2.)*0.0001;
        max_sum_term = exp(-threshold*threshold/2.)*0.0001;
        for (dz=u; max_sum_term/sum_term < 1000000. ; dz+=0.0001) {
            sum_term = exp(-dz*dz/2.)*0.0001;
            EN += sum_term;
            max_sum_term = (max_sum_term>sum_term) ? max_sum_term : sum_term;
        }
        EN *= (S/sqrt(2.*pi));

        /* ALTERNATE WAY TO CALCULATE E{N}, NOT USED.
        EN = (double) S * erffc((float)u/(float)sqrt(2.))/2.; */


        /* ******************************************************************** */
        /* CALCULATE W.                                                         */
        /* Equation 6.                            (-1/(2D))                     */
        /*                  W = product[Var{SPMi'}          ]                   */

        W = pow( (derivative_varX * derivative_varY * derivative_varZ),
                                                                (-1./(2.*D)));

        if (argc==13) W2 = pow(FWHMx*FWHMy*FWHMz,1./D)/sqrt(4.*log(2.));
        if (argc==11) W2 = (double) atof(argv[10]);


        /* ******************************************************************** */
        /* CALCULATE Em.                                                        */
        /* Here use equation 4,        S*[(2*pi)^(-(D+1)/2)]*u^(D-1)            */
        /*                      E{m} = -----------------------------            */
        /*                                 W^D   *    exp((u^2)/2)              */
        /* with D set to D=3.                                                   */

        Em = S * u * u / (4.*pi*pi*W*W*W*exp(u*u/2.));

        if ((argc==13) || (argc==11)) Em2 = S * u * u / (4.*pi*pi*W2*W2*W2*exp(u*u/2.));


        /* ******************************************************************** */
        /* CALCULATE beta.                                                      */
        /* Equation 12.            Gamma(D/2+1) * E{m}    (2/D)                 */
        /*               beta =  [----------------------]                       */
        /*                                 E{N}                                 */

        gamma = exp((double)gammln((float)(D/2.+1.)));
        beta = pow((gamma*Em/EN),(2./D));

        if ((argc==13) || (argc==11)) beta2 = pow((gamma*Em2/EN),(2./D));


        /* ******************************************************************** */
        /* PRINTOUT INTERMEDIATE VARIABLES.                                     */

        printf("\nINTERMEDIATE VARIABLES:\n");
        printf("\nmean of %s = %f\n",argv[2],mean);
        printf("variance of %s is %f\n",argv[2],var);
        printf("mean{SPMx'} = %f\n",derivative_meanX);
        printf("mean{SPMy'} = %f\n",derivative_meanY);
        printf("mean{SPMz'} = %f\n",derivative_meanZ);
        printf("Var{SPMx'} = %f\n",derivative_varX);
        printf("Var{SPMy'} = %f\n",derivative_varY);
        printf("Var{SPMz'} = %f\n",derivative_varZ);
        printf("z-score threshold u = %g\n",u);
        printf("probability threshold = %f\n",probability_threshold);
        printf("D = %g dimensions\n",D);
        printf("S = %g voxels\n",S);
        printf("Pi = %f\n",pi);
        printf("E{N} = %g\n",EN);
        printf("W = %f\n",W);
        printf("effective FWHM = W * sqrt(4*ln(2)) = %f\n",W*sqrt(4.*log(2)));
        printf("FWHMx = sqrt(4*ln(2)/Var{SPMx'}) = %f VOXEL LENGTHS\n",sqrt(4.*log(2)/derivative_varX));
        printf("FWHMy = sqrt(4*ln(2)/Var{SPMy'}) = %f VOXEL LENGTHS\n",sqrt(4.*log(2)/derivative_varY));
        printf("FWHMz = sqrt(4*ln(2)/Var{SPMz'}) = %f VOXEL LENGTHS\n",sqrt(4.*log(2)/derivative_varZ));
        printf("FWHMx = sqrt(4*ln(2)/Var{SPMx'}) = %f PHYSICAL UNITS\n",xvoxdim*sqrt(4.*log(2)/derivative_varX));
        printf("FWHMy = sqrt(4*ln(2)/Var{SPMy'}) = %f PHYSICAL UNITS\n",yvoxdim*sqrt(4.*log(2)/derivative_varY));
        printf("FWHMz = sqrt(4*ln(2)/Var{SPMz'}) = %f PHYSICAL UNITS\n",zvoxdim*sqrt(4.*log(2)/derivative_varZ));
        printf("E{m} = %f\n",Em);
        printf("Gamma[%f] = %g\n",D/2.+1.,gamma);
        printf("beta = %f\n\n\n",beta);

        fprintf(ftext,"\nINTERMEDIATE VARIABLES:\n");
        fprintf(ftext,"\nmean of %s = %f\n",argv[2],mean);
        fprintf(ftext,"variance of %s is %f\n",argv[2],var);
        fprintf(ftext,"mean{SPMx'} = %f\n",derivative_meanX);
        fprintf(ftext,"mean{SPMy'} = %f\n",derivative_meanY);
        fprintf(ftext,"mean{SPMz'} = %f\n",derivative_meanZ);
        fprintf(ftext,"Var{SPMx'} = %f\n",derivative_varX);
        fprintf(ftext,"Var{SPMy'} = %f\n",derivative_varY);
        fprintf(ftext,"Var{SPMz'} = %f\n",derivative_varZ);
        fprintf(ftext,"z-score threshold u = %g\n",u);
        fprintf(ftext,"probability threshold = %f\n",probability_threshold);
        fprintf(ftext,"D = %g dimensions\n",D);
        fprintf(ftext,"S = %g voxels\n",S);
        fprintf(ftext,"Pi = %f\n",pi);
        fprintf(ftext,"E{N} = %g\n",EN);
        fprintf(ftext,"W = %f\n",W);
        fprintf(ftext,"effective FWHM = W * sqrt(4*ln(2)) = %f\n",W*sqrt(4.*log(2)));
        fprintf(ftext,"FWHMx = sqrt(4*ln(2)/Var{SPMx'}) = %f VOXEL LENGTHS\n",sqrt(4.*log(2)/derivative_varX));
        fprintf(ftext,"FWHMy = sqrt(4*ln(2)/Var{SPMy'}) = %f VOXEL LENGTHS\n",sqrt(4.*log(2)/derivative_varY));
        fprintf(ftext,"FWHMz = sqrt(4*ln(2)/Var{SPMz'}) = %f VOXEL LENGTHS\n",sqrt(4.*log(2)/derivative_varZ));
        fprintf(ftext,"FWHMx = sqrt(4*ln(2)/Var{SPMx'}) = %f PHYSICAL UNITS\n",xvoxdim*sqrt(4.*log(2)/derivative_varX));
        fprintf(ftext,"FWHMy = sqrt(4*ln(2)/Var{SPMy'}) = %f PHYSICAL UNITS\n",yvoxdim*sqrt(4.*log(2)/derivative_varY));
        fprintf(ftext,"FWHMz = sqrt(4*ln(2)/Var{SPMz'}) = %f PHYSICAL UNITS\n",zvoxdim*sqrt(4.*log(2)/derivative_varZ));
        fprintf(ftext,"E{m} = %f\n",Em);
        fprintf(ftext,"Gamma[%f] = %g\n",D/2.+1.,gamma);
        fprintf(ftext,"beta = %f\n\n\n",beta);

        if (argc==13) {
            printf("CALCULATED intermediate variables, based on FWHM of (%g,%g,%g)\n",FWHMx,FWHMy,FWHMz);
            printf("in VOXEL LENGTHS in the three spatial dimensions:\n\n");
            printf("W = %f\n",W2);
            printf("E{m} = %f\n",Em2);
            printf("beta = %f\n\n\n",beta2);

            fprintf(ftext,"CALCULATED intermediate variables, based on FWHM of (%g,%g,%g)\n",FWHMx,FWHMy,FWHMz);
            fprintf(ftext,"in VOXEL LENGTHS in the three spatial dimensions:\n\n");
            fprintf(ftext,"W = %f\n",W2);
            fprintf(ftext,"E{m} = %f\n",Em2);
            fprintf(ftext,"beta = %f\n\n\n",beta2);

        }
        if (argc==11) {
            printf("CALCULATED intermediate variables, based on input value of W = %f :\n\n",W2);
            printf("E{m} = %f\n",Em2);
            printf("beta = %f\n\n\n",beta2);

            fprintf(ftext,"CALCULATED intermediate variables, based on input value of W = %f :\n\n",W2);
            fprintf(ftext,"E{m} = %f\n",Em2);
            fprintf(ftext,"beta = %f\n\n\n",beta2);
        }


    /* ************************************************************************ */
    /* WRITE HEADER.                                                            */

    printf("TABLE OF CLUSTERS:\n");
    printf("Cluster#    Size  Mean Value  meanX   meanY   meanZ");
    fprintf(ftext,"TABLE OF CLUSTERS:\n");
    fprintf(ftext,"Cluster#    Size  Mean Value  meanX   meanY   meanZ");
    if (spatially_normalized=='y') {
        printf("   Tal_x  Tal_y  Tal_z");
        fprintf(ftext,"   Tal_x  Tal_y  Tal_z");
    }
    printf("  P(nmax>=k)");
    fprintf(ftext,"  P(nmax>=k)");
    if ((argc==13) || (argc==11)) {
        printf("   P'(nmax>=k)");
        fprintf(ftext,"   P'(nmax>=k)");
    }
    printf("\n");
    fprintf(ftext,"\n");


    /* ************************************************************************ */
    /* INITIALIZE num_size[MAX_SIZE_CLUSTERS].                                  */

    for (clustersize_index=0; clustersize_index<num_vox; clustersize_index++)
        num_size[clustersize_index]=0;

    /* ************************************************************************ */
    /* CLEAR clustermask, intensity mask, & probmask.                           */

    for (vox_index = 0; vox_index<num_vox; vox_index++) {
        clustermask[vox_index] = 0;
        intensitymask[vox_index]=0.;
        probmask[vox_index]=0.;
    }

    num_cluster=0;
    for (vox_index=0; vox_index<num_suprathresh_vox; vox_index++) {
        xval=xc[vox_index];
        yval=yc[vox_index];
        zval=zc[vox_index];

        if ((image[(zval*ydim+yval)*xdim+xval]>=threshold) 
                && (clustermask[(zval*ydim+yval)*xdim+xval]!=1)) {

            /* **************************************************************** */
            /* SINCE DATA HAS BEEN PRE-PROCESSED, VOXELS IN image ARE KNOWN TO  */
            /* BE LOCAL MAXIMA.  PRINT THE VALUE & COORDINATES OF CURRENT       */
            /* MAXIMUM, AND PRINT COORDINATES OF CONTIGUOUS VOXELS OF THE SAME  */
            /* VALUE.                                                           */

            num_cluster++;


            /* **************************************************************** */
            /* "GROW" REGION AROUND "SEED" POINT (xval,yval,zval); RETURN 1's   */
            /* WHERE THE CORRESPONDING VOXEL IN image IS LESS THAN THRESHOLD    */
            /* AND IS CONTIGUOUS WITH THE VOXEL JUST OUTPUTED.                  */

            clustersize=0;
            xmean=0.;
            ymean=0.;
            zmean=0.;
            meanval=0.;
            growcluster(xval,yval,zval);
            Probability  = PnmaxGEk( beta ,(double)clustersize,D,Em);
            Probability2 = PnmaxGEk2(beta2,(double)clustersize,D,Em2);

            /* **************************************************************** */
            /* FILL CLUSTER MASK IMAGE.                                         */
            if ((argc==10) && (Probability <= probability_threshold))
                                                filloutimg(xval,yval,zval);
            if (((argc==13) || (argc==11)) && (Probability2 <= probability_threshold))
                                                filloutimg(xval,yval,zval);

            /* **************************************************************** */
            /* FILL PROBABILITY MASK IMAGE.                                     */
            if ((argc==10) && (Probability <= probability_threshold)) {
                if (Probability==0) Probability = pow(10.,-39.);
                fillprobimg(xval,yval,zval,(float)-log10(Probability));
            }
            if (((argc==13) || (argc==11)) && (Probability2 <= probability_threshold)) {
                if (Probability2==0) Probability2 = pow(10.,-39.);
                fillprobimg(xval,yval,zval,(float)-log10(Probability2));
            }

            /* **************************************************************** */
            /* FILL INTENSITY MASK IMAGE.                                       */
            if ((argc==10) && (Probability <= probability_threshold))
                                fillintensityimg(xval,yval,zval);
            if (((argc==13) || (argc==11)) && (Probability2 <= probability_threshold))
                                fillintensityimg(xval,yval,zval);

            /* **************************************************************** */
            /* CALCULATE MEAN XYZ COORDINATES.                                  */

            xmean/=(double)clustersize;
            ymean/=(double)clustersize;
            zmean/=(double)clustersize;
            meanval/=(double)clustersize;


            /* **************************************************************** */
            /* OUTPUT CLUSTERS FOUND.                                           */

            if (spatially_normalized=='y') {
                tx = ((xmean+1)-AX)*xvoxdim;
                ty = ((ymean+1)-AY)*yvoxdim;
                tz = ((zmean+1)-AZ)*zvoxdim;
            }

            printf("    %3d      %3d  %9g  %6.2f  %6.2f  %6.2f",
                        num_cluster,clustersize,meanval,xmean+1,ymean+1,zmean+1);
            fprintf(ftext,"    %3d      %3d  %9g  %6.2f  %6.2f  %6.2f",
                        num_cluster,clustersize,meanval,xmean+1,ymean+1,zmean+1);
            if (spatially_normalized=='y') {
                printf("  %+6.2f %+6.2f %+6.2f",tx,ty,tz);
                fprintf(ftext,"  %+6.2f %+6.2f %+6.2f",tx,ty,tz);
            }
            printf("  %11g",Probability);
            fprintf(ftext,"  %11g",Probability);

            if ((argc==13) || (argc==11)) {
                printf("  %11g",PnmaxGEk2(beta2,(double)clustersize,D,Em2));
                fprintf(ftext,"  %11g",Probability2);
            }

            printf("\n");
            fprintf(ftext,"\n");
                

            /* **************************************************************** */
            /* UPDATE num_size.                                                 */

            if ((argc==10) && (Probability <= probability_threshold))
                num_size[clustersize-1]++;
            if (((argc==13) || (argc==11)) && (Probability2 <= probability_threshold))
                num_size[clustersize-1]++;

        } /* END if */

    } /* END for LOOP */

    printf("\n\nSUMMARY TABLE:\n");
    printf("Size  Number of Clusters\n");
    fprintf(ftext,"\n\nSUMMARY TABLE:\n");
    fprintf(ftext,"Size  Number of Clusters\n");
    for (clustersize_index=0; clustersize_index<num_vox; clustersize_index++)
        if (num_size[clustersize_index]>0) {
            printf("%3d    %3d\n",clustersize_index+1,num_size[clustersize_index]);
            fprintf(ftext,
                   "%3d    %3d\n",clustersize_index+1,num_size[clustersize_index]);
        }


    /* ************************************************************************ */
    /* OUTPUT LOCAL MAXIMA TO TEXT FILE.                                        */

    printf("\n\nLOCAL MAXIMA:\n");
    printf("    Value    x   y   z");

    if (spatially_normalized=='y') {
        printf("    Tal_x    Tal_y    Tal_z");
        fprintf(ftext,"    Tal_x    Tal_y    Tal_z");
    }
    fprintf(ftext,"\n\nLOCAL MAXIMA:\n");
    fprintf(ftext,"    Value    x   y   z");

    printf("  Blob #\n");
    fprintf(ftext,"  Blob #\n");

    for (vox_index=0; vox_index<num_suprathresh_vox; vox_index++) {
        x = xc[vox_index];
        y = yc[vox_index];
        z = zc[vox_index];
        not_max = 0;
        for (xx = -1; xx <= 1; xx++)
            for (yy = -1; yy <= 1; yy++)
                for (zz = -1; zz <= 1; zz++)
                    if ((x+xx<=xdim-1) && (x+xx>=0)
                     && (y+yy<=ydim-1) && (y+yy>=0)
                     && (z+zz<=zdim-1) && (z+zz>=0)
                     && ((xx!=0) || (yy!=0) || (zz!=0))
                     && ((xx==0) || (yy==0) || (zz==0)))
                        if(image[((z+zz)*ydim+y+yy)*xdim+x+xx]
                                                > image[((z*ydim+y)*xdim)+x])
                            not_max = 1;
        if ((not_max!=1) && (outimg[((z*ydim+y)*xdim)+x]!=0)) {
            printf("%9.2f  %3d %3d %3d",image[((z*ydim+y)*xdim)+x],x+1,y+1,z+1);
            fprintf(ftext,"%9.2f  %3d %3d %3d",image[((z*ydim+y)*xdim)+x],x+1,y+1,z+1);
            if (spatially_normalized=='y') {
                tx = ((x+1)-AX)*xvoxdim;
                ty = ((y+1)-AY)*yvoxdim;
                tz = ((z+1)-AZ)*zvoxdim;
                printf("   %+6.2f   %+6.2f   %+6.2f",tx,ty,tz);
                fprintf(ftext,"   %+6.2f   %+6.2f   %+6.2f",tx,ty,tz);
            }
            printf("    %d\n",outimg[((z*ydim+y)*xdim)+x]);
            fprintf(ftext,"    %d\n",outimg[((z*ydim+y)*xdim)+x]);
        }
    }


    /* ************************************************************************ */
    /* WRITE OUTPUT CLUSTERMASK IMAGE.                                          */

    if ((fOUT = fopen(argv[4],"w")) == NULL) {
        printf ("blobs_nocorners: can't open %s\n",argv[4]);
        exit(9);
    }
    printf("\n\nWriting output image %s...\n",argv[4]);
    fprintf(ftext,"\n\nWriting output image %s...\n",argv[4]);
    num_wrote_out=fwrite(outimg,sizeof(unsigned short),num_vox,fOUT);
    if(num_wrote_out!=num_vox) {
        printf("blobs_nocorners: %d voxels wrote out; should have been %d.\n",
                                                        num_wrote_out,num_vox);

        exit(18);
    }
    fclose(fOUT);

    /* ************************************************************************ */
    /* WRITE OUTPUT INTENSITYMASK IMAGE.                                        */

    if ((fOUT = fopen(argv[5],"w")) == NULL) {
        printf ("blobs_nocorners: can't open %s\n",argv[5]);
        exit(9);
    }
    printf("\nWriting output image %s...\n",argv[5]);
    fprintf(ftext,"\n\nWriting output image %s...\n",argv[5]);
    num_wrote_out=fwrite(clusterintensity,sizeof(float),num_vox,fOUT);
    if(num_wrote_out!=num_vox) {
        printf("blobs_nocorners: %d voxels wrote out; should have been %d.\n",
                                                        num_wrote_out,num_vox);

        exit(19);
    }
    fclose(fOUT);

    /* ************************************************************************ */
    /* WRITE OUTPUT PROBABILITYMASK IMAGE.                                      */

    if ((fOUT = fopen(argv[6],"w")) == NULL) {
        printf ("blobs_nocorners: can't open %s\n",argv[6]);
        exit(9);
    }
    printf("\nWriting output image %s...\n",argv[6]);
    fprintf(ftext,"\n\nWriting output image %s...\n",argv[6]);
    num_wrote_out=fwrite(clusterprob,sizeof(float),num_vox,fOUT);
    if(num_wrote_out!=num_vox) {
        printf("blobs_nocorners: %d voxels wrote out; should have been %d.\n",
                                                        num_wrote_out,num_vox);

        exit(20);
    }
    fclose(fOUT);


    /* ************************************************************************ */
    /* CLOSE OUTPUT FILE, FREE MEMORY, EXIT PROGRAM.                            */
    fclose(ftext);
    printf("blobs_nocorners: releasing memory...\n");
    free(image);
    free(outimg);
    free(clustermask);
    free(clustermask2);
    free(clusterintensity);
    free(clusterprob);
    free(intensitymask);
    free(probmask);
    free(num_size);
    free(xc);
    free(yc);
    free(zc);

    return(0);
}

/* **************************************************************************** */
void growcluster(int x, int y, int z)
/* RECURSIVE FUNCTION TO "GROW" CLUSTERS.                                       */
/* ACCEPTS AS INPUT XYZ COORDINATES OF VOXEL IN *image AROUND WHICH             */
/* CONTIGUOUS VOXELS ABOVE THRESHOLD ARE TO BE FOUND.                           */
/* FINDS ALL SUPRA-THRESHOLD PIXELS IN image CONTIGUOUS WITH PIXEL(s) MASKED    */
/* OUT IN mask.  RETURNS "1" IN CORRESPONDING ENTRIES IN *clustermask.          */
/*                                                                              */
/* 11.07.95: JMM : MODIFIED TO NOT INCLUDE ONLY-SHARED-CORNERS CASES AS         */
/*                   CONTIGUOUS.                                                */
{
    int xx,yy,zz;
    int vox_index;

    vox_index = (z*ydim+y)*xdim+x;
    if ((image[vox_index] >= threshold) && (clustermask[vox_index] == 0)) {
        clustermask[vox_index]=1;
        clustersize++;
        xmean+=x;
        ymean+=y;
        zmean+=z;
        meanval+=image[vox_index];
        for (xx = -1; xx <= 1; xx++)
            for (yy = -1; yy <= 1; yy++)
                for (zz = -1; zz <= 1; zz++)
                    if ((x+xx<=xdim-1) && (x+xx>=0)
                     && (y+yy<=ydim-1) && (y+yy>=0)
                     && (z+zz<=zdim-1) && (z+zz>=0)
                     && ((xx!=0) || (yy!=0) || (zz!=0))
                     && ((xx==0) || (yy==0) || (zz==0)))
                        growcluster(x+xx,y+yy,z+zz);
    }
}

/* **************************************************************************** */
void filloutimg(int x, int y, int z)
/* RECURSIVE FUNCTION TO FILL IN OUTIMG WITH CLUSTERS.                          */
/*                                                                              */
/* 11.07.95: JMM : MODIFIED TO NOT INCLUDE ONLY-SHARED-CORNERS CASES AS         */
/*                   CONTIGUOUS.                                                */
{
    int xx,yy,zz;
    int vox_index;

    vox_index = (z*ydim+y)*xdim+x;
    if ((clustermask[vox_index] == 1) && (clustermask2[vox_index]==0)) {
        clustermask2[vox_index]=1;
        outimg[vox_index]=num_cluster;
        for (xx = -1; xx <= 1; xx++)
            for (yy = -1; yy <= 1; yy++)
                for (zz = -1; zz <= 1; zz++)
                    if ((x+xx<=xdim-1) && (x+xx>=0)
                     && (y+yy<=ydim-1) && (y+yy>=0)
                     && (z+zz<=zdim-1) && (z+zz>=0)
                     && ((xx!=0) || (yy!=0) || (zz!=0))
                     && ((xx==0) || (yy==0) || (zz==0)))
                        filloutimg(x+xx,y+yy,z+zz);
    }
}

/* **************************************************************************** */
void fillintensityimg(int x, int y, int z)
/* RECURSIVE FUNCTION TO FILL IN clusterintensity WITH CLUSTERS.                */
/*                                                                              */
/* 11.07.95: JMM : MODIFIED TO NOT INCLUDE ONLY-SHARED-CORNERS CASES AS         */
/*                   CONTIGUOUS.                                                */
{
    int xx,yy,zz;
    int vox_index;

    vox_index = (z*ydim+y)*xdim+x;
    if ((clustermask[vox_index] == 1) && (intensitymask[vox_index]==0)) {
        intensitymask[vox_index]=1;
        clusterintensity[vox_index]=image[vox_index];
        for (xx = -1; xx <= 1; xx++)
            for (yy = -1; yy <= 1; yy++)
                for (zz = -1; zz <= 1; zz++)
                    if ((x+xx<=xdim-1) && (x+xx>=0)
                     && (y+yy<=ydim-1) && (y+yy>=0)
                     && (z+zz<=zdim-1) && (z+zz>=0)
                     && ((xx!=0) || (yy!=0) || (zz!=0))
                     && ((xx==0) || (yy==0) || (zz==0)))
                        fillintensityimg(x+xx,y+yy,z+zz);
    }
}

/* **************************************************************************** */
void fillprobimg(int x, int y, int z, float fillval)
/* RECURSIVE FUNCTION TO FILL IN clusterprob WITH CLUSTERS.                     */
/*                                                                              */
/* 11.07.95: JMM : MODIFIED TO NOT INCLUDE ONLY-SHARED-CORNERS CASES AS         */
/*                   CONTIGUOUS.                                                */
{
    int xx,yy,zz;
    int vox_index;

    vox_index = (z*ydim+y)*xdim+x;
    if ((clustermask[vox_index] == 1) && (probmask[vox_index]==0)) {
        probmask[vox_index]=1;
        clusterprob[vox_index]=fillval;
        for (xx = -1; xx <= 1; xx++)
            for (yy = -1; yy <= 1; yy++)
                for (zz = -1; zz <= 1; zz++)
                    if ((x+xx<=xdim-1) && (x+xx>=0)
                     && (y+yy<=ydim-1) && (y+yy>=0)
                     && (z+zz<=zdim-1) && (z+zz>=0)
                     && ((xx!=0) || (yy!=0) || (zz!=0))
                     && ((xx==0) || (yy==0) || (zz==0)))
                        fillprobimg(x+xx,y+yy,z+zz,fillval);
    }
}

/* **************************************************************************** */
double PnGEx(double beta, double x, double D)
/* Probability that n is greater than or equal to x.                            */
/* See equation 11.                                                             */
/*                 P(n>=x) = exp(-beta*x^(2/D))                                 */
{
    return exp(-beta*pow(x,2./D));
}

/* **************************************************************************** */
double PnmaxGEk(double beta, double k, double D, double Em)
/* Probability that nmax is greater than or equal to k.                         */
/* See equation 14.                                                             */
/*                 P(nmax>=k) = 1 - exp(-E{m}*P(n>=k))                          */
{
    return (1.-exp(-Em*PnGEx(beta,k,D)));
}

/* **************************************************************************** */
double PnGEx2(double beta2, double x, double D)
/* Probability that n is greater than or equal to x.                            */
/* See equation 11.                                                             */
/*                 P(n>=x) = exp(-beta*x^(2/D))                                 */
{
    return exp(-beta2*pow(x,2./D));
}

/* **************************************************************************** */
double PnmaxGEk2(double beta2, double k, double D, double Em2)
/* Probability that nmax is greater than or equal to k.                         */
/* See equation 14.                                                             */
/*                 P(nmax>=k) = 1 - exp(-E{m}*P(n>=k))                          */
{
    return (1.-exp(-Em2*PnGEx2(beta2,k,D)));
}

