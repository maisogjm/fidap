#include <stdio.h>
#include <stdlib.h>
#include <float.h>      /* FOR FLT_MAX                                       */
#include "analyze6.h"
#include <math.h>

                        /* GLOBAL VARIABLES.                                 */
float *image;           /* Input image to be analyzed.                       */
int xdim, ydim, zdim;   /* Dimensions of image.                              */
double gamma_val;       /* Value of gamma function, used in evaluating beta. */
    double pi;          /* Constant containing value for Pi                  */
    double S;           /* Total volume of voxels in mask,                   */
                        /*    = num_vox_in_mask * vox_vol                    */
    double D=3.;        /* Number of dimensions.                             */
    double W;           /* W variable in Friston/Worsley's                   */
                        /* paper.                                            */
    double Em;          /* Expectation or mean of m.                         */
    double dz;          /* Partial of z, used in                             */
                        /* calculating EN.                                   */
    double sum_term;    /* Summation term in calculation of                  */
                        /* integral in formula for EN.                       */
    double max_sum_term;/* Maximum summation term.                           */
    double EN;          /* Expectation of N.                                 */
    double u;           /* Threshold for z-scores, in                        */
                        /* double-precision floating point.                  */
    double beta;        /* Derived parameter for calculating                 */
                        /* P(n=x)                                            */


    int k;              /* Number of voxels in a cluster.                    */

struct analyze_struct readhdr(char *argv[]);
float gammln(float xx);


/* **************************************************************************** */
main(int argc, char *argv[])

/*
Given an ANALYZE image file determines an effective smoothness.

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
05.24.95 : JMM L Added code to also output smoothness calculation as per
                 Forman SD, Cohen JD, Fitzgerald M, Eddy WF, Mintun MA,
                 Noll DC, "Improved Assessment of Significant Activation in
                 Functional Magnetic Resonance Imaging (fMRI): Use of a
                 Cluster-Size Threshold," MRM 33:636-647 (1995).

                 Also added code to allow program to read in 8-bit, 16-bit, OR
                 floating point images.


This program makes no assumption about image dimensions.
They are read in from the header file.
*/
{
    FILE *ftext;
    FILE *fIN;
    FILE *fOUT;

    struct analyze_struct img_info; /* Header information variable.     */
    int num_bits_per_pixel;         /* Bits per pixel in image.         */
    unsigned char *Bbuffer;         /* Buffer to hold 8-bit image as    */
                                    /* read in from disk.               */
    signed short int *Wbuffer;      /* Buffer to hold 16-bit image as   */
                                    /* read in from disk.               */
    float *Fbuffer;                 /* Buffer to hold floating point    */
                                    /* image as read in from disk.      */
    int num_read_in;                /* Number of voxels read in.        */
    int num_wrote_out;              /* Number of voxels written out.    */
    int *num_size;                  /* Number of clusters of a certain  */
                                    /* size.                            */
    int vox_index,                  /* Index into voxels.               */
        num_vox,                    /* Total number of voxels.          */
        cluster_index,              /* Index into clusters.             */
        num_suprathresh_vox;        /* Number of SUPRA-threshold voxels.*/
    unsigned char *mask;            /* Mask to be applied to input      */
                                    /* statistical image.               */
    int x, y, z;                    /* Indices for looping over 3       */
                                    /* spatial dimensions.              */
    int xx, yy, zz;                 /* Indices for looping over 3       */
                                    /* spatial dimensions.              */

                                    /* VARIABLES FOR DOING PARTICLE     */
                                    /* ANALYSIS.                        */
    int num_vox_in_mask;            /* Number of voxels in the mask.    */
    double xvoxdim, yvoxdim,        /* Voxel sizes along the 3          */
                        zvoxdim;    /*         dimensions.              */
    double vox_vol;                 /* Volume of one voxel.             */
    int num_derivatives;            /* Number of derivatives in one row */
                                    /* or column.                       */
    double derivative_sum;          /* Sum of derivatives.              */
    int total_num_derivatives;      /* Number of derivatives in sum     */
                                    /* along one axis.                  */
    double derivative_meanX;        /* Mean of derivatives along X axis.*/
    double derivative_meanY;        /* Mean of derivatives along Y axis.*/
    double derivative_meanZ;        /* Mean of derivatives along Z axis.*/
    double residual_derivative;     /* Difference from mean derivative. */
    double derivative_varX;         /* Variance of derivatives along X. */
    double derivative_varY;         /* Variance of derivatives along Y. */
    double derivative_varZ;         /* Variance of derivatives along Z. */

    int min;                        /* Minimum value of cost function.  */
    double mean, var;               /* Mean & Variance across all voxels*/

                                    /* ALTERNATE SMOOTHNESS CALCULATIONS*/
    double voxdim,                  /* Effective voxel dim in 3D.       */
        sx,sy,sz,                   /* SD of effective Gaussian filter  */
        FWHMx,FWHMy,FWHMz,          /* FWHM's of effective filter       */
        FWHM,                       /* Overall FWHM in 3D.              */
        W2;                         /* Effective W2.                    */


    /* ************************************************************************ */
    /* CHECK COMMAND LINE ARGUMENTS                                             */

    if (argc != 5) {
        printf("vardir2: Need name of input .hdr file, input image file to estimate\n");
        printf("         smoothness on, input 8-bit mask image, and output textfile name.\n");
        exit(1);
    }


    /* ************************************************************************ */
    /* OPEN OUTPUT TEXT FILE.                                                   */

    if ((ftext = fopen(argv[4],"w")) == NULL) {
        printf ("blobs: can't open %s\n",argv[4]);
        exit(2);
    }

    /* ************************************************************************ */
    /*  READ HEADER FILE.                                                       */

    img_info = readhdr(argv+1);

    xdim = img_info.dime.dim[1];
    ydim = img_info.dime.dim[2];
    zdim = img_info.dime.dim[3];
    num_bits_per_pixel=img_info.dime.bitpix;
    num_vox = xdim*ydim*zdim;
    xvoxdim = (double) img_info.dime.pixdim[1];
    yvoxdim = (double) img_info.dime.pixdim[2];
    zvoxdim = (double) img_info.dime.pixdim[3];
    vox_vol = xvoxdim * yvoxdim * zvoxdim;
    printf("Image dimensions are %d x %d x %d.\n",xdim,ydim,zdim);
    printf("Voxel size is %g x %g x %g.\n",xvoxdim,yvoxdim,zvoxdim);
    printf("Volume of one voxel = %g\n",vox_vol);
    printf("Bits per pixel = %d\n\n",num_bits_per_pixel);

    fprintf(ftext,"Image dimensions are %d x %d x %d.\n",xdim,ydim,zdim);
    fprintf(ftext,"Voxel size is %g x %g x %g.\n",xvoxdim,yvoxdim,zvoxdim);
    fprintf(ftext,"Volume of one voxel = %g\n",vox_vol);
    fprintf(ftext,"Bits per pixel = %d\n\n",num_bits_per_pixel);


    /* ************************************************************************ */
    /* ALLOCATE MEMORY FOR BUFFER.                                              */

    if (num_bits_per_pixel==8) {
        Bbuffer = (unsigned char *) malloc(sizeof(unsigned char)*num_vox);
        if (Bbuffer==NULL) {
            printf("convolve3D: unable to malloc for Bbuffer.\n");
            exit(11);
        }
    }
    else if (num_bits_per_pixel==16) {
        Wbuffer = (signed short int *) malloc(sizeof(signed short int)*num_vox);
        if (Wbuffer==NULL) {
            printf("convolve3D: unable to malloc for Wbuffer.\n");
            exit(12);
        }
    }
    else if (num_bits_per_pixel==32) {
        Fbuffer = (float *) malloc(sizeof(float)*num_vox);
        if (Fbuffer==NULL) {
            printf("convolve3D: unable to malloc for Fbuffer.\n");
            exit(13);
        }
    }

    /* ************************************************************************ */
    /* READ IN Z-SCORE IMAGE.                                                   */

    image = (float *) malloc(sizeof(float)*num_vox);
    if (image==NULL) {
        printf("vardir2: unable to malloc for image.\n");
        exit(3);
    }
    printf("  Loading input image %s...\n",argv[2]);
    fprintf(ftext,"  Loading input image %s...\n",argv[2]);
    if ((fIN = fopen(argv[2],"r")) == NULL) {
        printf ("vardir2: can't open %s\n",argv[2]);
        exit(4);
    }
    if (num_bits_per_pixel==8)
        num_read_in=fread(Bbuffer,sizeof(unsigned char),num_vox,fIN);
    else if (num_bits_per_pixel==16)
        num_read_in=fread(Wbuffer,sizeof(signed short int),num_vox,fIN);
    else if (num_bits_per_pixel==32)
        num_read_in=fread(Fbuffer,sizeof(float),num_vox,fIN);
    fclose(fIN);
    if(num_read_in!=num_vox) {
        printf("vardir2: %d voxels read in; should have been %d.\n",
                                                num_read_in,num_vox);
        exit(5);
    }
    for (vox_index=0; vox_index<num_vox; vox_index++)
        if (num_bits_per_pixel==8)
            image[vox_index]=(float)Bbuffer[vox_index];
        else if (num_bits_per_pixel==16)
            image[vox_index]=(float)Wbuffer[vox_index];
        else if (num_bits_per_pixel==32)
            image[vox_index]=Fbuffer[vox_index];

    if (num_bits_per_pixel==8)
        free(Bbuffer);
    else if (num_bits_per_pixel==16)
        free(Wbuffer);
    else if (num_bits_per_pixel==32)
        free(Fbuffer);

    /* ************************************************************************ */
    /* READ IN mask, COUNT NUMBER OF VOXELS INCLUDED IN MASK.                   */

    mask      = (unsigned char *) malloc(sizeof(unsigned char)*num_vox);
    if (mask==NULL) {
        printf("vardir2: unable to malloc for mask.\n");
        exit(6);
    }
    printf("  Loading mask image %s...\n",argv[3]);
    fprintf(ftext,"  Loading mask image %s...\n",argv[3]);
    if ((fIN = fopen(argv[3],"r")) == NULL) {
        printf ("vardir2: can't open %s\n",argv[3]);
        exit(7);
    }
    num_read_in=fread(mask,sizeof(unsigned char),num_vox,fIN);
    fclose(fIN);
    if(num_read_in!=num_vox) {
        printf("vardir2: %d voxels read in; should have been %d.\n",
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

    printf("masking %s with %s...\n",argv[2],argv[3]);
    fprintf(ftext,"masking %s with %s...\n",argv[2],argv[3]);
    for (vox_index=0; vox_index<num_vox; vox_index++)
        if (mask[vox_index]!=1) image[vox_index]=0.;


    /* ************************************************************************ */
    /* CALCULATE MEAN OF INPUT IMAGE.                                           */
    mean = 0.;
    for (vox_index=0; vox_index<num_vox; vox_index++)
        if (mask[vox_index]!=0)
            mean += (double) image[vox_index];

    mean /= (double)num_vox_in_mask;
    printf("\nmean of %s = %g\n",argv[2],(float)mean);
    fprintf(ftext,"\nmean of %s = %g\n",argv[2],(float)mean);


    /* ************************************************************************ */
    /* CALCULATE VARIANCE OF INPUT IMAGE.                                       */
    var = 0.;
    for (vox_index=0; vox_index<num_vox; vox_index++)
        if (mask[vox_index]!=0)
            var += (((double)image[vox_index]-mean)*((double)image[vox_index]-mean));

    var /= (double) (num_vox_in_mask-1);
    printf("variance of %s is %g\n",argv[2],(double)var);
    fprintf(ftext,"variance of %s is %g\n",argv[2],(double)var);


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
        /* CALCULATE W.                                                         */
        /* Equation 6.                            (-1/(2D))                     */
        /*                  W = product[Var{SPMi'}          ]                   */

        W = pow( (derivative_varX * derivative_varY * derivative_varZ),
                                                                (-1./(2.*D)));



        /* ******************************************************************** */
        /* CALCULATE Gamma((D/2)+1) AND PI.                                     */

        gamma_val = exp((double)gammln((float)(D/2.+1.)));
        pi = 2.*acos(0.);


        /* ******************************************************************** */
        /* ALTERNATE CALCULATION OF SMOOTHNESS.                                 */

        voxdim = pow(xvoxdim*yvoxdim*zvoxdim,1./3.);
	if ( ((1 - (derivative_varX / (2.*var)))>0. ) &&
	            ((-1./(4. * log( 1 - (derivative_varX / (2.*var)) )))>0.) )
            sx = sqrt(-1./(4. * log( 1 - (derivative_varX / (2.*var)) )));
        else
            sx = 0.;

	if ( ((1 - (derivative_varY / (2.*var)))>0. ) &&
	            ((-1./(4. * log( 1 - (derivative_varY / (2.*var)) )))>0.) )
            sy = sqrt(-1./(4. * log( 1 - (derivative_varY / (2.*var)) )));
        else
            sy = 0.;

	if ( ((1 - (derivative_varZ / (2.*var)))>0. ) &&
	            ((-1./(4. * log( 1 - (derivative_varZ / (2.*var)) )))>0.) )
            sz = sqrt(-1./(4. * log( 1 - (derivative_varZ / (2.*var)) )));
        else
            sz = 0.;

        FWHMx  = sx*sqrt(8.*log(2.));
        FWHMy  = sy*sqrt(8.*log(2.));
        FWHMz  = sz*sqrt(8.*log(2.));
        FWHM   = pow(FWHMx*FWHMy*FWHMz,1./3.);
        W2     = FWHM/sqrt(4.*log(2.));
        

        /* ******************************************************************** */
        /* PRINTOUT INTERMEDIATE VARIABLES.                                     */

        printf("\nINTERMEDIATE VARIABLES:\n");
        printf("mean{SPMx'} = %g\n",derivative_meanX);
        printf("mean{SPMy'} = %g\n",derivative_meanY);
        printf("mean{SPMz'} = %g\n",derivative_meanZ);
        printf("Var{SPMx'}  = %g\n",derivative_varX);
        printf("Var{SPMy'}  = %g\n",derivative_varY);
        printf("Var{SPMz'}  = %g\n",derivative_varZ);
        printf("Geometric mean of Var{SPMi'}                  = %g\n",
                pow(derivative_varX*derivative_varY*derivative_varZ,1./3.));
        printf("Geometric mean of Var{SPMi'} / Image Variance = %g\n",
                pow(derivative_varX*derivative_varY*derivative_varZ,1./3.)/var);
        printf("D  = %g dimensions\n",D);
        printf("S  = %g voxels\n\n",S);
        /* printf("Pi = %g\n\n",pi); */

        printf("*********************************************************************\n");
	printf("ESTIMATE OF SPATIAL SMOOTHNESS ASSUMING UNIT VARIANCE:\n");
        printf("W  = %g\n",W);
        printf("effective FWHM = W * sqrt(4*ln(2)) = %g VOXEL LENGTHS\n",W*sqrt(4.*log(2)));
        printf("FWHMx = sqrt(4*ln(2)/Var{SPMx'})   = %g VOXEL LENGTHS\n",sqrt(4.*log(2)/derivative_varX));
        printf("FWHMy = sqrt(4*ln(2)/Var{SPMy'})   = %g VOXEL LENGTHS\n",sqrt(4.*log(2)/derivative_varY));
        printf("FWHMz = sqrt(4*ln(2)/Var{SPMz'})   = %g VOXEL LENGTHS\n",sqrt(4.*log(2)/derivative_varZ));
        printf("FWHMx = sqrt(4*ln(2)/Var{SPMx'})   = %g PHYSICAL UNITS\n",xvoxdim*sqrt(4.*log(2)/derivative_varX));
        printf("FWHMy = sqrt(4*ln(2)/Var{SPMy'})   = %g PHYSICAL UNITS\n",yvoxdim*sqrt(4.*log(2)/derivative_varY));
        printf("FWHMz = sqrt(4*ln(2)/Var{SPMz'})   = %g PHYSICAL UNITS\n",zvoxdim*sqrt(4.*log(2)/derivative_varZ));
        printf("Gamma[%g] = %g\n\n",D/2.+1.,gamma_val);

        printf("*********************************************************************\n");
        printf("SMOOTHNESS WITH CORRECTION FOR NON-UNIT VARIANCE:\n");
        printf("W  = %g\n",W*sqrt(var));
        printf("effective FWHM = W * sqrt(4*ln(2)) = %g VOXEL LENGTHS\n",W*sqrt(var*4.*log(2)));
        printf("FWHMx = sqrt(4*ln(2)/Var{SPMx'})   = %g VOXEL LENGTHS\n",sqrt(var*4.*log(2)/derivative_varX));
        printf("FWHMy = sqrt(4*ln(2)/Var{SPMy'})   = %g VOXEL LENGTHS\n",sqrt(var*4.*log(2)/derivative_varY));
        printf("FWHMz = sqrt(4*ln(2)/Var{SPMz'})   = %g VOXEL LENGTHS\n",sqrt(var*4.*log(2)/derivative_varZ));
        printf("FWHMx = sqrt(4*ln(2)/Var{SPMx'})   = %g PHYSICAL UNITS\n",xvoxdim*sqrt(var*4.*log(2)/derivative_varX));
        printf("FWHMy = sqrt(4*ln(2)/Var{SPMy'})   = %g PHYSICAL UNITS\n",yvoxdim*sqrt(var*4.*log(2)/derivative_varY));
        printf("FWHMz = sqrt(4*ln(2)/Var{SPMz'})   = %g PHYSICAL UNITS\n\n",zvoxdim*sqrt(var*4.*log(2)/derivative_varZ));

        printf("*********************************************************************\n");
        printf("ALTERNATE ESTIMATE OF SMOOTHNESS (as per Forman et al.):\n");
        printf("Effective standard deviation in X, in VOXEL LENGTHS  = %g\n",sx);
        printf("Effective standard deviation in Y, in VOXEL LENGTHS  = %g\n",sy);
        printf("Effective standard deviation in Z, in VOXEL LENGTHS  = %g\n",sz);
        printf("Effective FWHM in X, in VOXEL LENGTHS                = %g\n",FWHMx);
        printf("Effective FWHM in Y, in VOXEL LENGTHS                = %g\n",FWHMy);
        printf("Effective FWHM in Z, in VOXEL LENGTHS                = %g\n",FWHMz);
        printf("Effective overall FWHM in VOXEL LENGTHS              = %g\n",FWHM);
        printf("Effective overall W in VOXEL LENGTHS                 = %g\n\n",W2);
        printf("Effective standard deviation in X, in PHYSICAL UNITS = %g\n",xvoxdim*sx);
        printf("Effective standard deviation in Y, in PHYSICAL UNITS = %g\n",yvoxdim*sy);
        printf("Effective standard deviation in Z, in PHYSICAL UNITS = %g\n",zvoxdim*sz);
        printf("Effective FWHM in X, in PHYSICAL UNITS               = %g\n",xvoxdim*FWHMx);
        printf("Effective FWHM in Y, in PHYSICAL UNITS               = %g\n",yvoxdim*FWHMy);
        printf("Effective FWHM in Z, in PHYSICAL UNITS               = %g\n",zvoxdim*FWHMz);
        printf("Effective overall FWHM in PHYSICAL UNITS             = %g\n",voxdim*FWHM);

        fprintf(ftext,"\nINTERMEDIATE VARIABLES:\n");
        fprintf(ftext,"mean{SPMx'} = %g\n",derivative_meanX);
        fprintf(ftext,"mean{SPMy'} = %g\n",derivative_meanY);
        fprintf(ftext,"mean{SPMz'} = %g\n",derivative_meanZ);
        fprintf(ftext,"Var{SPMx'}  = %g\n",derivative_varX);
        fprintf(ftext,"Var{SPMy'}  = %g\n",derivative_varY);
        fprintf(ftext,"Var{SPMz'}  = %g\n",derivative_varZ);
        fprintf(ftext,"Geometric mean of Var{SPMi'}                  = %g\n",
                pow(derivative_varX*derivative_varY*derivative_varZ,1./3.));
        fprintf(ftext,"Geometric mean of Var{SPMi'} / Image Variance = %g\n",
                pow(derivative_varX*derivative_varY*derivative_varZ,1./3.)/var);
        fprintf(ftext,"D  = %g dimensions\n",D);
        fprintf(ftext,"S  = %g voxels\n\n",S);
        /* fprintf(ftext,"Pi = %g\n\n",pi); */

        fprintf(ftext,"*********************************************************************\n");
	fprintf(ftext,"ESTIMATE OF SPATIAL SMOOTHNESS ASSUMING UNIT VARIANCE:\n");
        fprintf(ftext,"W  = %g\n",W);
        fprintf(ftext,"effective FWHM = W * sqrt(4*ln(2)) = %g VOXEL LENGTHS\n",W*sqrt(4.*log(2)));
        fprintf(ftext,"FWHMx = sqrt(4*ln(2)/Var{SPMx'})   = %g VOXEL LENGTHS\n",sqrt(4.*log(2)/derivative_varX));
        fprintf(ftext,"FWHMy = sqrt(4*ln(2)/Var{SPMy'})   = %g VOXEL LENGTHS\n",sqrt(4.*log(2)/derivative_varY));
        fprintf(ftext,"FWHMz = sqrt(4*ln(2)/Var{SPMz'})   = %g VOXEL LENGTHS\n",sqrt(4.*log(2)/derivative_varZ));
        fprintf(ftext,"FWHMx = sqrt(4*ln(2)/Var{SPMx'})   = %g PHYSICAL UNITS\n",xvoxdim*sqrt(4.*log(2)/derivative_varX));
        fprintf(ftext,"FWHMy = sqrt(4*ln(2)/Var{SPMy'})   = %g PHYSICAL UNITS\n",yvoxdim*sqrt(4.*log(2)/derivative_varY));
        fprintf(ftext,"FWHMz = sqrt(4*ln(2)/Var{SPMz'})   = %g PHYSICAL UNITS\n",zvoxdim*sqrt(4.*log(2)/derivative_varZ));
        fprintf(ftext,"Gamma[%g] = %g\n\n",D/2.+1.,gamma_val);

        fprintf(ftext,"*********************************************************************\n");
        fprintf(ftext,"SMOOTHNESS WITH CORRECTION FOR NON-UNIT VARIANCE:\n");
        fprintf(ftext,"W  = %g\n",W*sqrt(var));
        fprintf(ftext,"effective FWHM = W * sqrt(4*ln(2)) = %g VOXEL LENGTHS\n",W*sqrt(var*4.*log(2)));
        fprintf(ftext,"FWHMx = sqrt(4*ln(2)/Var{SPMx'})   = %g VOXEL LENGTHS\n",sqrt(var*4.*log(2)/derivative_varX));
        fprintf(ftext,"FWHMy = sqrt(4*ln(2)/Var{SPMy'})   = %g VOXEL LENGTHS\n",sqrt(var*4.*log(2)/derivative_varY));
        fprintf(ftext,"FWHMz = sqrt(4*ln(2)/Var{SPMz'})   = %g VOXEL LENGTHS\n",sqrt(var*4.*log(2)/derivative_varZ));
        fprintf(ftext,"FWHMx = sqrt(4*ln(2)/Var{SPMx'})   = %g PHYSICAL UNITS\n",xvoxdim*sqrt(var*4.*log(2)/derivative_varX));
        fprintf(ftext,"FWHMy = sqrt(4*ln(2)/Var{SPMy'})   = %g PHYSICAL UNITS\n",yvoxdim*sqrt(var*4.*log(2)/derivative_varY));
        fprintf(ftext,"FWHMz = sqrt(4*ln(2)/Var{SPMz'})   = %g PHYSICAL UNITS\n\n",zvoxdim*sqrt(var*4.*log(2)/derivative_varZ));

        fprintf(ftext,"*********************************************************************\n");
        fprintf(ftext,"ALTERNATE ESTIMATE OF SMOOTHNESS (as per Forman et al.):\n");
        fprintf(ftext,"Effective standard deviation in X, in VOXEL LENGTHS  = %g\n",sx);
        fprintf(ftext,"Effective standard deviation in Y, in VOXEL LENGTHS  = %g\n",sy);
        fprintf(ftext,"Effective standard deviation in Z, in VOXEL LENGTHS  = %g\n",sz);
        fprintf(ftext,"Effective FWHM in X, in VOXEL LENGTHS                = %g\n",FWHMx);
        fprintf(ftext,"Effective FWHM in Y, in VOXEL LENGTHS                = %g\n",FWHMy);
        fprintf(ftext,"Effective FWHM in Z, in VOXEL LENGTHS                = %g\n",FWHMz);
        fprintf(ftext,"Effective overall FWHM in VOXEL LENGTHS              = %g\n",FWHM);
        fprintf(ftext,"Effective overall W in VOXEL LENGTHS                 = %g\n\n",W2);
        fprintf(ftext,"Effective standard deviation in X, in PHYSICAL UNITS = %g\n",xvoxdim*sx);
        fprintf(ftext,"Effective standard deviation in Y, in PHYSICAL UNITS = %g\n",yvoxdim*sy);
        fprintf(ftext,"Effective standard deviation in Z, in PHYSICAL UNITS = %g\n",zvoxdim*sz);
        fprintf(ftext,"Effective FWHM in X, in PHYSICAL UNITS               = %g\n",xvoxdim*FWHMx);
        fprintf(ftext,"Effective FWHM in Y, in PHYSICAL UNITS               = %g\n",yvoxdim*FWHMy);
        fprintf(ftext,"Effective FWHM in Z, in PHYSICAL UNITS               = %g\n",zvoxdim*FWHMz);
        fprintf(ftext,"Effective overall FWHM in PHYSICAL UNITS             = %g\n",voxdim*FWHM);



    /* ************************************************************************ */
    /* FREE MEMORY, EXIT PROGRAM.                                               */

    free(image);
    fclose(ftext);

    return(0);
}
