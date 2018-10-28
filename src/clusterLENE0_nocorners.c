#define MAX_SIZE_CLUSTERS 1000
#include <stdio.h>
#include <stdlib.h>
#include "analyze6.h"

                        /* GLOBAL VARIABLES.                                    */
int *clustermask;       /* Mask used to fill in clusters in recursive routine.  */
int *clustermask2;      /* Mask used to fill in outimg in recursive routine.    */
float *image;           /* Input image to be analyzed.                          */
unsigned short *outimg; /* Output image of clusters.                            */
int xdim, ydim, zdim;   /* Dimensions of image.                                 */
int num_cluster;        /* Cumulative number of clusters found.                 */
int clustersize;        /* Size of a cluster.  Initialize to zero before        */
                        /* invoking growcluster.                                */
float threshold;        /* Upper threshold for probabilities.                   */
float xmean, ymean,     /* Mean coordinate values and intensity value for       */
        zmean, meanval; /* cluster.                                             */

struct  analyze_struct readhdr(char *argv[]);
void growcluster(int x, int y, int z);
void filloutimg(int x, int y, int z);
unsigned char *unsigned_char_malloc(int num_items);
int read_unsigned_char_data(char *filename, unsigned char *array, int num_items);

/* ************************************************************************** */
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
01.25.96 : JMM : Modified from clusterGE to generate clusterGE_nocorners.
                 Program doesn't include corners in connectivity rule.

Given a floating pointing ANALYZE image file and a threshold, finds
clusters of voxels LESS THAN OR EQUAL TO a threshold.
A weighted average coordinate is generated and outputed as well.

This program makes no assumption about image dimensions.
They are read in from the header file.
*/
{
    FILE *ftext1, *ftext2;
    FILE *fIN;
    FILE *fOUT;

    struct analyze_struct img_info; /* Header information variable.     */
    int num_read_in;                /* Number of voxels read in.        */
    int num_wrote_out;              /* Number of voxels written out.    */
    int num_size[MAX_SIZE_CLUSTERS];/* Number of clusters of a certain  */
                                    /* size.                            */
    int vox_index,                  /* Index into voxels.               */
        num_vox,                    /* Total number of voxels.          */
        cluster_index,              /* Index into clusters.             */
        num_infrathresh_vox;        /* Number of INFRA-threshold voxels.*/
    int clustersize_index;          /* Index into sizes of clusters.    */
    int size_threshold;             /* Minimum size of clusters to be   */
                                    /* kept (optional argument)         */
    int xval,yval,zval;             /* Coordinates of a INFRA-threshold */
                                    /* voxel.                           */
    unsigned char *mask;            /* 8-bit mask.                      */
    int *xc, *yc, *zc;              /* List of coordinates of INFRA-    */
                                    /* threshold voxels.                */


    /* ****************************************************************** */
    /* CHECK COMMAND LINE ARGUMENTS*/
    if ((argc != 8) && (argc!=9)) {
        printf("clusterLENE0_nocorners: Need name of input .hdr file, input image file,\n");
        printf("output image file, output text file #1, output text file #2,\n");
        printf("threshold for output, and 8-bit mask.\n");
        printf("Clusters of voxels <= threshold will be found.\n");
        printf("You can include an optional eighth argument which is a size threshold.\n");
        printf("Only clusters of this size or greater will be appear in the output image file.\n");
        exit(1);
    }

    /* ****************************************************************** */
    /*  READ HEADER FILE.                                                 */
    img_info = readhdr(argv+1);

    xdim = img_info.dime.dim[1];
    ydim = img_info.dime.dim[2];
    zdim = img_info.dime.dim[3];
    num_vox = xdim*ydim*zdim;

    /* ****************************************************************** */
    /* READ IN IMAGE.                                                     */
    image      = (float *) malloc(sizeof(float)*num_vox);
    if (image==NULL) {
        printf("clusterLENE0_nocorners: unable to malloc for image.\n");
        exit(2);
    }
    printf("Loading %s...\n",argv[2]);
    if ((fIN = fopen(argv[2],"r")) == NULL) {
        printf ("clusterLENE0_nocorners: can't open %s\n",argv[2]);
        exit(2);
    }
    num_read_in=fread(image,sizeof(float),num_vox,fIN);
    fclose(fIN);
    if(num_read_in!=num_vox) {
        printf("clusterLENE0_nocorners: %d voxels read in; should have been %d.\n",
                                                num_read_in,num_vox);
        exit(3);
    }


    /* ****************************************************************** */
    /* MASK FLOATING POINT MAP FOR COLORATION.                            */
    mask = unsigned_char_malloc(num_vox);
    printf("Reading in mask from %s...\n",argv[7]);
    read_unsigned_char_data(argv[7],mask,num_vox);
    for (vox_index=0; vox_index<num_vox; vox_index++)
	if (mask[vox_index]==0) image[vox_index]=0.;

    /* ****************************************************************** */
    /* OPEN OUTPUT IMAGE AND TEXT FILES.                                  */
    if ((fOUT = fopen(*(argv+3),"w")) == NULL) {
        printf ("clusterLENE0_nocorners: can't open %s\n",*(argv+3));
        exit(2);
    }
    if ((ftext1 = fopen(*(argv+4),"w")) == NULL) {
        printf ("clusterLENE0_nocorners: can't open %s\n",*(argv+4));
        exit(2);
    }
    if ((ftext2 = fopen(*(argv+5),"a")) == NULL) {
        printf ("clusterLENE0_nocorners: can't open %s\n",*(argv+5));
        exit(2);
    }
 
    /* ****************************************************************** */
    /* ALLOCATE MEMORY FOR clustermask, clustermask2, AND outimg.         */
    clustermask = (int *) malloc(sizeof(int)*num_vox);
    if (clustermask==NULL) {
        printf("clusterLENE0_nocorners: unable to malloc for clustermask.\n");
        exit(4);
    }
    clustermask2 = (int *) malloc(sizeof(int)*num_vox);
    if (clustermask2==NULL) {
        printf("clusterLENE0_nocorners: unable to malloc for clustermask2.\n");
        exit(4);
    }
    outimg = (unsigned short *) malloc(sizeof(unsigned short)*num_vox);
    if (outimg==NULL) {
        printf("clusterLENE0_nocorners: unable to malloc for outimg.\n");
        exit(4);
    }

    /* ****************************************************************** */
    /* ZERO OUT outimg AND clustermask2.                                  */
    for (vox_index=0; vox_index<num_vox; vox_index++) {
        outimg[vox_index]=0;
        clustermask2[vox_index]=0;
    }

    /* ****************************************************************** */
    /* ACCEPT threshold AS COMMAND LINE ARGUMENT.                         */
    threshold = (float) atof(argv[6]);

    printf("Threshold   = %g\n",threshold);
    if (argc==9) {
        size_threshold= (int) atoi(argv[8]);
        printf("Size threshold = %d\n",size_threshold);
    }


    /******************************************************************** */
    /* SET ALL NEGATIVE VALUES TO 1.0                                     */
/*    for (vox_index=0; vox_index<num_vox; vox_index++)
        if (image[vox_index]<0)
            image[vox_index]=1.;
*/

    /* ****************************************************************** */
    /* COUNT NUMBER OF VOXELS LESS THAN/EQUAL TO threshold AND != 0.0     */
    num_infrathresh_vox = 0;
    for (vox_index=0; vox_index<num_vox; vox_index++)
        if ((image[vox_index] <= threshold) && (image[vox_index]!=0.))
            num_infrathresh_vox++;

    xc = (int *) malloc(num_infrathresh_vox*sizeof(int));
    yc = (int *) malloc(num_infrathresh_vox*sizeof(int));
    zc = (int *) malloc(num_infrathresh_vox*sizeof(int));

    num_infrathresh_vox=0;
    for (zval = 0; zval <= zdim-1; zval++)
        for (yval = 0; yval <= ydim-1; yval++)
            for (xval = 0; xval <= xdim-1; xval++) {
                if ((image[(zval*ydim+yval)*xdim+xval] <= threshold)
                 && (image[(zval*ydim+yval)*xdim+xval]!=0.0)) {
                    xc[num_infrathresh_vox]=xval;
                    yc[num_infrathresh_vox]=yval;
                    zc[num_infrathresh_vox]=zval;
                    num_infrathresh_vox++;
                }
            }
    printf("There are %d voxels less than/equal to Threshold = %g.\n",num_infrathresh_vox,threshold);
    printf("(and not equal to zero)\n");

    /* ****************************************************************** */
    /* WRITE HEADER                                                       */
    printf("Cluster#   Size     Value      x       y       z\n");
    fprintf(ftext1,"# Cluster#   Size     Value      x       y       z\n");


    /* ****************************************************************** */
    /* INITIALIZE num_size[MAX_SIZE_CLUSTERS]                             */
    for (clustersize_index=0; clustersize_index<MAX_SIZE_CLUSTERS; clustersize_index++)
        num_size[clustersize_index]=0;

    /* ****************************************************************** */
    /* CLEAR clustermask MATRIX.                                          */
    for (vox_index = 0; vox_index<num_vox; vox_index++)
        clustermask[vox_index] = 0;

    num_cluster=0;
    for (vox_index=0; vox_index<num_infrathresh_vox; vox_index++) {
        xval=xc[vox_index];
        yval=yc[vox_index];
        zval=zc[vox_index];

        if ((image[(zval*ydim+yval)*xdim+xval]<=threshold) 
                && (image[(zval*ydim+yval)*xdim+xval]!=0.)
                && (clustermask[(zval*ydim+yval)*xdim+xval]!=1)) {
        /* SINCE DATA HAS BEEN PRE-PROCESSED, VOXELS IN image ARE         */
        /* KNOWN TO BE LOCAL MAXIMA.  PRINT THE VALUE & COORDINATES OF    */
        /* CURRENT MAXIMUM, AND PRINT COORDINATES OF CONTIGUOUS VOXELS OF */
        /* THE SAME VALUE.                                                */
        num_cluster++;

        /* "GROW" REGION AROUND "SEED" POINT (xval,yval,zval); RETURN 1's */
        /* WHERE THE CORRESPONDING VOXEL IN image IS LESS THAN            */
        /* THRESHOLD AND IS CONTIGUOUS WITH THE VOXEL JUST OUTPUTED.      */
        clustersize=0;
        xmean=0.;
        ymean=0.;
        zmean=0.;
        meanval=0.;
        growcluster(xval,yval,zval);
        if (clustersize>MAX_SIZE_CLUSTERS) {
            printf("Warning: just found a cluster of size %d, greater than %d\n",
                                                                clustersize,MAX_SIZE_CLUSTERS);
            fprintf(ftext1," #Warning: just found a cluster of size %d, greater than %d\n",
                                                                clustersize,MAX_SIZE_CLUSTERS);
        }
        if (argc==8)
            filloutimg(xval,yval,zval);
        else if ((argc==9) && (clustersize>=size_threshold))
            filloutimg(xval,yval,zval);

        /* CALCULATE MEAN XYZ COORDINATES.                                */
        xmean/=(float)clustersize;
        ymean/=(float)clustersize;
        zmean/=(float)clustersize;
        meanval/=(float)clustersize;

        /* OUTPUT CLUSTERS FOUND.                                         */
        printf("    %3d      %3d  %9g  %6.2f  %6.2f  %6.2f %4d %4d %4d\n",
                num_cluster,clustersize,meanval,xmean+1,ymean+1,zmean+1,
                                                xval+1,yval+1,zval+1);
        fprintf(ftext1,"#    %3d      %3d  %9g  %6.2f  %6.2f  %6.2f %4d %4d %4d\n",
                num_cluster,clustersize,meanval,xmean+1,ymean+1,zmean+1,
                                                xval+1,yval+1,zval+1);
                
        /* UPDATE num_size.                                                  */
        num_size[clustersize-1]++;

        } /* END if */

    } /* END for LOOP */

    printf("\n\n Summary Table:\n");
    printf("Size  Number of Clusters\n");
    for (clustersize_index=0; clustersize_index<MAX_SIZE_CLUSTERS; clustersize_index++) {
        if (num_size[clustersize_index]!=0) {
            printf("%3d    %3d\n",clustersize_index+1,num_size[clustersize_index]);
            fprintf(ftext2, "%3d    %3d\n",clustersize_index+1,num_size[clustersize_index]);
        }
    }

    /* ****************************************************************** */
    /* WRITE OUTPUT IMAGE.                                                */
    num_wrote_out=fwrite(outimg,sizeof(unsigned short),num_vox,fOUT);
    if(num_wrote_out!=num_vox) {
        printf("clusterLENE0_nocorners: %d voxels wrote out; should have been %d.\n",
                                                        num_wrote_out,num_vox);

        exit(10);
    }

    /* ****************************************************************** */
    /* CLOSE OUTPUT FILE, FREE MEMORY, EXIT PROGRAM.                      */
    fclose(ftext1);
    fclose(ftext2);
    fclose(fOUT);
    free(image);
    free(outimg);
    free(clustermask);
    free(clustermask2);
    free(xc);
    free(yc);
    free(zc);

    return(0);
}

/* ************************************************************************** */
void growcluster(int x, int y, int z)
/* RECURSIVE FUNCTION TO "GROW" CLUSTERS.                                     */
/* ACCEPTS AS INPUT XYZ COORDINATES OF VOXEL IN *image AROUND WHICH           */
/* CONTIGUOUS VOXELS BELOW THRESHOLD ARE TO BE FOUND.                         */
/* FINDS ALL SUB-THREHOLD PIXELS IN image CONTIGUOUS WITH PIXEL(s) MASKED     */
/* OUT IN mask.  RETURNS "1" IN CORRESPONDING ENTRIES IN *clustermask.        */
{
    int xx,yy,zz;
    int vox_index;

    vox_index = (z*ydim+y)*xdim+x;
    if ((image[vox_index] <= threshold)
     && (image[vox_index] != 0.)
     && (clustermask[vox_index] == 0)) {
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

/* ************************************************************************** */
void filloutimg(int x, int y, int z)
/* RECURSIVE FUNCTION TO FILL IN OUTIMG WITH CLUSTERS.                        */
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

