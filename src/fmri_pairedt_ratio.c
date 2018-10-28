#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "analyze6.h"


		       /* FUNCTION DECLARATIONS.				*/
float *float_malloc(int num_items);
float ***float_malloc3(int num_items1, int num_items2, int num_items3);
unsigned char *unsigned_char_malloc(int num_items);
unsigned short *unsigned_short_malloc(int num_items);
char ***char_malloc3(int num_items1, int num_items2, int num_items3);
int read_unsigned_short_int_data(char *filename,
					unsigned short int *array, int num_items);

struct  analyze_struct readhdr(char *argv[]);
int float_fwrite(char *filename, float *array, int num_items);

/* **************************************************************************** */
/* Program reads in names of images from a text file, then reads in those	*/
/* images.  Images are scaled, and the global mean of each image is also scaled	*/
/* to some weighted average of the global means of the input images.		*/
/* Parametric images are calculated by voxel-by-voxel curve-fitting, and they	*/
/* are output to disk as floating point ANALYZE .img files.                     */
/*										*/
/* HISTORY: 02.11.93 : fitcurv2voxA written by Jose' Ma. Maisog, M.D.		*/
/* 	    02.12.93 : fitcurv2voxA modified by JMM to make fitcurv2voxB	*/
/*		       fitcurv2voxB normalizes the global mean of the input	*/
/*		       images to the global mean of the resting state image,	*/
/*		       or to some weighted average of the global means of the	*/
/*		       resting state images if there are more than one.  This	*/
/*		       weighted average is specified by weights given in the	*/
/*		       input text file.						*/
/*	    03.02.93 : added code to fitcurv2voxB to output in addition a	*/
/*		       MATLAB .mat file containing the variables scale_factor,	*/
/*		       NUMRTR and glob_mean, because these variables will be	*/
/*		       needed for image analysis in MATLAB			*/
/*		       (see read_param_filesA.m).				*/
/*	    03.03.93 : changed X-offset parameter into a steepnesS parameter,	*/
/*		       because an X-offset parameter is redundant -- the K	*/
/*		       parameter controls X-offset.  Needed a parameter to	*/
/*		       describe how steep the sigmoid is, so changed X-offset	*/
/*		       to a "slope" parameter.  Note that this steepness	*/
/*		       parameter is NOT a slope in the sense of a change in Y	*/
/*		       per change in X.						*/
/*		       Also, note that the function being minimized is not	*/
/*		       the squared error per se but the squared error plus 1.	*/
/*		       This is because the Nelder-Mead Simplex search algorithm	*/
/*		       as currently implemented relies on termination of the	*/
/*		       search by calculating an intermediate value rtol, but	*/
/*		       the way this rtol is calculated doesn't work for a cost	*/
/*		       function which is minimized at a value of zero; and the	*/
/*		       squared error function is minimized at zero.  However,	*/
/*		       minimizing the squared error plus 1 will be equivalent	*/
/*		       to minimizing the squared error itself.			*/
/*		       When the squared error image is written to disk as a	*/
/*		       MATLAB .mat data file, the image is decremented by 1 so	*/
/*		       the what is saved to disk is the true squared error, and	*/
/*		       not the squared error plus 1.				*/
/*	    03.15.93 : JMM: Changed ALL floats to doubles.  All calculations	*/
/*		       now done in double.					*/
/*	    10.06.93 : JMM : modified fitcurv2voxB.c to generate		*/
/*		       fmri_pairedt_ratio.c fmri_pairedt_ratio reads in 16-bit	*/
/*		       images from two groups, calculates a t-map, and outputs	*/
/*		       the t-map in a MATLAB .mat data file.  This .mat data	*/
/*		       file can be converted into an 8-bit ANALYZE image using	*/
/*		       mtz.							*/
/*          04.25.96 : JMM : Modified fmri_pairedt_ratio to output floating     */
/*                     point ANALYZE format .img files rather than MATLAB .mat  */
/*                     data files.  Changed all doubles back to floats.         */

main(int argc, char *argv[])
{
			    /* READING IN RAW DATA.				*/
    struct analyze_struct   /* Header information variable.			*/
		img_info;   /* 							*/
    FILE *fimgin;           /* Pointer for first input image file.		*/
    char choice;	    /* For keyboard input to yes/no question.		*/
    char dummy[80];	    /* Dummy pointer for counting number of entries	*/
			    /* in a text file.					*/
    char **header_name;	    /* Pointer to pointer to name of header file.	*/
    int grp_index;	    /* Index indicating group 1 vs group 2.		*/
    int n;		    /* Number of images in groups 1 and 2.		*/
    int n_index;	    /* Index indicating which image in group.		*/
    char ***img_name;	    /* Pointer to pointers to input file names.		*/
    int xdim,ydim,zdim;     /* Image dimensions.				*/
    int num_vox;	    /* Number of voxels in an image = xdim*ydim*zdim.	*/
    unsigned short *temp;   /* Temporary buffer for reading in input images.	*/
    float ***imgin;	    /* Pointer to pointers to input image data in	*/
			    /* floating point form.				*/

			    /* INTENSITY SCALING & GLOBAL MEAN NORMALIZATION	*/
    int vox_index;	    /* Index into imgin[*][*][], extracting a specific	*/
			    /* voxel value.					*/
    FILE *fintrmdt;         /* Pointer for output intermediate values .txt file.*/
    FILE *fmean;	    /* Pointer for output mean file.			*/
    float **mean;	    /* Pointer to mean images.				*/
    float *Pdiff;	    /* Pointer to difference of means images.	        */
    float *outmean;         /* floating point version of  mean[2], for writing  */
			    /* to disk.						*/
    float max_vox_val;	    /* Maximum voxel value in the mean image.		*/
    float min_vox_val;	    /* Minimum voxel value in the mean image.		*/
    float float_thresh;     /* Temporary variable for reading in threshold.     */
    float threshold;	    /* If the images are to be normalized to some	*/
			    /* global mean intensity, this is a threshold used	*/
			    /* to generate a mask from the average image.	*/
    float *mask;	    /* Pointer to mask image.				*/
    unsigned char *zeros;   /* Mask set to 0 whereever a zero voxel in any of   */
                            /* the original scans is found.                     */
    int num_mask_vox;	    /* Number of non-zero voxels in *mask.		*/
    float **glob_mean;	    /* Global mean intensity for images.		*/
    float float_c;	    /* Temporary variable for reading in *c.		*/
    float *c;		    /* Weighting factor for normalizing global means.	*/
    float NUMRTR;	    /* Numerator used for normalizing global means.	*/
    float SUM;		    /* Sum used for normalizing global means.		*/

			    /* t-IMAGE AND OUTPUT TO DISK.			*/
    float diff;	    /* Difference between pairs of voxels.		*/
    float avg_diff;	    /* Average difference between pairs of voxels.	*/
    float sd2;		    /* Variance of differences.				*/
    float *timg;	    /* Pointer to output image data, ANALYZE format.	*/
    int x,y,z;		    /* Indices into ANALYZE format image arrays.	*/

/* ---------------------------------------------------------------------------- */
/*  CHECK TO SEE IF CORRECT NUMBER OF ARGUMENTS WERE PASSED.			*/

    system("clear");
    if (argc < 5) {
	printf("fmri_pairedt_ratio: need header file, input textfile,\n");
	printf("                    output .img file name for t image,\n");
	printf("                    AND output .img file name for difference of means.\n\n");

	printf("The input textfile should contain the pathnames\n");
	printf("(either absolute or relative) of the group 1 input images,\n") ;
	printf("and the corresponding paired input group 2 textfile,\n");
	printf("with two paired filenames per line.\n");
	printf("To correct for differences in global mean intensity, the\n");
	printf("following multiplicative scaling will be done on an image-\n");
	printf("by-image basis.  A global mean GBM(i)  will be calculated for\n");
	printf("each image IMG(i).  Then each image IMG(i) will be multiplied by\n\n");

	printf("                  NUMRTR/GBM(i)\n\n");

	printf("where NUMRTR is a weighted average of all the GBM's.\n");
	printf("Let c(i) be the weight assigned to GBM(i).  Let\n\n");

	printf("           SUM = abs(c(1))+abs(c(2))+...\n\n");

	printf("for all the c(i)`s.  Then\n\n");

	printf("         NUMRTR = [c(1)*GBM(1)+c(2)*GBM(2)+...]/SUM.\n\n");

	printf("For example, if you want the numerator to be the average of\n");
	printf("GBM(1) and GBM(6) set\n\n");

	printf("          c(1) = 1.0\n");
	printf("          c(2) = 0.0\n");
	printf("          c(3) = 0.0\n");
	printf("          c(4) = 0.0\n");
	printf("          c(5) = 0.0\n");
	printf("          c(6) = 1.0\n\n");

	printf("The intermediate variable file intrmPT.txt will contain the\n");
        printf("intermediate variables c, NUMRTR, glob_mean, and threshold\n");
        printf("(the number used to generate the mask for calculation of the\n");
        printf("global means).\n\n");

	exit(1);
    }

/* ---------------------------------------------------------------------------- */
/*  READ HEADER FILE; ASSUME IMAGES ARE SAME SIZE, SO USE SAME HEADER FILE	*/
/*  FOR ALL IMAGES.								*/

    img_info = readhdr(argv+1);

    xdim = img_info.dime.dim[1];
    ydim = img_info.dime.dim[2];
    zdim = img_info.dime.dim[3];

    num_vox = xdim*ydim*zdim;

/* ---------------------------------------------------------------------------- */
/*  COUNT NUMBER OF ENTRIES IN INPUT TEXT FILE. 				*/

    if ((fimgin = fopen(argv[2],"r")) == NULL) {
	printf ("fmri_pairedt_ratio: can't open %s\n",argv[1]);
	exit(2);
    }
    n = 0;
    fscanf(fimgin,"%s",dummy);
    while(fscanf(fimgin,"%s %s %f",dummy,dummy,&float_c)!=EOF)
	n++;
    fclose(fimgin);
    printf("There are %d paired images.\n",n);

/* ---------------------------------------------------------------------------- */
/*  READ IN FILENAMES AND GLOBAL MEAN WEIGHTS IN INPUT TEXT FILE.		*/

    c = (float *) float_malloc(n);
    img_name=char_malloc3(2, n, 80); 

    fimgin = fopen(argv[2],"r");
    for (n_index=0; n_index<=n-1; n_index++) {
	fscanf(fimgin,"%s %s %f",
				img_name[0][n_index],img_name[1][n_index],&float_c);
	c[n_index] = (float) float_c;
    }
    fclose(fimgin);


/* ---------------------------------------------------------------------------- */
/*  SHOW USER INFORMATION READ IN SO FAR.					*/

    printf("\nThese are the %d paired images and corresponding weighting factors:\n",
										n);
    for (n_index=0; n_index<=n-1; n_index++)
	printf("  %20s  %20s %f\n", img_name[0][n_index],img_name[1][n_index],
									c[n_index]);

    printf("\nImages have dimensions:\n");
    printf("  xdim = %d\n  ydim = %d\n  zdim = %d\n",xdim,ydim,zdim);

    printf("\nThe t-test map will be written to the floating point file %s.\n", argv[3]);
    printf("The map of difference in means will be written to the floating point file %s.\n", argv[4]);

    printf("\nIt is assumed that the input images are registered.  A mean\n");
    printf("image will be calculated.  A mask will be made from this mean\n");
    printf("image by thresholding it.  You will need to specify the\n");
    printf("number T below.  All voxels in the mean image which have pixel\n");
    printf("intensity T %% the maximum voxel intensity in the mean image or\n");
    printf("greater will be taken to represent gray matter.\n");
    printf("This mask will be used to calculate the global mean\n");
    printf("intensities of the individual images.\n");
    printf("\nENTER the threshold T for the mask: ");
    scanf("%f",&float_thresh);
    threshold = (float) float_thresh;
    printf("\nThreshold = %3.2f %%.\n",threshold);

    printf("Intermediate variables will be written to the textfile intrmPT.txt.\n");

    printf("\nDoes the above information look okay <y/n>? ");
    scanf("%s",&choice);
    while((choice!='y') && (choice!='n')) {
	printf("Type 'y' or 'n': ");
	scanf("%s",&choice);
    }
    if (choice=='n') exit(6);

/* ---------------------------------------------------------------------------- */
/*  READ IN GROUP 1 AND GROUP 2 IMAGES.  ALLOCATE MEMORY FOR OUTPUT IMAGE.	*/

    temp  = unsigned_short_malloc(num_vox);
    zeros = unsigned_char_malloc(num_vox);
    imgin  = float_malloc3(2, n, num_vox);

    /* INITIALIZE zeros[] */
    for (vox_index=0; vox_index<num_vox; vox_index++)
        zeros[vox_index]=0;

    for (n_index=0; n_index<=n-1; n_index++) {
	printf("Reading in image from %s....\n",img_name[0][n_index]);
        read_unsigned_short_int_data(img_name[0][n_index], temp, num_vox);

	for (vox_index=0; vox_index<num_vox; vox_index++) {
	    imgin[0][n_index][vox_index]=(float)temp[vox_index];
            zeros[vox_index]=(temp[vox_index]==0)?zeros[vox_index]:1;
        }

	printf("Reading in image from %s....\n",img_name[1][n_index]);
        read_unsigned_short_int_data(img_name[1][n_index], temp, num_vox);

	for (vox_index=0; vox_index<num_vox; vox_index++) {
	    imgin[1][n_index][vox_index]=(float)temp[vox_index];
            zeros[vox_index]=(temp[vox_index]==0)?zeros[vox_index]:1;
        }
    }

    timg = (float *) malloc(num_vox*sizeof(float));
    if (timg==NULL) {
	printf("fmri_pairedt_ratio: error malloc'ing for timg.\n");
	exit(10);
    }


    /* ------------------------------------------------------------------------ */
    /* ALLOCATE MEMORY FOR mean, glob_mean, AND outmean.                	*/

    mean = (float **) malloc(3*sizeof(float *));
    mean[0] = (float *) malloc(num_vox*sizeof(float));
    mean[1] = (float *) malloc(num_vox*sizeof(float));
    mean[2] = (float *) malloc(num_vox*sizeof(float));
    if (mean==NULL || mean[0]==NULL || mean[1]==NULL || mean[2]==NULL) {
	printf("fmri_pairedt_ratio: error malloc'ing for mean.\n");
	exit(11);
    }

    glob_mean = (float **) malloc(2*sizeof(float *));
    glob_mean[0] = (float *) malloc(n*sizeof(float));
    glob_mean[1] = (float *) malloc(n*sizeof(float));
    if (glob_mean==NULL || glob_mean[0]==NULL || glob_mean[1]==NULL) {
	printf("fmri_pairedt_ratio: error malloc'ing for glob_mean.\n");
	exit(12);
    }
    Pdiff = (float *) malloc(num_vox*sizeof(float));
    if (Pdiff == NULL) {
	printf("fmri_pairedt_ratio: failed to malloc for Pdiff.\n");
	exit(5);
    }
    outmean = (float *) malloc(num_vox*sizeof(float));
    if (outmean==NULL) {
	printf("fmri_pairedt_ratio: error malloc'ing for outmean.\n");
	exit(13);
    }


    /* ------------------------------------------------------------------------ */
    /* CALCULATE MEAN IMAGES.  mean[0] AND mean[1] ARE THE MEANS FOR THE TWO	*/
    /* GROUPS, TO BE USED IN CALCULATING THE t MAPS,  mean[2] IS THE MEAN OF	*/
    /* ALL THE IMAGES, WHICH IS EQUAL TO THE MEAN OF mean[0] AND mean[1].	*/
    /*										*/
    /* FOR EACH GROUP....							*/

    for (grp_index=0; grp_index<=1; grp_index++) {

	/* -------------------------------------------------------------------- */
	/* INITIALIZE MEAN IMAGE.						*/

	printf("Initializing group %d mean image...\n",grp_index+1);
	for (vox_index=0; vox_index<=num_vox-1; vox_index++)
	    mean[grp_index][vox_index] = 0.;

	/* -------------------------------------------------------------------- */
	/* CALCULATE MEAN IMAGE.						*/
	/* This mean image is used to generate a mask, which in turn is		*/
	/* used to calculate global means.  Later, the mean is recalculated	*/
	/* using the normalized values, and then used to calculate the t-test.	*/

	printf("  Calculating mean image...\n");
	for (n_index=0; n_index<=n-1; n_index++)
	    for (vox_index=0; vox_index<=num_vox-1; vox_index++)
		mean[grp_index][vox_index] += imgin[grp_index][n_index][vox_index];

	for (vox_index=0; vox_index<=num_vox-1; vox_index++)
	    mean[grp_index][vox_index] /= ((float) n);

    }

    printf("Making overall mean image.\n");
    for (vox_index=0; vox_index<=num_vox-1; vox_index++)
	mean[2][vox_index] = (mean[0][vox_index] + mean[1][vox_index])/2.;

    /* ------------------------------------------------------------------------ */
    /* MAKE MASK FROM mean[2].							*/

    printf("Making mask image.\n");

	/* -------------------------------------------------------------------- */
	/* DETERMINE MAXIMUM AND MINIMUM PIXEL VALUES OF MEAN IMAGE.		*/

	max_vox_val = 0.;
	min_vox_val = 100000.;
	for (vox_index=0; vox_index<=num_vox-1; vox_index++) {
	    max_vox_val = (mean[2][vox_index] > max_vox_val)
					? mean[2][vox_index] : max_vox_val;
	    min_vox_val = (mean[2][vox_index] < min_vox_val)
					? mean[2][vox_index] : min_vox_val;
	}

	/* -------------------------------------------------------------------- */
	/* THRESHOLD MEAN IMAGE USING threshold, MAKING *mask.			*/

	mask = (float *) malloc(num_vox*sizeof(float));
	if (mask==NULL) {
	    printf("fmri_pairedt_ratio: error malloc'ing for mask.\n");
	    exit(13);
	}
	for (vox_index=0; vox_index<=num_vox-1; vox_index++)
	    mask[vox_index] =
		(mean[2][vox_index] >= (max_vox_val*threshold/100.)) ? 1. : 0.;

	/* -------------------------------------------------------------------- */
	/* DELETE FROM MASK ANY VOXELS WHICH ARE SET TO ZERO IN ANY ORIGINAL    */
        /* OF THE INPUT fMRI SCANS.                                             */

    	for (vox_index=0; vox_index<=num_vox-1; vox_index++)
             if (zeros[vox_index]==0) mask[vox_index]=0.;


	/* -------------------------------------------------------------------- */
	/* MASK OUT mean[2].							*/

	for (vox_index=0; vox_index<=num_vox-1; vox_index++)
	    mean[2][vox_index] *= mask[vox_index];


	/* -------------------------------------------------------------------- */
	/*  WRITE outmean TO DISK.						*/

        float_fwrite("meanPT.img",mean[2],num_vox);

	/* -------------------------------------------------------------------- */
	/* COUNT NUMBER OF NON-ZERO VOXELS IN *mask.				*/

	num_mask_vox = 0.;
	for (vox_index=0; vox_index<=num_vox-1; vox_index++)
	    if (mask[vox_index] != 0.)
		num_mask_vox++;


    /* ------------------------------------------------------------------------ */
    /* NOW NORMALIZE IMAGES.  FOR EACH GROUP....				*/

    for (grp_index=0; grp_index<=1; grp_index++) {

	printf("\nNormalizing group %d images....\n",grp_index+1);

	/* -------------------------------------------------------------------- */
	/* CALCULATE GLOBAL MEAN OF EACH INPUT IMAGE USING mask.		*/

	for (n_index=0; n_index<=n-1; n_index++) {
	    glob_mean[grp_index][n_index]=0.;
	    for (vox_index=0; vox_index<=num_vox-1; vox_index++)
		glob_mean[grp_index][n_index] +=
			imgin[grp_index][n_index][vox_index] * mask[vox_index];
	    glob_mean[grp_index][n_index] /= ((float) num_mask_vox);
	    printf("  global mean of %s = %f\n", img_name[grp_index][n_index],
						glob_mean[grp_index][n_index]);
	}

	/* -------------------------------------------------------------------- */
	/* CALCULATE SUM.							*/

	SUM = 0;
	for (n_index=0; n_index<=n-1; n_index++)
	    SUM+=((int) (fabs(c[n_index])));

	/* -------------------------------------------------------------------- */
	/* CALCULATE NUMRTR.							*/

	NUMRTR = 0.;
	for (n_index=0; n_index<=n-1; n_index++)
	    NUMRTR += (c[n_index] * glob_mean[grp_index][n_index]);

	NUMRTR /= SUM;
	printf("  Group # %d NUMRTR = %f\n", grp_index+1, NUMRTR);
	    
	/* -------------------------------------------------------------------- */
	/* NORMALIZE INDIVIDUAL IMAGES.  MULTIPLY EACH VOXEL BY			*/
	/* NUMRTR/glob_mean[i].							*/
	/* ALSO MASK OUT IMAGES, USING mask.					*/

	printf("  Normalizing group # %d image global means...\n", grp_index+1);
	for (n_index=0; n_index<=n-1; n_index++)
	    for (vox_index=0; vox_index<=num_vox-1; vox_index++)
		imgin[grp_index][n_index][vox_index] *=
			(mask[vox_index]*NUMRTR/glob_mean[grp_index][n_index]);

    }

    /* ------------------------------------------------------------------------ */
    /* SAVE NUMRTR, glob_mean, & scale_factor TO OUTPUT FILE.			*/
    /* ALLOCATE MEMORY FOR xi.							*/

    if ((fintrmdt = fopen("intrmPT.txt","w")) == NULL) {
	printf ("fmri_pairedt_ratio: can't open intrmPT.txt.\n");
	exit(15);
    }

    fprintf(fintrmdt,"NUMRTR     = %g\n\n\n",NUMRTR);
    fprintf(fintrmdt,"threshold  = %g\n\n\n",threshold);

    fprintf(fintrmdt,"Global means:\n");
    for (n_index=0; n_index<=n-1; n_index++) {
        for (grp_index=0; grp_index<=1; grp_index++)
            fprintf(fintrmdt,"%g ",glob_mean[grp_index][n_index]);
        fprintf(fintrmdt,"\n");
    }
    fprintf(fintrmdt,"\n\n");


    fprintf(fintrmdt,"Global means:\n");
    for (n_index=0; n_index<=n-1; n_index++)
        fprintf(fintrmdt,"%g ",c[n_index]);

    fclose(fintrmdt);


    /* ------------------------------------------------------------------------ */
    /* RECALCULATE mean[][], USING NEWLY NORMALIZED VALUES.			*/
    /* FOR EACH GROUP....							*/

    for (grp_index=0; grp_index<=1; grp_index++) {

	/* -------------------------------------------------------------------- */
	/* INITIALIZE MEAN IMAGE.						*/

	for (vox_index=0; vox_index<=num_vox-1; vox_index++)
	    mean[grp_index][vox_index] = 0.;

	/* -------------------------------------------------------------------- */
	/* RECALCULATE MEAN IMAGE.						*/

	printf("Calculating normalized mean image for group %d...\n",grp_index+1);
	for (n_index=0; n_index<=n-1; n_index++)
	    for (vox_index=0; vox_index<=num_vox-1; vox_index++)
		mean[grp_index][vox_index] += imgin[grp_index][n_index][vox_index];

	for (vox_index=0; vox_index<=num_vox-1; vox_index++)
	    mean[grp_index][vox_index] /= ((float) n);
    }

/* ---------------------------------------------------------------------------- */
/*  CALCULATE DIFFERENCE OF MEAN IMAGES.                                        */

    printf("Calculating difference of mean images...\n");
    for (vox_index=0; vox_index<=num_vox-1; vox_index++)
	Pdiff[vox_index]=(mask[vox_index]==1.)
		? mean[1][vox_index]-mean[0][vox_index] : 0.;

    /* ------------------------------------------------------------------------ */
    /*  AND THEN SAVE TO DISK.                                                  */

    printf("Writing output to %s...\n",argv[4]);
    float_fwrite(argv[4],Pdiff,num_vox);

    /* ------------------------------------------------------------------------ */
    /* NOW CALCULATE t VOXEL-BY-VOXEL.						*/

    printf("Now calculating t images....\n");

    for (vox_index=0; vox_index<=num_vox-1; vox_index++) {
	avg_diff = 0.;
	for (n_index=0; n_index<=n-1; n_index++)
	    avg_diff +=
		((float) imgin[1][n_index][vox_index]
					- (float) imgin[0][n_index][vox_index]);
	avg_diff /= (float) n;

	sd2 = 0.;
	for (n_index=0; n_index<=n-1; n_index++) {
	    diff = (float) imgin[1][n_index][vox_index] 
					- (float) imgin[0][n_index][vox_index];
	    sd2 += ((diff - avg_diff) * (diff - avg_diff));
	}
	sd2 /= (float) (n-1);
	if (sd2==0.)
	    timg[vox_index] = 0;
	else
	    timg[vox_index] = avg_diff / (sqrt(sd2/(float)n));
/*	if (isnan(timg[vox_index]) || isinf(timg[vox_index])) {
	    printf("NaN or Inf generated!\n");
	    printf("n         = %d\n",n);
	    printf("vox_index = %d\n",vox_index);
	    for (n_index=0; n_index<=n-1; n_index++)
		printf("  %d %d %d\n",n_index+1,imgin[0][n_index][vox_index],
							imgin[1][n_index][vox_index]);
	    printf("avg_diff  = %f\n",avg_diff);
	    printf("sd2       = %f\n",sd2);
	    timg[vox_index] = 30.;
	}
*/
    }


/* ---------------------------------------------------------------------------- */
/*  WRITE OUT T-TEST MAP TO OUTPUT FILE.                                        */

    printf("Writing output to %s...\n",argv[3]);
    float_fwrite(argv[3],timg,num_vox);

/* ---------------------------------------------------------------------------- */
/*  FREE MEMORY.								*/

    printf("Releasing memory...\n");
    free(header_name);
    free(c);
    free(img_name);
    free(temp);
    free(imgin);
    free(timg);
    free(mean);
    free(Pdiff);
    free(glob_mean);
    free(outmean);
    free(mask);
}
