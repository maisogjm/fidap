#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "analyze6.h"

		       /* FUNCTION DECLARATIONS.				*/
struct  analyze_struct readhdr(char *argv[]);
float factorial(int input);

/* **************************************************************************** */
/* Program fmri_df.c							*/
/* ANSI C Code by Jose' Ma. Maisog.						*/
/*                                                                              */

main(int argc, char *argv[])
{
			    /* INPUT IMAGES.					*/
    FILE *ftxt;             /* Pointer for input text file containing filenames	*/
    char **header_name;	    /* Pointer to pointer to name of header file.	*/
    char **img_name;	    /* Pointer to pointers to input image file names.	*/
    char choice;	    /* For keyboard input to yes/no question.		*/
    char dummy[80];	    /* Dummy pointer for counting number of entries	*/
			    /* in a text file.					*/
    struct analyze_struct img_info; /* Header information variable.		*/
    int     xdim,ydim,zdim; /* Image dimensions.				*/
    int x,y,z;		    /* Indices into image arrays.			*/
    int num_vox;	    /* Number of voxels in an image.			*/
    int     vox_index;	    /* Index into imgin[][*], extracting a specific	*/
			    /* voxel.						*/
    FILE *fimgin;           /* Pointer for input image file.			*/
    unsigned		    /* Pointer to temporary short int array, for	*/
	short int *temp;    /* reading in 16-bit fMRI data.			*/
			    /* byte data.					*/
    int num_read_in;	    /* Number of objects read in in fread operation.	*/
    float **imgin;	    /* Pointer to input image data.			*/
    float FWHM;		    /* FWHM of smoothing kernel used to smooth time	*/
			    /* series.						*/
    float SD;		    /* Standard deviation of gaussian smoothing kernel,	*/
			    /* calculated from FWHM.				*/
    float pi;		    /* Value of Pi.					*/
    float *gaussian;	    /* Filter kernel used to smooth time series.	*/
    int kernel_index;	    /* Index used to smooth time series with kernel.	*/
    float *smooth_series;   /* Smoothed version of time series.			*/
    float *dx_dt;	    /* Derivatives, first differences of time series.	*/
    int j,k;		    /* Indices for generating the vector xp at various	*/
			    /* lags, following the same names used in Friston's	*/
			    /* MATLAB code.					*/
    float imgin_mean;	    /* Mean of input image data.			*/
    float sum_imgin_mean;   /* Sum of means of time series.			*/
    float mean_imgin_mean;  /* Mean of means of time series.			*/
    float imgin_var;        /* Variance of input image data.			*/
    float sum_imgin_var;    /* Sum of variances of time series.			*/
    float mean_imgin_var;   /* Mean of variances of time series.		*/
    int num_dat_pts;	    /* The number of abscissa values.			*/
    int dat_index;	    /* Index into vector of abscissa values		*/
    int dat_index2;	    /* Index into vector of abscissa values		*/
			    /* Used for convolving smoothing kernel with fMR	*/
			    /* time series.					*/

			    /* MEASURES OF IMAGE SMOOTHNESS.			*/
    float derivative_sum;   /* Sum of derivatives at a particular voxel.	*/
    float derivative_mean;  /* Mean of derivatives at a voxel.			*/
    float		    /* Derivatives minus derivative_mean.		*/
	residual_derivative;
    float mean_dir;	    /* Mean of derivatives for a voxel.			*/
    float var_dir;	    /* Variance of derivatives for a voxel.		*/
    float sum_var_dir;	    /* Sum of variances of derivatives across voxels.	*/
    float mean_var_dir;	    /* Mean of variances of derivatives across voxels.	*/
    float W;		    /* Smoothness factor of time series.		*/
    float sum_W;	    /* Sum of W's across voxels.			*/
    float mean_W;	    /* Mean of W's across voxels.			*/
    float effective_FWHM;   /* Effective FWHM of time series.			*/
    float sum_effective_FWHM;/* Sum of effective FWHM of time series.		*/
    float mean_effective_FWHM;/* Sum of effective FWHM of time series.		*/

    float v;		    /* Effective degrees of freedom.			*/
    float sum_v;	    /* Sum of effective degrees of freedom.		*/
    float mean_v;	    /* Mean of effective degrees of freedom.		*/
    int num_nonzero_vardir; /* Number of voxels for which effective df was	*/
			    /* actually calculated.  Effective df is undefined	*/
			    /* when Var{dx/dt} is 0, since that term appears in	*/
			    /* the denominator.					*/

			    /* MASK						*/
    FILE *fMASK;	    /* Pointer for input mask .img file.		*/
    unsigned char *mask;    /* Pointer to mask image.				*/
    int num_mask_vox;	    /* Number of non-zero voxels in *mask.		*/



/* ---------------------------------------------------------------------------- */
/*  CHECK TO SEE IF CORRECT NUMBER OF ARGUMENTS WERE PASSED.			*/

    if (argc != 4) {
	printf("fmri_df: need input textfile of file names, name of\n");
	printf("               8-bit mask, and FWHM of smoothing kernel.\n\n");

	printf("The first line in the input text file should be the name of\n");
	printf("the header file to be used with the input images.  The rest of\n") ;
	printf("the input text file should contain the names of the input .img\n") ;
	printf("files, with one filename/number per line.\n");
	printf("If you wish NO temporal smoothing to be done, type in a\n");
	printf("zero \"0\" for the FWHM.\n\n");

	exit(1);
    }


/* ---------------------------------------------------------------------------- */
/*  COUNT NUMBER OF DATA POINTS IN INPUT TEXT FILE.				*/

    if ((ftxt = fopen(argv[1],"r")) == NULL) {
	printf ("fmri_df: can't open %s\n",*(argv+1));
	exit(2);
    }
    num_dat_pts = 0;
    fscanf(ftxt,"%s",dummy);
    while(fscanf(ftxt,"%s",dummy)!=EOF)
	num_dat_pts++;


/* ---------------------------------------------------------------------------- */
/*  READ IN NAMES OF INPUT IMAGE FILES.						*/

    header_name = (char **) malloc(sizeof(char *));
    if (header_name==NULL) {
	printf("fmri_df: error malloc'ing for header_name.\n");
	exit(3);
    }
    *header_name = (char *) malloc(80*sizeof(char));
    if (*header_name==NULL) {
	printf("fmri_df: error malloc'ing for *header_name.\n");
	exit(4);
    }
    img_name = (char **) malloc(num_dat_pts*sizeof(char *));
    if (img_name==NULL) {
	printf("fmri_df: error malloc'ing for img_name.\n");
	exit(5);
    }

    rewind(ftxt);
    fscanf(ftxt,"%s",*header_name);
    for (dat_index=0; dat_index<=num_dat_pts-1; dat_index++) {
	img_name[dat_index] = (char *) malloc(80*sizeof(char));
	if (img_name[dat_index]==NULL) {
	    printf("fmri_df: error malloc'ing for img_name[%d].\n",dat_index);
	    exit(6);
	}
	fscanf(ftxt,"%s %f",img_name[dat_index]);
    }
    fclose(ftxt);


/* ---------------------------------------------------------------------------- */
/*  READ HEADER FILE; ASSUME IMAGES ARE SAME SIZE, SO USE SAME HEADER FILE	*/
/*  FOR ALL IMAGES.								*/

    img_info = readhdr(header_name);

    xdim = img_info.dime.dim[1];
    ydim = img_info.dime.dim[2];
    zdim = img_info.dime.dim[3];

    num_vox = xdim*ydim*zdim;


/* ---------------------------------------------------------------------------- */
/*  SHOW USER INFORMATION READ IN SO FAR.					*/

    printf("There are %d input images.\n", num_dat_pts);
    printf("\nImages have dimensions:\n");
    printf("  xdim = %d\n  ydim = %d\n  zdim = %d\n\n",xdim,ydim,zdim);

/*    printf("Here are the %d images:\n",num_dat_pts);
    for (dat_index=0; dat_index<=num_dat_pts-1; dat_index++)
	printf("%s\n",img_name[dat_index]);

    printf("\nMask will be read in from 8-bit image %s.\n",argv[2]);

    printf("\nDoes the above information look okay <y/n>? ");
    scanf("%s",&choice);
    while((choice!='y') && (choice!='n')) {
	printf("Type 'y' or 'n': ");
	scanf("%s",&choice);
    }
    if (choice=='n') exit(4);
*/

/* ---------------------------------------------------------------------------- */
/*  ALLOCATE MEMORY.								*/

    temp = (unsigned short int *) malloc(num_vox*sizeof(unsigned short int));
    if (temp==NULL) {
	printf("fmri_df: unable to malloc for temp.\n");
	exit(10);
    }
    imgin = (float **) malloc(num_vox*sizeof(float *));
    if (imgin==NULL) {
	printf("fmri_df: unable to malloc for imgin.\n");
	exit(11);
    }
    for (vox_index=0; vox_index<=num_vox-1; vox_index++) {
	imgin[vox_index]=(float *)malloc(num_dat_pts*sizeof(float));
	if (imgin[vox_index]==NULL) {
	    printf("fmri_df: unable to malloc for imgin[%d].\n",vox_index);
	    exit(12);
	}
    }
    gaussian = (float *) malloc(num_dat_pts*sizeof(float));
    if (gaussian==NULL) {
	printf("fmri_df: unable to malloc for gaussian.\n");
	exit(11);
    }
    smooth_series = (float *) malloc(num_dat_pts*sizeof(float));
    if (smooth_series==NULL) {
	printf("fmri_df: unable to malloc for smooth_series.\n");
	exit(11);
    }
    dx_dt = (float *) malloc((num_dat_pts-1)*sizeof(float));
    if (dx_dt==NULL) {
	printf("fmri_df: unable to malloc for dx_dt.\n");
	exit(11);
    }


/* ---------------------------------------------------------------------------- */
/*  READ IN IMAGES.								*/

    for (dat_index=0; dat_index<num_dat_pts; dat_index++) {
	if ((fimgin = fopen(img_name[dat_index],"r")) == NULL) {
	    printf("fmri_df: can't open %s\n",img_name[dat_index]);
	    exit(6);
	}
	printf("Reading in data from %s...\n",img_name[dat_index]);
	num_read_in=fread(temp,sizeof(unsigned short int),num_vox,fimgin);
	fclose(fimgin);
	if (num_read_in<num_vox) {
	    printf("Read in only %d objects, not %d.\n",num_read_in,num_vox);
	    exit(7);
	}
	for (vox_index=0; vox_index<num_vox; vox_index++)
	    imgin[vox_index][dat_index]=(float)temp[vox_index];
    }

/* ---------------------------------------------------------------------------- */
/* READ IN 8-BIT MASK.								*/

    mask = (unsigned char *) malloc(num_vox*sizeof(unsigned char));
    if (mask == NULL) {
	printf("fmri_df: malloc failed for mask.\n");
	exit(15);
    }
    if ((fMASK = fopen(argv[2],"r")) == NULL) {
	printf("fmri_df: can't open %s\n",argv[2]);
	exit(13);
    }
    printf("Reading in mask from %s...\n",argv[2]);
    num_read_in=fread(mask,sizeof(unsigned char),num_vox,fMASK);
    fclose(fMASK);
    if(num_read_in!=num_vox) {
	printf("fmri_df: %d voxels read in; should have been %d.\n",
						num_read_in,num_vox);
	exit(16);
    }
    num_mask_vox=0;
    for (vox_index=0; vox_index<num_vox; vox_index++)
	if (mask[vox_index]!=0)
	    num_mask_vox++;
    printf("\nN = %d\n",num_dat_pts);

    printf("Number of voxels in mask = %d.\n",num_mask_vox);


/* ---------------------------------------------------------------------------- */
/* READ IN FWHM.  CALCULATE SD OF GAUSSIAN SMOOTHING KERNEL, WHICH IS RELATED	*/
/* TO FWHM.									*/

    FWHM=atof(argv[3]);
    SD = FWHM/(2.*(float)sqrt(2.*log(2.)));
    pi = 2.*acos(0.);
    if (FWHM!=0.) {
	printf("FWHM = %f 'voxels'.\n",FWHM);
	printf("SD = %f 'voxels'.\n",SD);
    }
    else printf("No temporal smoothing will be done.\n");

    printf("PI = %f\n",pi);


/* ---------------------------------------------------------------------------- */
/* CALCULATE GAUSSIAN SMOOTHING KERNEL.						*/
/* Time series data may be smoothed to suppress thermal noise.			*/

    if (FWHM!=0.) {
	for (dat_index=0; dat_index<num_dat_pts; dat_index++) {
	    kernel_index= (float)(num_dat_pts/2)
					-
			(float) (fabs((double)(dat_index-(num_dat_pts/2))));
	    gaussian[dat_index]=(1/(SD*(float)sqrt((double)(2.*pi))))
						*
			exp(-((float)(kernel_index*kernel_index))/(2*SD*SD));
	}
    }


/* ---------------------------------------------------------------------------- */
/* CALCULATE Var{dx/dt}, VOXEL-BY-VOXEL.					*/
/* Calculate smoothness factor on all voxels masked in by mask.  Then will	*/
/* calculate the mean smoothness factor over all these voxels.			*/

    printf("\nNow calculating mean values for Var{x(t)} and Var{dx/dt}...\n");
    num_nonzero_vardir=0;
    sum_imgin_mean=0.;
    sum_imgin_var=0.;
    sum_var_dir=0.;
    sum_W=0.;
    sum_effective_FWHM=0.;
    sum_v=0.;
    for (vox_index=0; vox_index<num_vox; vox_index++) {

    /* ------------------------------------------------------------------------ */
    /* IF THIS VOXEL IS IN MASK, CALCULATE Var{dx/dt} ON IT.			*/

	if (mask[vox_index]==1) {

	/* -------------------------------------------------------------------- */
	/* SMOOTH OUT fMR SERIES WITH GAUSSIAN SMOOTHING KERNEL IF FWHM IS	*/
	/* NON-ZERO.  IF FWHM==0, THEN NO SMOOTHING IS TO BE DONE.		*/
	/* Perform convolution in time domain.  Could have done it in frequency	*/
	/* domain, but that would require doing FFT's.  In his paper, Friston	*/
	/* states that this is done to suppress thermal noise.			*/

	    if (FWHM!=0.) {
		for (dat_index=0; dat_index<num_dat_pts; dat_index++) {
		    smooth_series[dat_index]=0.;
		    for (dat_index2=0; dat_index2<num_dat_pts; dat_index2++) {
			kernel_index = dat_index2 - dat_index;
			if (kernel_index<0) kernel_index += num_dat_pts;
			smooth_series[dat_index] +=
		    		imgin[vox_index][dat_index2] * gaussian[kernel_index];
		    }
		}
	    }
	    else {
		for (dat_index=0; dat_index<num_dat_pts; dat_index++) {
		    smooth_series[dat_index]=imgin[vox_index][dat_index];
		    if (isnan(smooth_series[dat_index]))
			printf("smooth_series[%d] = imgin[%d] = NaN!\n",
							dat_index,dat_index);
		}
	    }
/*	    printf("Time series data:\n");
	    for (dat_index=0; dat_index<num_dat_pts; dat_index++)
		printf("%f\n",smooth_series[dat_index]);
*/

	/* -------------------------------------------------------------------- */
	/* CALCULATE Var{x(t)}.							*/

	imgin_mean=0.;
	for (dat_index=0; dat_index<num_dat_pts; dat_index++)
	    imgin_mean+=smooth_series[dat_index];
	imgin_mean/=(float)num_dat_pts;
	for (dat_index=0; dat_index<num_dat_pts; dat_index++)
	    smooth_series[dat_index]-=imgin_mean;
	sum_imgin_mean+=imgin_mean;
	imgin_var=0.;
	for (dat_index=0; dat_index<num_dat_pts; dat_index++)
	    imgin_var+=(smooth_series[dat_index]*smooth_series[dat_index]);
	imgin_var/=(float)(num_dat_pts-1);
	sum_imgin_var+=imgin_var;

	/* -------------------------------------------------------------------- */
	/* CALCULATE dx_dt[*], VECTOR OF FIRST DIFFERENCES.			*/

	    for (dat_index=0; dat_index<num_dat_pts-1; dat_index++)
		dx_dt[dat_index]=smooth_series[dat_index+1]-smooth_series[dat_index];

	/* -------------------------------------------------------------------- */
	/* CALCULATE Var{dx/dt}.						*/

	    mean_dir=0.;
	    for (dat_index=0; dat_index<num_dat_pts-1; dat_index++)
		mean_dir+=dx_dt[dat_index];
	    mean_dir/=(float)(num_dat_pts-1);
	    var_dir=0.;
	    for (dat_index=0; dat_index<num_dat_pts-1; dat_index++)
		var_dir+=(dx_dt[dat_index]-mean_dir)*(dx_dt[dat_index]-mean_dir);
	    var_dir/=(num_dat_pts-2);
	    sum_var_dir+=var_dir;

	/* -------------------------------------------------------------------- */
	/* CALCULATE W, EFFECTIVE FWHM, AND EFFECTIVE df.			*/

	    if (var_dir!=0.) {
		W = sqrt(imgin_var)/(float)sqrt((double)var_dir);
		effective_FWHM = W * sqrt(4.*log(2.));
		v = ((float)num_dat_pts)/effective_FWHM;
		sum_W+=W;
		sum_effective_FWHM+=effective_FWHM;
		sum_v+=v;
		num_nonzero_vardir++;
	    }
	}
    }


    /* ------------------------------------------------------------------------ */
    /* CALCULATE MEAN OF VARIANCES OF DERIVATIVES ACROSS VOXELS.		*/
    /* ALSO MEAN OF VARIANCES OF TIME SERIES ACROSS VOXELS.			*/

    printf("\nNumber of voxels for which effective df was calculated = %d\n",
							num_nonzero_vardir);
    mean_imgin_mean = sum_imgin_mean/(float)num_nonzero_vardir;;
    mean_imgin_var = sum_imgin_var/(float)num_nonzero_vardir;
    mean_var_dir = sum_var_dir/(float)num_nonzero_vardir;
    printf("\nMean value of E{x(t)}    = E{E{x(t)}}     = %.2f\n",mean_imgin_mean);
    printf("Mean value of Var{x(t)}  = E{Var{x(t)}}   = %.2f\n",mean_imgin_var);
    printf("Mean value of Var{dx/dt} = E{Var{dx/dt}}  = %.2f\n",mean_var_dir);


/* ---------------------------------------------------------------------------- */
/* CALCULATE EFFECTIVE DEGREES OF FREEDOM.					*/

    mean_W=sum_W/(float)num_nonzero_vardir;
    mean_effective_FWHM=sum_effective_FWHM/(float)num_nonzero_vardir;;
    mean_v=sum_v/(float)num_nonzero_vardir;
    printf("\nMEAN W = sqrt(E{Var{x(t)}}/sqrt(E{Var{dx/dt}}) = %.2f\n",mean_W);
    printf("MEAN Effective FWHM = W * sqrt(4*ln(2))        = %.2f scans\n",mean_effective_FWHM);
    printf("\nMEAN Effective degrees of freedom              = %.2f\n",mean_v);


/* ---------------------------------------------------------------------------- */
/*  FREE MEMORY.								*/

    printf("fmri_df: releasing memory...\n");
    free(header_name);
    free(img_name);
    free(temp);
    free(imgin);
    free(gaussian);
    free(smooth_series);
    free(dx_dt);
    free(mask);

}
