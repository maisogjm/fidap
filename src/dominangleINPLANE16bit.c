#include <stdio.h>	    /* For writing to files and	standard output.    */
#include <stdlib.h>	    /* For allocating memory with malloc.	    */
#include <time.h>	    /* For determining elapsed time.		    */
#include <math.h>	    /* For use of floor	function in generating	    */
			    /* the final corrected image.		    */
#include "analyze6.h"       /* For reading header information.		    */

		       /* EXTERNAL VARIABLES				      */
		       /* Not explicitly passed to the function F(*parm).     */
		       /* This was done to preserve the generality of the     */
		       /* search algorithm search(), which calls F.	      */
unsigned short int *mat1; /* Input image1 TO WHICH input image2 will be aligned. */
unsigned short int *mat2; /* Input image2 to align to input image1.		 */
unsigned short int *newmat; /* Result of aligning image2 to image1.		 */
int   *inewmat;	       /* Integer copy of *newmat			      */
int    xdim,ydim,zdim; /* X, Y, & Z dimensions of BOTH image1 and image2.     */
double xpixdim;        /* X pixel dimension size.			      */
double ypixdim;        /* Y pixel dimension size.			      */
double zpixdim;        /* Z pixel dimension size.			      */
double *finalparam;     /* Pointer to seven movement parameters.		      */
double zscalefactor;   /* Multiplicative scale factor for z pixel dimension.  */
int    step;	       /* Undersampling factor.  Applied to X & Y axes, so    */
		       /* step of 3 yields undersampling of 3*3=9	      */
double *rotpt;	       /* Point about which image2 is rotated.  In this       */
		       /* implementation, this is just the center of the      */
		       /* image.					      */
int xlo,xhi,ylo,yhi;   /* Boundaries of	portion	of image to be fitted.	The   */
		       /* whole	image is NOT fitted, because large portions   */
		       /* of the image contain little or no information; it   */
		       /* is hoped that	this will save time.		      */

/*int   *history;*/	       /* Image of history voxels.  Keeps track of which */
		       /* voxels in *newmat were fully sampled in *mat2.      */
/*int   *newhist; */	       /* Result of moving *history.		      */
/* int     tally_history(int xdim, int ydim, int zdim, double zscalefactor,
		int xlo, int xhi, int ylo, int yhi, double *param, double *rotpt);
*/


double   F(double *parm);
void    centergrav16bit(unsigned short int *mat, int xdim, int ydim, int zdim, double *center);
void    find_lims16bit(unsigned short int *mat, int xdim, int ydim, int zdim,
		int *xlo, int *xhi, int *ylo, int *yhi);
int     amoebaD(double **p,double *y,int ndim,double ftol,double (*F) (double *),
								    int *nfunk);
void    moveimg16bitD(unsigned short int *mat, int xdim, int ydim, int zdim, double zscalefactor,
		  double *param, double *rotpt, unsigned short *newmat);
double   dotprod16bit(unsigned short int *mat1, unsigned short int *mat2, int xdim, int ydim, int zdim, int step,
	       int xlo, int xhi, int ylo, int yhi);

struct  analyze_struct readhdr(char *argv[]);
double **matrix(int a, int b, int c, int d);
double  *vector(int a, int b);
void    free_vector(double *v,int nl);
void    free_matrix(double **m, int nrl, int nrh, int ncl, int nch);

/* **************************************************************************** */
/*	This program calls C subroutines and registers two images.		*/
/*	In this	implementation,	zero crossovers	are maximized; the		*/
/*	Nelder-Mead Simplex search algorithm is	used to	search the		*/
/*	parameter space.  Of note are:						*/
/*										*/
/*		(1) The	use of the distance between the	centers	of		*/
/*		    gravity of the two images as the initial estimate		*/
/*		    of the translation artifact					*/
/*		(2) The	two-step search	employed.  First the search is		*/
/*		    done with a	set of 8 initial parameter settings,		*/
/*		    and	with an	initial	undersampling factor.  Then		*/
/*		    it is done with another set	of 8 initial parameter		*/
/*		    settings,	one being the best parameter settings		*/
/*		    found in the first search, with no undersampling.		*/
/*		(3) The	use of a trilinear interpolation routine to		*/
/*		    generate a translated/rotated image.			*/
/*		(4) The	function F, included in	this file.  This is		*/
/*		    defined as the negative number of crossovers;		*/
/*		    minimizing this function will maximize the zero		*/
/*		    crossovers (most search algorithms seem to be		*/
/*		    written in terms of	minimizing a function).			*/
/*		(5) Use	of the 1-offset	matrix/vector memory allocaters		*/
/*		    matrix() and vector(), for compatibility with code for	*/
/*		    the	Nelder-Mead Simplex method.  This makes	dereferencing	*/
/*		    of pointer *param in moveimg7.c a little trickier.		*/
/*		(6) Although a multiplicative scale factor (parameter 7) is	*/
/*		    also optimized, the	final output image is generated	from	*/
/*		    image 2 with a scaling factor of 1.0.  If the user wishes	*/
/*		    to rescale this image, he/she should multiply each voxel	*/
/*		    image being	registered by the scaling factor.		*/
/*										*/
/*	COMPILE: gcc dominangleINPLANE16bit.c imgio.c linalgebra.c	moveimg16bitD.c		*/
/*		 tally_history.c amoebaD.c dotprod16bit.c centergrav16bit.c		*/
/*		 find_lims16bit.c -o	/afs/nih/nia/joem/bin/dominangleINPLANE16bit -lm -O	*/


main(int argc, char *argv[])
{
	struct analyze_struct img_info;	/* Header information variable		*/
	FILE *fp;		/* Pointer for output text file			*/
	FILE *fimgin;		/* Pointer for input image files.		*/
	FILE *fimgout;		/* Pointer for output image files.		*/

	double **param;		/* Pointer to parameters			*/
        double tol;             /* Tolerance for search algorithm		*/
	int numparam,i,j;	/* Number of parameters and two index variables */
	int num_vox;		/* Number of voxels in an image.		*/
	double *y;		/* Function F evaluated	at each	point in	*/
				/* the simplex					*/
        int *nfunk;		/* Number of times F is evaluated		*/
	double *center;		/* Center of gravity temporary variable		*/
	double cx1,cy1,cz1;	/* Center of gravity of image 1			*/
	double cx2,cy2,cz2;	/* Center of gravity of image 2			*/
	int lowest;		/* Index to entry in **param which has lowest	*/
				/* cost function value				*/
	int vol_restrict;	/* Flag set to nonzero value if it is desired   */
				/* to restrict the area of the image resampled  */
				/* to only that portion encompassing all voxels */
				/* of intensity 1/3 or greater the maximum      */
				/* intensity value                              */


	time_t timeval;		/* Numeric representation for time		*/
	time_t timeval2;	/* Numeric representation for time		*/
	char timestr[26];	/* String representation for time		*/


/* ---------------------------------------------------------------------------- */
/*	TIME STAMP								*/

	time(&timeval);
	printf("Time stamp:\n");
	system("date");

/* ---------------------------------------------------------------------------- */
/*	CHECK TO SEE IF	CORRECT	NUMBER OF ARGUMENTS WERE PASSED.		*/

	if (argc < 9) {
	    printf("dominangleINPLANE16bit: need input headerfile, input imagefile1,\n");
	    printf("            input imagefile2, output imagefile names,\n");
	    printf("            output textfile name, tol,\n");
	    printf("            sampling factor, volume restriction flag.\n");
	    exit(1);
	}

/* ---------------------------------------------------------------------------- */
/*	OPEN OUTPUT TEXT FILE.							*/

	if ((fp = fopen(*(argv+5),"w")) == NULL) {
	    printf ("dominangleINPLANE16bit: can't open %s\n",*(argv+5));
	    exit(4);
	}

/* ---------------------------------------------------------------------------- */
/*	READ HEADER FILE; ASSUME IMAGES	ARE SAME SIZE, SO USE SAME HEADER FILE 	*/
/*	FOR BOTH IMAGES.  CHECK TO MAKE SURE PIXEL DIMENSIONS LOOK OKAY.	*/
/*	CALCULATE zscalefactor, A MULTIPLICATE SCALE FACTOR WHICH WILL BE	*/
/*	APPLIED TO THE Z COORDINATE TO ALLOW RESAMPLING WITH NON-CUBIC VOXELS.	*/

	img_info = readhdr(argv+1);

	if (img_info.hk.session_error != 0) {
	    printf("dominangleINPLANE16bit: error reading header.\n");
	    exit(2);
	}

	xdim = img_info.dime.dim[1];
	ydim = img_info.dime.dim[2];
	zdim = img_info.dime.dim[3];
	num_vox=xdim*ydim*zdim;

	xpixdim = img_info.dime.pixdim[1];
	ypixdim = img_info.dime.pixdim[2];
	zpixdim = img_info.dime.pixdim[3];

	printf("X pixel dimension = %f\n",xpixdim);
	printf("Y pixel dimension = %f\n",ypixdim);
	printf("Z pixel dimension = %f\n",zpixdim);
	fprintf(fp,"X pixel dimension = %f\n",xpixdim);
	fprintf(fp,"Y pixel dimension = %f\n",ypixdim);
	fprintf(fp,"Z pixel dimension = %f\n",zpixdim);

	if ((xpixdim<=0) || (ypixdim<=0) || (zpixdim<=0)
	    || (xpixdim!=ypixdim)) {
	    printf("dominangleINPLANE16bit: something wrong with pixel dimensions.\n");
	    printf("            Perhaps xpixdim is not equal to ypixdim.\n");
	    exit(3);
	}

	zscalefactor = (double) zpixdim / (double) xpixdim;
	printf("Z scale factor = zpixdim/xpixdim = %f\n",zscalefactor);
	fprintf(fp,"Z scale factor = zpixdim/xpixdim = %f\n",zscalefactor);

/* ---------------------------------------------------------------------------- */
/*	SET PARAMETERS,	ALLOCATE MEMORY	FOR POINTERS.  Each integer		*/
/*	datum takes up FOUR bytes.  Each doubling point datum takes up		*/
/*	FOUR bytes.								*/
/*										*/
/*	There are 7 degrees of freedom:	1 for each translation and one		*/
/*	for each rotation about	the three axes,	and pixel intensity		*/
/*	scaling.								*/
/*										*/

	mat1 = (unsigned short int *) malloc(num_vox*sizeof(unsigned short int));
	mat2 = (unsigned short int *) malloc(num_vox*sizeof(unsigned short int));

	numparam = 3;	/* Number of parameters/degrees	of freedom.		*/

	newmat = (unsigned short int *)  malloc(sizeof(unsigned short int)*num_vox);
	nfunk = (int *) malloc(4);
	param  = matrix(1,4,1,3);
	finalparam = vector(1,7);
	y = vector(1,4);
	rotpt  = (double *) malloc(4*3);
	center = (double *) malloc(4*3);

	tol    = atof(*(argv+6));
	step   = atoi(*(argv+7));

        printf("\nRtol = %f\n",tol);
        printf("Sampling factor = %d\n\n",step);
	printf("Will register %s to %s\n",*(argv+3),*(argv+2));
	fprintf(fp,"Will register %s to %s\n",*(argv+3),*(argv+2));

/* ---------------------------------------------------------------------------- */
/*	READ IN	IMAGES.								*/

	if ((fimgin = fopen(argv[2],"r")) == NULL) {
	    printf ("dominangleINPLANE16bit: can't open %s\n",*(argv+2));
	    exit(2);
	}
	fread(mat1,sizeof(signed short),num_vox,fimgin);
	fclose(fimgin);

	if ((fimgin = fopen(argv[3],"r")) == NULL) {
	    printf ("dominangleINPLANE16bit: can't open %s\n",*(argv+3));
	    exit(2);
	}
	fread(mat2,sizeof(signed short),num_vox,fimgin);
	fclose(fimgin);

/* ---------------------------------------------------------------------------- */
/*	MAKE HISTORY IMAGE.							*/

/*	history = (int *) malloc(4*xdim*ydim*zdim);
	newhist = (int *)  malloc(4*xdim*ydim*zdim);

	for (i=0; i<=xdim*ydim*zdim-1; i++)
	    *(history+i)=1;
*/
/* ---------------------------------------------------------------------------- */
/*	FIND CENTERS OF	GRAVITY	FOR EACH INPUT IMAGE.				*/
/*	Set initial estimates for translation parameters according		*/
/*	to the difference between the centers of gravity.			*/
/*	Set the	rotation point to be the center	of the image.			*/
/*	The z-coordinate of the rotation point will NOT have to be rescaled	*/
/*	with zscalefactor (calculated above) because it is used in the NON-	*/
/*	isometric voxel space.							*/

	printf("dominangleINPLANE16bit: finding centers of gravity....\n");
	centergrav16bit(mat1, xdim, ydim, zdim, center);
	cx1 = *(center+0);
	cy1 = *(center+1);
	cz1 = *(center+2);

	centergrav16bit(mat2, xdim, ydim, zdim, center);
	cx2 = *(center+0);
	cy2 = *(center+1);
	cz2 = *(center+2);

        *(rotpt+0) = (double) (xdim-1)/2.;
        *(rotpt+1) = (double) (ydim-1)/2.;
        *(rotpt+2) = (double) (zdim-1)/2.;

	printf("Delta xyz for centers of gravity = %f %f %f\n",
	       cx1-cx2, cy1-cy2, cz1-cz2);
	fprintf(fp,"Delta xyz for centers of gravity = %f %f %f\n",
	       cx1-cx2, cy1-cy2, cz1-cz2);

/* ---------------------------------------------------------------------------- */
/*	FIND xlo, xhi, ylo, yhi.						*/
/*	The subroutine	find_lims16bit will find the	lowest and highest x and y	*/
/*	coordinates for points which are one third the value of the maximum	*/
/*	pixel value or greater.  These limits will constrain the moveimg16bitD and	*/
/*	dotprod16bit subroutines to only a	portion	of each	slice, thus hopefully	*/
/*	saving time.								*/

	vol_restrict = atoi(*(argv+8));

	if (vol_restrict == 0) {
	    xlo = 0;
	    xhi = xdim-1;
	    ylo = 0;
	    yhi = ydim-1;
	} else
	    find_lims16bit(mat1,xdim,ydim,zdim,&xlo,&xhi,&ylo,&yhi);

	printf("\nBoundaries of image actually processed:\n");
	printf("xlo = %d\n",xlo);
	printf("xhi = %d\n",xhi);
	printf("ylo = %d\n",ylo);
	printf("yhi = %d\n",yhi);

	fprintf(fp,"\nBoundaries of image actually processed:\n");
	fprintf(fp,"xlo = %d\n",xlo);
	fprintf(fp,"xhi = %d\n",xhi);
	fprintf(fp,"ylo = %d\n",ylo);
	fprintf(fp,"yhi = %d\n",yhi);

/* ---------------------------------------------------------------------------- */
/*	DO INITIAL SEARCH.							*/

	/* ---------------------------------------------------- */
	/* PUT IN INITIAL POINTS IN param[][].			*/

	    param[1][1]= cx1 - cx2 + (0.1 * ((double) xdim));
	    param[1][2]= cy1 - cy2;
	    param[1][3]= 0.;

	    param[2][1]= cx1 - cx2;
	    param[2][2]= cy1 - cy2 + (0.1 * ((double) ydim));
	    param[2][3]= 0.;

	    param[3][1]= cx1 - cx2;
	    param[3][2]= cy1 - cy2;
	    param[3][3]= 5.;

	    param[4][1]= cx1 - cx2;
	    param[4][2]= cy1 - cy2;
	    param[4][3]= 0.;

/*	    for (i=1;i<=4;i++) {
		for (j=1;j<=3;j++)
		    printf("%8.4f ",param[i][j]);
		printf("\n");
	    }
*/

	/* ---------------------------------------------------- */
	/* PUT IN INITIAL FUNCTION VALUES IN y[].		*/

	printf("\ndominangleINPLANE16bit: making y[] matrix...\n");
	printf("\nInitial vertices of simplex, and corresponding cost function values:\n");
	printf(" X-Trans  Y-Trans   Z-Rot    Angle\n");
	for (i=1;i<=4;i++) {
	    y[i]=F(param[i]);
	    for (j=1;j<=3;j++)
		printf("%8.4f ",param[i][j]);
	    printf(": %8.4f\n",y[i]);
	}


	/* ---------------------------------------------------- */
	/* CALL	AMOEBA.						*/

	printf("\ndominangleINPLANE16bit: doing initial search (may take a while)....\n");

        lowest=amoebaD(param,y,numparam,tol,(double (*)(double *)) F,nfunk);

        printf("dominangleINPLANE16bit: finished initial search.\n");

/* -------------------------------------------------------------------- */
/*	SHOW INTERMEDIATE RESULTS.					*/

        printf("\nIntermediate results:\n");
        printf("Rtol = %f\n",tol);
        printf("Sampling factor = %d\n",step);
        printf("  X Trans = %f\n",param[lowest][1]);
        printf("  Y Trans = %f\n",param[lowest][2]);
        printf("  Z Rot   = %f\n",param[lowest][3]);
	printf("Angle = %f\n",(int) y[lowest]);

        fprintf(fp,"\nIntermediate results:\n");
        fprintf(fp,"Rtol = %f\n",tol);
        fprintf(fp,"Sampling factor = %d\n",step);
        fprintf(fp,"  X Trans = %f\n",param[lowest][1]);
        fprintf(fp,"  Y Trans = %f\n",param[lowest][2]);
        fprintf(fp,"  Z Rot   = %f\n",param[lowest][3]);
	fprintf(fp,"Angle = %f\n",(int) y[lowest]);

	strcpy(timestr, ctime(&timeval));
	timestr[19] = '\0';
	printf("Thus far, program has run from %s to", timestr+11);
	fprintf(fp,"Thus far, program has run from %s to", timestr+11);

	time(&timeval2);
	strcpy(timestr, ctime(&timeval2));
	timestr[19] = '\0';
	printf(" %s\n", timestr+11);
	fprintf(fp," %s\n", timestr+11);

/* ---------------------------------------------------------------------------- */
/*	DO SECOND SEARCH.							*/
/*	Do second search using best fit	from first search as a starting	point.	*/

	/* ---------------------------------------------------- */
	/* Set sampling	factor to 1.				*/

	step=1;


	/* ---------------------------------------------------- */
	/* PUT IN INITIAL POINTS IN param[][].			*/
	/* Set param[4]	to the param[lowest] resulting from	*/
	/* initial search.  The other vertices will differ from */
	/* param[lowest] by some small amount, e.g. 1, at each	*/
	/* degree of freedom.					*/

	for (i=1;i<=4;i++)
	    for (j=1;j<=3;j++)
		param[i][j]=param[lowest][j];

	param[1][1]-=(0.1 * ((double) xdim));
	param[2][2]-=(0.1 * ((double) ydim));
	param[3][3]-=5.;

	/* ---------------------------------------------------- */
	/* PUT IN INITIAL FUNCTION VALUES IN y[].		*/

	printf("\ndominangleINPLANE16bit: reinitializing simplex for second search (may take a while)...\n");
	printf("             (will use FULL sampling for second search.)\n");
	printf("\ndominangleINPLANE16bit: making y[] matrix...\n");
	printf("\nInitial vertices of simplex, and corresponding cost function values:\n");
	printf(" X-Trans  Y-Trans  Z-Rot    Angle\n");
	for (i=1;i<=4;i++) {
	    y[i]=F(param[i]);
	    for (j=1;j<=3;j++)
		printf("%8.4f ",param[i][j]);
	    printf(": %8.4f\n",y[i]);
	}


	/* ---------------------------------------------------- */
	/* CALL	AMOEBA.						*/

        printf("\ndominangleINPLANE16bit: doing second search (may take a while)....\n");

        lowest=amoebaD(param,y,numparam,tol,(double (*)(double *)) F,nfunk);

        printf("dominangleINPLANE16bit: finished second search.\n");

/* -------------------------------------------------------------------- */
/*	SHOW FINAL RESULTS.						*/

        printf("\nFinal results:\n");
        printf("Rtol = %f\n",tol);
        printf("Sampling factor = %d\n",step);
        printf("  X Trans = %f\n",param[lowest][1]);
        printf("  Y Trans = %f\n",param[lowest][2]);
        printf("  Z Rot   = %f\n",param[lowest][3]);
	printf("Angle = %f\n",(int) y[lowest]);

        fprintf(fp,"\nFinal results:\n");
        fprintf(fp,"Rtol = %f\n",tol);
        fprintf(fp,"Sampling factor = %d\n",step);
        fprintf(fp,"  X Trans = %f\n",param[lowest][1]);
        fprintf(fp,"  Y Trans = %f\n",param[lowest][2]);
        fprintf(fp,"  Z Rot   = %f\n",param[lowest][3]);
	fprintf(fp,"Angle = %f\n",(int) y[lowest]);

/* ---------------------------------------------------------------------------- */
/*	GENERATE CORRECTED IMAGE AND WRITE TO DISK.				*/

	printf("\ndominangleINPLANE16bit: resampling image %s...\n",*(argv+3));
        finalparam[1]=param[lowest][1];
        finalparam[2]=param[lowest][2];
        finalparam[3]=0.;
        finalparam[4]=0.;
        finalparam[5]=0.;
        finalparam[6]=param[lowest][3];
        finalparam[7]=1.;

	moveimg16bitD(mat2, xdim, ydim, zdim, zscalefactor, finalparam, rotpt, newmat);

	printf("Writing image to %s.\n",*(argv+4));
	fprintf(fp,"Writing image to %s.\n",*(argv+4));
	if ((fimgout = fopen(argv[4],"w")) == NULL) {
	    printf("dominangleINPLANE16bit: error writing to %s.\n",argv[4]);
	    exit(7);
	}
	fwrite(newmat,sizeof(unsigned short int),num_vox,fimgout);
	fclose(fimgout);

/* ---------------------------------------------------------------------------- */
/*	Print out time elapsed.							*/

	strcpy(timestr, ctime(&timeval));
	timestr[19] = '\0';
	printf("Program ran from %s to", timestr+11);
	fprintf(fp,"Program ran from %s to", timestr+11);

	time(&timeval);
	strcpy(timestr, ctime(&timeval));
	timestr[19] = '\0';
	printf(" %s\n", timestr+11);
	fprintf(fp," %s\n", timestr+11);

	fclose(fp);

/* ---------------------------------------------------------------------------- */
/*	TIME STAMP								*/

	printf("Time stamp:\n");
	system("date");

/* ---------------------------------------------------------------------------- */
/*	DEALLOCATE MEMORY							*/
	printf("dominangleINPLANE16bit: deallocating memory...\n");
	free(mat1);
	free(mat2);
	free(newmat);
	free(rotpt);
	free(center);
	free(nfunk);
	free_vector(y,1);
	free_vector(finalparam,1);
	free_matrix(param,1,8,1,7);

	exit(0);
}

/* **************************************************************************** */
double F(double *parm)
{
	int count;	/* Number of fully sampled voxels in *newhist.		*/
	int index;

	/* ---------------------------------------------------- */
	/* Move	image 2, then calculate	ANGLE resultant image	*/
	/* and image 1.						*/
    
	finalparam[1]=parm[1];
	finalparam[2]=parm[2];
	finalparam[3]=0.;
	finalparam[4]=0.;
	finalparam[5]=0.;
	finalparam[6]=parm[3];
	finalparam[7]=1.;

/*	printf("Moving image with these parameters:\n");
	for (index=1; index<=7; index++)
	    printf("param[%d] = %f\n",index,finalparam[index]);
*/
	moveimg16bitD(mat2, xdim, ydim, zdim, zscalefactor, finalparam, rotpt, newmat);

/*	count=tally_history(xdim, ydim, zdim, zscalefactor, xlo, xhi, ylo, yhi,
								parm, rotpt);
	if (count==0)
	    count=1;

	printf("count = %d\n",count);
*/
	return (double) acos(dotprod16bit(mat1,newmat,xdim,ydim,zdim,step,xlo,xhi,ylo,yhi));

}
