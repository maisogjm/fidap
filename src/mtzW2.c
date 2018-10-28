/*
 * mtz_noscal - A C-language program modified from the program testls.c
 *	 (written by J.N. Little) which came in the MATLAB toolbox.
 *       mtz reads in a MATLAB .mat file and writes out an
 *	 ANALYZE format .img file.  Formerly named mtz.c
 *
 * MODIFIED DEC-11-1991 : JMM : Reads in MAT-files, outputs ANALYZE format
 *				.img files.
 * MODIFIED APR-17-1992 : JMM : Modified to do NO rescaling; assumes that
 *                              data in .mat file is already scaled from
 *                              0 to 255.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int loadmat2(FILE *fp, long int *type, char *pname, long int *mrows,
	long int *ncols,
	long int *imagf,
	double **preal,
	double **pimag);
double max(double *spmtmat, int xdim, int ydim, int zdim);
double min(double *spmtmat, int xdim, int ydim, int zdim);
main(int argc, char *argv[])
{
	FILE *fp;/*, *fopen();*/
	char name[20];
	long int type, mrows, ncols, imagf;
	double *xr, *xi;
	int xdim,ydim,zdim,loop;
	signed short int *wordp,Wmin,Wmax;
	double maxval,minval;
	int x,y,z,mx,my;
	double scalefactor;

/*
 * Take the "b" out of the fopen on systems that do not do newline
 * translation, for example Unix and Macintosh systems.  Leave
 * the "b" in for MS-DOS and VAX/VMS C compilers.  The "rfm=var"
 * is needed for VAX/VMS C compilers to make the file type the
 * same as when MATLAB saves a .MAT file (it does this for
 * performance reasons).
 */

	/* CHECK COMMAND LINE ARGUMENTS					*/
	if (argc < 7) {
	    printf("mtzW2: Need name of input .mat file,\n");
	    printf("      name of output .img file,\n");
	    printf("      xdim, ydim, zdim,\n");
	    printf("      AND scalefactor.\n");

	    exit(1);
	}


	/* OPEN DATA FILE						*/
        if ((fp = fopen(*(argv+1),"r")) == NULL)
        {
                printf ("mtzW2: can't open %s\n",*(argv+1));
                exit (2);
        }


	/* LOAD DATA FROM DISK INTO DOUBLE ARRAY			*/
	if (loadmat2(fp, &type, name, &mrows, &ncols, &imagf, &xr, &xi)) {
		printf("\nmtzW2: read error\n");
		exit(3);
	}
	fclose(fp);


	/* SHOW DIMENSIONS						*/
	printf("type  = %d\n",type);
	printf("mrows = %d\n",mrows);
	printf("ncols = %d\n",ncols);

	xdim = atoi(*(argv+3));
	ydim = atoi(*(argv+4));
	zdim = atoi(*(argv+5));
	printf("xdim = %d\n",xdim);
	printf("ydim = %d\n",ydim);
	printf("zdim = %d\n",zdim);
        scalefactor = (double) atof(argv[6]);
        printf("scalefactor = %f\n",scalefactor);

	printf("\nmrows * ncols = %d\n",mrows*ncols);
	printf("\nxdim * ydim * zdim = %d\n",xdim*ydim*zdim);


	/* CHECK DIMENSIONS						*/
	if (mrows*ncols != xdim*ydim*zdim) {
	    printf("Error in dimensions; mrows * ncols does not equal\n");
	    printf("xdim * ydim * zdim!  Cannot do mtz.\n");
	    exit(4);
	}

	minval = min(xr, xdim, ydim, zdim);
	maxval = max(xr, xdim, ydim, zdim);
	printf("Floating point min is %f.\n",minval);
	printf("Floating point max is %f.\n",maxval);

        /* TRUNCATE DATA, IF IT EXCEEDS MAX OR FALLS BELOW MIN. */
        for (loop = 0; loop <= xdim*ydim*zdim-1; loop++) {
            if (xr[loop]>32767.) xr[loop]=32767.;
            if (xr[loop]<-32768.) xr[loop]=-32768.;
        }

	/* TRANSFER SCALED DATA TO 16-BIT ARRAY			*/
	wordp = (signed short int *)
			malloc(sizeof(signed short int)*xdim*ydim*zdim);
	for (z = 0; z <= zdim-1; z++)
	    for (y = 0; y <= ydim-1; y++)
		for (x = 0; x <= xdim-1; x++) {
		    mx = x + z*xdim;
		    my = y;
		    *(wordp + (z*ydim+y)*xdim +x) =
		(signed short int) floor((*(xr + mx + my*xdim*zdim)*scalefactor)+0.5);
	}

        /*
	Wmin=32767;
	Wmax=0;
	for (loop = 0; loop <= xdim*ydim*zdim-1; loop++) {
	    Wmin=(wordp[loop] < Wmin) ? wordp[loop] : Wmin;
	    Wmax=(wordp[loop] > Wmax) ? wordp[loop] : Wmax;
	}
	printf("Pre-truncation Maximum of wordp is %d.\n",Wmax);
	printf("Pre-truncation Minimum of wordp is %d.\n",Wmin);
	*/

	Wmin=32767;
	Wmax=0;
	for (loop = 0; loop <= xdim*ydim*zdim-1; loop++) {
	    Wmin=(wordp[loop] < Wmin) ? wordp[loop] : Wmin;
	    Wmax=(wordp[loop] > Wmax) ? wordp[loop] : Wmax;
	}
	printf("Maximum of wordp is %d.\n",Wmax);
	printf("Minimum of wordp is %d.\n",Wmin);

	/* WRITE INTEGER ARRAY TO DISK					*/
        if ((fp = fopen(argv[2],"w")) == NULL) {
            printf("mtzW2: error writing to %s.\n",argv[2]);
            exit(7);
        }
        fwrite(wordp,sizeof(signed short int),xdim*ydim*zdim,fp);
        fclose(fp);

	free(xr);
	free(wordp);
	if (imagf) {
		free(xi);
	}
}

 
/* ************************************************************************** */
int writimg2(char *argv[], int *mat, int xdim, int ydim, int zdim)

/*	Writes image to disk given filename in *argv[], data to		*/
/*	write starting at address pointed to by mat, and dimensions	*/
/*	xdim, ydim, and zdim.						*/
/*	DEC-11-91: JMM : Modified from existing code for writimg.c to   */
/*	handle conversion from MATLAB matrix format to ANALYZE 3 coor-	*/
/*	dinate format.							*/
{
	FILE *fp;
	int x,y,z;
	int mx,my;

	if ((fp=fopen(*argv,"w"))==0) {
	    printf ("writeimg2: unable to create: %s\n",*argv);
	    return 1;
	}

	printf("writimg2: writing image to %s...\n",*argv);
	for (z = 0; z <= zdim-1; z++)
	for (y = 0; y <= ydim-1; y++)
	for (x = 0; x <= xdim-1; x++) {
		    mx = x + z*xdim;
		    my = y;
 		    fputc((char) *(mat + mx + my*xdim*zdim), fp);
	}
	fclose(fp);
	return 0;
}


/* ************************************************************************** */
double max(double *spmtmat, int xdim, int ydim, int zdim)
/* FINDS MAXIMUM OF MATRIX spmtmat[][][]				      */
{
	int x,y,z;
	double max;

	max = *spmtmat;

	for (x = 0; x <= xdim-1; x++)
	    for (y = 0; y <= ydim-1; y++)
		for (z = 0; z <= zdim-1; z++)
		    max   = (*(spmtmat+(x*ydim+y)*zdim+z) > max) ?
			     *(spmtmat+(x*ydim+y)*zdim+z) : max;

	return max;
}


/* ************************************************************************** */
double min(double *spmtmat, int xdim, int ydim, int zdim)
/* FINDS MINIMUM OF MATRIX spmtmat[][][]				      */
{
        int x,y,z;
        double min;

        min = *spmtmat;

        for (x = 0; x <= xdim-1; x++)
            for (y = 0; y <= ydim-1; y++)
                for (z = 0; z <= zdim-1; z++)
                    min   = (*(spmtmat+(x*ydim+y)*zdim+z) < min) ?
                             *(spmtmat+(x*ydim+y)*zdim+z) : min;

        return min;
}
