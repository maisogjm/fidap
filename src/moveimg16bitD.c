#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void   mkrotmat(double x_theta, double y_theta, double z_theta,
     						double rotmat[3][3]);
void multmatvect(double mat[3][3], double in_vect[3], double outvect[3]);

/* **************************************************************************** */
void moveimg16bitD(unsigned short int *mat, int xdim, int ydim, int zdim, double zscalefactor,
				double *param, double *rotpt, unsigned short int *newmat)

/*    Takes as input a pointer mat to a matrix with dimensions xdim, ydim,	*/
/*    and zdim; translation, rotation, and scaling parameters in *param;	*/
/*    and the coordinates of the point in the image about which to do		*/
/*    rotations, in *rotpt; and a pointer to allocated space in which to	*/
/*    write results.  Returns the result of applying the transformation		*/
/*    specified by *param on *mat, in the space pointed to by newmat. Note	*/
/*    that, unlike readimg, moveimg does NOT allocate memory for the new	*/
/*    matrix (i.e., newmat); this must be done in the calling program with	*/
/*    the usual call to malloc.  (If moveimg WERE TO allocate memory as		*/
/*    readimg does, and if there were repeated calls to moveimg in the		*/
/*    calling program, then each call to moveimg after the first call would	*/
/*    have to be preceded by a call to the stdlib.h function free, to		*/
/*    deallocate memory.  If the programmer were to forget this, memory		*/
/*    might get used up, the program would stop execution, and a nonspec-	*/
/*    ific/noninformative error message such as "Segmentation fault" might	*/
/*    be generated, which may be potentially difficult to debug.)		*/

/*    07.02.92 : JMM : altered to contain code for interpolpix2.c, so that	*/
/*		       repeated calls to an outside function can be avoided	*/
/*		       (potentially speeding things up)				*/
/*    07.09.92 : JMM : changed pointer param offsets to reflect 1-offset	*/
/*		       data; to be used with domaxzero3a2.c, amoeba.c		*/
/*    07.15.92 : JMM : added input arguments step,xlo,xhi,ylo.yhi		*/
/*    07.23.92 : JMM : DELETED input arguments step,xlo,xhi,ylo,yhi to		*/
/*		       make code for generation of final registered image	*/
/*    09.02.92 : JMM : added code to take into account pixel dimensions		*/
/*		       when voxels are NOT isometric.  Thus, this moveimg	*/
/*		       will now work on 15-slice O-15 PET data without		*/
/*		       having to interpolate up to 46 slices to get cubic	*/
/*		       voxels.							*/
/*    12.16.92 : JMM : type cast argument of floor to double; for use with	*/
/*		       domoveimg3.c						*/
/*    03.11.93 : JMM : modified moveimg9.c to generate moveimg11.c		*/
/*		       moveimg11 uses an incremental technique to resample	*/
/*		       images, resulting in a savings in computation and	*/
/*		       increase in speed.  Also, rewrote code to handle new way	*/
/*		       of reading in images into memory, using fread.  Reading	*/
/*		       in images using fread places the voxels with x varying	*/
/*		       the fastest, y varying next as fast, and z varying	*/
/*		       slowest.  Also, changed all floating point numbers to	*/
/*		       doubles.  Thus, this code is to be used with		*/
/*		       linalgebraD.c, the double version of linalgebra.c	*/
/*		       This required a major rewriting of the code.		*/
/*		       Also, removed "sample" macro, inserting sampling		*/
/*		       formula directly into code.				*/ 
/*    05.12.93 : JMM : modified moveimg11.c to generate moveimg16bit.c		*/
/*		       code works on 16-bit data now.				*/
/*    07.20.93 : JMM : modified moveimg16bit.c to generate moveimg16bitD.c	*/
/*		       changed input argument param from float to double.	*/

{
    double  xtran,ytran,ztran;		/* Translations.			*/
    double x_theta,y_theta,z_theta;	/* Rotations.				*/
    double  scale;			/* Multiplicative scaling factor.	*/
    double  rotptx,rotpty,rotptz;	/* Point about which to rotate.		*/
    double inv_mat[3][3];		/* Inverse rotation matrix.		*/
    double  invect[3], outvect[3];	/* Input and output coordinates to be	*/
					/* transformed by inv_mat.		*/
    int    x,y,z;			/* Indices into output image array.	*/

    register double dxx, dyx, dzx;	/* Partial increments for unit x step.	*/
    register double dxy, dyy, dzy;	/* Partial increments for unit y step.  */
    register double dxz, dyz, dzz;	/* Partial increments for unit z step.  */

    register double basexx, baseyx, basezx;
    register double basexy, baseyy, basezy;
    register double basexz, baseyz, basezz;

    double u,v,w;			/* Indices into input image array.	*/
    double vol0, vol1, vol2, vol3,	/* Partial volumes.			*/
	vol4, vol5, vol6, vol7;
    double val0, val1, val2, val3,	/* Voxel values of partial volumes.	*/
	val4, val5, val6, val7;
    register double flu, flv, flw,	/* Floating point representations of	*/ 
    flup1, flvp1, flwp1;		/* u, v, w and u+1, v+1, w+1.		*/
 
    xtran   = (double) *(param+1);		/* NOTE ONE-OFFSET ARRAY.  DOES NOT	*/
    ytran   = (double) *(param+2);		/* BEGIN DEREFERENCING AT ZERO!		*/
    ztran   = (double) *(param+3);
    x_theta = (double) *(param+4);
    y_theta = (double) *(param+5);
    z_theta = (double) *(param+6);
    scale   = 1.;
    rotptx  = *(rotpt+0);
    rotpty  = *(rotpt+1);
    rotptz  = *(rotpt+2);
/*    printf("\nmoveimg16bit: moving image with these parameters:\n");
    printf("  X Trans = %f\n",xtran);
    printf("  Y Trans = %f\n",ytran);
    printf("  Z Trans = %f\n",ztran);
    printf("  X Rot   = %f\n",x_theta);
    printf("  Y Rot   = %f\n",y_theta);
    printf("  Z Rot   = %f\n",z_theta);
    printf("  Scale   = %f\n",scale);
*/
    mkrotmat(-x_theta,-y_theta,-z_theta,inv_mat);

    /* ********************************************************	*/
    /* Calculate dxx, dyx, dzx.					*/

	invect[0] = 1.;
	invect[1] = 0.;

	/* Map Z coordinate to isometric voxel space.		*/

	invect[2] = 0. * zscalefactor;
	multmatvect(inv_mat,invect,outvect);

	/* Now map moved z coordinate back to			*/
	/* NON-isometric voxel space.				*/

	outvect[2] /= zscalefactor;

	dxx = outvect[0];
	dyx = outvect[1];
	dzx = outvect[2];

    /* ********************************************************	*/
    /* Calculate dxy, dyy, dzy.					*/

	invect[0] = 0.;
	invect[1] = 1.;

	/* Map Z coordinate to isometric voxel space.		*/

	invect[2] = (0.) * zscalefactor;
	multmatvect(inv_mat,invect,outvect);
	/* Now map moved z coordinate back to			*/
	/* NON-isometric voxel space.				*/

	outvect[2] /= zscalefactor;

	dxy = outvect[0];
	dyy = outvect[1];
	dzy = outvect[2];

    /* ********************************************************	*/
    /* Calculate dxz, dyz, dzz.					*/

	invect[0] = 0.;
	invect[1] = 0.;

	/* Map Z coordinate to isometric voxel space.		*/

	invect[2] = (1.) * zscalefactor;
	multmatvect(inv_mat,invect,outvect);

	/* Now map moved z coordinate back to			*/
	/* NON-isometric voxel space.				*/

	outvect[2] /= zscalefactor;

	dxz = outvect[0];
	dyz = outvect[1];
	dzz = outvect[2];

    /* ********************************************************	*/
    /* Calculate where (-1,-1,-1) of output image maps to. 	*/

	invect[0] = (double) (-1.-rotptx);
	invect[1] = (double) (-1.-rotpty);

	/* Map Z coordinate to isometric voxel space.		*/

	invect[2] = ((double) (-1.-rotptz)) * zscalefactor;
	multmatvect(inv_mat,invect,outvect);

	/* Now map moved z coordinate back to			*/
	/* NON-isometric voxel space.				*/

	outvect[2] /= zscalefactor;

	basexz = outvect[0]+rotptx-xtran;
	baseyz = outvect[1]+rotpty-ytran;
	basezz = outvect[2]+rotptz-ztran;

/*	printf("outvect[2] = %f\n",outvect[2]);
	printf("rotptz     = %f\n",rotptz);
	printf("ztran      = %f\n",ztran);

	printf("(-1,-1,-1) maps to (%f, %f, %f).\n",
					basexz, baseyz, basezz);
*/
    /* ********************************************************	*/
    /* Now resample image, stepping through with increments	*/
    /* determined above, using (-1,-1,-1) as the starting	*/
    /* point.							*/
/*
    printf("dxx = %f;  dyx = %f;  dzx = %f\n",dxx,dyx,dzx);
    printf("dxy = %f;  dyy = %f;  dzy = %f\n",dxy,dyy,dzy);
    printf("dxz = %f;  dyz = %f;  dzz = %f\n",dxz,dyz,dzz);
*/
    for (z = 0; z <= zdim-1; z++) {

	basexz+=dxz;
	baseyz+=dyz;
	basezz+=dzz;


	basexy=basexz;
	baseyy=baseyz;
	basezy=basezz;

        for (y = 0; y <= ydim-1; y++) {

	    basexy+=dxy;
	    baseyy+=dyy;
	    basezy+=dzy;


	    basexx=basexy;
	    baseyx=baseyy;
	    basezx=basezy;

	    for (x = 0; x <= xdim-1; x++) {

		basexx+=dxx;
		baseyx+=dyx;
		basezx+=dzx;

		u=basexx;
		v=baseyx;
		w=basezx;
/*
		if (((int)u!=x) || ((int)v!=y) || ((int)w!=z))
		    printf("%d %d %d ---> %6.2f %6.2f %6.2f\n",x,y,z,u,v,w);
*/
		flu = floor((double)u); flup1 = flu+1.;
		flv = floor((double)v); flvp1 = flv+1.;
		flw = floor((double)w); flwp1 = flw+1.;

		vol0 =    (flup1-u) *    (flvp1-v) *    (flwp1-w);
		vol1 =    (flup1-u) *    (flvp1-v) *      (w-flw);
		vol2 =    (flup1-u) *      (v-flv) *    (flwp1-w);
		vol3 =    (flup1-u) *      (v-flv) *      (w-flw);
		vol4 =    (u-flu)   *    (flvp1-v) *    (flwp1-w);
		vol5 =    (u-flu)   *    (flvp1-v) *      (w-flw);
		vol6 =    (u-flu)   *      (v-flv) *    (flwp1-w);
		vol7 =    (u-flu)   *      (v-flv) *      (w-flw);

		val0 = ((int)flu>=0 && (int)flu<xdim && (int)flv>=0 && (int)flv<ydim && (int)flw>=0
		     && (int)flw<zdim)   ? (double) mat[((int)flw*ydim+(int)flv)*xdim+(int)flu] : 0.;

		val1 = ((int)flu>=0 && (int)flu<xdim && (int)flv>=0 && (int)flv<ydim && (int)flwp1>=0
		     && (int)flwp1<zdim) ? (double) mat[((int)flwp1*ydim+(int)flv)*xdim+(int)flu] : 0.;

		val2 = ((int)flu>=0 && (int)flu<xdim && (int)flvp1>=0 && (int)flvp1<ydim && (int)flw>=0
		     && (int)flw<zdim)   ? (double) mat[((int)flw*ydim+(int)flvp1)*xdim+(int)flu] : 0.;

		val3 = ((int)flu>=0 && (int)flu<xdim && (int)flvp1>=0 && (int)flvp1<ydim && (int)flwp1>=0
		     && (int)flwp1<zdim) ? (double) mat[((int)flwp1*ydim+(int)flvp1)*xdim+(int)flu] : 0.;

		val4 = ((int)flup1>=0 && (int)flup1<xdim && (int)flv>=0 && (int)flv<ydim && (int)flw>=0
		     && (int)flw<zdim)   ? (double) mat[((int)flw*ydim+(int)flv)*xdim+(int)flup1] : 0.;

		val5 = ((int)flup1>=0 && (int)flup1<xdim && (int)flv>=0 && (int)flv<ydim && (int)flwp1>=0
		     && (int)flwp1<zdim) ? (double) mat[((int)flwp1*ydim+(int)flv)*xdim+(int)flup1] : 0.;

		val6 = ((int)flup1>=0 && (int)flup1<xdim && (int)flvp1>=0 && (int)flvp1<ydim && (int)flw>=0
		     && (int)flw<zdim)   ? (double) mat[((int)flw*ydim+(int)flvp1)*xdim+(int)flup1] : 0.;

		val7 = ((int)flup1>=0 && (int)flup1<xdim && (int)flvp1>=0 && (int)flvp1<ydim && (int)flwp1>=0
		     && (int)flwp1<zdim) ? (double) mat[((int)flwp1*ydim+(int)flvp1)*xdim+(int)flup1] : 0.;

		newmat[(z*ydim+y)*xdim+x] =
					(unsigned short int) (scale * (val0*vol0 + val1*vol1 + val2*vol2 + val3*vol3
							+ val4*vol4 + val5*vol5 + val6*vol6 + val7*vol7));
	    }
        }
    }
}

