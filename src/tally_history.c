#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define sample(xdim, ydim, zdim, x, y, z) ((x)>=0 && (x)<(xdim) && (y)>=0 \
	&& (y)<(ydim) && (z)>=0 && \
	(z)<(zdim)) ? 1:0

void   mkrotmat(double x_theta, double y_theta, double z_theta,
 							float rotmat[3][3]);
/* ************************************************************************** */
int tally_history(int xdim, int ydim, int zdim, float zscalefactor,
              int xlo, int xhi, int ylo, int yhi,
	      float *param, float *rotpt)

/*	Takes as input pointer history to a matrix with dimensions xdim, ydim,*/
/*	and zdim; translation, rotation, and scaling parameters in *param;    */
/*	and the coordinates of the point in the image about which to do	      */
/*	rotations, in *rotpt; and a pointer to allocated space in which	to    */
/*	write results.  Returns the result of applying the transformation     */
/*	specified by *param on *history, in the space pointed to by newhist.  */
/*	Unlike readimg, moveimg does NOT allocate memory for the new	      */
/*	matrix (i.e., newhist); this must be done in the calling program with */
/*	the usual call to malloc.  (If moveimg WERE TO allocate memory as     */
/*	readimg does, and if there were repeated calls to moveimg in the      */
/*	calling program, then each call to moveimg after the first call would */
/*	have to be preceded by a call to the stdlib.h function free, to	      */
/*	deallocate memory.  If the programmer were to forget this, memory     */
/*	might get used up, the program would stop execution, and a nonspec-   */
/*	ific/noninformative error message such as "Segmentation fault" might  */
/*	be generated, which may be potentially difficult to debug.)    	      */

/*	07.02.92 : JMM : altered to contain code for interpolpix2.c, so that  */
/*		         repeated calls to an outside function can be avoided */
/*		         (potentially speeding things up)		      */
/*	07.09.92 : JMM : changed pointer param offsets to reflect 1-offset    */
/*		         data; to be used with domaxzero3a2.c, amoeba.c	      */
/*	07.15.92 : JMM : added input arguments step,xlo,xhi,ylo.yhi	      */
/*	09.02.92 : JMM : added code to take into account pixel dimensions     */
/*			 when voxels are NOT isometric.  Thus, this moveimg   */
/*			 will now work on 15-slice O-15 PET data without      */
/*			 having to interpolate up to 46 slices to get cubic   */
/*			 voxels.					      */
/*	09.10.92 : JMM : added code to skip doing lines in y direction        */
/*			 if step=1, because then all the necessary resampling */
/*			 would already have been done in doing the lines in   */
/*			 the x direction.				      */
/*	02.17.93 : JMM : Altered from moveimg7.c to create tally_history.c.   */
/*			 Routine now generates binary image which has 1's     */
/*			 where image is fully sampled from original image,    */
/*			 and 0's where image is NOT fully sampled from	      */
/*			 original image.				      */
/*	02.18.93 : JMM : Altered tally_history so that it never really creates*/
/*			 binary image, but still returns correct count of     */
/*			 fully sampled voxels.				      */
{
	float  xtran,ytran,ztran;
	double x_theta,y_theta,z_theta;
	float  scale;
	float  rotptx,rotpty,rotptz;
	float  rot_mat[3][3],inv_mat[3][3];
	float  invect[3], outvect[3];
	int    x,y,z;
        int i,j;
	int count;

	float u,v,w;
	float vol0, vol1, vol2, vol3, vol4, vol5, vol6, vol7;
	float val0, val1, val2, val3, val4, val5, val6, val7;
	register float flu, flv, flw, flup1, flvp1, flwp1;
 
	xtran   = *(param+1);
	ytran   = *(param+2);
	ztran   = *(param+3);
	x_theta = (double) *(param+4);
	y_theta = (double) *(param+5);
	z_theta = (double) *(param+6);
	scale   = *(param+7);
	rotptx  = *(rotpt+0);
	rotpty  = *(rotpt+1);
	rotptz  = *(rotpt+2);

	mkrotmat(-x_theta,-y_theta,-z_theta,inv_mat);

	/* Map Z coordinate to isometric voxel space   */
	count = 0;
	for (z = 0; z <= zdim-1; z++) {
	    invect[2] = ((float) (z-rotptz)) * zscalefactor;

	    for (y = ylo; y <= yhi; y++) {
		invect[1] = (float) (y-rotpty);

		for (x = xlo; x <= xhi; x++) {
		    invect[0] = (float) (x-rotptx);

		    for (i = 0; i <= 2; i++) {
			outvect[i] = 0.;
			for (j = 0; j <= 2; j++)
			    outvect[i] += inv_mat[i][j] * invect[j];
			}

		    /* Now map moved z coordinate back to */
		    /* NON-isometric voxel space          */

		    outvect[2] /= zscalefactor;

		    u = outvect[0]+rotptx-xtran;
		    v = outvect[1]+rotpty-ytran;
		    w = outvect[2]+rotptz-ztran;

		    flu = floor(u); flup1 = flu+1.;
		    flv = floor(v); flvp1 = flv+1.;
		    flw = floor(w); flwp1 = flw+1.;

		    val0 = (float) sample(xdim,ydim,zdim,
				    (int)flu,   (int)flv,   (int)flw  );
		    val1 = (float) sample(xdim,ydim,zdim,
				    (int)flu,   (int)flv,   (int)flwp1);
		    val2 = (float) sample(xdim,ydim,zdim,
				    (int)flu,   (int)flvp1,(int)flw  );
		    val3 = (float) sample(xdim,ydim,zdim,
				    (int)flu,   (int)flvp1,(int)flwp1);
		    val4 = (float) sample(xdim,ydim,zdim,
				    (int)flup1,(int)flv,   (int)flw  );
		    val5 = (float) sample(xdim,ydim,zdim,
				    (int)flup1,(int)flv,   (int)flwp1);
		    val6 = (float) sample(xdim,ydim,zdim,
				    (int)flup1,(int)flvp1,(int)flw  );
		    val7 = (float) sample(xdim,ydim,zdim,
				    (int)flup1,(int)flvp1,(int)flwp1);

		    if ((val0==1) && (val1==1) && (val2==1) && (val3==1)
			&& (val4==1) && (val5==1) && (val6==1) && (val7==1))
			count++;
		}
	    }
	}
	return count;
}

