#include <math.h>

double dotprod16bit(unsigned short int *mat1, unsigned short int *mat2, int xdim, int ydim, int zdim, int step,
	       int xlo, int xhi, int ylo, int yhi)
{
        double sum_of_prods, mag1, mag2;
        int index;

/*	Measure of how alike two images are.  Method: generate a third image
	which is the pixel-by-pixel difference between the two input images.
	Count the total number of zero crossovers in the third image along the
	X- and Y-axes.  This number is returned as the zerosum.  Zero
	crossovers can be either a change from positive to negative number,
	or a change from negative to positive number.

	07.20.93 : JMM : Modified to handle 16-bit images.		*/

/*printf("Entered dotprod16bit...\n");
printf("xdim = %d\n",xdim);
printf("ydim = %d\n",ydim);
printf("zdim = %d\n",zdim);
printf("step = %d\n",step);
printf("xlo  = %d\n",xlo);
printf("xhi  = %d\n",xhi);
printf("ylo  = %d\n",ylo);
printf("yhi  = %d\n",yhi);
*/

        sum_of_prods = mag1 = mag2 = 0.;

        for (index = 0; index <= xdim*ydim*zdim-1; index+=step) {
            /* Calculate magnitude of mat1 represented as a vector */
            mag1 += ((double) *(mat1+index)) * ((double) *(mat1+index));

            /* Calculate magnitude of mat2 represented as a vector */
            mag2 += ((double) *(mat2+index)) * ((double) *(mat2+index));

            /* Calculate dot product between mat1 and mat2 */
            sum_of_prods += ((double) *(mat1+index)) * ((double) *(mat2+index));
        }

        mag1 = sqrt(mag1);
        mag2 = sqrt(mag2);

        /* Normalize the result to value between 0 and 1 */
        /* Equivalent to normalizing mat1 and mat2       */
        /* BEFORE calculating the dot product, but       */
        /* doing the division by this large number AFTER */
        /* calculating sum_of_prods reduces roundoff     */
        /* error.                                        */
        sum_of_prods /= (mag1*mag2);

        return sum_of_prods;

}
