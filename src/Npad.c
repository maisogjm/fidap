#include "math.h"

int main(int argc, char *argv[])
{
    /* 09.06.96 : Jose' Ma. Maisog, M.D.
       Takes as input the number of slices in a volume to be smoothed, and the
       FWHM of the smoothing to be done along the Z-dimension.  Returns a
       suggested number of blank "zero" slices to pad with.                     */

    int zdim;    /* Number of slices in volume to be smoothed.                  */
    float Nwrap; /* Minimum number of blank "zero" slices required to avoid     */
		 /* significant wrap-around effects.                            */
    float b;     /* Nearest integer power of 2 >= a.                            */
    float log2;  /* log base 2 of a, rounded to next lowest integer.            */
    float FWHMz; /* FWHM of the smoothing to be done along the Z-dimension.     */

    /* ------------------------------------------------------------------------ */
    /* CHECK INPUT ARGUMENTS.                                                   */

    if (argc<3) {
        printf("Npad: Need zdim, and FWHM of the smoothing to be done along\n");
	printf("      the Z-dimension.  Express FWHM in terms of VOXEL LENGTHS,\n");
	printf("      not mm.  For example, if you wish to smooth with FWHM\n");
	printf("      of 10 mm, and slice thickness was 5 mm, this is a FWHM\n");
	printf("      of 2 VOXEL LENGTHS.\n\n");
	printf("      Will return suggested number of blank \"zero\" slices to pad with.\n");
        exit(1);
    }


    /* ------------------------------------------------------------------------ */
    /* OUTPUT SUGGESTED NUMBER OF BLANK SLICES.                                 */

    zdim  = atoi(argv[1]);
    FWHMz = atoi(argv[2]);
    Nwrap = zdim+(3*FWHMz);

    log2 = (int)(log((double)Nwrap)/log(2.));
    if ((int)pow(2.,log2)==Nwrap)
        b=Nwrap;
    else
        b=(float)pow(2,log2+1);
    printf("%d\n",(int)b-zdim);
}
