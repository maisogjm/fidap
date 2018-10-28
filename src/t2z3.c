#include <stdio.h>
#include <stdlib.h>
#include <math.h>
float betai(float a, float b, float x);
float *float_malloc(int num_vox);
int read_float_data(char *filename, float *array, int num_items);
int float_fwrite(char *filename, float *array, int num_items);

main(int argc, char *argv[])
{
    /* 10.07.96 : JM Maisog
        Converts an t-test image to a zscore image.
	Z = sqrt(df*log(1+t^2/df) * (1-(1/(2*df))))
        Reference: Federighi, E.T. (1959) "Extended tables of the percentage points of
        Student's t-distribution" J. Am. Stat. Assoc., 54, 683-688.             */

    FILE *fIN;           /* Pointer to input image file.                        */
    float *stats_map,    /* t-test                                              */
	  *Z_score_map,  /* Z-score map.                                        */
              df;        /* degrees of freedom.                                 */
    int num_voxX4;       /* Number of bytes in image.                           */
    int num_vox;         /* Number of voxels in image.                          */
    int vox_index;       /* Index into voxels in image.                         */




/* ---------------------------------------------------------------------------- */
/*  CHECK TO SEE IF CORRECT NUMBER OF ARGUMENTS WERE PASSED.                    */

    if(argc<4) {
        printf("t2z: need input t-test image, degrees of freedom, ");
        printf("and output filename.\n");
        exit(1);
    }

/* ---------------------------------------------------------------------------- */
/*  READ IN DEGREES OF FREEDOM FROM COMMAND LINE.                               */

    df     = atof(argv[2]);

/* ---------------------------------------------------------------------------- */
/* OPEN INPUT FILE.                                                             */
    if ((fIN=fopen(argv[1],"r")) == NULL) {
        printf("t2z: can't open %s\n",argv[1]);
        exit(2);
    }

/* ---------------------------------------------------------------------------- */
/* COUNT NUMBER OF BYTES IN INPUT FILE.                                         */
    num_voxX4 = 0;
    while (fgetc(fIN)!=EOF)
        num_voxX4++;
    fclose(fIN);

    num_vox = num_voxX4/4;

/* ---------------------------------------------------------------------------- */
/* ASSIGN MEMORY.                                                               */
    stats_map   = float_malloc(num_vox);
    Z_score_map = float_malloc(num_vox);

/* ---------------------------------------------------------------------------- */
/* READ IN INPUT IMAGE FILE.                                                    */
    read_float_data(argv[1], stats_map, num_vox);

/* ---------------------------------------------------------------------------- */
/* CONVERT INPUT IMAGE DATA TO PROBABILITIES, AND THEN TO ZSCORES.              */
    
    for (vox_index=0; vox_index<num_vox; vox_index++) {
        Z_score_map[vox_index]=
                sqrt(df*log(1.+stats_map[vox_index]*stats_map[vox_index]/df)
                        * (1.-(1./(2.*df))));
	if (stats_map[vox_index]<0.) Z_score_map[vox_index]*=-1.;
    }

/* ---------------------------------------------------------------------------- */
/* WRITE TO OUTPUT FILE.                                                        */
    float_fwrite(argv[3], Z_score_map, num_vox);
    free(stats_map);
    free(Z_score_map);
}
