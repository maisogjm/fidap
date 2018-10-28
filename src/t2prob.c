#include <stdio.h>
#include <stdlib.h>
float betai(float a, float b, float x);

main(int argc, char *argv[])
{
    /* 03.14.94 : JM Maisog
	Converts an t-test image to a probability image.			*/
    FILE *fIN;		/* Pointer to input image file.				*/
    FILE *fOUT;		/* Pointer to output image file.			*/
    float *imgin,	/* F-test						*/
	      df;	/* degrees of freedom.					*/
    float *imgout;	/* Pointer to output probability image.			*/
    int num_voxX4;	/* Number of bytes in image.				*/
    int num_vox;	/* Number of voxels in image.				*/
    int vox_index;	/* Index into voxels in image.				*/
    int read_status;	/* Check on read status of fread.			*/
    int write_status;	/* Check on write status of fwrite.			*/


/* ---------------------------------------------------------------------------- */
/*  CHECK TO SEE IF CORRECT NUMBER OF ARGUMENTS WERE PASSED.                    */

    if(argc<4) {
	printf("t2prob: need input t-test image, degrees of freedom, ");
	printf("and output filename.\n");
	exit(1);
    }

/* ---------------------------------------------------------------------------- */
/*  READ IN COMMAND LINE ARGUMENTS.						*/

    df     = atof(argv[2]);

/* ---------------------------------------------------------------------------- */
/* OPEN INPUT FILE.								*/
    if ((fIN=fopen(argv[1],"r")) == NULL) {
	printf("t2prob: can't open %s\n",argv[1]);
	exit(2);
    }

/* ---------------------------------------------------------------------------- */
/* COUNT NUMBER OF BYTES IN INPUT FILE.						*/
    num_voxX4 = 0;
    while (fgetc(fIN)!=EOF)
	num_voxX4++;

    num_vox = num_voxX4/4;
    printf("t2prob: Number of voxels = %d\n",num_vox);

/* ---------------------------------------------------------------------------- */
/* ASSIGN MEMORY.								*/
    imgin = (float *) malloc(num_vox*sizeof(float));
    if (imgin == NULL) {
	printf("t2prob: failed to malloc for imgin.\n");
	exit(4);
    }
    imgout = (float *) malloc(num_vox*sizeof(float));
    if (imgout == NULL) {
    printf("t2prob: failed to malloc for imgout.\n");
    exit(5);
    }

/* ---------------------------------------------------------------------------- */
/* REWIND INPUT IMAGE FILE AND READ IN DATA.					*/
    rewind(fIN);
    read_status=fread(imgin, sizeof(float), num_vox, fIN);
    fclose(fIN);
    if (read_status != num_vox) {
	printf("%d voxels read in.  Should have been %d voxels.\n",read_status,
									num_vox);
	exit(8);
    }

/* ---------------------------------------------------------------------------- */
/* CONVERT INPUT IMAGE DATA TO PROBABILITIES.					*/
    
    for (vox_index=0; vox_index<num_vox; vox_index++)
	imgout[vox_index]=betai(df/2.,0.5,(df/(df+(imgin[vox_index]*imgin[vox_index]))));

/* ---------------------------------------------------------------------------- */
/* WRITE TO OUTPUT FILE.							*/
    if ((fOUT=fopen(argv[3],"w")) == NULL) {
	printf("t2prob: can't open %s\n",argv[2]);
	exit(3);
    }
    write_status=fwrite(imgout, sizeof(float), num_vox, fOUT);
    fclose(fOUT);
    if (write_status != num_vox) {
	printf("%d voxels written.  Should have been %d voxels.\n",write_status,
									num_vox);
	exit(10);
    }
    fclose(fOUT);

}
