#include <stdio.h>
#include <stdlib.h>
#include <math.h>

main (int argc, char *argv[])
{
    /* 03.04.94 : JM Maisog
	Converts a single-precision floating point flat file to
	signed 16-bit flat file, with scalefactor user specifies. */



    FILE *fIN, *fOUT;
    float *imgin;
    float scalefactor;
    signed short int *imgout;
    int num_voxX4;
    int num_vox;
    int vox_index;
    int read_status;
    int write_status;


    /* **************************************************************** */
    /* VERIFY NUMBER OF COMMAND LINE ARGUMENTS.				*/
    if (argc < 4) {
	printf("fp2word2: need input file, output file, and scalefactor.\n");
	exit(1);
    }

    /* **************************************************************** */
    /* DETERMINE SCALEFACTOR.						*/
    scalefactor = atof(argv[3]);
    printf("scalefactor = %f\n",scalefactor);

    /* **************************************************************** */
    /* OPEN FILES FOR INPUT AND OUTPUT.					*/
    if ((fIN=fopen(argv[1],"r")) == NULL) {
	printf("fp2word2: can't open %s\n",argv[1]);
	exit(2);
    }
    if ((fOUT=fopen(argv[2],"w")) == NULL) {
	printf("fp2word2: can't open %s\n",argv[2]);
	exit(3);
    }

    /* **************************************************************** */
    /* COUNT NUMBER OF BYTES IN INPUT FILE				*/
    num_voxX4 = 0;
    while (fgetc(fIN)!=EOF)
	num_voxX4++;

    num_vox = num_voxX4/4;
    printf("fp2word2: Number of voxels = %d\n",num_vox);

    /* **************************************************************** */
    /* ASSIGN MEMORY.							*/
    imgin = (float *) malloc(num_vox*sizeof(float));
    if (imgin == NULL) {
	printf("fp2word2: failed to malloc for imgin.\n");
	exit(4);
    }

    imgout = (signed short int *) malloc(num_vox*sizeof(signed short int));
    if (imgout == NULL) {
	printf("fp2word2: failed to malloc for imgout.\n");
	exit(5);
    }

    /* **************************************************************** */
    /* CLOSE, THEN RE-OPEN INPUT FILE, AND  READ IN FP DATA		*/
    fclose(fIN);
    fIN=fopen(argv[1],"r");
    read_status=fread(imgin, sizeof(float), num_vox, fIN);
    fclose(fIN);
    if (read_status != num_vox) {
        printf("%d voxels read in.  Should have been %d voxels.\n",read_status,
                                                                num_vox);
        exit(8);
    }


    /* **************************************************************** */
    /* CONVERT FP DATA TO SIGNED 16 BIT					*/

    for (vox_index=0; vox_index<num_vox; vox_index++) {
	imgin[vox_index] *= scalefactor;
	if (imgin[vox_index] >  32767.) imgin[vox_index] =  32767;
	if (imgin[vox_index] < -32768.) imgin[vox_index] = -32768;
	imgout[vox_index] = (signed short int) floor ((double) imgin[vox_index]+0.5);
    }

    /* **************************************************************** */
    /* WRITE OUT SIGNED 16 BIT TO OUTPUT FILE				*/
    write_status=fwrite(imgout, sizeof(signed short int), num_vox, fOUT);
    fclose(fOUT);
    if (write_status != num_vox) {
        printf("%d voxels written.  Should have been %d voxels.\n",write_status,
                                                                num_vox);
        exit(10);
    }


    /* CLOSE FILES							*/
/*    fclose(fIN); Joe - 25 October 2002 - commenting out extra fclose()
    fclose(fOUT);  Joe - 25 October 2002 - commenting out extra fclose() */
    exit(0);
}
