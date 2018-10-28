#include <stdio.h>
#include <stdlib.h>
#include <math.h>

main (int argc, char *argv[])
{
    /* 09.08.93 : JM Maisog
	Converts a single-precision floating point flat file to
	unsigned 8-bit flat file.  RESCALES VALUES TO RANGE OF
	0 TO 255.						*/



    FILE *fIN, *fOUT;
    float *imgin, min, max;
    unsigned char *imgout;
    int num_voxX4;
    int num_vox;
    int vox_index;
    int read_status;
    int write_status;


    /* **************************************************************** */
    /* VERIFY NUMBER OF COMMAND LINE ARGUMENTS.				*/
    if (argc < 3) {
	printf("fp2byte: need input file and output file.\n");
	exit(1);
    }

    /* **************************************************************** */
    /* OPEN FILES FOR INPUT AND OUTPUT.					*/
    if ((fIN=fopen(argv[1],"r")) == NULL) {
	printf("fp2byte: can't open %s\n",argv[1]);
	exit(2);
    }
    if ((fOUT=fopen(argv[2],"w")) == NULL) {
	printf("fp2byte: can't open %s\n",argv[2]);
	exit(3);
    }

    /* **************************************************************** */
    /* COUNT NUMBER OF BYTES IN INPUT FILE				*/
    num_voxX4 = 0;
    while (fgetc(fIN)!=EOF)
	num_voxX4++;

    num_vox = num_voxX4/4;
    printf("fp2byte: Number of voxels = %d\n",num_vox);

    /* **************************************************************** */
    /* ASSIGN MEMORY.							*/
    imgin = (float *) malloc(num_vox*sizeof(float));
    if (imgin == NULL) {
	printf("fp2byte: failed to malloc for imgin.\n");
	exit(4);
    }

    imgout = (unsigned char *) malloc(num_vox*sizeof(unsigned char));
    if (imgout == NULL) {
	printf("fp2byte: failed to malloc for imgout.\n");
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
    /* DETERMINE MIN AND MAX OF FP DATA					*/
    min = 10000.;
    max = -10000.;
    for (vox_index=0; vox_index<num_vox; vox_index++) {
	min =  (min<imgin[vox_index]) ? min : imgin[vox_index];
	max = (max>imgin[vox_index]) ? max : imgin[vox_index];
    }
    printf("fp2byte: Minimum of %g will be mapped to 0.\n",min);
    printf("fp2byte: Maximum of %g will be mapped to 255.\n",max);

    /* **************************************************************** */
    /* CONVERT FP DATA TO UNSIGNED 8 BIT				*/

    for (vox_index=0; vox_index<num_vox; vox_index++)
	imgout[vox_index] = (unsigned char)
		floor((double) (((imgin[vox_index] - min)/(max-min)*255.0)+0.5));

    /* **************************************************************** */
    /* WRITE OUT UNSIGNED 8 BIT TO OUTPUT FILE				*/
    write_status=fwrite(imgout, sizeof(unsigned char), num_vox, fOUT);
    fclose(fOUT);
    if (write_status != num_vox) {
        printf("%d voxels written.  Should have been %d voxels.\n",write_status,
                                                                num_vox);
        exit(10);
    }
}
