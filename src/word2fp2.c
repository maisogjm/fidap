#include <stdio.h>
#include <stdlib.h>

main (int argc, char *argv[])
{
    /* 09.08.93 : JM Maisog
	Converts a signed 16-bit flat file to
	floating point flat file.				*/



    FILE *fIN, *fOUT;
    signed short int *imgin, min, max;
    float scaling_factor;
    float *imgout;
    int num_voxX2;
    int num_vox;
    int vox_index;
    int read_status;
    int write_status;


    /* **************************************************************** */
    /* VERIFY NUMBER OF COMMAND LINE ARGUMENTS.				*/
    if (argc != 4) {
	printf("word2fp2: need input file, output file, and scaling factor.\n");
	exit(1);
    }

    /* **************************************************************** */
    /* OPEN FILES FOR INPUT AND OUTPUT.					*/
    if ((fIN=fopen(argv[1],"r")) == NULL) {
	printf("word2fp2: can't open %s\n",argv[1]);
	exit(2);
    }
    if ((fOUT=fopen(argv[2],"w")) == NULL) {
	printf("word2fp2: can't open %s\n",argv[2]);
	exit(3);
    }

    /* **************************************************************** */
    /* READ SCALING FACTOR OFF OF COMMAND LINE.				*/
    scaling_factor = atof(argv[3]);
    printf("Output will be scaled by factor of %f\n",scaling_factor);

    /* **************************************************************** */
    /* COUNT NUMBER OF BYTES IN INPUT FILE				*/
    num_voxX2 = 0;
    while (fgetc(fIN)!=EOF)
	num_voxX2++;

    num_vox = num_voxX2/2;
    printf("word2fp2: Number of voxels = %d\n",num_vox);

    /* **************************************************************** */
    /* ASSIGN MEMORY.							*/
    imgin = (signed short int *) malloc(num_vox*sizeof(signed short int));
    if (imgin == NULL) {
	printf("word2fp2: failed to malloc for imgin.\n");
	exit(4);
    }

    imgout = (float *) malloc(num_vox*sizeof(float));
    if (imgout == NULL) {
	printf("word2fp2: failed to malloc for imgout.\n");
	exit(5);
    }

    /* **************************************************************** */
    /* REWIND INPUT DATA FILE AND READ IN DATA.				*/

    rewind(fIN);
    read_status=fread(imgin, sizeof(signed short int), num_vox, fIN);
    fclose(fIN);
    if (read_status != num_vox) {
        printf("%d voxels read in.  Should have been %d voxels.\n",
						read_status, num_vox);
        exit(8);
    }


    /* **************************************************************** */
    /* CONVERT WORD DATA TO FLOATING POINT.				*/

    for (vox_index=0; vox_index<num_vox; vox_index++) {
	imgout[vox_index] = scaling_factor *(float) imgin[vox_index];
    }

    /* **************************************************************** */
    /* WRITE OUT FLOATING POINT IMAGE TO OUTPUT FILE.			*/

    write_status=fwrite(imgout, sizeof(float), num_vox, fOUT);
    fclose(fOUT);
    if (write_status != num_vox) {
        printf("%d voxels written.  Should have been %d voxels.\n",write_status,
                                                                num_vox);
        exit(10);
    }


    /* **************************************************************** */
    /* CLOSE FILES							*/
/*    fclose(fIN);
    fclose(fOUT); */
    exit(0);
}
