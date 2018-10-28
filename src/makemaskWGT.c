#include <stdio.h>
#include <stdlib.h>

main (int argc, char *argv[])
{
    /* 11.14.94 : JM Maisog
	Makes a mask out of a unsigned short image.		*/



    FILE *fIN, *fOUT;
    unsigned short *imgin, threshold;
    unsigned char *imgout;
    int num_voxX2;
    int num_vox;
    int vox_index;
    int read_status;
    int write_status;


    /* **************************************************************** */
    /* VERIFY NUMBER OF COMMAND LINE ARGUMENTS.				*/
    if (argc != 4) {
	printf("makemaskWGT: need input file, output file, and threshold.\n");
	printf("An 8-bit mask will be made from the unsigned short input.\n");
	printf("All unsigned short values greater than the\n");
	printf("threshold will be set to 1 in the mask.\n");
	exit(1);
    }

    /* **************************************************************** */
    /* OPEN FILES FOR INPUT AND OUTPUT.					*/
    if ((fIN=fopen(argv[1],"r")) == NULL) {
	printf("makemaskWGT: can't open %s\n",argv[1]);
	exit(2);
    }
    if ((fOUT=fopen(argv[2],"w")) == NULL) {
	printf("makemaskWGT: can't open %s\n",argv[2]);
	exit(3);
    }

    /* **************************************************************** */
    /* ASSIGN VALUE TO threshold.					*/
    threshold = atof(argv[3]);
    printf("threshold = %f\n",threshold);

    /* **************************************************************** */
    /* COUNT NUMBER OF BYTES IN INPUT FILE				*/
    num_voxX2 = 0;
    while (fgetc(fIN)!=EOF)
	num_voxX2++;

    num_vox = num_voxX2/2;
    printf("makemaskWGT: Number of voxels = %d\n",num_vox);

    /* **************************************************************** */
    /* ASSIGN MEMORY.							*/
    imgin = (unsigned short *) malloc(num_vox*sizeof(unsigned short));
    if (imgin == NULL) {
	printf("makemaskWGT: failed to malloc for imgin.\n");
	exit(4);
    }

    imgout = (unsigned char *) malloc(num_vox*sizeof(unsigned char));
    if (imgout == NULL) {
	printf("makemaskWGT: failed to malloc for imgout.\n");
	exit(5);
    }

    /* **************************************************************** */
    /* READ IN FP DATA.							*/
    rewind(fIN);
    read_status=fread(imgin, sizeof(unsigned short), num_vox, fIN);
    fclose(fIN);
    if (read_status != num_vox) {
        printf("%d voxels read in.  Should have been %d voxels.\n",read_status,
                                                                num_vox);
        exit(8);
    }

    /* **************************************************************** */
    /* MASK FP DATA TO UNSIGNED 8 BIT.					*/

    for (vox_index=0; vox_index<num_vox; vox_index++)
	imgout[vox_index] = (imgin[vox_index]>threshold) ? 1 : 0;

    /* **************************************************************** */
    /* WRITE OUT UNSIGNED 8 BIT TO OUTPUT FILE				*/
    write_status=fwrite(imgout, sizeof(unsigned char), num_vox, fOUT);
    fclose(fOUT);
    if (write_status != num_vox) {
        printf("%d voxels written.  Should have been %d voxels.\n",write_status,
                                                                num_vox);
        exit(10);
    }


    /* CLOSE FILES							*/
    fclose(fIN);
    fclose(fOUT);
}
