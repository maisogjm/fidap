#include <stdio.h>
#include <stdlib.h>

int main (int argc, char *argv[])
{
    /* 03.04.99 : JM Maisog
	Makes a 16-bit mask out of a unsigned short image		*/

    FILE *fIN, *fOUT;
    unsigned char *imgin, threshold;
    unsigned short *imgout;
    int num_voxX2;
    int num_vox;
    int vox_index;
    int read_status;
    int write_status;

    /* **************************************************************** */
    /* VERIFY NUMBER OF COMMAND LINE ARGUMENTS.				*/
    if (argc != 4) {
	printf("makemaskBGT2: need input file, output file, and threshold.\n");
	printf("A 16-bit mask will be made from the unsigned char input.\n");
	printf("All unsigned char values greater than the\n");
	printf("threshold will be set to 1 in the mask.\n");
	exit(1);
    }

    /* **************************************************************** */
    /* OPEN FILE FOR INPUT.						*/
    if ((fIN=fopen(argv[1],"r")) == NULL) {
	printf("makemaskBGT2: can't open %s\n",argv[1]);
	exit(2);
    }

    /* **************************************************************** */
    /* ASSIGN VALUE TO threshold.					*/
    threshold = (unsigned char)atoi(argv[3]);
    printf("threshold = %d\n",threshold);

    /* **************************************************************** */
    /* COUNT NUMBER OF BYTES IN INPUT FILE				*/
    num_vox = 0;
    while (fgetc(fIN)!=EOF)
	num_vox++;

    printf("makemaskBGT2: Number of voxels = %d\n",num_vox);

    /* **************************************************************** */
    /* ASSIGN MEMORY.							*/
    imgin = (unsigned char *) malloc(num_vox*sizeof(unsigned char));
    if (imgin == NULL) {
	printf("makemaskBGT2: failed to malloc for imgin.\n");
	exit(4);
    }

    imgout = (unsigned short *) malloc(num_vox*sizeof(unsigned short));
    if (imgout == NULL) {
	printf("makemaskBGT2: failed to malloc for imgout.\n");
	exit(5);
    }

    /* **************************************************************** */
    /* READ IN 8-BIT DATA.						*/
    rewind(fIN);
    read_status=fread(imgin, sizeof(unsigned char), num_vox, fIN);
    fclose(fIN);
    if (read_status != num_vox) {
        printf("%d voxels read in.  Should have been %d voxels.\n",
				read_status,num_vox);
        exit(8);
    }

    /* **************************************************************** */
    /* MASK 8-bit DATA TO UNSIGNED 16 BIT.					*/
    for (vox_index=0; vox_index<num_vox; vox_index++)
	imgout[vox_index] = (imgin[vox_index]>threshold) ? 1 : 0;

    /* **************************************************************** */
    /* WRITE OUT UNSIGNED 16 BIT TO OUTPUT FILE				*/
    if ((fOUT=fopen(argv[2],"w")) == NULL) {
	printf("makemaskBGT2: can't open %s\n",argv[2]);
	exit(3);
    }
    write_status=fwrite(imgout, sizeof(unsigned short), num_vox, fOUT);
    fclose(fOUT);
    if (write_status != num_vox) {
        printf("%d voxels written.  Should have been %d voxels.\n",
				write_status, num_vox);
        exit(10);
    }
    return 0;
}
