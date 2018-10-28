#include <stdio.h>
#include <stdlib.h>

main (int argc, char *argv[])
{
    /* 12.03.93 : JM Maisog : Subtracts floating point images.	*/

    FILE *fIN, *fOUT;
    float *imgin, *imgout;
    float subtrahend;
    int num_voxX4;
    int num_vox;
    int vox_index;
    int read_status;
    int write_status;


    /* **************************************************************** */
    /* VERIFY NUMBER OF COMMAND LINE ARGUMENTS.				*/
    if (argc < 4) {
	printf("subtractFP2: need input file, number which to subtract,\n");
	printf("             and output file.\n");
	exit(1);
    }

    /* **************************************************************** */
    /* COUNT NUMBER OF BYTES IN FIRST INPUT FILE			*/
    printf("subtractFP2: opening file %s...\n",argv[1]);
    if ((fIN=fopen(argv[1],"r")) == NULL) {
	printf("subtractFP2: can't open %s\n",argv[1]);
	exit(2);
    }
    num_voxX4 = 0;
    while (fgetc(fIN)!=EOF)
	num_voxX4++;

    fclose(fIN);
    num_vox = num_voxX4/4;
    printf("subtractFP2: Number of voxels = %d\n",num_vox);


    /* **************************************************************** */
    /* ASSIGN MEMORY.							*/
    imgin = (float *) malloc(num_vox*sizeof(float));
    if (imgin == NULL) {
	printf("subtractFP2: failed to malloc for imgin.\n");
	exit(4);
    }
    imgout = (float *) malloc(num_vox*sizeof(float));
    if (imgout == NULL) {
	printf("subtractFP2: failed to malloc for imgout.\n");
	exit(5);
    }

    /* **************************************************************** */
    /* READ IN VOXELS FROM FIRST INPUT FILE.				*/
    if ((fIN=fopen(argv[1],"r")) == NULL) {
	printf("subtractFP2: unable to open %s...\n",argv[1]);
	exit(6);
    }
    printf("subtractFP2: reading in voxels from %s...\n",argv[1]);
    read_status=fread(imgin, sizeof(float), num_vox, fIN);
    fclose(fIN);
    if (read_status != num_vox) {
	printf("%d voxels read in.  Should have been %d voxels.\n",
						read_status, num_vox);
	exit(8);
    }

    /* **************************************************************** */
    /* PUT DIFFERENCE IN *imgout.					*/
    subtrahend=atof(argv[2]);
    printf("Subtracting %f from %s...\n",subtrahend, argv[1]);
    for (vox_index=0; vox_index<num_vox; vox_index++)
	imgout[vox_index] = imgin[vox_index] - subtrahend;


    /* **************************************************************** */
    /* WRITE TO OUTPUT FILE, AND CLOSE FILES.				*/
    printf("subtractFP2: opening file %s...\n",argv[3]);
    if ((fOUT=fopen(argv[3],"w")) == NULL) {
	printf("subtractFP2: can't open %s\n",argv[3]);
	exit(3);
    }
    printf("subtractFP2: writing to %s...\n",argv[3]);
    write_status=fwrite(imgout, sizeof(float), num_vox, fOUT);
    fclose(fOUT);
    if (write_status != num_vox) {
        printf("%d voxels written.  Should have been %d voxels.\n",
						write_status,num_vox);
        exit(10);
    }
}
