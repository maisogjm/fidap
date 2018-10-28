#include <stdio.h>
#include <stdlib.h>

main(int argc, char *argv[])
{
    /* 06.15.94 : JM Maisog : Divides floating point image by scalar. */

    FILE *fIN, *fOUT;
    float *imgin, *imgout;
    float divisor;
    int num_voxX4;
    int num_vox;
    int vox_index;
    int read_status;
    int write_status;

    /* ****************************************************************	*/
    /* CHECK COMMAND LINE ARGUMENTS.					*/
    if (argc!=4) {
	printf("divideFP2: need input file, divisor, and output file.\n");
	exit(1);
    }

    /* **************************************************************** */
    /* COUNT NUMBER OF BYTES IN INPUT FILE.				*/
    printf("divideFP2: opening file %s...\n",argv[1]);
    if ((fIN=fopen(argv[1],"r"))==NULL) {
	printf("divideFP2: can't open %s\n",argv[1]);
	exit(2);
    }
    num_voxX4=0;
    while (fgetc(fIN)!=EOF)
	num_voxX4++;

    num_vox=num_voxX4/4;
    printf("divideFP2: number of voxels = %d\n",num_vox);

    /* **************************************************************** */
    /* ASSIGN MEMORY.							*/
    imgin=(float *)malloc(num_vox*sizeof(float));
    if (imgin==NULL) {
	printf("divideFP2: failed to malloc for imgin.\n");
	exit(4);
    }
    imgout=(float *)malloc(num_vox*sizeof(float));
    if (imgout==NULL) {
	printf("divideFP2: failed to malloc for imgout.\n");
	exit(5);
    }

    divisor=atof(argv[2]);
    printf("divideFP2: divisor = %f\n",divisor);
    if (divisor==0) {
	printf("divideFP2: can't divide by 0.\n");
	exit(6);
    }

    /* **************************************************************** */
    /* READ IN VOXELS FROM INPUT FILE.					*/
    rewind(fIN);
    printf("divideFP2: reading in voxels from %s\n",argv[1]);
    read_status=fread(imgin,sizeof(float),num_vox,fIN);
    fclose(fIN);
    if (read_status!=num_vox) {
	printf("divideFP2: %d voxels read in.  Should have been %d voxels.\n",
					read_status, num_vox);
	exit(8);
    }

    /* **************************************************************** */
    /* PUT DIVIDEND IN *imgout.						*/
    printf("divideFP2: dividing %s by %f...\n",argv[1],divisor);
    for (vox_index=0; vox_index<num_vox; vox_index++)
	imgout[vox_index]=imgin[vox_index]/divisor;


    /* **************************************************************** */
    /* WRITE TO OUTPUT FILE, AND CLOSE FILES.				*/
    printf("divideFP2: opening file %s...\n",argv[3]);
    if ((fOUT=fopen(argv[3],"w"))==NULL) {
	printf("divideFP2: can't open %s\n",argv[3]);
	exit(3);
    }
    printf("divideFP2: writing to %s...\n",argv[3]);
    write_status=fwrite(imgout,sizeof(float),num_vox,fOUT);
    fclose(fOUT);
    if (write_status!=num_vox) {
	printf("divideFP2: %d voxels written.  Should have been %d voxels.\n",
					write_status,num_vox);
	exit(10);
    }
    free(imgin);
    free(imgout);
}
