#include <stdio.h>
#include <stdlib.h>
#include <math.h>

main (int argc, char *argv[])
{
    /* 12.03.93 : JM Maisog : Multiplies floating point images.	*/

    FILE *fIN, *fOUT;
    float *imgin, *product;
    float multiplicand;
    int num_voxX4;
    int num_vox;
    int vox_index;
    int read_status;
    int write_status;


    /* **************************************************************** */
    /* VERIFY NUMBER OF COMMAND LINE ARGUMENTS.				*/
    if (argc != 4) {
	printf("multFP2: need input file, number by which to multiply and output file.\n");
	exit(1);
    }

    /* **************************************************************** */
    /* OPEN FILES FOR INPUT AND OUTPUT.					*/
    printf("multFP2: opening file %s...\n",argv[1]);
    if ((fIN=fopen(argv[1],"r")) == NULL) {
	printf("multFP2: can't open %s\n",argv[1]);
	exit(2);
    }

    /* **************************************************************** */
    /* COUNT NUMBER OF BYTES IN FIRST INPUT FILE			*/
    /* THEN READ THEM IN.						*/

    num_voxX4 = 0;
    while (fgetc(fIN)!=EOF)
	num_voxX4++;

    num_vox = num_voxX4/4;
    printf("multFP2: Number of voxels = %d\n",num_vox);


    /* **************************************************************** */
    /* ASSIGN MEMORY, AND READ IN DATA.					*/
    imgin = (float *) malloc(num_vox*sizeof(float));
    if (imgin == NULL) {
	printf("multFP2: failed to malloc for imgin.\n");
	exit(4);
    }
    rewind(fIN);
    read_status=fread(imgin, sizeof(float), num_vox, fIN);
    fclose(fIN);

    product = (float *) malloc(num_vox*sizeof(float));
    if (product == NULL) {
	printf("multFP2: failed to malloc for product.\n");
	exit(5);
    }

    multiplicand = atof(argv[2]);

    /* **************************************************************** */
    /* MAKE PRODUCT.							*/
    printf("multFP2: multiplying data from %s by %f...\n",argv[1],multiplicand);
    for (vox_index=0; vox_index<num_vox; vox_index++)
	product[vox_index] = imgin[vox_index] * multiplicand;


    /* **************************************************************** */
    /* WRITE OUT UNSIGNED 8 BIT TO OUTPUT FILE, AND CLOSE FILES.	*/
    printf("multFP2: opening file %s...\n",argv[argc-1]);
    if ((fOUT=fopen(argv[argc-1],"w")) == NULL) {
	printf("multFP2: can't open %s\n",argv[argc-1]);
	exit(3);
    }
    printf("multFP2: writing to %s...\n",argv[argc-1]);
    write_status=fwrite(product, sizeof(float), num_vox, fOUT);
    fclose(fOUT);
    if (write_status != num_vox) {
        printf("%d voxels written.  Should have been %d voxels.\n",
						write_status,num_vox);
        exit(10);
    }
}
