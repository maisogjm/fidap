#include <stdio.h>
#include <stdlib.h>

main (int argc, char *argv[])
{
    /* 12.23.93 : JM Maisog
        Converts a unsigned 8-bit flat file to single-precision
        floating point flat file. */



    FILE *fIN, *fOUT;
    float *imgout;
    unsigned char *imgin;
    int num_vox;
    int vox_index;
    int read_status;
    int write_status;


    /* **************************************************************** */
    /* VERIFY NUMBER OF COMMAND LINE ARGUMENTS.                         */
    if (argc < 3) {
        printf("byte2fp: need input file and output file.\n");
        exit(1);
    }

    /* **************************************************************** */
    /* OPEN INPUT FILE.                                                 */
    if ((fIN=fopen(argv[1],"r")) == NULL) {
        printf("byte2fp: can't open %s\n",argv[1]);
        exit(2);
    }

    /* **************************************************************** */
    /* COUNT NUMBER OF BYTES IN INPUT FILE                              */
    num_vox = 0;
    while (fgetc(fIN)!=EOF)
        num_vox++;

    printf("byte2fp: Number of voxels = %d\n",num_vox);

    /* **************************************************************** */
    /* ASSIGN MEMORY.                                                   */
    imgin = (unsigned char *) malloc(num_vox*sizeof(unsigned char));
    if (imgin == NULL) {
        printf("byte2fp: failed to malloc for imgin.\n");
        exit(4);
    }

    imgout = (float *) malloc(num_vox*sizeof(float));
    if (imgout == NULL) {
        printf("byte2fp: failed to malloc for imgout.\n");
        exit(5);
    }

    /* **************************************************************** */
    /* CLOSE, THEN RE-OPEN INPUT FILE, AND  READ IN FP DATA             */
    fclose(fIN);
    fIN=fopen(argv[1],"r");
    read_status=fread(imgin, sizeof(unsigned char), num_vox, fIN);
    fclose(fIN);
    if (read_status != num_vox) {
        printf("%d voxels read in.  Should have been %d voxels.\n",read_status,
                                                                num_vox);
        exit(8);
    }

    /* **************************************************************** */
    /* CONVERT UNSIGNED 8-BIT DATA TO FP.                               */

    for (vox_index=0; vox_index<num_vox; vox_index++)
        imgout[vox_index] = (float) imgin[vox_index];

    /* **************************************************************** */
    /* WRITE OUT UNSIGNED 8 BIT TO OUTPUT FILE                          */

    if ((fOUT=fopen(argv[2],"w")) == NULL) {
        printf("byte2fp: can't open %s\n",argv[2]);
        exit(3);
    }
    write_status=fwrite(imgout, sizeof(float), num_vox, fOUT);
    fclose(fOUT);
    if (write_status != num_vox) {
        printf("%d voxels written.  Should have been %d voxels.\n",write_status,
                                                                num_vox);
        exit(10);
    }
}
