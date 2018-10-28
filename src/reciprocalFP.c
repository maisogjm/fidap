#include <stdio.h>
#include <stdlib.h>
#include <math.h>

main (int argc, char *argv[])
{
    /* 01.19.96 : JM Maisog : Returns voxel-wise reciprocals.           */

    FILE *fIN, *fOUT;
    float *imgin, *imgout;
    int cmd_lin_arg;
    int num_voxX4;
    int num_vox;
    int vox_index;
    int read_status;
    int write_status;


    /* **************************************************************** */
    /* VERIFY NUMBER OF COMMAND LINE ARGUMENTS.                         */
    if (argc < 3) {
        printf("reciprocalFP: need input file and output file.\n");
        exit(1);
    }

    /* **************************************************************** */
    /* COUNT NUMBER OF BYTES IN FIRST INPUT FILE                        */
    /* printf("reciprocalFP: opening file %s...\n",argv[1]); */
    if ((fIN=fopen(argv[1],"r")) == NULL) {
        printf("reciprocalFP: can't open %s\n",argv[1]);
        exit(2);
    }
    num_voxX4 = 0;
    while (fgetc(fIN)!=EOF)
        num_voxX4++;

    fclose(fIN);
    num_vox = num_voxX4/4;
    /* printf("reciprocalFP: Number of voxels = %d\n",num_vox); */


    /* **************************************************************** */
    /* ASSIGN MEMORY.                                                   */
    imgin = (float *) malloc(num_vox*sizeof(float));
    if (imgin == NULL) {
        printf("reciprocalFP: failed to malloc for imgin.\n");
        exit(4);
    }
    imgout = (float *) malloc(num_vox*sizeof(float));
    if (imgout == NULL) {
        printf("reciprocalFP: failed to malloc for imgout.\n");
        exit(5);
    }


    /* **************************************************************** */
    /* READ IN VOXELS FROM INPUT FILE.                                  */
    if ((fIN=fopen(argv[1],"r")) == NULL) {
        printf("reciprocalFP: unable to open %s...\n",argv[1]);
        exit(6);
    }
    /* printf("reciprocalFP: reading in voxels from %s...\n",argv[1]); */
    read_status=fread(imgin, sizeof(float), num_vox, fIN);
    fclose(fIN);
    if (read_status != num_vox) {
        printf("%d voxels read in.  Should have been %d voxels.\n",
                                                read_status, num_vox);
        exit(8);
    }

    /* **************************************************************** */
    /* PUT RECIPROCAL IN *imgout.                                       */
    printf("If any voxel in %s is zero, the output voxel will be set to zero.\n",
                                                                    argv[2]);
    for (vox_index=0; vox_index<num_vox; vox_index++) {
        if (fabs(imgin[vox_index])<pow(10.,-32.))
            imgout[vox_index] = 0.;
        else
            imgout[vox_index] = 1./imgin[vox_index];
    }


    /* **************************************************************** */
    /* WRITE TO OUTPUT FILE, AND CLOSE FILES.                           */
    /* printf("reciprocalFP: opening file %s...\n",argv[2]); */
    if ((fOUT=fopen(argv[2],"w")) == NULL) {
        printf("reciprocalFP: can't open %s\n",argv[2]);
        exit(3);
    }
    printf("reciprocalFP: writing to %s...\n",argv[2]);
    write_status=fwrite(imgout, sizeof(float), num_vox, fOUT);
    fclose(fOUT);
    if (write_status != num_vox) {
        printf("%d voxels written.  Should have been %d voxels.\n",
                                                write_status,num_vox);
        exit(10);
    }
    fclose(fOUT);
}
