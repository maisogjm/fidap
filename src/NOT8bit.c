#include <stdio.h>
#include <stdlib.h>
#include <math.h>

main (int argc, char *argv[])
{
    /* 10.02.95 : JM Maisog
        Performs NOT operation on voxels in an 8-bit mask.    */


    FILE *fIN;
    unsigned char *imgin;
    int num_vox;
    int vox_index;
    int read_status,write_status;


    /* **************************************************************** */
    /* VERIFY NUMBER OF COMMAND LINE ARGUMENTS.                         */
    if (argc!=3) {
        printf("NOT8bit: need input 8-bit mask file, name of output file.\n");
        printf("         Performs NOT operation on voxels.\n");
        exit(1);
    }

    /* **************************************************************** */
    /* OPEN FILES FOR INPUT AND OUTPUT.                                 */
    if ((fIN=fopen(argv[1],"r")) == NULL) {
        printf("NOT8bit: can't open %s\n",argv[1]);
        exit(2);
    }

    /* **************************************************************** */
    /* COUNT NUMBER OF BYTES IN INPUT FILE                              */
    num_vox = 0;
    while (fgetc(fIN)!=EOF)
        num_vox++;


    /* **************************************************************** */
    /* ASSIGN MEMORY.                                                   */
    imgin = (unsigned char *) malloc(num_vox*sizeof(unsigned char));
    if (imgin == NULL) {
        printf("NOT8bit: failed to malloc for imgin.\n");
        exit(4);
    }


    /* **************************************************************** */
    /* READ IN FP DATA                                                  */
    rewind(fIN);
    read_status=fread(imgin, sizeof(unsigned char), num_vox, fIN);
    fclose(fIN);
    if (read_status != num_vox) {
        printf("%d voxels read in.  Should have been %d voxels.\n",read_status,
                                                                num_vox);
        exit(8);
    }

    /* **************************************************************** */
    /* PERFORM NOT OPERATION.                                           */

    for (vox_index=0; vox_index<num_vox; vox_index++)
        if (imgin[vox_index]!=0) 
            imgin[vox_index]=0;
        else 
            imgin[vox_index]=1;

    if ((fIN=fopen(argv[2],"w")) == NULL) {
        printf("NOT8bit: can't open %s\n",argv[2]);
        exit(9);
    }
    write_status=fwrite(imgin,sizeof(unsigned char),num_vox,fIN);
    if (write_status != num_vox) {
        printf("%d voxels written.  Should have been %d voxels.\n",write_status,
                                                                num_vox);
        exit(10);
    }

    fclose(fIN);
}
