#include <stdio.h>
#include <stdlib.h>

main (int argc, char *argv[])
{
    /* 05.24.95 : JM Maisog
        Counts number of non-zero voxels in an 8-bit file. */


    FILE *fIN;
    unsigned char *imgin;
    int num_vox;
    int vox_index;
    int read_status;
    int num_nonzero;


    /* **************************************************************** */
    /* VERIFY NUMBER OF COMMAND LINE ARGUMENTS.                         */
    if (argc!=2) {
        printf("countmask: need input 8-bit file.\n");
        printf("           Returns number of non-zero voxels.\n");
        exit(1);
    }

    /* **************************************************************** */
    /* OPEN FILES FOR INPUT AND OUTPUT.                                 */
    if ((fIN=fopen(argv[1],"r")) == NULL) {
        printf("countmask: can't open %s\n",argv[1]);
        exit(2);
    }

    /* **************************************************************** */
    /* COUNT NUMBER OF BYTES IN INPUT FILE                              */
    num_vox = 0;
    while (fgetc(fIN)!=EOF)
        num_vox++;

/*    printf("countmask: Number of voxels = %d\n",num_vox);
*/

    /* **************************************************************** */
    /* ASSIGN MEMORY.                                                   */
    imgin = (unsigned char *) malloc(num_vox*sizeof(unsigned char));
    if (imgin == NULL) {
        printf("countmask: failed to malloc for imgin.\n");
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
    /* COUNT NON-ZERO VOXELS.                                           */

    num_nonzero=0;
    for (vox_index=0; vox_index<num_vox; vox_index++)
        if (imgin[vox_index]!=0) num_nonzero++;

     printf("%d\n",num_nonzero);
}
