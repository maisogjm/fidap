#include <stdio.h>

main (int argc, char *argv[])
{
    /* 10.27.95 : JM Maisog, M.D. : Counts number of bytes in a file. */

    FILE *fIN;
    int num_bytes;

    /* ************************************************************** */
    /* VERIFY NUMBER OF COMMAND LINE ARGUMENTS.                       */
    if (argc!=2) {
        printf("countbytes: need input file.\n");
        printf("            Returns number of bytes.\n");
        exit(1);
    }

    /* ************************************************************** */
    /* OPEN INPUT FILE.                                               */
    if ((fIN=fopen(argv[1],"r")) == NULL) {
        printf("countbytes: can't open %s\n",argv[1]);
        exit(2);
    }

    /* ************************************************************** */
    /* COUNT NUMBER OF BYTES IN INPUT FILE                            */
    num_bytes = 0;
    while (fgetc(fIN)!=EOF)
        num_bytes++;
    fclose(fIN);
    printf("%d\n",num_bytes);
}
