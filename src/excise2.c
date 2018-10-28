#include <stdio.h>
#include <stdlib.h>

unsigned char *unsigned_char_malloc(int num_items);

main (int argc, char *argv[])
{
    /* 04.09.96 : JM Maisog, M.D.
        Like excise, skips a specified number of bytes from an input
        file, and then passes a specified number of bytes to an output
        file, APPENDING to that file if it already exist.  Uses fseek.  */

    FILE *fIN, *fOUT;
    unsigned char *c;
    int count,fset;

    if (argc < 5) {
        printf("excise2: need input file, output file, number of\n");
        printf("         bytes to skip, number of bytes to keep.\n");
        exit(1);
    }

    /* OPEN FILES FOR INPUT AND OUTPUT.                                 */
    if ((fIN=fopen(argv[1],"r")) == NULL) {
        printf("excise2: can't open %s\n",argv[1]);
        exit(2);
    }
    if ((fOUT=fopen(argv[2],"a")) == NULL) {
        printf("excise2: can't open %s\n",argv[2]);
        exit(3);
    }
    c=unsigned_char_malloc(atoi(argv[4]));

    /* SKIP NUMBER OF BYTES SPECIFIED.                                  */
    fset = fseek(fIN,atoi(argv[3]),0);
    fread(c,sizeof(unsigned char),atoi(argv[4]),fIN);
    fclose(fIN);
    fwrite(c,sizeof(unsigned char),atoi(argv[4]),fOUT);
    fclose(fOUT);
}
