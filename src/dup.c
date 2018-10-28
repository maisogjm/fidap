#include <stdio.h>
#include <stdlib.h>
int *int_malloc(int num_items);

/* Program to duplicate a file num_dups times.                          */

main(int argc, char *argv[]) {
    FILE *fIN, *fOUT; /* Input and output file pointers.                */
    int num_dups;     /* Number of duplications.                        */
    int num_char;     /* Number of characters in input file.            */
    int dup_index;    /* Index over number of duplications.             */
    int c;            /* Buffer for reading in a character.             */
    int *char_array;  /* Character array to be duplicated.              */
    int char_index;   /* Index for looping over characters.             */

    /* **************************************************************** */
    /* VERIFY NUMBER OF COMMAND LINE ARGUMENTS.                         */
    if (argc!=4) {
        printf("dup: need input textfile, output textfile, number of duplications.\n");
        exit(1);
    }
    num_dups=atoi(argv[3]);

    /* **************************************************************** */
    /* OPEN INPUT FILE.                                                 */
    if ((fIN=fopen(argv[1],"r")) == NULL) {
        printf("dup: can't open %s\n",argv[1]);
        exit(2);
    }

    /* **************************************************************** */
    /* COUNT NUMBER OF BYTES IN INPUT FILE                              */
    num_char = 0;
    while (fgetc(fIN)!=EOF)
        num_char++;

    /* **************************************************************** */
    /* ALLOCATE MEMORY FOR CHARACTER DATA, READ IT IN.                  */
    char_array=int_malloc(num_char);
    rewind(fIN);
    for (char_index=0; char_index<num_char; char_index++) {
        c=fgetc(fIN);
        char_array[char_index]=c;
    }
    fclose(fIN);

    /* **************************************************************** */
    /* OPEN OUTPUT FILE, DUPLICATE INPUT FILE num_dups TIMES.           */
    if ((fOUT=fopen(argv[2],"w")) == NULL) {
        printf("dup: can't open %s\n",argv[2]);
        exit(3);
    }
    for (dup_index=0; dup_index<num_dups; dup_index++) {
        for (char_index=0; char_index<num_char; char_index++)
            fputc(char_array[char_index],fOUT);
    }
    fclose(fOUT);
}
