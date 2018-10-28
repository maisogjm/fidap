#include <stdio.h>
#include <stdlib.h>

main (int argc, char *argv[])
{
    /* 01.24.94 : JM Maisog
	Writes out a "mac" header.					*/



    FILE *fOUT;
    signed short int *header,max,min,xdim,ydim;
    int write_status;


    /* **************************************************************** */
    /* VERIFY NUMBER OF COMMAND LINE ARGUMENTS.				*/
    if (argc < 6) {
	printf("make_mac_header: need output file, data max, data min,\n");
	printf("                 xdim, AND ydim.\n");
	exit(1);
    }

    /* **************************************************************** */
    /* OPEN FILE FOR OUTPUT.						*/
    if ((fOUT=fopen(argv[1],"w")) == NULL) {
	printf("make_mac_header: can't open %s\n",argv[2]);
	exit(3);
    }

    /* **************************************************************** */
    /* ASSIGN MEMORY.							*/
    header = (signed short int *) malloc(4*sizeof(signed short int));
    if (header == NULL) {
	printf("make_mac_header: failed to malloc for header.\n");
	exit(5);
    }

    /* **************************************************************** */
    /* FILL IN HEADER INFORMATION.					*/
    header[0]=(signed short int) atoi(argv[2]);
    header[1]=(signed short int) atoi(argv[3]);
    header[2]=(signed short int) atoi(argv[4]);
    header[3]=(signed short int) atoi(argv[5]);

    /* **************************************************************** */
    /* WRITE OUT UNSIGNED 8 BIT TO OUTPUT FILE				*/

    write_status=fwrite(header, sizeof(signed short int), 4, fOUT);
    fclose(fOUT);
    if (write_status != 4) {
        printf("%d header elements written.  Should have been 4.\n");
        exit(10);
    }


    /* **************************************************************** */
    /* CLOSE FILES							*/
    fclose(fOUT);
}
