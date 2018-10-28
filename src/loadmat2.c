#include <stdio.h>
#include <stdlib.h>

typedef struct {
     long int type;   /* type */
     long int mrows;  /* row dimension */
     long int ncols;  /* column dimension */
     long int imagf;  /* flag indicating imag part */
     long int namlen; /* name length (including NULL) */
} Fmatrix;

int loadmat2(FILE *fp, long int *type, char *pname, long int *mrows,
		long int *ncols, long int *imagf, double **preal, double **pimag)
{
    Fmatrix x;
    long int mn, namlen;
    int readstatus;


    /* Get Fmatrix structure from file */
    if ((readstatus=fread((char *)&x, sizeof(Fmatrix), 1, fp)) != 1) {
	printf("loadmat2: fread = %d\n",readstatus);
	return 1;
    }

    *type = x.type;
    *mrows = x.mrows;
    *ncols = x.ncols;
    *imagf = x.imagf;
    namlen = x.namlen;
    mn = x.mrows * x.ncols;

    /* Get matrix name from file */
    if (fread(pname, sizeof(char), (int) namlen, fp) != namlen) {
	printf("loadmat2: bad namelength.\n");
	return 2;
    }

    /* Get Real part of matrix from file */
    if (!(*preal = (double *)malloc((int) mn*sizeof(double)))) {
	printf("\nError: Variable too big to load\n");
	printf("\nError: Variable too big to load\n");
	printf("mn = %d\n",mn);
	printf("mrows = %d\n",*mrows);
	printf("ncols = %d\n",*ncols);
	return 3;
    }
    
    if (fread(*preal, sizeof(double), (int) mn, fp) != mn) {
	printf("loadmat2: bad dimensions.\n");
	free(*preal);
	return 4;
    }

    /* Get Imag part of matrix from file, if it exists */
    if (x.imagf) {
/*	printf("Imaginary data!\n"); */
	if (!(*pimag = (double *)malloc((int) mn*sizeof(double)))) {
	    printf("\nError: Variable too big to load\n");
	    printf("mn = %d\n",mn);
	    printf("mrows = %d\n",*mrows);
	    printf("ncols = %d\n",*ncols);
	    free(*preal);
	    return 5;
	}
	if (fread(*pimag, sizeof(double), (int) mn, fp) != mn) {
	    printf("loadmat2: bad dimensions.\n");
	    free(*pimag);
	    free(*preal);
	    return 6;
	}
    }

    return 0;
}
