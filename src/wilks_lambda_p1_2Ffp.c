#include <stdio.h>
#include <stdlib.h>

main(int argc, char *argv[])
{
    /* 01.20.96 : JM Maisog
        Converts a Wilks' Lambda map to an F-test map, using an approximate
        transformation.  See Rencher AC, Methods of Multivariate Analysis,
        New York: John Wiley & Sons, 1995, eq. (6.13), p. 182.                  */

    int num_voxX4;            /* Number of bytes in image.                      */
    int num_vox;              /* Number of voxels in image.                     */
    int vox_index;            /* Index into voxels in image.                    */

    FILE *fIN;                /* Pointer to input image file.                   */
    int read_status;          /* Check on read status of fread.                 */
    float *wilks_lambda;      /* Wilks' Lambda.                                 */
    float p, nu_H, nu_E;      /* Parameters of Wilks' Lambda.                   */

    float df_ratio;           /* df2/df1                                        */
    float *F_test;            /* Pointer to output probability image.           */
    FILE *fOUT;               /* Pointer to output image file.                  */
    int write_status;         /* Check on write status of fwrite.               */


/* ---------------------------------------------------------------------------- */
/*  CHECK TO SEE IF CORRECT NUMBER OF ARGUMENTS WERE PASSED.                    */

    if (argc!=5) {
        printf("wilks_lambda_p1_2Ffp: need input Wilks' Lambda map, nu-H, nu-E,\n");
        printf("                  and output filename.\n");

        printf("Will convert a Wilk's lambda to an F-test map using an EXACT\n");
        printf("Wilks' Lambda-to-F-test transform.  This program assumes\n");
        printf("that p = 1.\n");
        exit(1);
    }

/* ---------------------------------------------------------------------------- */
/*  READ IN COMMAND LINE ARGUMENTS.                                             */

    p    = 1.;
    nu_H = atof(argv[2]);
    nu_E = atof(argv[3]);
    printf("p    = %f\n",p);
    printf("nu_H = %f\n",nu_H);
    printf("nu_E = %f\n",nu_E);
    

/* ---------------------------------------------------------------------------- */
/* OPEN INPUT FILE.                                                             */

    if ((fIN=fopen(argv[1],"r")) == NULL) {
        printf("wilks_lambda_p1_2Ffp: can't open %s\n",argv[1]);
        exit(2);
    }

/* ---------------------------------------------------------------------------- */
/* COUNT NUMBER OF BYTES IN INPUT FILE.                                         */

    num_voxX4 = 0;
    while (fgetc(fIN)!=EOF)
        num_voxX4++;

    num_vox = num_voxX4/4;
    /* printf("wilks_lambda_p1_2Ffp: number of voxels = %d\n",num_vox); */

/* ---------------------------------------------------------------------------- */
/* ASSIGN MEMORY.                                                               */

    wilks_lambda = (float *) malloc(num_vox*sizeof(float));
    if (wilks_lambda == NULL) {
        printf("wilks_lambda_p1_2Ffp: failed to malloc for wilks_lambda.\n");
        exit(4);
    }
    F_test = (float *) malloc(num_vox*sizeof(float));
    if (F_test == NULL) {
    printf("wilks_lambda_p1_2Ffp: failed to malloc for F_test.\n");
    exit(5);
    }

/* ---------------------------------------------------------------------------- */
/* REWIND INPUT IMAGE FILE AND READ IN DATA.                                    */

    /* Wilks' Lambda Map. */
    rewind(fIN);
    read_status=fread(wilks_lambda, sizeof(float), num_vox, fIN);
    fclose(fIN);
    if (read_status != num_vox) {
        printf("wilks_lambda_p1_2Ffp: %d voxels read in from %s.  Should have been %d voxels.\n",
                                        read_status,argv[1],num_vox);
        exit(8);
    }

/* ---------------------------------------------------------------------------- */
/* COMPUTE df_ratio.                                                            */

    df_ratio = nu_E/nu_H;
    
/* ---------------------------------------------------------------------------- */
/* CONVERT WILKS' LAMBDA TO EXACT F-TEST.                                       */
    
    printf("wilks_lambda_p1_2Ffp: converting Wilks' Lambda's to F-tests...\n");
    for (vox_index=0; vox_index<num_vox; vox_index++) {
        if (wilks_lambda[vox_index]!=0.) {
            F_test[vox_index] = ((1.-wilks_lambda[vox_index])
                                        /wilks_lambda[vox_index])*df_ratio;
        }
        else F_test[vox_index]=0.;
    }

/* ---------------------------------------------------------------------------- */
/* WRITE TO OUTPUT FILE.                                                        */

    if ((fOUT=fopen(argv[4],"w")) == NULL) {
        printf("wilks_lambda_p1_2Ffp: can't open %s\n",argv[4]);
        exit(3);
    }
    write_status=fwrite(F_test, sizeof(float), num_vox, fOUT);
    fclose(fOUT);
    if (write_status != num_vox) {
        printf("wilks_lambda_p1_2Ffp: %d voxels written.  Should have been %d voxels.\n",write_status,
                                                                        num_vox);
        exit(10);
    }

    free(wilks_lambda);
    free(F_test);
}

