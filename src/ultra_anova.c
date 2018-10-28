#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "analyze6.h"


                       /* FUNCTION DECLARATIONS.                                */
struct  analyze_struct readhdr(char *argv[]);
float betai(float a, float b, float x);

/* **************************************************************************** */
/* PROGRAM ultra_anova.c                                                        */
/*                                                                              */
/* ANSI C code completed 03/07/94 by Jose' Ma. Maisog, M.D.                     */
/*                                                                              */
/* Uses C code from Numerical Recipes in C.                                     */
/*                                                                              */
/* Building 10, Rm 4C110                                                        */
/* NIH/NIMH/LPP/SFBI                                                            */
/* 9000 Rockville Pike                                                          */
/* Bethesda, Maryland  20892                                                    */
/*                                                                              */
/* Tel: (301) 402-0416                                                          */
/* Fax: (301) 402-0921                                                          */
/* email: joem@alw.nih.gov                                                      */
/*                                                                              */
/* Program reads in names of images from a text file, then reads in those       */
/* images.  Global mean of each image is adjusted to some weighted average of   */
/* the global means of the input images.  User is asked to input parameters     */
/* determining the experimental design, and a pixel-by-pixel ANOVA is then      */
/* computed.  Floating point images are generated.  A text file containing the  */
/* ANOVA for one voxel is also generated.                                       */
/*                                                                              */
/* HISTORY: 11.24.93 : JMM : started modifying fmri_pairedt_ratio.c to generate */
/*                       ultra_anova.c.  ultra_anova reads in 16-bit images and */
/*                       performs an N-WAY ANOVA.                               */
/*            01.19.94 : JMM : added code to ultra_anova to output in addition  */
/*                       the error SS image.                                    */
/*            03.04.94 : JMM : final debugging of ultra_anova.                  */
/*            03.07.94 : JMM : final touches on output.                         */
/*          11.23.94 : JMM : modified to output mean images of effects and of   */
/*                     interactions.                                            */
/*          03.14.95 : JMM : modified to read in an 8-bit mask, rather than     */
/*                     creating it on the fly.                                  */

main(int argc, char *argv[])
{
                            /* READING IN RAW DATA.                             */
    struct analyze_struct   /* Header information variable.                     */
                img_info;   /*                                                  */
    FILE *fimgin;           /* Pointer for first input image file.              */
    char choice;            /* For keyboard input to yes/no question.           */
    char dummy[80];         /* Dummy pointer for counting number of entries     */
                            /* in a text file.                                  */
    int num_read_in;        /* Number of voxels read into memory from disk.     */
    int total_n;            /* Total number of images.                          */
    int n_index;            /* Index indicating which image in group.           */
                            /* Range: [0,total_n-1]                             */
    char **img_name;        /* Pointer to pointers to input file names.         */
    int num_vox;            /* Number of voxels in an image = xdim*ydim*zdim.   */
    int vox_index;          /* Index into imgin[*][*][], extracting a specific  */
                            /* voxel value.  Range: [0,vox_index-1]             */
    int xdim,ydim,zdim;     /* Image dimensions.                                */
    signed short *temp;     /* Temporary buffer for reading in images.          */
    float **imgin;          /* Pointer to pointers to input image data in       */
                            /* floating point form.                             */

                            /* VOXEL OF INTEREST.                               */
    int x_int, y_int,       /* X, Y, and Z coordinates of a voxel of interest.  */
                z_int;      /* The ANOVA for this voxel of interest will be     */
                            /* output to a text file.                           */
    int vox_index_int;      /* The voxel index for the voxel of interest,       */
                            /* equal to (z_int*ydim+y_int)*xdim+x_int, after    */
                            /* x_int, y_int, and z_int have been decremented by */
                            /* 1 (pointer arithmetic starts counting at zero).  */

                            /* EXPERIMENTAL DESIGN.                             */
    int num_fact;           /* Number of factors in the experiment.             */
    int fact_index;         /* Index into factors.  Range: [0,num_fact-1]       */
    int fact_index2;        /* Secondary index into factors.                    */
                            /* Range: [0,num_fact-1]                            */
    int *num_level;         /* Number of levels each factor has.                */
    int lvl_index;          /* Index into levels.                               */
                            /* Range: [0,num_level[fact_index] ]                */
    int *lvl_index2;        /* Secondary index into levels.                     */
                            /* Range: [0,num_level[fact_index] ]                */
    int **lvl_equiv;        /* Gives equivalent factor levels for an image      */
                            /* index number n_index, given factor specified     */
                            /* by fact_index.                                   */
    int num_interactions;   /* There are (2^num_fact)-1 possible interactions   */
                            /* plus main effects.  This variable is that number.*/
    int interaction_index;  /* Index for looping over the different             */
                            /* interaction classes.  There is no interaction    */
                            /* associated with interaction_index = 0.           */
                            /* Range: [1,num_interactions]                      */
    int interaction_index2; /* Secondary index for looping over the different   */
                            /* interactions.  There is no interaction           */
                            /* associated with interaction_index = 0.           */
                            /* Range: [1,num_interactions]                      */
    int *num_types;         /* Number of types within a class of interactions.  */
                            /* For main effects, this is merely the number of   */
                            /* levels of a factor, and thus this variable is    */
                            /* closely related to the variable *num_level.      */
                            /* *num_types will need to be  calculated for       */
                            /* interactions.                                    */
    int type_index;         /* Index over the number of types within a class of */
                            /* interactions.                                    */
                            /* Range: [0,num_types[interaction_index] ]         */
    int ***type_equiv1;     /* Given an interaction index specified by          */
                            /* interaction_index, a lower-order interaction     */
                            /* specified by interaction_index2, and a           */
                            /* type/level within the interaction class          */
                            /* interaction_index, returns the corresponding     */
                            /* type/level of the interaction specified by       */
                            /* interaction_index2.                              */
    int **type_equiv2;      /* Given an interaction class specified by          */
                            /* interaction_index and an image index number      */
                            /* specified by n_index, returns the corresponding  */
                            /* type within that interaction class.  Closely     */
                            /* related to lvl_equiv.                            */
    int *num_fact_in_interaction;
                            /* Number of factors involved in interaction        */
                            /* specified by interaction_index.                  */
    int *df;                /* Degrees of freedom for a class of interactions.  */
    int **interaction_equiv;/* Given an interaction specified by                */
                            /* interaction_index and a factor specified by      */
                            /* fact_index, this array is set to 0 or 1, a 1     */
                            /* indicating that factor fact_index is involved in */
                            /* calculating interaction interaction_index, a     */
                            /* 0 indicating that it isn't.  Thus, if            */
                            /* interaction_equiv[2][3] is set to 1, then factor */
                            /* 3 is involved in interaction 2.                  */
    int **multiples;        /* Multiples in converting level equivalents to     */
                            /* interaction class types.                         */
    int num_cells;          /* Number of cells.                                 */
    int num_obs;            /* Number of observations per cell.                 */
    int remainder;          /* Remainder after performing modulus operations    */
                            /* on n_index, translating n_index to the           */
                            /* corresponding levels of the factors.             */
    int *level;             /* Corresponding level of factor, given n_index.    */
    int *major_inc;         /* Major increment for factors.                     */
    int *num_cell_per_level;/* Number of cells a factor has per level.          */
                            /* This is merely total_n/num_level[*].             */
    int offset;             /* Offset used to dereference names and images.     */
                            /* Necessary because of variable number of factors  */
                            /* and levels.                                      */

                            /* INTENSITY SCALING & GLOBAL MEAN NORMALIZATION    */
    FILE *fMASK;            /* Pointer for input mask image.                    */
    float prob_thresh;      /* Probability threshold for significance.          */
    unsigned char *mask;    /* Pointer to mask image.                           */
    int num_mask_vox;       /* Number of non-zero voxels in *mask.              */
    float *glob_mean;       /* Global mean intensity for images.                */
    float float_c;          /* Temporary variable for reading in *c.            */
    float *c;               /* Weighting factor for normalizing global means.   */
    float NUMRTR;           /* Numerator used for normalizing global means.     */
    float SUM;              /* Sum used for normalizing global means.           */

                            /* FOR OUTPUT TO TEXT FILE.                         */
    char outtext_name[80];  /* Name of output text file.                        */
    FILE *fOUTTEXT;         /* Pointer to output text file.                     */

                            /* CALCULATION OF F-IMAGE AND OUTPUT TO DISK.       */
    float grandmean;        /* Grand mean, used in calculating sums of squares. */
    float **interaction_mean;
                            /* First used to store interaction sums; then       */
                            /* adjusted to contain interaction means by dividing*/
                            /* by the number of cells in each interaction       */
                            /* The interaction terms for the 4th level of       */
                            /* interaction class 2 will be kept in              */
                            /* interaction[2][3].  Classes will be 1-offset,    */
                            /* whereas levels are zero-offset.  Thus            */
                            /* interaction[0] is unused.                        */
    float ***interaction_mean_store;
                            /* Added so that interaction mean images could be   */
                            /* output as well.                                  */
    float **interaction_effects;
                            /* Generated by subtracting lower order interaction */
                            /* means and the grandmean from the factor mean of  */
                            /* the interaction effect in question.  Interaction */
                            /* classes will be 1-offset, because there is no    */
                            /* interaction associated with                      */
                            /* interaction_index = 0; whereas levels are zero-  */
                            /* offset, e.g., the interaction effects for the    */
                            /* 4th level of interaction class 2 will be kept in */
                            /* interaction[2][3].                               */
    int num_factors_index;  /* Index into the number of factors involved in a   */
                            /* particular interaction.                          */
                            /* Range: [0,                                       */
                            /*    num_fact_in_interaction[interaction_index]-1] */
    int type;               /* Given an interaction class specified by          */
                            /* interaction_index, and given a data point        */
                            /* specified by n_index = [0,total_n-1], used to    */
                            /* map what type the data point is mapped to.       */
    int subset_flag;        /* Binary flag indicating whether or not a lower    */
                            /* order interaction's factors are a subset of a    */
                            /* higher order interaction's factors.              */
    int neg_one_exponent;   /* Exponent for -1, to determine the appropriate    */
                            /* sign for lower order interactions when factoring */
                            /* out their contributions to higher order          */
                            /* interactions.                                    */
    int exp_index;          /* Loop index for determining the appropriate sign  */
                            /* for lower order interactions.                    */
    float sign;             /* Sign for lower order interactions, +1 or -1.     */
    float *residual_error;  /* Residuals for interactions.  Could have used     */
                            /* the unused space in interaction_effects[0], but  */
                            /* for clarity's sake shall allocate a few          */
                            /* additional bytes for this term.                  */
    float *MSE;             /* MS's of the residual errors.                     */
    float *SS;              /* Sums-of-squares, one per interaction.  SS is     */
                            /* divided by the appropriate df, to yield MEAN     */
                            /* SQUARE ERRORS, which will still be stored in SS. */
                            /* SS for total SS is stored in MSE.                */
    float **F;              /* F images, one per factor.                        */
    float v1,v2;            /* Degrees of freedom for F test.                   */
    float Prob;             /* Probability, given F and degrees of freedom.     */
    float *probability;     /* Pointer to floating point output for probability */
                            /* map output.                                      */
    int x,y,z;              /* Indices into ANALYZE format image arrays.        */
    char Fout[30];          /* Pointer to name of output t map file.            */
    FILE *fOUT;             /* Pointer for output parameter .mat file.          */

/* ---------------------------------------------------------------------------- */
/*  CHECK TO SEE IF CORRECT NUMBER OF ARGUMENTS WERE PASSED.                    */

    if (argc < 7) {
        printf("ultra_anova: need header file, input textfile, name of 8-bit mask\n");
        printf("             image file, and X, Y,& Z coordinates of a voxel of interest.\n\n");

        printf("  e.g., ultra_anova test.hdr text.txt mask_10.img 35 4 5\n\n");

        printf("ultra_anova will perform an N-WAY ANOVA on a\n");
        printf("group of images, on a voxel-by-voxel basis.\n");
        printf("Any number of factors and levels are allowed, but\n");
        printf("each level within a factor must contain the same\n");
        printf("number of elements, and there must be no missing cells.\n");

        printf("\nThe input textfile should contain the pathnames\n");
        printf("(either absolute or relative) of the input images.\n") ;
        printf("Each name should be followed by a space and then by\n");
        printf("a weighting factor, described below.\n\n");
        printf("The ordering of the names in the list is *critical*\n");
        printf("and should reflect the experimental design.\n");
        printf("Suppose that there were a total of N images, and\n");
        printf("that the first factor had N1 levels.  Then the first\n");
        printf("N/N1 names are those images in the first level of\n");
        printf("the first factor, the next N/N1 are the images in\n");
        printf("the second level of the first factor, and so on.\n");
        printf("*Within* each group of N/N1 images, the names are\n");
        printf("grouped in a similar manner for the second factor,\n");
        printf("and so on.\n\n");

        printf("\nTo correct for differences in global mean intensity,\n");
        printf("the following multiplicative scaling will be done on\n");
        printf("an image-by-image basis.  A global mean GBM(i) will\n");
        printf("be calculated for each image IMG(i).  Then each image\n");
        printf("IMG(i) will be multiplied by\n\n");

        printf("                  NUMRTR/GBM(i)\n\n");

        printf("where NUMRTR is a weighted average of all the GBM's.\n");
        printf("Let c(i) be the weight assigned to GBM(i).  Let\n\n");

        printf("           SUM = abs(c(1))+abs(c(2))+...\n\n");

        printf("for all the c(i)`s.  Then\n\n");

        printf("         NUMRTR = [c(1)*GBM(1)+c(2)*GBM(2)+...]/SUM.\n\n");

        printf("For example, if you want the numerator to be the\n");
        printf("average of GBM(1) and GBM(6), set\n\n");

        printf("          c(1) = 1.0\n");
        printf("          c(2) = 0.0\n");
        printf("          c(3) = 0.0\n");
        printf("          c(4) = 0.0\n");
        printf("          c(5) = 0.0\n");
        printf("          c(6) = 1.0\n\n");

        printf("An 8-bit image of 0's and 1's must also be specifed.\n");
        printf("Voxels set to 1 in the mask will be counted as within\n");
        printf("brain, and calculations will be restricted to those voxels.\n\n");

        printf("The output text file will contain the\n");
        printf("intermediate variables c, NUMRTR, glob_mean, and\n");
        printf("threshold (the number used to generate the mask for\n");
        printf("calculation of the global means).  It will also\n");
        printf("contain the ANOVA for the voxel of interest.\n");

        exit(1);
    }

/* ---------------------------------------------------------------------------- */
/*  READ HEADER FILE; ASSUME IMAGES ARE SAME SIZE, SO USE SAME HEADER FILE      */
/*  FOR ALL IMAGES.                                                             */

    img_info = readhdr(argv+1);

    xdim = img_info.dime.dim[1];
    ydim = img_info.dime.dim[2];
    zdim = img_info.dime.dim[3];

    num_vox = xdim*ydim*zdim;

/* ---------------------------------------------------------------------------- */
/* DETERMINE VOXEL OF INTEREST COORDINATES AND NAME OF OUTPUT TEXT FILE.        */
    x_int = atoi(argv[4]);
    y_int = atoi(argv[5]);
    z_int = atoi(argv[6]);

    vox_index_int = ((z_int-1)*ydim+(y_int-1))*xdim+(x_int-1);

    sprintf(outtext_name,"anova_%d_%d_%d.txt\0",x_int,y_int,z_int);

    if ((fOUTTEXT = fopen(outtext_name,"w")) == NULL) {
        printf ("ultra_anova: can't open %s.\n",outtext_name);
        exit(2);
    }
    fprintf(fOUTTEXT,"ANOVA FOR VOXEL AT COORDINATES [%d,%d,%d]\n\n\n",
                                                        x_int,y_int,z_int);
    fprintf(fOUTTEXT,"Voxel index of interest = %d\n",vox_index_int);
    if ((vox_index_int<0) || (vox_index_int)>=num_vox) {
        printf("ultra_anova: bad coordinates of interest, given image\n");
        printf("             dimensions in %s.\n",argv[1]);
        fclose(fOUTTEXT);
        exit(3);
    }

/* ---------------------------------------------------------------------------- */
/*  COUNT NUMBER OF ENTRIES IN INPUT TEXT FILE.                                 */

    if ((fimgin = fopen(argv[2],"r")) == NULL) {
        printf ("ultra_anova: can't open %s\n",argv[2]);
        exit(4);
    }
    total_n = 1;
    fscanf(fimgin,"%s",dummy);
    while(fscanf(fimgin,"%s %s %f",dummy,dummy,&float_c)!=EOF)
        total_n++;
    printf("There are %d images.\n",total_n);

/* ---------------------------------------------------------------------------- */
/*  READ IN FILENAMES AND GLOBAL MEAN WEIGHTS IN INPUT TEXT FILE.               */

    c = (float *) malloc(total_n*sizeof(float));
    if (c==NULL) {
        printf("ultra_anova: error malloc'ing for c.\n");
        exit(5);
    }
    img_name = (char **) malloc(total_n*sizeof(char *));
    if (img_name==NULL) {
        printf("ultra_anova: error malloc'ing for img_name.\n");
        exit(6);
    }

    rewind(fimgin);
    for (n_index=0; n_index<=total_n-1; n_index++) {
        img_name[n_index] = (char *) malloc(80*sizeof(char));
        if (img_name[n_index]==NULL) {
            printf("ultra_anova: error malloc'ing for img_name.\n");
            exit(7);
        }
        fscanf(fimgin,"%s %f", img_name[n_index],c+n_index);
    }
    fclose(fimgin);

/* ---------------------------------------------------------------------------- */
/*  QUERY USER FOR FACTORS AND LEVELS.                                                */

    printf("How many factors does this experiment have? ");
    scanf("%d",&num_fact);
    num_level = (int *) malloc((num_fact+1)*sizeof(int));
    num_cell_per_level  = (int *) malloc(num_fact*sizeof(int));
    if ((num_level==NULL) || (num_cell_per_level==NULL)) {
        printf("ultra_anova: error malloc'ing for num_level or num_cell_per_level.\n");
        exit(8);
    }
    num_cells = 1;
    for (fact_index=0; fact_index<num_fact; fact_index++) {
        printf("How many levels does factor #%d have? ", fact_index+1);
        scanf("%d",num_level+fact_index);
        num_cell_per_level[fact_index] = total_n / num_level[fact_index];
        printf("  factor # %d has %d cells per level.\n",fact_index+1,
                                        num_cell_per_level[fact_index]);
        num_cells *= num_level[fact_index];
    }
    num_level[num_fact]=total_n;
    num_obs = total_n / num_cells;
    printf("\nThere are a total of %d cells in this study, with %d observations per cell.\n",
                                        num_cells,num_obs);

    /* ------------------------------------------------------------------------ */
    /* CALCULATE MAJOR MINOR INCREMENTS FOR FACTORS.                            */
    /* Major increment increments to the same level of the factor, but the next */
    /* level of the next-greatest factor.                                       */

    major_inc = (int *) malloc(num_fact*sizeof(int));
    if (major_inc==NULL) {
        printf("ultra_anova: error malloc'ing for major_inc.\n");
        exit(9);
    }
    for (fact_index=num_fact-1; fact_index>=0; fact_index--) {
        if (fact_index==num_fact-1)
            major_inc[fact_index] = num_obs;
        else
            major_inc[fact_index] = num_level[fact_index+1] * major_inc[fact_index+1];
    }

/* ---------------------------------------------------------------------------- */
/* CALCULATE CLASSES OF INTERACTIONS, AND NUMBER OF TYPES WITHIN EACH CLASS.    */
/* Here determine the experimental design based on information input so far.    */

    /* ------------------------------------------------------------------------ */
    /* CALCULATE NUMBER OF INTERACTIONS.                                        */
    num_interactions = 1;
    for (fact_index=1; fact_index<=num_fact; fact_index++)
        num_interactions *= 2;
    num_interactions--;
    printf("\nNumber of possible interactions plus main effects = %d.\n",
        num_interactions);

    printf("\nDoes the above information look okay <y/n>? ");
    scanf("%s",&choice);
    while((choice!='y') && (choice!='n')) {
    printf("\n  Type 'y' or 'n': ");
        scanf("%s",&choice);
    }
    if (choice=='n') exit(10);


    /* ------------------------------------------------------------------------ */
    /* ALLOCATE MEMORY FOR interaction_equiv.                                        */
    interaction_equiv = (int **) malloc(num_interactions*sizeof(int*));
    if (interaction_equiv==NULL) {
        printf("ultra_anova: error malloc'ing for interaction_equiv.\n");
        exit(11);
    }
    for (interaction_index=1; interaction_index<=num_interactions;
                                                                interaction_index++) {
        interaction_equiv[interaction_index]=
                        (int *) malloc(num_fact*sizeof(int));
        if (interaction_equiv[interaction_index]==NULL) {
            printf("ultra_anova: error malloc'ing for interaction_equiv.\n");
            exit(12);
        }
    }

    /* ------------------------------------------------------------------------ */
    /* CALCULATE interaction_equiv'S                                            */
    for (interaction_index=1; interaction_index<=num_interactions;
                                                                interaction_index++)
        for (fact_index=0; fact_index<num_fact; fact_index++)
            interaction_equiv[interaction_index][fact_index]
                                = (interaction_index >> fact_index) & 1;

    /* ------------------------------------------------------------------------ */
    /* CALCULATE NUMBER OF TYPES WITHIN EACH CLASS OF INTERACTION.              */
    /* ALSO, CALCULATE NUMBER OF DEGREES OF FREEDOM ASSOCIATED WITH EACH        */
    /* INTERACTION.                                                             */
    num_types = (int *) malloc(num_interactions*sizeof(int));
    if (num_types==NULL) {
        printf("ultra_anova: error malloc'ing for num_types.\n");
        exit(13);
    }
    df = (int *) malloc(num_interactions*sizeof(int));
    if (df==NULL) {
        printf("ultra_anova: error malloc'ing for df.\n");
        exit(14);
    }
    for (interaction_index=1; interaction_index<=num_interactions;
                                                                interaction_index++) {
        num_types[interaction_index]=1;
        for (fact_index=0; fact_index<num_fact; fact_index++)
            if (interaction_equiv[interaction_index][fact_index]==1)
                num_types[interaction_index] *= num_level[fact_index];

        df[interaction_index]=1;
        for (fact_index=0; fact_index<num_fact; fact_index++)
            if (interaction_equiv[interaction_index][fact_index]==1)
                df[interaction_index] *= (num_level[fact_index]-1);
    }

    /* ------------------------------------------------------------------------ */
    /* CALCULATE df FOR RESIDUALS.  THEN SHOW USER INTERACTIONS AND df's.       */
    df[0] = total_n-1;
    for (interaction_index=1; interaction_index<=num_interactions;
                                                                interaction_index++)
        df[0]-=df[interaction_index];

    fprintf(fOUTTEXT,"Interaction table:\n\n");
    fprintf(fOUTTEXT,"  #  Factor On/Off\n");
    for (interaction_index=1; interaction_index<=num_interactions;
                                                                interaction_index++) {
        fprintf(fOUTTEXT,"  %d : ",interaction_index);
        for (fact_index=num_fact-1; fact_index>=0; fact_index--)
            fprintf(fOUTTEXT," %d",interaction_equiv[interaction_index][fact_index]);
        fprintf(fOUTTEXT," : Number of types = %d",num_types[interaction_index]);
        fprintf(fOUTTEXT," : df = %d\n",df[interaction_index]);
    }
    fprintf(fOUTTEXT,"\n");
    fprintf(fOUTTEXT,"df of MSE is %d.\n\n",df[0]);


    /* ------------------------------------------------------------------------ */
    /* CALCULATE NUMBER OF FACTORS INVOLVED IN EACH INTERACTION.                */
    /* Count number of 1's in interaction_equiv[interaction_index][].           */
    /* This is the number of factors involved.                                  */
    num_fact_in_interaction = (int *) malloc(num_interactions*sizeof(int *));
    if (num_fact_in_interaction==NULL) {
        printf("ultra_anova: error malloc'ing for num_fact_in_interaction.\n");
        exit(15);
    }

    for (interaction_index=1; interaction_index<=num_interactions; interaction_index++) {
        num_fact_in_interaction[interaction_index]=0;
        for (fact_index=num_fact-1; fact_index>=0; fact_index--)
            if (interaction_equiv[interaction_index][fact_index]==1)
                num_fact_in_interaction[interaction_index]++;

        fprintf(fOUTTEXT,"  There are %d factors in interaction %d.\n",
                        num_fact_in_interaction[interaction_index],interaction_index);
    }

    /* ------------------------------------------------------------------------ */
    /* CALCULATE LEVEL EQUIVALENCES.                                            */
    /* lvl_equiv maps n_index to a level of a given factor, specified by        */
    /* fact_index.                                                              */

    lvl_index2 = (int *) malloc(num_fact*sizeof(int));
    if ((lvl_index2==NULL)) {
        printf("ultra_anova: error malloc'ing for lvl_index2.\n");
        exit(16);
    }
    lvl_equiv = (int **) malloc(total_n*sizeof(int*));
    if ((lvl_equiv==NULL)) {
        printf("ultra_anova: error malloc'ing for lvl_equiv.\n");
        exit(17);
    }
    for (n_index=0; n_index<total_n; n_index++) {
        lvl_equiv[n_index] = (int *) malloc(num_fact*sizeof(int));
        if (lvl_equiv[n_index]==NULL) {
            printf("ultra_anova: error malloc'ing for lvl_equiv.\n");
            exit(18);
        }
    }

    fprintf(fOUTTEXT,"\n\nThese are the number of levels each of the %d factors has:\n",num_fact);
    for (fact_index=0; fact_index<num_fact; fact_index++) {
        fprintf(fOUTTEXT,"\n  Factor #%d has %d levels, %d cells per level.\n",
                        fact_index+1,num_level[fact_index],num_cell_per_level[fact_index]);
        for (lvl_index=0; lvl_index<num_level[fact_index]; lvl_index++) {
            fprintf(fOUTTEXT,"    Level #%d:\n",lvl_index+1);

            for (n_index=0; n_index<total_n; n_index++) {
                remainder = n_index;
                for (fact_index2=0; fact_index2<num_fact; fact_index2++) {
                    lvl_index2[fact_index2] = remainder / major_inc[fact_index2];
                    remainder = remainder % major_inc[fact_index2];
                }
                if (lvl_index2[fact_index]==lvl_index) {
                    lvl_equiv[n_index][fact_index] = lvl_index;
                    fprintf(fOUTTEXT,"      %s\n",img_name[n_index]);
                }
            }
        }
    }
    fprintf(fOUTTEXT,"\nCheck on level equivalence:\n");
    for (n_index=0; n_index<total_n; n_index++) {
        fprintf(fOUTTEXT,"      %d : ",n_index+1);
        for (fact_index2=0; fact_index2<num_fact; fact_index2++)
            fprintf(fOUTTEXT,"%d ",lvl_equiv[n_index][fact_index2]+1);
        fprintf(fOUTTEXT,"\n");
    }


/* ---------------------------------------------------------------------------- */
/*  SHOW USER INFORMATION READ IN SO FAR.                                       */
/*  Allow user to quit program if incorrect information has been input.         */

    printf("\nThese are the %d images and weighting factors:\n", total_n);
    for (n_index=0; n_index<=total_n-1; n_index++)
        printf("  %20s  %f\n", img_name[n_index], c[n_index]);
    printf("\nDoes the above information look okay <y/n>? ");
    scanf("%s",&choice);
    while((choice!='y') && (choice!='n')) {
        printf("\n  Type 'y' or 'n': ");
        scanf("%s",&choice);
    }
    if (choice=='n') exit(19);

    printf("\nImages have dimensions:\n");
    printf("  xdim = %d\n  ydim = %d\n  zdim = %d\n",xdim,ydim,zdim);

    printf("\nThe F-image files will be stored in files named\n");
    printf("F?_d1_d2_img, where ? is the number of the\n");
    printf("factor/interaction in a convenient binary\n");
    printf("representation, and d1 and d2 are the numerator\n");
    printf("and denominator degrees of freedom.  They are\n");
    printf("floating point flat files which can be treated as\n");
    printf("floating point ANALYZE images.\n");

/*  printf("\nIt is assumed that the input images are registered.  A mean\n");
    printf("image will be calculated.  A mask will be made from this mean\n");
    printf("image by thresholding it.  Below, you will need to specify the\n");
    printf("number T.  All voxels in the mean image which have pixel\n");
    printf("intensity T %% the maximum voxel intensity in the mean image or\n");
    printf("greater will be set to 1 in the mask.\n");
    printf("E.g., if you want a 10 %% threshold, type in <10>.\n");
    printf("This mask will be used to calculate the global mean\n");
    printf("intensities of the individual images.\n");
    printf("\nENTER the threshold T for the mask: ");
    scanf("%f",&threshold);
    printf("\nThreshold = %3.2f %%.\n",threshold);
*/
    printf("\nPlease enter probability threshold (e.g., 0.05): ");
    scanf("%f",&prob_thresh);
    printf("\nProbability threshold = %f.\n",prob_thresh);

    printf("\nANOVA for the voxel of interest at [%d,%d,%d] will be\n",
                                                        x_int,y_int,z_int);
    printf("written to the text file %s.\n",outtext_name);

    printf("\nDoes the above information look okay <y/n>? ");
    scanf("%s",&choice);
    while((choice!='y') && (choice!='n')) {
        printf("\n  Type 'y' or 'n': ");
        scanf("%s",&choice);
    }
    if (choice=='n') exit(20);
    printf("\n");

/* ---------------------------------------------------------------------------- */
/*  READ IN IMAGES.                                                             */

    temp = (signed short *) malloc(num_vox*sizeof(unsigned short));
    imgin = (float **) malloc(total_n*sizeof(float *));

    if (temp==NULL) {
        printf("ultra_anova: error malloc'ing for temp.\n");
        exit(21);
    }
    if (imgin==NULL) {
        printf("ultra_anova: error malloc'ing for imgin.\n");
        exit(21);
    }

    for (n_index=0; n_index<total_n; n_index++) {
        imgin[n_index] = (float *) malloc(num_vox*sizeof(float));
        if (imgin[n_index]==NULL) {
            printf("ultra_anova: error malloc'ing for imgin.\n");
            exit(22);
        }
        printf("ultra_anova: reading in image from %s....\n",img_name[n_index]);
        if ((fimgin=fopen(img_name[n_index],"r")) == NULL) {
            printf("ultra_anova: can't open %s\n",img_name[n_index]);
            exit(23);
        }
        num_read_in=fread(temp, sizeof(signed short), num_vox, fimgin);
        fclose(fimgin);
        if(num_read_in!=num_vox) {
            printf("fmri_lagNdisp4: %d voxels read in; should have been %d.\n",
                                                    num_read_in,num_vox);
            exit(30);
        }
        for (vox_index=0; vox_index<num_vox; vox_index++)
            imgin[n_index][vox_index]=(float)temp[vox_index];
    }

/* ---------------------------------------------------------------------------- */
/*  ALLOCATE MEMORY FOR OUTPUT IMAGES.                                          */

    /* ------------------------------------------------------------------------ */
    /* ALLOCATE MEMORY FOR glob_mean.                                           */

    glob_mean = (float *) malloc(total_n*sizeof(float));
    if (glob_mean==NULL) {
        printf("ultra_anova: error malloc'ing for glob_mean.\n");
        exit(25);
    }


    /* ------------------------------------------------------------------------ */
    /* READ IN MASK FROM DISK.  THEN COUNT NUMBER OF NON-ZERO VOXELS IN *mask.  */

    mask = (unsigned char *) malloc(num_vox*sizeof(unsigned char));
    if (mask==NULL) {
	printf("ultra_anova: error malloc'ing for mask.\n");
	exit(26);
    }
    if ((fMASK = fopen(argv[3],"r")) == NULL) {
        printf("fmri_lagNdisp4: can't open %s\n",argv[3]);
        exit(29);
    }
    printf("ultra_anova: reading in mask from %s...\n",argv[3]);
    num_read_in=fread(mask,sizeof(unsigned char),num_vox,fMASK);
    fclose(fMASK);
    if(num_read_in!=num_vox) {
        printf("fmri_lagNdisp4: %d voxels read in; should have been %d.\n",
                                                num_read_in,num_vox);
        exit(30);
    }

    num_mask_vox = 0.;
    for (vox_index=0; vox_index<num_vox; vox_index++)
        if (mask[vox_index] != 0.)
            num_mask_vox++;

    /* ------------------------------------------------------------------------ */
    /* NOW NORMALIZE IMAGES.                                                    */

        /* -------------------------------------------------------------------- */
        /* CALCULATE GLOBAL MEAN OF EACH INPUT IMAGE USING mask.                */
        printf("ultra_anova: calculating global means...\n");
        for (n_index=0; n_index<total_n; n_index++) {
            glob_mean[n_index]=0.;
            for (vox_index=0; vox_index<num_vox; vox_index++)
                glob_mean[n_index] +=
                        imgin[n_index][vox_index] * mask[vox_index];
            glob_mean[n_index] /= ((float) num_mask_vox);
            printf("  global mean of %s = %f\n", img_name[n_index],
                                                glob_mean[n_index]);
        }

        /* -------------------------------------------------------------------- */
        /* CALCULATE SUM.                                                       */

        SUM = 0;
        for (n_index=0; n_index<total_n; n_index++)
            SUM+=((float) (abs(c[n_index])));

        /* -------------------------------------------------------------------- */
        /* CALCULATE NUMRTR.                                                    */

        NUMRTR = 0.;
        for (n_index=0; n_index<total_n; n_index++)
            NUMRTR += (c[n_index] * glob_mean[n_index]);

        NUMRTR /= SUM;
            
        /* -------------------------------------------------------------------- */
        /* NORMALIZE INDIVIDUAL IMAGES.  MULTIPLY EACH VOXEL BY                 */
        /* NUMRTR/glob_mean[i].                                                 */
        /* ALSO MASK OUT IMAGES, USING mask.                                    */

        printf("ultra_anova: normalizing image means...\n");
        for (n_index=0; n_index<total_n; n_index++)
            for (vox_index=0; vox_index<num_vox; vox_index++)
                imgin[n_index][vox_index] *=
                        ((float)mask[vox_index]*NUMRTR/glob_mean[n_index]);

    /* ------------------------------------------------------------------------ */
    /* WRITE NUMRTR, glob_mean, GLOBAL MEAN WEIGHTS, THRESHOLD FOR MASK, AND        */
    /* PROBABILITY THRESHOLD TO OUTPUT TEXT FILE.                               */

    fprintf(fOUTTEXT,"\nNUMRTR = %f\n",NUMRTR);
    fprintf(fOUTTEXT,"\nglobal means and weights for global means:\n");
    for (n_index=0; n_index<total_n; n_index++)
        fprintf(fOUTTEXT,"  %s %f %f\n",img_name[n_index],glob_mean[n_index],
                                                                        c[n_index]);

    fprintf(fOUTTEXT,"\nprobability threshold = %f\n",prob_thresh);
    fprintf(fOUTTEXT,"\nnumber of voxels included in mask = %d\n",num_mask_vox);
    fprintf(fOUTTEXT,"\nnumber of voxels NOT included mask = %d\n\n",num_vox-num_mask_vox);

    fprintf(fOUTTEXT,"Number of voxels in mask on per slice basis:\n");
    for (z=0; z<zdim; z++) {
        num_mask_vox=0;
        for (y=0; y<ydim; y++)
            for(x=0; x<xdim; x++)
                if(mask[(z*ydim+y)*xdim+x]==1.)
                    num_mask_vox++;
        fprintf(fOUTTEXT,"  slice # %2d : %d\n",z+1,num_mask_vox);
    }
    fprintf(fOUTTEXT,"\n\n");


    /* ------------------------------------------------------------------------ */
    /* ALLOCATE MEMORY FOR multiples, grandmean, INTERACTION MEANS & EFFECTS,   */
    /* TYPE EQUIVALENTS, SUMS-OF-SQUARES, AND F IMAGES.                         */
        printf("ultra_anova: allocating memory for ANOVA calculations...\n");

        interaction_mean = (float **) malloc((num_interactions+1)*sizeof(float *));
        if (interaction_mean==NULL) {
            printf("ultra_anova: error malloc'ing for interaction_mean.\n");
            exit(29);
        }
        for (interaction_index=1; interaction_index<=num_interactions; interaction_index++) {
            interaction_mean[interaction_index] =
                (float *) malloc(num_types[interaction_index]*sizeof(float));
            if (interaction_mean[interaction_index]==NULL) {
                printf("ultra_anova: error malloc'ing for interaction_mean.\n");
                exit(30);
            }
        }
        interaction_mean_store = (float ***) malloc((num_interactions+1)*sizeof(float **));
        if (interaction_mean_store==NULL) {
            printf("ultra_anova: error malloc'ing for interaction_mean_store.\n");
            exit(29);
        }
        for (interaction_index=1; interaction_index<=num_interactions; interaction_index++) {
            interaction_mean_store[interaction_index] =
                (float **) malloc(num_types[interaction_index]*sizeof(float *));
            if (interaction_mean_store[interaction_index]==NULL) {
                printf("ultra_anova: error malloc'ing for interaction_mean_store.\n");
                exit(30);
            }
            for (type_index=0; type_index<num_types[interaction_index]; type_index++) {
                interaction_mean_store[interaction_index][type_index]=
                    (float *) malloc(num_vox*sizeof(float));
                if (interaction_mean_store[interaction_index][type_index]==NULL) {
                    printf("ultra_anova: error malloc'ing for interaction_mean_store.\n");
                    exit(31);
                }
                for (vox_index=0; vox_index<num_vox; vox_index++)
                    interaction_mean_store[interaction_index][type_index][vox_index]=0.;
            }
        }
        interaction_effects =
                (float **) malloc((num_interactions+1)*sizeof(float *));
        if (interaction_effects==NULL) {
            printf("ultra_anova: error malloc'ing for interaction.\n");
            exit(31);
        }
        for (interaction_index=1; interaction_index<=num_interactions;
                                                        interaction_index++) {
            interaction_effects[interaction_index] =
                (float *) malloc(num_types[interaction_index]*sizeof(float));
            if (interaction_effects[interaction_index]==NULL) {
                printf("ultra_anova: error malloc'ing for interaction_effects.\n");
                exit(32);
            }
        }

        /* -------------------------------------------------------------------- */
        /* type_equiv1                                                          */

        type_equiv1 = (int ***) malloc((num_interactions+1)*sizeof(int **));
        if (type_equiv1==NULL) {
            printf("ultra_anova: error malloc'ing for type_equiv1.\n");
            exit(33);
        }
        for (interaction_index=1;interaction_index<=num_interactions;
                                                        interaction_index++) {
            type_equiv1[interaction_index] =
                (int **) malloc((num_interactions+1)*sizeof(int *));
            if (type_equiv1[interaction_index]==NULL) {
                printf("ultra_anova: error malloc'ing for type_equiv1.\n");
                exit(34);
            }
            for (interaction_index2=1; interaction_index2<=num_interactions;
                                                        interaction_index2++) {
                type_equiv1[interaction_index][interaction_index2] =
                    (int *) malloc(num_types[interaction_index]*sizeof(int));
                if (type_equiv1[interaction_index][interaction_index2]==NULL) {
                    printf("ultra_anova: error malloc'ing for type_equiv1.\n");
                    exit(35);
                }
            }
        }

        /* -------------------------------------------------------------------- */
        /* type_equiv2                                                          */
        type_equiv2 = (int **) malloc((num_interactions+1)*sizeof(int *));
        if (type_equiv2==NULL) {
            printf("ultra_anova: error malloc'ing for type_equiv2.\n");
            exit(36);
        }
        for (interaction_index=1;interaction_index<=num_interactions;
                                                        interaction_index++) {
            type_equiv2[interaction_index] = (int *) malloc(total_n*sizeof(int));
            if (type_equiv2[interaction_index]==NULL) {
                printf("ultra_anova: error malloc'ing for type_equiv2.\n");
                exit(37);
            }
        }

        residual_error = (float *) malloc(total_n*sizeof(float));
        if (residual_error==NULL) {
            printf("ultra_anova: error malloc'ing for residual_error.\n");
            exit(38);
        }
        MSE = (float *) malloc(num_vox*sizeof(float));
        if (MSE==NULL) {
            printf("ultra_anova: error malloc'ing for MSE.\n");
            exit(39);
        }
        SS = (float *) malloc((num_interactions+1)*sizeof(float));
        if (SS==NULL) {
            printf("ultra_anova: error malloc'ing for SS.\n");
            exit(40);
        }
        F = (float **) malloc((num_interactions+1)*sizeof(float *));
        if (F==NULL) {
            printf("ultra_anova: error malloc'ing for F.\n");
            exit(41);
        }
        for (interaction_index=1; interaction_index<=num_interactions; interaction_index++) {
            F[interaction_index] = (float *) malloc(num_vox*sizeof(float));
            if (F[interaction_index]==NULL){
                printf("ultra_anova: error malloc'ing for F.\n");
                exit(42);
            }
        }


        /* -------------------------------------------------------------------- */
        /* INITIALIZE INTERACTION TERMS.                                        */
        printf("ultra_anova: initializing interaction terms...\n");
        for (interaction_index=1;interaction_index<=num_interactions;
                                                        interaction_index++)
            for (type_index=0;type_index<num_types[interaction_index];type_index++)
                interaction_mean[interaction_index][type_index] = 0;


        /* -------------------------------------------------------------------- */
        /* CALCULATE multiples, USED IN CONVERTING LEVEL EQUIVALENTS TO FACTOR  */
        /* CLASS TYPES.                                                         */
        /*  Needed to map an image index n_index to a type within an            */
        /* interaction class, which is done below in the section entitled       */
        /*           SUM VOXEL VALUES INTO interaction_mean[][].....            */
        /* Not difficult to calculate, but conceptually tricky.                 */
        multiples = (int **) malloc((num_interactions+1)*sizeof(int));
        if (multiples==NULL) {
            printf("ultra_anova: error malloc'ing for multiples.\n");
            exit(43);
        }
        for (interaction_index=1;interaction_index<=num_interactions;
                                                        interaction_index++) {
            multiples[interaction_index] = (int *) malloc(num_fact*sizeof(int));
            if (multiples[interaction_index]==NULL) {
                printf("ultra_anova: error malloc'ing for multiples.\n");
                exit(44);
            }
        }
        for (interaction_index=1;interaction_index<=num_interactions;
                                                        interaction_index++) {
            fprintf(fOUTTEXT,"Interaction #%d : ",interaction_index);
            num_factors_index=num_fact_in_interaction[interaction_index]-1;
            multiples[interaction_index][num_factors_index]=1;
            for (fact_index=num_fact-1; fact_index>=0; fact_index--)
                if (interaction_equiv[interaction_index][fact_index]==1) {
                    num_factors_index--;
                    multiples[interaction_index][num_factors_index] =
                        multiples[interaction_index][num_factors_index+1]
                                                * num_level[fact_index];
                }

            for (num_factors_index=0;
                num_factors_index<num_fact_in_interaction[interaction_index];
                                                        num_factors_index++)
                fprintf(fOUTTEXT,"  multiples[%d] = %d",
                                num_factors_index,
                                multiples[interaction_index][num_factors_index]);
            fprintf(fOUTTEXT,"\n\n");
        }

        /* -------------------------------------------------------------------- */
        /* CALCULATE type_equiv2.                                               */
        /* Perhaps the most difficult part of this programming effort.  Given   */
        /* n_index, the index number of the image, and interaction_index, the   */
        /* class of interaction, map the image index number to a type within    */
        /* that interaction class.  Note the use of multiples[].                */
        for (interaction_index=1;interaction_index<=num_interactions;
                                                        interaction_index++) {
            fprintf(fOUTTEXT,"Interaction #%d:\n",interaction_index);
            for (n_index=0; n_index<total_n; n_index++) {
                type = 0;
                num_factors_index=0;
                for (fact_index=0; fact_index<num_fact; fact_index++) {
                    if (interaction_equiv[interaction_index][fact_index]==1) {
                        type += lvl_equiv[n_index][fact_index]
                                * multiples[interaction_index][num_factors_index];
                        num_factors_index++;
                    }
                }
                fprintf(fOUTTEXT,"  image %d, type %d\n", n_index+1,type+1);
                type_equiv2[interaction_index][n_index]=type;
            }
            fprintf(fOUTTEXT,"\n");
        }

            
        /* -------------------------------------------------------------------- */
        /* CALCULATE type_equiv1.                                               */
        /* Given an interaction specified by interaction_index, a lower-order   */
        /* interaction specified by interaction_index2, and a type/level within */
        /* interaction_index specified by type_index, returns the               */
        /* corresponding type/level of interaction_index2.  If there is no      */
        /* corresponding type/level, returns a zero.                            */
        /* This is used to calculate interaction class effects, after having    */
        /* calculated interaction class means.                                  */

        /* FIRST INITIALIZE type_equiv1.                                        */
        /* Set to -1.                                                           */
        for (interaction_index=1;interaction_index<=num_interactions;
                                                        interaction_index++)
            for (interaction_index2=1; interaction_index2<=num_interactions;
                                                        interaction_index2++)
                for (type_index=0; type_index<num_types[interaction_index];
                                                                type_index++)
                    type_equiv1[interaction_index][interaction_index2][type_index]=-1;

        fprintf(fOUTTEXT,"\nCheck mapping between interaction classes:\n");
        fprintf(fOUTTEXT,"Format: [interaction1,type1] --> [interaction2,type2]\n\n");
        for (interaction_index=1;interaction_index<=num_interactions;
                                                        interaction_index++) {
            for (interaction_index2=1; interaction_index2<=num_interactions;
                                                        interaction_index2++) {

                /* ------------------------------------------------------------ */
                /* SET subset_flag TO 1.                                        */
                /* If subset_flag is still 1 after the next two sections, the   */
                /* current interaction specified by interaction_index2 will be  */
                /* factored out of the interaction effect specified by          */
                /* interaction_index.  The next two sections determine whether  */
                /* or not interaction_index2 is indeed a subset of              */
                /* interaction_index.                                           */
                subset_flag = 1;

                /* --------------------------------------------------------     */
                /* Must be a subset of higher order interaction.                */
                for (fact_index=0; fact_index<num_fact; fact_index++)
                    if ((interaction_equiv[interaction_index2][fact_index]==1)
                        && (interaction_equiv[interaction_index][fact_index]==0))
                            subset_flag = 0;

                /* ------------------------------------------------------------ */
                /* But also must not be the same as the higher order            */
                /* interaction.  Check this by comparing the number of factors  */
                /* involved in the lower order interaction to the number of     */
                /* factors involved in the higher order interaction.  There     */
                /* should be fewer factors in the lower order interaction than  */
                /* in the higher order interation.                              */
                if (num_fact_in_interaction[interaction_index2]
                                >= num_fact_in_interaction[interaction_index])
                    subset_flag = 0;

                /* --------------------------------------------------------     */
                /* If the lower order interaction has met these criteria, then  */
                /* determine appopriate sign and add (in effect, subtract, if   */
                /* negative sign) to factor effect.  A difficult part of        */
                /* the program, almost identical to the problem dealt with in   */
                /* the section below entitled                                   */
                /*           SUM VOXEL VALUES INTO interaction_mean[][].....    */
                /* Given that the lower order interaction specified by          */
                /* interaction_index2 is indeed to be subtracted from/added to  */
                /* the higher order interaction specified by                    */
                /* interaction_index, subtract the appropriate types within     */
                /* interaction_index2 from the appropriate types within         */
                /* interaction_index.  In other words, given interaction_index, */
                /* interaction_index2, and type_index, where type_index         */
                /* indicates one of the types of the interaction specified by   */
                /* interaction_index, map that type to one of the types of the  */
                /* interaction specified by interaction_index2.  After this is  */
                /* determined, one can add/subtract the lower order interaction */
                /* types from the appropriate higher order interaction term     */
                /* types.                                                       */
                if (subset_flag==1) {
                    fprintf(fOUTTEXT,"\n---------------------------------\n");
                    fprintf(fOUTTEXT,"Interaction %d is a subset of interaction %d\n",
                                                interaction_index2,interaction_index);

                    /* -------------------------------------------------------- */
                    /* LOOP OVER INTERACTION TYPES AND MAP FROM HIGHER ORDER    */
                    /* TYPE TO LOWER ORDER TYPE.                                */
                    /* This is done by first mapping the higher order           */
                    /* interaction class and type to an n_index, where n_index  */
                    /* is in the range [0,total_n-1], and then mapping that     */
                    /* n_index to the corresponding lower order interaction     */
                    /* type.                                                    */
                    for (type_index=0; type_index<num_types[interaction_index];
                                                                type_index++) {

                        /* ---------------------------------------------------- */
                        /* FIRST, FIND n_index WHICH MAPS TO INTERACTION CLASS  */
                        /* interaction_index, TYPE type_index.                  */
                        n_index = 0;
                        type = -1;
                        while (type!=type_index) {
                            n_index++;
                            type = type_equiv2[interaction_index][n_index];
                        }

                        /* ---------------------------------------------------- */
                        /* NOW MAP n_index TO APPROPRIATE TYPE IN INTERACTION   */
                        /* CLASS interaction_index2.                            */
                        type_equiv1[interaction_index][interaction_index2][type_index] =
                                        type_equiv2[interaction_index2][n_index];
                        fprintf(fOUTTEXT,"[%d,%d] --> [%d,%d]\n\n",interaction_index,type_index+1,
                            interaction_index2,
                            type_equiv1[interaction_index][interaction_index2][type_index]+1);

                    } /* END LOWER ORDER INTERACTION type_index LOOP */
                } /* END subset_flag IF STATEMENT */
            } /* END interaction_index2 LOOP */
        } /* END interaction_index LOOP */

/* ---------------------------------------------------------------------------- */
/*  CALCULATE F IMAGES VOXEL-BY-VOXEL.                                          */

    /* ------------------------------------------------------------------------ */
    /* LOOP OVER VOXELS.                                                        */
    printf("ultra_anova: Now calculating F's voxel-by-voxel...\n");
    for (vox_index=0; vox_index<num_vox; vox_index++) {

        /* -------------------------------------------------------------------- */
        /* CALCULATE grandmean FOR CURRENT VOXEL.                               */
        grandmean = 0;
        for (n_index=0; n_index<total_n; n_index++)
            grandmean += (float) imgin[n_index][vox_index];
        grandmean /= ((float)total_n);



        /* -------------------------------------------------------------------- */
        /* LOOP OVER INTERACTIONS, CALCULATE INTERACTION MEANS.                 */
        /* AFTER interaction means are calculated, loop over interactions again */
        /* to calculate interaction effects.  Can't calculate effects without   */
        /* first having calculated means.                                       */
        for (interaction_index=1;interaction_index<=num_interactions;
                                                        interaction_index++) {

            /* ---------------------------------------------------------------- */
            /* INITIALIZE interaction_mean[][].                                 */
            for (type_index=0;type_index<num_types[interaction_index];
                                                                type_index++)
                interaction_mean[interaction_index][type_index]=0.;

            /* ---------------------------------------------------------------- */
            /* SUM VOXEL VALUES INTO interaction_mean[][]....                   */
            /* Use type_equiv2, calculated above.                               */
            /* Given n_index, the index number of the image, and                */
            /* interaction_index, the class of interaction, map the image index */
            /* number to a type within that interaction class.  After this is   */
            /* determined, one can add that image's voxel value to the          */
            /* appropriate interaction term type.                               */
            for (n_index=0; n_index<total_n; n_index++) {
                type = type_equiv2[interaction_index][n_index];
                interaction_mean[interaction_index][type] += imgin[n_index][vox_index];
            }


            /* ---------------------------------------------------------------- */
            /* THE FOLLOWING SECTION IS COMMENTED OUT.  IT WOULD WRITE OUT THE  */
            /* INTERACTION SUMS TO THE TEXT FILE.  USER CAN INSTEAD LOOK AT THE */
            /* INTERACTION MEANS.  MAY DECOMMENT THIS SECTION FOR DEBUGGING.    */
            /*
            if (vox_index==vox_index_int)
                for (type_index=0;type_index<num_types[interaction_index];type_index++)
                    fprintf(fOUTTEXT,"  interaction_sum %d = %f\n",
                        type_index,interaction_mean[interaction_index][type_index]);
            */

            /* ---------------------------------------------------------------- */
            /* DIVIDE interaction_mean[][] BY total_n/num_types[] TO GET        */
            /* INTERACTION MEANS.                                               */
            /* total_n/num_types[] is the number of elements within a type.     */
            for (type_index=0;type_index<num_types[interaction_index];
                                                                type_index++) {
                interaction_mean[interaction_index][type_index] 
                                        /= (total_n/num_types[interaction_index]);
                interaction_mean_store[interaction_index][type_index][vox_index]=
                        interaction_mean[interaction_index][type_index];
            }

        } /* END INTERACTION LOOP */



        /* -------------------------------------------------------------------- */
        /* FILL IN interaction_effects WITH INTERACTION MEANS NOW STORED IN     */
        /* interaction[][].                                                     */
        for (interaction_index=1;interaction_index<=num_interactions;
                                                        interaction_index++)
            for (type_index=0;type_index<num_types[interaction_index];
                                                                type_index++)
                interaction_effects[interaction_index][type_index]
                                = interaction_mean[interaction_index][type_index];

        /* -------------------------------------------------------------------- */
        /* LOOP OVER INTERACTIONS, CALCULATE INTERACTION EFFECTS.               */
        /* Now that interaction means are known, calculate effects.  This is    */
        /* done by subtracting/adding lower order interaction terms and the     */
        /* grand mean from/to the higher order interaction means.               */
        /* The grand mean will be added/subtracted after the lower order        */
        /* interaction effects.                                                 */

        for (interaction_index=1;interaction_index<=num_interactions;
                                                        interaction_index++) {

            /* ---------------------------------------------------------------- */
            /* NOW SUBTRACT LOWER ORDER INTERACTION EFFECTS WITH APPROPRIATE    */
            /* SIGN.  OPEN LOWER ORDER INTERACTION LOOP.                        */
            /* Subtract only if the lower order interaction's set of involved   */
            /* factors is a subset of the higher order interaction's set, as    */
            /* long as the two sets are not identical (in which case the two    */
            /* interactions are the same interaction!)  Here, interaction_index */
            /* will indicate the higher order interaction, and                  */
            /* interaction_index2 will indicate the lower order interaction.    */

            for (interaction_index2=1; interaction_index2<=num_interactions;
                                                        interaction_index2++) {

                /* ------------------------------------------------------------ */
                /* SET subset_flag TO 1.                                        */
                /* If subset_flag is still 1 after the next two sections, the   */
                /* current interaction specified by interaction_index2 will be  */
                /* factored out of the interaction effect specified by          */
                /* interaction_index.                                           */
                subset_flag = 1;

                /* --------------------------------------------------------     */
                /* Must be a subset of higher order interaction.                */
                for (fact_index=0; fact_index<num_fact; fact_index++)
                    if ((interaction_equiv[interaction_index2][fact_index]==1)
                        && (interaction_equiv[interaction_index][fact_index]==0))
                            subset_flag = 0;

                /* ------------------------------------------------------------ */
                /* But also must not be the same as the higher order            */
                /* interaction.  Check this by comparing the number of factors  */
                /* involved in the lower order interaction to the number of     */
                /* factors involved in the higher order interaction.  There     */
                /* should be fewer factors in the lower order interaction than  */
                /* in the higher order interation.                              */
                if (num_fact_in_interaction[interaction_index2]
                                >= num_fact_in_interaction[interaction_index])
                    subset_flag = 0;

                /* --------------------------------------------------------     */
                /* If the lower order interaction has met these criteria, then  */
                /* determine appopriate sign and add (in effect, subtract, if   */
                /* negative sign) to factor effect.  Again, a difficult part of */
                /* the program, almost identical to the problem dealt with in   */
                /* the above section entitled                                   */
                /*           SUM VOXEL VALUES INTO interaction_mean[][].....    */
                /* Given that the lower order interaction specified by          */
                /* interaction_index2 is indeed to be subtracted from/added to  */
                /* the higher order interaction specified by                    */
                /* interaction_index, subtract the appropriate types within     */
                /* interaction_index2 from the appropriate types within         */
                /* interaction_index.  In other words, given interaction_index, */
                /* interaction_index2, and type_index, where type_index         */
                /* indicates one of the types of the interaction specified by   */
                /* interaction_index2, map that type to one of the types of the */
                /* interaction specified by interaction_index.  After this is   */
                /* determined, one can add/subtract the lower order interaction */
                /* types from the appropriate higher order interaction term     */
                /* types.                                                       */
                if (subset_flag==1) {
                    /* -------------------------------------------------------- */
                    /* DETERMINE SIGN.                                          */
                    neg_one_exponent=num_fact_in_interaction[interaction_index]
                                - num_fact_in_interaction[interaction_index2];
                        
                    sign = 1.;
                    for (exp_index=1; exp_index<=neg_one_exponent; exp_index++)
                        sign *= -1.;

                    /* -------------------------------------------------------- */
                    /* LOOP OVER INTERACTION TYPES AND ADD/SUBTRACT LOWER ORDER */
                    /* TYPE MEANS TO/FROM HIGHER ORDER TYPE.                    */
                    /* This is done by simply reading entries from the array    */
                    /* type_equiv1[][][], which was filled with data above.     */
                    /* Again, type_equiv maps a type within the higher order    */
                    /* interaction to a type in the lower order interaction.    */
                    /* Addition/subtraction is determined by the sign.          */
                    for (type_index=0;type_index<num_types[interaction_index];
                                                                type_index++) {
                        type =
                            type_equiv1[interaction_index][interaction_index2][type_index];
                        interaction_effects[interaction_index][type_index] +=
                            sign * interaction_mean[interaction_index2][type];
                    }
                } /* END subset_flag IF STATEMENT */
                    
            } /* END interaction_index2 LOOP */


            /* ---------------------------------------------------------------- */
            /* NOW ADD grandmean WITH APPROPRIATE SIGN.                         */

                /* ------------------------------------------------------------ */
                /* DETERMINE SIGN.                                              */
                neg_one_exponent=num_fact_in_interaction[interaction_index];
                sign = 1.;
                for (exp_index=1; exp_index<=neg_one_exponent; exp_index++)
                    sign *= -1.;

                /* ------------------------------------------------------------ */
                /* LOOP OVER INTERACTION TYPES AND ADD/SUBTRACT grandmean       */
                for (type_index=0;type_index<num_types[interaction_index];
                                                                type_index++)
                    interaction_effects[interaction_index][type_index] +=
                        sign * grandmean;
                

        } /* END interaction_index LOOP */




        /* -------------------------------------------------------------------- */
        /* CALCULATE RESIDUAL ERRORS.                                           */
        /* Now that interaction effects are known, calculate residual errors.   */
        /* Each data point, i.e., each input image, will have associated with   */
        /* it a residual error; therefore the indexing of residual_error[] runs */
        /* from 0 to total_n-1.  Calculate residual errors by initializing to   */
        /* the original data value and subtracting from that the grand mean and */
        /* all the interaction effects.                                         */
        for (n_index=0; n_index<total_n; n_index++) {

            /* ---------------------------------------------------------------- */
            /* INITIALIZE TO ORIGINAL DATA MINUS GRAND MEAN.                    */
            residual_error[n_index]=imgin[n_index][vox_index]-grandmean;

            /* ---------------------------------------------------------------- */
            /* THEN SUBTRACT ALL THE INTERACTION EFFECTS.                       */
            /* Loop over interaction classes, determine appropriate type        */
            /* within a class, and subtract from value in residual_error[].     */

            for (interaction_index=1;interaction_index<=num_interactions;
                                                        interaction_index++) {
                type = type_equiv2[interaction_index][n_index];
                residual_error[n_index]-=interaction_effects[interaction_index][type];
            }
        } /* END n_index LOOP */


        /* -------------------------------------------------------------------- */
        /* INITIALIZE SS's AND MSE[vox_index].                                  */
        for (interaction_index=1;interaction_index<=num_interactions;
                                                        interaction_index++)
            SS[interaction_index] = 0.;

        MSE[vox_index] = 0.;

        /* -------------------------------------------------------------------- */
        /* CALCULATE SS's.                                                      */
        /* Need to multiply afterwards by total_n/num_types[] because within an */
        /* interaction type, there are more than one image.                     */
        /* Note that indexing for SS begins at 0 (zero-offset), where as the    */
        /* first index of interaction_effects begins at 1 (one-offset).         */
        for (interaction_index=1; interaction_index<=num_interactions;
                                                        interaction_index++) {
            for (type_index=0; type_index<num_types[interaction_index];
                                                                type_index++)
                SS[interaction_index] +=
                        (interaction_effects[interaction_index][type_index]
                                                *
                         interaction_effects[interaction_index][type_index]);
            SS[interaction_index] *= (total_n/num_types[interaction_index]);
        }

        /* -------------------------------------------------------------------- */
        /* CALCULATE SSE FROM residual_error[].  STORE IN MSE[vox_index] FOR    */
        /* NOW.  LATER, WILL DIVIDE BY df[0] TO GET MSE.                        */
        for (n_index=0; n_index<total_n; n_index++)
            MSE[vox_index] += residual_error[n_index]*residual_error[n_index];

        /* -------------------------------------------------------------------- */
        /* IF CURRENT VOXEL IS VOXEL OF INTEREST, OUTPUT ANALYSIS SO FAR TO     */
        /* OUTPUT TEXT FILE.                                                    */
        if (vox_index==vox_index_int) {
            fprintf(fOUTTEXT,"\n\nOriginal voxel values:\n");
            for (n_index=0; n_index<total_n; n_index++)
                fprintf(fOUTTEXT,"  Image # %d = %f\n",n_index+1,imgin[n_index][vox_index]);
            fprintf(fOUTTEXT,"\n\n  grandmean = %f\n",grandmean);
            fprintf(fOUTTEXT,"\n\nInteraction means:\n");
            for (interaction_index=1;interaction_index<=num_interactions;
                                                        interaction_index++)
                for (type_index=0;type_index<num_types[interaction_index];
                                                                type_index++)
                    fprintf(fOUTTEXT,"  interaction mean [%d][%d] = %f\n",
                        interaction_index,type_index+1,
                                interaction_mean[interaction_index][type_index]);
            fprintf(fOUTTEXT,"\n\nInteraction effects:\n");
            for (interaction_index=1;interaction_index<=num_interactions;
                                                        interaction_index++)
                for (type_index=0;type_index<num_types[interaction_index];
                                                                type_index++)
                    fprintf(fOUTTEXT,"  interaction effect [%d][%d] = %f\n",
                        interaction_index,type_index+1,
                        interaction_effects[interaction_index][type_index]);
            fprintf(fOUTTEXT,"\n\nResidual errors:\n");
            for (n_index=0; n_index<total_n; n_index++)
                fprintf(fOUTTEXT,"  res_err[%d] = %f\n",n_index+1,
                                                residual_error[n_index]);


            fprintf(fOUTTEXT,"\n\nSS's:\n");
            for (interaction_index=1; interaction_index<=num_interactions;
                                                        interaction_index++)
                fprintf(fOUTTEXT,"  SS[%d] = %f, df = %d\n",interaction_index,
                                SS[interaction_index],df[interaction_index]);
            fprintf(fOUTTEXT,"\n\n  SSE = %f, df = %d\n",MSE[vox_index],df[0]);

        }

        /* -------------------------------------------------------------------- */
        /* CALCULATE MS's.                                                      */
        for (interaction_index=1; interaction_index<=num_interactions;
                                                        interaction_index++)
            SS[interaction_index] /= df[interaction_index];


        /* -------------------------------------------------------------------- */
        /* DIVIDE SSE (STORED IN MSE[vox_index]) BY df[0] TO GET MSE.           */

        MSE[vox_index] /= df[0];


        /* -------------------------------------------------------------------- */
        /* CALCULATE F's.                                                       */
        for (interaction_index=1; interaction_index<=num_interactions;
                                                        interaction_index++)
            if (MSE[vox_index]!=0)
                F[interaction_index][vox_index]=
                        SS[interaction_index]*(float)mask[vox_index]/MSE[vox_index];
            else
                F[interaction_index][vox_index]=0.;


        /* -------------------------------------------------------------------- */
        /* IF CURRENT VOXEL IS VOXEL OF INTEREST, OUTPUT THE REST OF THE        */
        /* ANALYSIS TO OUTPUT TEXT FILE, AND CLOSE OUTPUT TEXT FILE.            */
        /* ALSO OUTPUT PROBABILITY FOR EACH F-TEST.                             */
        if (vox_index==vox_index_int) {
            fprintf(fOUTTEXT,"\n\nMS's:\n");
            for (interaction_index=1; interaction_index<=num_interactions;
                                                        interaction_index++)
                fprintf(fOUTTEXT,"  MS[%d] = %f, df = %d\n",
                        interaction_index,SS[interaction_index],df[interaction_index]);

            fprintf(fOUTTEXT,"\n\n  MSE = %f, df = %d\n",MSE[vox_index],df[0]);
            fprintf(fOUTTEXT,"\n\nF's:\n");
            printf("\n  F tests for voxel[%d,%d,%d]:\n",x_int,y_int,z_int);
            for (interaction_index=1; interaction_index<=num_interactions;
                                                        interaction_index++) {
                v1 = (float) df[interaction_index];
                v2 = (float) df[0];
                Prob = (float)betai(v2/2.,v1/2,
                        (v2/(v2+(v1*F[interaction_index][vox_index]))));
                fprintf(fOUTTEXT,"  F[%3d] = %9f, df1 = %3d, df2 = %3d, P = %9f\n",interaction_index,
                        F[interaction_index][vox_index], df[interaction_index],df[0],Prob);
                printf("    F[%3d] = %9f, df1 = %3d, df2 = %3d, P = %9f\n\n",interaction_index,
                        F[interaction_index][vox_index], df[interaction_index],df[0],Prob);
            }
            fclose(fOUTTEXT);

        }

    } /* END VOXEL LOOP */


/* ---------------------------------------------------------------------------- */
/*  WRITE OUT F IMAGES TO OUTPUT FILES.                                         */
    for (interaction_index=1; interaction_index<=num_interactions;
                                                        interaction_index++) {

        /* -------------------------------------------------------------------- */
        /* FIRST, GENERATE OUTPUT FILE NAME.                                    */
        /* CONVENTION:         F.a.b.c_v1_v2.img                                */
        /* WHERE a, b, c, etc. ARE 1's AND 0's INDICATING WHICH FACTORS ARE     */
        /* INVOLVED IN THAT INTERACTION, AND v1 AND v2 ARE THE DEGREES OF       */
        /* FREEDOM FOR THE NUMERATOR AND DENOMINATOR OF THE F-TEST.             */
        /* IN THIS EXAMPLE, c IS THE FIRST FACTOR, b IS THE SECOND FACTOR, AND  */
        /* a IS THE LAST FACTOR, BUT IF YOU HAD FOUR FACTORS, THERE WOULD THEN  */
        /* BE FOUR 1's AND 0's INSTEAD OF MERELY THREE.                         */
        sprintf(Fout,"F");
        for (fact_index=num_fact-1; fact_index>=0; fact_index--)
            sprintf(Fout,"%s%d",Fout,
                                interaction_equiv[interaction_index][fact_index]);
        
        sprintf(Fout,"%s_%d_%d.img\0",Fout,df[interaction_index],df[0]);


        /* -------------------------------------------------------------------- */
        /*  AND THEN SAVE TO DISK.                                              */

        if ((fOUT = fopen(Fout,"w")) == NULL) {
            printf ("ultra_anova: can't open %s\n",Fout);
            exit(45);
        }
        printf("ultra_anova: writing output to %s...\n",Fout);
        fwrite(F[interaction_index],sizeof(float),num_vox,fOUT);
        fclose(fOUT);
    }

/* ---------------------------------------------------------------------------- */
/*  WRITE OUT INTERACTION/LEVEL MEAN IMAGES TO OUTPUT FILES.                    */
    for (interaction_index=1; interaction_index<=num_interactions;
                                                        interaction_index++) {
        for (type_index=0;type_index<num_types[interaction_index];
                                                                type_index++) {

        /* -------------------------------------------------------------------- */
        /* FIRST, GENERATE OUTPUT FILE NAME.                                    */
        /* CONVENTION:         interaction_a.b.c_level_d.img                    */
        /* WHERE a, b, c, etc. ARE 1's AND 0's INDICATING WHICH FACTORS ARE     */
        /* INVOLVED IN THAT INTERACTION, AND d IS THE LEVEL OF THAT INTERACTION.*/
        /* IN THIS EXAMPLE, c IS THE FIRST FACTOR, b IS THE SECOND FACTOR, AND  */
        /* a IS THE LAST FACTOR, BUT IF YOU HAD FOUR FACTORS, THERE WOULD THEN  */
        /* BE FOUR 1's AND 0's INSTEAD OF MERELY THREE.                         */
        sprintf(Fout,"interaction_");
        for (fact_index=num_fact-1; fact_index>=0; fact_index--)
            sprintf(Fout,"%s%d",Fout,
                                interaction_equiv[interaction_index][fact_index]);
        
        sprintf(Fout,"%s_level_%d.img\0",Fout,type_index+1);


        /* -------------------------------------------------------------------- */
        /*  AND THEN SAVE TO DISK.                                              */

        if ((fOUT = fopen(Fout,"w")) == NULL) {
            printf ("ultra_anova: can't open %s\n",Fout);
            exit(45);
        }
        printf("ultra_anova: writing output to %s...\n",Fout);
        fwrite(interaction_mean_store[interaction_index][type_index],sizeof(float),num_vox,fOUT);
        fclose(fOUT);
        }
    }

/* ---------------------------------------------------------------------------- */
/*  WRITE OUT TOTAL SS IMAGE.                                                   */

    if ((fOUT = fopen("MSE.img","w")) == NULL) {
        printf ("ultra_anova: can't open MSE.img.\n");
        exit(46);
    }
    printf("ultra_anova: writing output to MSE.img...\n");
    fwrite(MSE,sizeof(float),num_vox,fOUT);
    fclose(fOUT);


/* ---------------------------------------------------------------------------- */
/*  WRITE OUT PROBABILITY IMAGES.                                               */

    probability = (float *) malloc(num_vox*sizeof(float));
    if (probability==NULL) {
        printf("ultra_anova: error malloc'ing for probability.\n");
        exit(47);
    }

    v2 = (float) df[0];
    for (interaction_index=1; interaction_index<=num_interactions;
                                                        interaction_index++) {
        v1 = (float) df[interaction_index];

        /* -------------------------------------------------------------------- */
        /* CONVERT F TO PROBABILITY.                                            */

        for (vox_index=0; vox_index<num_vox; vox_index++) {
            Prob = (float)betai(v2/2.,v1/2,
                        (v2/(v2+(v1*F[interaction_index][vox_index]))));
            if (Prob<=0.)
                probability[vox_index] = 0.;
            else
                probability[vox_index] = (float) (-log10((double)Prob));

        }

        /* -------------------------------------------------------------------- */
        /* INSERT THRESHOLD SQUARE INTO LOWER LEFT CORNER OF F IMAGE.           */

        for (x = 0; x < 5; x++)
            for (y = 0; y < 5; y++)
                for (z = 0; z < zdim; z++)
                    probability[(z*ydim+y)*xdim+x] =
                                        (float) (-log10((double)prob_thresh));


        /* -------------------------------------------------------------------- */
        /*  AND THEN SAVE TO DISK.                                              */

        sprintf(Fout,"pF");
        for (fact_index=num_fact-1; fact_index>=0; fact_index--)
            sprintf(Fout,"%s%d",Fout,
                                interaction_equiv[interaction_index][fact_index]);
        
        sprintf(Fout,"%s_%d_%d.img\0",Fout,df[interaction_index],df[0]);
        if ((fOUT = fopen(Fout,"w")) == NULL) {
            printf ("ultra_anova: can't open %s\n",Fout);
            exit(48);
        }
        printf("ultra_anova: writing output to %s...\n",Fout);
        fwrite(probability,sizeof(float),num_vox,fOUT);
        fclose(fOUT);
    }


/* ---------------------------------------------------------------------------- */
/*  FREE MEMORY.                                                                */

    printf("\nultra_anova: releasing memory...\n");
    free(c);
    free(img_name);
    free(num_level);
    free(num_cell_per_level);
    free(major_inc);
    free(interaction_equiv);
    free(num_types);
    free(df);
    free(num_fact_in_interaction);
    free(lvl_index2);
    free(lvl_equiv);
    free(temp);
    free(imgin);
    free(glob_mean);
    free(mask);
    free(interaction_mean);
    free(type_equiv1);
    free(type_equiv2);
    free(residual_error);
    free(MSE);
    free(SS);
    free(F);
    free(multiples);
    free(probability);
}
