#include "analyze6.h"
#include <stdlib.h>
#include <stdio.h>

/* ************************************************************************** */
struct analyze_struct readhdr(char *argv[])

/*        Returns header of type struct analyze_struct, given file            */
/*        *argv[] from which to read.                                         */
{
        struct analyze_struct img_info;
        FILE *fp;

/* -------------------------------------------------------------------- */
/*        Open header file *argv[]                                      */

        if ((fp = fopen(*argv,"r")) == NULL) {
            printf("readhdr: can't open %s\n",*argv);
            img_info.hk.session_error = 1;
            exit(1);
        }

/* -------------------------------------------------------------------- */
/*        Read header from file.                                        */

         if(fread(&img_info, sizeof(struct analyze_struct), 1, fp) != 1){
            printf("readhdr: error reading %s\n",*argv);
            img_info.hk.session_error = 2;
            exit(2);
        }
        fclose(fp);

        img_info.hk.session_error = 0;
        return img_info;
}

/* ************************************************************************** */
int *readimg(char *argv[], int xdim, int ydim, int zdim)

/*        Returns pointer to image read into memory from *argv[], having      */
/*        dimensions xdim, ydim, and zdim.  Returns a NULL pointer if         */
/*        data file cannot be opened.                                         */
{
        FILE *fp;
        int *mat;
        int x,y,z;

/* -------------------------------------------------------------------- */
/*        Allocate memory for pointer.                                  */
/*        Integer data take up FOUR bytes each.                         */

        mat   = (int *)  malloc(4*xdim*ydim*zdim);

/* -------------------------------------------------------------------- */
/*        Read in image.  Note that data is read into an integer array  */
/*        pointed to by *mat, and that data is read in as UNSIGNED CHAR */
/*        before being converted to data type INT.  This gives byte data*/
/*        a range of 0 to 255.                                          */

        if ((fp = fopen(*argv,"r")) == NULL) {
            printf("readimg: can't open %s\n",*argv);
            return NULL;
        }

        printf("readimg: reading in image from %s...\n",*argv);
        for (z = 0; z <= zdim-1; z++)
            for (y = 0; y <= ydim-1; y++)
                for (x = 0; x <= xdim-1; x++)
                    *(mat + x*ydim*zdim + y*zdim + z) =
                                                    (int) (unsigned char) fgetc(fp);
        fclose(fp);

        return mat;
}
 
/* ************************************************************************** */
/*int writimg(char *argv[], int *mat, int xdim, int ydim, int zdim) */

/*        Writes image to disk given filename in *argv[], data to             */
/*        write starting at address pointed to by mat, and dimensions         */
/*        xdim, ydim, and zdim.                                              
{
        FILE *fp;
        int x,y,z;

        if ((fp=fopen(*argv,"w"))==0) {
            printf ("writeimg: unable to create: %s\n",*argv);
            return 1;
        }

        printf("writimg: writing image to %s...\n",*argv);
        for (z = 0; z <= zdim-1; z++)
            for (y = 0; y <= ydim-1; y++)
                for (x = 0; x <= xdim-1; x++)
                     fputc((char) *(mat + x*ydim*zdim + y*zdim + z), fp);
        fclose(fp);
        return 0;
}
 */
/* ************************************************************************** */
double *readDimg(char *argv[], int xdim, int ydim, int zdim)

/*      Returns pointer to double image read into memory from *argv[],        */
/*        having dimensions xdim, ydim, and zdim.  Returns a NULL pointer     */
/*        if data file cannot be opened.                                      */
{
        FILE *fp;
        double *mat;
        int x,y,z;

/* -------------------------------------------------------------------- */
/*      Allocate memory for pointer.                                    */
/*      Double data take up EIGHT bytes each.                           */

        mat   = (double *)  malloc(8*xdim*ydim*zdim);

/* -------------------------------------------------------------------- */
/*      Read in image.  Note that data is read into a DOUBLE array      */
/*      pointed to by *mat, and that data is read in as UNSIGNED CHAR   */
/*      before being converted to data type DOUBLE.                     */

        if ((fp = fopen(*argv,"r")) == NULL) {
            printf("readimg: can't open %s\n",*argv);
            return NULL;
        }

        printf("readDimg: reading in image from %s...\n",*argv);
        for (z = 0; z <= zdim-1; z++)
            for (y = 0; y <= ydim-1; y++)
                for (x = 0; x <= xdim-1; x++)
                    *(mat + x*ydim*zdim + y*zdim + z) =
                                        (double) (unsigned char) fgetc(fp);
        fclose(fp);

        return mat;
}

/* ************************************************************************** */
int writimg2(char *argv[], int *mat, int numpix)

/*        Writes image to disk given filename in *argv[], data to             */
/*        write starting at address pointed to by mat, and number of          */
/*        pixels numpix.                                                      */
{
        FILE *fp;
        int x;

        if ((fp=fopen(*argv,"w"))==0) {
            printf ("writeimg2: unable to create: %s\n",*argv);
            return 1;
        }

        printf("writimg2: writing image to %s...\n",*argv);
        for (x = 0; x <= numpix-1; x++) {
        printf("%d ",*(mat+x));
             fputc((char) *(mat + x), fp);
            }
        fclose(fp);
        return 0;
}

/* ************************************************************************** */
int unsigned_char_fwrite(char *filename, unsigned char *array, int num_items)

/*        Writes unsigned char data to disk given filename in *argv[], data to  */
/*        write starting at address pointed to by array, and number of voxels   */
/*        num_items.                                                            */
{
        FILE *fp;
        int x,num_wrote_out;

        if ((fp=fopen(filename,"w"))==NULL) {
            printf ("unsigned_char_fwrite: unable to create: %s\n",filename);
            return 1;
        }
        printf("unsigned_char_fwrite: writing unsigned char data to %s...\n",filename);
        num_wrote_out=fwrite(array,sizeof(unsigned char),num_items,fp);
        if (num_wrote_out!=num_items) {
            printf("unsigned_char_fwrite: %d voxels written out.  Should have been %d.\n",
                                                        num_wrote_out,num_items);
            return 2;
        }
        fclose(fp);
        return 0;
}

/* ************************************************************************** */
int float_fwrite(char *filename, float *array, int num_items)

/*        Writes floating point data to disk given filename in *argv[], data to */
/*        write starting at address pointed to by array, and number of voxels   */
/*        num_items.                                                            */
{
        FILE *fp;
        int x,num_wrote_out;

        if ((fp=fopen(filename,"w"))==NULL) {
            printf ("float_fwrite: unable to create: %s\n",filename);
            return 1;
        }
        printf("float_fwrite: writing floating point data to %s...\n",filename);
        num_wrote_out=fwrite(array,sizeof(float),num_items,fp);
        if (num_wrote_out!=num_items) {
            printf("float_fwrite: %d voxels written out.  Should have been %d.\n",
                                                        num_wrote_out,num_items);
            return 2;
        }
        fclose(fp);
        return 0;
}

/* ************************************************************************** */
int *readimgquiet(char *argv[], int xdim, int ydim, int zdim)

/*        Returns pointer to image read into memory from *argv[], having      */
/*        dimensions xdim, ydim, and zdim.  Returns a NULL pointer if         */
/*        data file cannot be opened.                                        */
{
        FILE *fp;
        int *mat;
        int x,y,z;

/* -------------------------------------------------------------------- */
/*        Allocate memory for pointer.                                  */
/*        Integer data take up FOUR bytes each.                         */

        mat   = (int *)  malloc(4*xdim*ydim*zdim);

/* -------------------------------------------------------------------- */
/*        Read in image.  Note that data is read into an integer array  */
/*        pointed to by *mat, and that data is read in as UNSIGNED CHAR */
/*        before being converted to data type INT.  This gives byte data*/
/*        a range of 0 to 255.                                          */

        if ((fp = fopen(*argv,"r")) == NULL) {
            printf("readimg: can't open %s\n",*argv);
            return NULL;
        }

        for (z = 0; z <= zdim-1; z++)
            for (y = 0; y <= ydim-1; y++)
                for (x = 0; x <= xdim-1; x++)
                    *(mat + x*ydim*zdim + y*zdim + z) =
                                                    (int) (unsigned char) fgetc(fp);
        fclose(fp);

        return mat;
}
 
/* ************************************************************************** */
int *int_malloc(int num_items)
{
    int *array;
    array=(int *)malloc(num_items*sizeof(int));
    if (array==NULL) {
        printf("int_malloc: unable to malloc %d bytes.\n",num_items);
        exit(1);
    }
    return array;
}

/* ************************************************************************** */
char *char_malloc(int num_items)
{
    char *array;
    array=(char *)malloc(num_items*sizeof(char));
    if (array==NULL) {
        printf("char_malloc: unable to malloc %d bytes.\n",num_items);
        exit(1);
    }
    return array;
}

/* ************************************************************************** */
char **char_malloc2(int num_items1, int num_items2)
{
    char **array;
    int index;

    array = (char **) malloc(num_items1*sizeof(char *));
    if (array==NULL) {
        printf("char_malloc2: error malloc'ing for array.\n");
        exit(12);
    }
    for (index=0; index<num_items1; index++)
        array[index]=char_malloc(num_items2);

    return array;
}

/* ************************************************************************** */
char ***char_malloc3(int num_items1, int num_items2, int num_items3)
{
    char ***array;
    int index,count;;

    array = (char ***) malloc(num_items1*sizeof(char **));
    if (array==NULL) {
        printf("char_malloc3: error malloc'ing for array.\n");
        exit(12);
    }
/*    printf("       "); */
    count=0;
    for (index=0; index<num_items1; index++) {
        count++;
/*        printf("\b\b\b\b\b\b\b%5.1f %%",(float)count*100./(float)num_items1); */
        array[index]=char_malloc2(num_items2,num_items3);
    }
/*    printf("\n"); */

    return array;
}

/* ************************************************************************** */
unsigned char *unsigned_char_malloc(int num_items)
{
    unsigned char *array;
    array=(unsigned char *)malloc(num_items*sizeof(unsigned char));
    if (array==NULL) {
        printf("unsigned_char_malloc: unable to malloc %d bytes.\n",num_items);
        exit(1);
    }
    return array;
}

/* ************************************************************************** */
float *float_malloc(int num_items)
{
    float *array;
    array=(float *)malloc(num_items*sizeof(float));
    if (array==NULL) {
        printf("float_malloc: unable to malloc %d bytes.\n",num_items);
        exit(1);
    }
    return array;
}

/* ************************************************************************** */
float **float_malloc2(int num_items1, int num_items2)
{
    float **array;
    int index;

    array = (float **) malloc(num_items1*sizeof(float *));
    if (array==NULL) {
        printf("float_malloc2: error malloc'ing for array.\n");
        exit(12);
    }
    for (index=0; index<num_items1; index++)
        array[index]=float_malloc(num_items2);

    return array;
}

/* ************************************************************************** */
float ***float_malloc3(int num_items1, int num_items2, int num_items3)
{
    float ***array;
    int index,count;

    array = (float ***) malloc(num_items1*sizeof(float **));
    if (array==NULL) {
        printf("float_malloc3: error malloc'ing for array.\n");
        exit(12);
    }
/*    printf("       "); */
    count=0;
    for (index=0; index<num_items1; index++) {
        count++;
/*        printf("\b\b\b\b\b\b\b%5.1f %%",(float)count*100./(float)num_items1); */
        array[index]=float_malloc2(num_items2,num_items3);
    }
/*    printf("\n"); */

    return array;
}

/* ************************************************************************** */
double *double_malloc(int num_items)
{
    double *array;
    array=(double *)malloc(num_items*sizeof(double));
    if (array==NULL) {
        printf("double_malloc: unable to malloc %d bytes.\n",num_items);
        exit(1);
    }
    return array;
}

/* ************************************************************************** */
double **double_malloc2(int num_items1, int num_items2)
{
    double **array;
    int index;

    array = (double **) malloc(num_items1*sizeof(double *));
    if (array==NULL) {
        printf("double_malloc2: error malloc'ing for array.\n");
        exit(12);
    }
    for (index=0; index<num_items1; index++)
        array[index]=double_malloc(num_items2);

    return array;
}

/* ************************************************************************** */
double ***double_malloc3(int num_items1, int num_items2, int num_items3)
{
    double ***array;
    int index,count;

    array = (double ***) malloc(num_items1*sizeof(double **));
    if (array==NULL) {
        printf("double_malloc3: error malloc'ing for array.\n");
        exit(12);
    }
/*    printf("       "); */
    count=0;
    for (index=0; index<num_items1; index++) {
        count++;
/*        printf("\b\b\b\b\b\b\b%5.1f %%",(double)count*100./(double)num_items1); */
        array[index]=double_malloc2(num_items2,num_items3);
    }
/*    printf("\n"); */

    return array;
}

/* ************************************************************************** */
unsigned short *unsigned_short_malloc(int num_items)
{
    unsigned short *array;
    array=(unsigned short *)malloc(num_items*sizeof(unsigned short));
    if (array==NULL) {
        printf("unsigned_short_malloc: unable to malloc %d bytes.\n",num_items);
        exit(1);
    }
    return array;
}

/* ************************************************************************** */
unsigned short **unsigned_short_malloc2(int num_items1, int num_items2)
{
    unsigned short **array;
    int index;

    array = (unsigned short **) malloc(num_items1*sizeof(unsigned short *));
    if (array==NULL) {
        printf("unsigned_short_malloc2: error malloc'ing for array.\n");
        exit(12);
    }
    for (index=0; index<num_items1; index++)
        array[index]=unsigned_short_malloc(num_items2);

    return array;
}

/* ************************************************************************** */
unsigned short ***unsigned_short_malloc3(int num_items1, int num_items2, int num_items3)
{
    unsigned short ***array;
    int index,count;

    array = (unsigned short ***) malloc(num_items1*sizeof(unsigned short **));
    if (array==NULL) {
        printf("unsigned_short_malloc3: error malloc'ing for array.\n");
        exit(12);
    }
/*    printf("       "); */
    count=0;
    for (index=0; index<num_items1; index++) {
        count++;
/*        printf("\b\b\b\b\b\b\b%5.1f %%",(float)count*100./(float)num_items1); */
        array[index]=unsigned_short_malloc2(num_items2,num_items3);
    }
/*    printf("\n"); */

    return array;
}

/* ************************************************************************** */
int unsigned_short_fwrite(char *filename, unsigned short *array, int num_items)

/*        Writes unsigned short data to disk given filename in *argv[], data to */
/*        write starting at address pointed to by array, and number of voxels   */
/*        num_items.                                                            */
{
        FILE *fp;
        int x,num_wrote_out;

        if ((fp=fopen(filename,"w"))==NULL) {
            printf ("unsigned_short_fwrite: unable to create: %s\n",filename);
            return 1;
        }
        printf("unsigned_short_fwrite: writing unsigned short data to %s...\n",filename);
        num_wrote_out=fwrite(array,sizeof(unsigned short),num_items,fp);
        if (num_wrote_out!=num_items) {
            printf("unsigned_short_fwrite: %d voxels written out.  Should have been %d.\n",
                                                        num_wrote_out,num_items);
            return 2;
        }
        fclose(fp);
        return 0;
}


/* ************************************************************************** */
signed short *signed_short_malloc(int num_items)
{
    signed short *array;
    array=(signed short *)malloc(num_items*sizeof(signed short));
    if (array==NULL) {
        printf("signed_short_malloc: unable to malloc %d bytes.\n",num_items);
        exit(1);
    }
    return array;
}

/* ************************************************************************** */
signed short **signed_short_malloc2(int num_items1, int num_items2)
{
    signed short **array;
    int index;

    array = (signed short **) malloc(num_items1*sizeof(signed short *));
    if (array==NULL) {
        printf("signed_short_malloc2: error malloc'ing for array.\n");
        exit(12);
    }
    for (index=0; index<num_items1; index++)
        array[index]=signed_short_malloc(num_items2);

    return array;
}

/* ************************************************************************** */
signed short ***signed_short_malloc3(int num_items1, int num_items2, int num_items3)
{
    signed short ***array;
    int index,count;

    array = (signed short ***) malloc(num_items1*sizeof(signed short **));
    if (array==NULL) {
        printf("signed_short_malloc3: error malloc'ing for array.\n");
        exit(12);
    }
/*    printf("       "); */
    count=0;
    for (index=0; index<num_items1; index++) {
        count++;
/*        printf("\b\b\b\b\b\b\b%5.1f %%",(float)count*100./(float)num_items1); */
        array[index]=signed_short_malloc2(num_items2,num_items3);
    }
/*    printf("\n"); */

    return array;
}

/* ************************************************************************** */
int signed_short_fwrite(char *filename, signed short *array, int num_items)

/*        Writes signed short data to disk given filename in *argv[], data to */
/*        write starting at address pointed to by array, and number of voxels   */
/*        num_items.                                                            */
{
        FILE *fp;
        int x,num_wrote_out;

        if ((fp=fopen(filename,"w"))==NULL) {
            printf ("signed_short_fwrite: unable to create: %s\n",filename);
            return 1;
        }
        printf("signed_short_fwrite: writing signed short data to %s...\n",filename);
        num_wrote_out=fwrite(array,sizeof(signed short),num_items,fp);
        if (num_wrote_out!=num_items) {
            printf("signed_short_fwrite: %d voxels written out.  Should have been %d.\n",
                                                        num_wrote_out,num_items);
            return 2;
        }
        fclose(fp);
        return 0;
}


/* ************************************************************************** */
int write_double_matrix2textfile(char *filename, double **array, int num_rows,
                                                                int num_cols)
{
        FILE *fp;
        int row_index,col_index,num_wrote_out;

        if ((fp=fopen(filename,"w"))==NULL) {
            printf ("write_double_matrix2textfile: unable to create: %s\n",
                                                                filename);
            return 1;
        }
/*        printf("write_double_matrix2textfile: writing floating point data to %s...\n",filename); */
        for(row_index=0;row_index<num_rows;row_index++) {
            for(col_index=0;col_index<num_cols;col_index++)
                fprintf(fp,"%f ",array[row_index][col_index]);
            fprintf(fp,"\n");
        }

        fclose(fp);
        return 0;
}

/* ************************************************************************** */
int write_float_matrix2textfile(char *filename, float **array, int num_rows,
                                                                int num_cols)
{
        FILE *fp;
        int row_index,col_index,num_wrote_out;

        if ((fp=fopen(filename,"w"))==NULL) {
            printf ("write_float_matrix2textfile: unable to create: %s\n",
                                                                filename);
            return 1;
        }
/*        printf("write_float_matrix2textfile: writing floating point data to %s...\n",filename); */
        for(row_index=0;row_index<num_rows;row_index++) {
            for(col_index=0;col_index<num_cols;col_index++)
                fprintf(fp,"%f ",array[row_index][col_index]);
            fprintf(fp,"\n");
        }

        fclose(fp);
        return 0;
}

/* ************************************************************************** */
int write_NRCdouble_matrix2textfile(char *filename, double **array, int num_rows,
                                                                int num_cols)
{
        FILE *fp;
        int row_index,col_index,num_wrote_out;

        if ((fp=fopen(filename,"w"))==NULL) {
            printf ("write_NRCdouble_matrix2textfile: unable to create: %s\n",
                                                                filename);
            return 1;
        }
/*        printf("write_NRCdouble_matrix2textfile: writing floating point data to %s...\n",filename); */
        for(row_index=1;row_index<=num_rows;row_index++) {
            for(col_index=1;col_index<=num_cols;col_index++)
                fprintf(fp,"%f ",array[row_index][col_index]);
            fprintf(fp,"\n");
        }

        fclose(fp);
        return 0;
}

/* ************************************************************************** */
int read_unsigned_char_data(char *filename, unsigned char *array, int num_items)
{
        FILE *fp;
        int num_objs_read_in;

        if ((fp = fopen(filename,"r")) == NULL) {
            printf("read_unsigned_char_data: can't open %s\n",filename);
            exit(1);
        }

        num_objs_read_in=fread(array,sizeof(unsigned char),num_items,fp);
        fclose(fp);

        if(num_objs_read_in!=num_items) {
            printf("read_unsigned_char_data: %d items read in; should have been %d.\n",
                                                num_objs_read_in,num_items);
            exit(2);
        }
        return 0;
}
 
/* ************************************************************************** */
int read_unsigned_short_int_data(char *filename, unsigned short int *array, int num_items)
{
        FILE *fp;
        int num_objs_read_in;

        if ((fp = fopen(filename,"r")) == NULL) {
            printf("read_unsigned_short_int_data: can't open %s\n",filename);
            exit(1);
        }

        num_objs_read_in=fread(array,sizeof(unsigned short int),num_items,fp);
        fclose(fp);

        if(num_objs_read_in!=num_items) {
            printf("read_unsigned_short_int_data: %d items read in; should have been %d.\n",
                                                num_objs_read_in,num_items);
            exit(2);
        }
        return 0;
}
 
/* ************************************************************************** */
int read_signed_short_int_data(char *filename, signed short int *array, int num_items)
{
        FILE *fp;
        int num_objs_read_in;

        if ((fp = fopen(filename,"r")) == NULL) {
            printf("read_signed_short_int_data: can't open %s\n",filename);
            exit(1);
        }

        num_objs_read_in=fread(array,sizeof(signed short int),num_items,fp);
        fclose(fp);

        if(num_objs_read_in!=num_items) {
            printf("read_signed_short_int_data: %d items read in; should have been %d.\n",
                                                num_objs_read_in,num_items);
            exit(2);
        }
        return 0;
}
 
/* ************************************************************************** */
int read_float_data(char *filename, float *array, int num_items)
{
        FILE *fp;
        int num_objs_read_in;

        if ((fp = fopen(filename,"r")) == NULL) {
            printf("read_float_data: can't open %s\n",filename);
            exit(1);
        }

        num_objs_read_in=fread(array,sizeof(float),num_items,fp);
        fclose(fp);

        if(num_objs_read_in!=num_items) {
            printf("read_float_data: %d items read in; should have been %d.\n",
                                                num_objs_read_in,num_items);
            exit(2);
        }
        return 0;
}
 
/* ************************************************************************** */
void mult_NRCdouble_matrices(double **input1, double **input2, int num_rows1,
        int num_cols1, int num_cols2, double **output)
/* Multiplies two Numerical Recipes in C matrices.                            */
{
        int index1,index2,index3;
        for(index1=1;index1<=num_rows1;index1++)
            for(index2=1;index2<=num_cols2;index2++) {
                output[index1][index2]=0.;
                for(index3=1;index3<=num_cols1;index3++)
                    output[index1][index2]
                        +=input1[index1][index3]*input2[index3][index2];
            }

}
 
/* ************************************************************************** */
void mult_NRCdouble_matrices_tr2(double **input1, double **input2, int num_rows1,
        int num_cols1, int num_rows2, double **output)
/* Same as mult_NRCdouble_matrices, only the second input matrix is           */
/* transposed before matrix multiplication.                                   */
{
        int index1,index2,index3;
        for(index1=1;index1<=num_rows1;index1++)
            for(index2=1;index2<=num_rows2;index2++) {
                output[index1][index2]=0.;
                for(index3=1;index3<=num_cols1;index3++)
                    output[index1][index2]
                        +=input1[index1][index3]*input2[index2][index3];
            }

}
 
/* ************************************************************************** */
void mult_NRCdouble_matrices_tr1(double **input1, double **input2, int num_rows1,
        int num_cols1, int num_cols2, double **output)
/* Same as mult_NRCdouble_matrices, only the first input matrix is            */
/* transposed before matrix multiplication.                                   */
{
        int index1,index2,index3;
        for(index1=1;index1<=num_cols1;index1++)
            for(index2=1;index2<=num_cols2;index2++) {
                output[index1][index2]=0.;
                for(index3=1;index3<=num_rows1;index3++)
                    output[index1][index2]
                        +=input1[index3][index1]*input2[index3][index2];
            }

}
 
/* ************************************************************************** */
void mult_NRCdouble_matrices_tr_both(double **input1, double **input2, int num_rows1,
        int num_cols1, int num_rows2, double **output)
/* Same as mult_NRCdouble_matrices, only BOTH input matrices are              */
/* transposed before matrix multiplication.                                   */
{
        int index1,index2,index3;
        for(index1=1;index1<=num_cols1;index1++)
            for(index2=1;index2<=num_rows2;index2++) {
                output[index1][index2]=0.;
                for(index3=1;index3<=num_rows1;index3++)
                    output[index1][index2]
                        +=input1[index3][index1]*input2[index2][index3];
            }

}
 
