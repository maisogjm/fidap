#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "analyze6.h"

struct analyze_struct readhdr(char *argv[]);

unsigned char *unsigned_char_malloc(int num_items);
int *int_malloc(int num_items);
float *float_malloc(int num_items);
unsigned short *unsigned_short_malloc(int num_items);

int read_unsigned_char_data(char *filename, unsigned char *array, int num_items);
int read_unsigned_short_int_data(char *filename, unsigned short int *array, int num_items);
int read_float_data(char *filename, float *array, int num_items);

int main(int argc, char *argv[])
{
                                    /* ***************************************** */
                                    /* FILE I/O.                                 */
    FILE *fp;                       /* General purpose file pointer.             */
    struct analyze_struct img_info; /* Header information variable.              */
    int tot_num_scans;              /* Total number of scans to be smoothed.     */

                                    /* ***************************************** */
                                    /* INPUT SCAN DATA.                          */
    int cmd_lin_arg;                /* Index for looping over command line args. */
    int num_vox;                    /* Number of voxels in image volume.         */
    int vox_index1;                 /* Index for looping over voxels.            */
    int vox_index2;                 /* Another index for looping over voxels.    */

    int num_bits_per_pixel;         /* Bits per pixel in image.                  */
    unsigned char *Bbuffer;         /* Buffer to hold 8-bit image as read in     */
                                    /* from disk.                                */
    unsigned short int *Wbuffer;    /* Buffer to hold 16-bit image as read in    */
                                    /* from disk.                                */
    float *Fbuffer;                 /* Buffer to hold floating point image as    */
                                    /* read in from disk.                        */
    float *template;                /* Template image volume.                    */
    int x,y,z;                      /* 3D coordinates.                           */
    int xdim, ydim, zdim;           /* dimensions of image data.                 */
    float mean, sum_of_squares;     /* Mean and sum of squares of image data.    */
    float var;                      /* Variance of image data.                   */
    float mag_template;             /* Magnitude of template image volume.       */
                                    /* Think of it as a vector in num_mask_vox-  */
                                    /* dimensional space.                        */
    float mag_test_image;           /* Magnitude of test image volume.           */
                                    /* Think of it as a vector in num_mask_vox-  */
                                    /* dimensional space.                        */
    float correlation;              /* Correlation between template and test.    */
    float max,min;                  /* Max and min values.                       */
    char outStr[1024];

                             /* ************************************************ */
                             /* MASK.                                            */
    unsigned char *mask;     /* Pointer to mask image.                           */
    int num_mask_vox;        /* Number of non-zero voxels in *mask.              */
    int count;               /* Number of voxels processed thus far.             */
    int *mask_voxel;         /* Voxel indices of nonzero mask voxels.            */
                             /* mask_voxel[N] is the index of the (N+1)th        */
                             /* nonzero voxel in the mask.                       */


    /* ---------------------------------------------------------------- */
    /* CHECK INPUT ARGUMENTS.  Make sure correct number of command      */
    /* line arguments are present.                                      */

    if (argc<5) {
        printf("correlate_images: need header file, name of template .img file,\n");
        printf("                  8-bit mask, output text file, then list of .img\n");
        printf("                  files to correlate with template file.\n");
        exit(1);
    }


    /* ---------------------------------------------------------------- */
    /* READ IN HEADER INFO.                                             */
    /* In particular, get the 3D dimensions from the ANALYZE header.    */
    /* readhdr is a routine to read the ANALYZE header, and returns the */
    /* C structure stored in the ANALYZE header.                        */

    img_info = readhdr(argv+1);

    xdim = img_info.dime.dim[1];
    ydim = img_info.dime.dim[2];
    zdim = img_info.dime.dim[3];
    num_bits_per_pixel=img_info.dime.bitpix;
    num_vox = xdim*ydim*zdim;

    /* ---------------------------------------------------------------- */
    /* READ IN MASK.                                                    */

    mask = unsigned_char_malloc(num_vox);
    read_unsigned_char_data(argv[3],mask,num_vox);

    num_mask_vox=0;
    for (vox_index1=0; vox_index1<num_vox; vox_index1++)
        if (mask[vox_index1]!=0)
            num_mask_vox++;

    mask_voxel=int_malloc(num_mask_vox);
    vox_index2=0;
    for (vox_index1=0; vox_index1<num_vox; vox_index1++)
        if (mask[vox_index1]!=0) {
            mask_voxel[vox_index2]=vox_index1;
            vox_index2++;
        }

    /* ---------------------------------------------------------------- */
    /* READ IN TEMPLATE IMAGE VOLUME.                                   */

    template=float_malloc(num_vox);
    if (num_bits_per_pixel==8) {
        Bbuffer=unsigned_char_malloc(num_vox);
        read_unsigned_char_data(argv[2], Bbuffer, num_vox);
        for (vox_index1=0; vox_index1<num_vox; vox_index1++)
            template[vox_index1]=(float)Bbuffer[vox_index1];
    }
    else if (num_bits_per_pixel==16) {
        Wbuffer=unsigned_short_malloc(num_vox);
        read_unsigned_short_int_data(argv[2], Wbuffer, num_vox);
        for (vox_index1=0; vox_index1<num_vox; vox_index1++)
            template[vox_index1]=(float)Wbuffer[vox_index1];
    }
    else if (num_bits_per_pixel==32) {
        read_float_data(argv[2], template, num_vox);
    }
    else {
        printf("Sorry, can't handle %d bits per pixel.\n",num_bits_per_pixel);
        exit(2);
    }

    /* ---------------------------------------------------------------- */
    /* SET MEAN OF TEMPLATE TO ZERO, THEN CALCULATE ITS MAGNITUDE.      */
    /* AC POWER = TOTAL POWER - DC POWER.                               */

    if ((fp=fopen(argv[4],"w"))==NULL) {
        printf("Unable to open %s.\n",argv[4]);
        exit(2);
    }
    mean=0.;
    for (vox_index2=0; vox_index2<num_mask_vox; vox_index2++)
        mean+=template[mask_voxel[vox_index2]];
    mean/=(float)num_mask_vox;

    sum_of_squares=0.;
    for (vox_index2=0; vox_index2<num_mask_vox; vox_index2++) {
        template[mask_voxel[vox_index2]]-=mean;
        sum_of_squares+=(template[mask_voxel[vox_index2]]*template[mask_voxel[vox_index2]]);
    }
    var=sum_of_squares/(float)num_mask_vox;
    if (var==0.) {
        printf("%s seems to be a blank zero image volume!\n",argv[2]);
        fprintf(fp,"%s seems to be a blank zero image volume!\n",argv[2]);
        mag_template=0.;
    }
    else {
        sum_of_squares=0.;
        for (vox_index2=0; vox_index2<num_mask_vox; vox_index2++) {
            template[mask_voxel[vox_index2]]/=(float)sqrt((double)var);
            sum_of_squares+=(template[mask_voxel[vox_index2]]*template[mask_voxel[vox_index2]]);
        }
        mag_template=(float)sqrt((double)sum_of_squares);
    }

    /* **************************************************************** */
    /* LOOP OVER ALL INPUT FILES....                                    */


    printf("The correlation between %s and...\n",argv[2]);
    fprintf(fp,"The correlation between %s and...\n",argv[2]);

    Fbuffer=float_malloc(num_vox);
    for (cmd_lin_arg = 5; cmd_lin_arg<=argc-1; cmd_lin_arg++) {

        /* ------------------------------------------------------------ */
        /* LOAD NEXT SCAN INTO MEMORY, THEN COMPRESS IT.                */

        if (num_bits_per_pixel==8) {
            read_unsigned_char_data(argv[cmd_lin_arg], Bbuffer, num_vox);
            for (vox_index1=0; vox_index1<num_vox; vox_index1++)
                Fbuffer[vox_index1]=(float)Bbuffer[vox_index1];
        }
        else if (num_bits_per_pixel==16) {
            read_unsigned_short_int_data(argv[cmd_lin_arg], Wbuffer, num_vox);
            for (vox_index1=0; vox_index1<num_vox; vox_index1++)
                Fbuffer[vox_index1]=(float)Wbuffer[vox_index1];
        }
        else if (num_bits_per_pixel==32) {
            read_float_data(argv[cmd_lin_arg], Fbuffer, num_vox);
        }

        /* ------------------------------------------------------------ */
        /* FIND MAX AND MIN VALUES.                                     */

        max=Fbuffer[mask_voxel[0]];
        min=Fbuffer[mask_voxel[0]];
        for (vox_index2=1; vox_index2<num_mask_vox; vox_index2++) {
            max=(Fbuffer[mask_voxel[vox_index2]] > max)
                                        ? Fbuffer[mask_voxel[vox_index2]] : max;
            min=(Fbuffer[mask_voxel[vox_index2]] < min)
                                        ? Fbuffer[mask_voxel[vox_index2]] : min;
        }

        /* ------------------------------------------------------------ */
        /* SET MEAN OF TEMPLATE TO ZERO, THEN CALCULATE ITS MAGNITUDE.  */
        /* AC POWER = TOTAL POWER - DC POWER.                           */

        mean=0.;
        for (vox_index2=0; vox_index2<num_mask_vox; vox_index2++)
            mean+=Fbuffer[mask_voxel[vox_index2]];
        mean/=(float)num_mask_vox;

        sum_of_squares=0.;
        for (vox_index2=0; vox_index2<num_mask_vox; vox_index2++) {
            Fbuffer[mask_voxel[vox_index2]]-=mean;
            sum_of_squares+=(Fbuffer[mask_voxel[vox_index2]]*Fbuffer[mask_voxel[vox_index2]]);
        }
        var=sum_of_squares/(float)num_mask_vox;

        if (var==0.) {
            printf("%s seems to be a blank zero image volume!\n",argv[cmd_lin_arg]);
            fprintf(fp,"%s seems to be a blank zero image volume!\n",argv[cmd_lin_arg]);
            mag_test_image=0.;
        }
        else {
            sum_of_squares=0.;
            for (vox_index2=0; vox_index2<num_mask_vox; vox_index2++) {
                Fbuffer[mask_voxel[vox_index2]]/=(float)sqrt((double)var);
                sum_of_squares+=(Fbuffer[mask_voxel[vox_index2]]*Fbuffer[mask_voxel[vox_index2]]);
            }
            mag_test_image=(float)sqrt((double)sum_of_squares);
        }

        /* ------------------------------------------------------------ */
        /* CALCULATE CORRELATION BETWEEN TEMPLATE AND TEST IMAGE.       */

        correlation=0.;
        for (vox_index2=0; vox_index2<num_mask_vox; vox_index2++)
            correlation+=
                Fbuffer[mask_voxel[vox_index2]]*template[mask_voxel[vox_index2]];

        if (mag_template!=0.)   correlation/=mag_template;
        if (mag_test_image!=0.) correlation/=mag_test_image;

        sprintf(outStr,"  %s is %f , max = %g , min = %g , mean = %g , std = %g\n",argv[cmd_lin_arg],correlation,max,min,mean,(float)sqrt((double)var));
        printf(outStr);
        fprintf(fp,outStr);
#if 0
        printf("  %s is %f (max = %g, min = %g, mean = %g, std = %g) \n",argv[cmd_lin_arg],correlation,max,min,mean,(float)sqrt((double)var));
        fprintf(fp,"  %s is %f (max = %g, min = %g, mean = %g, std = %g) \n",argv[cmd_lin_arg],correlation,max,min,mean,(float)sqrt((double)var));
#endif
    }
    fclose(fp);

    /* ---------------------------------------------------------------- */
    /* DEALLOCATE MEMORY.                                               */

    free(Bbuffer);
    free(Wbuffer);
    free(Fbuffer);
    free(template);
    free(mask);
    return 0;
}
