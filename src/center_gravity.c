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
    int start, end;                 /* Starting and ending slices to be included */
                                    /* in average.                               */
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

    float xx,yy,zz,x_av,y_av,z_av,denom;


    /* ---------------------------------------------------------------- */
    /* CHECK INPUT ARGUMENTS.  Make sure correct number of command      */
    /* line arguments are present.                                      */

    if (argc<6) {
        printf("center_gravity: need header file, output text file, starting\n");
        printf("                slice, ending slice, and list of .img files.\n");
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

    start=atoi(argv[3]);
    end=atoi(argv[4]);
    if (start<1) {
        printf("Starting slice must be 1 or greater.\n");
        exit(2);
    }
    if (start>zdim) {
        printf("Starting slice must be less than total number of slices.\n");
        exit(3);
    }
    if (start<1) {
        printf("Ending slice must be 1 or greater.\n");
        exit(4);
    }
    if (end>zdim) {
        printf("Ending slice must be less than total number of slices.\n");
        exit(5);
    }
    if (start>end) {
        printf("Ending slice must be greater than starting slice.\n");
        exit(6);
    }

    /* **************************************************************** */
    /* LOOP OVER ALL INPUT FILES....                                    */

    if ((fp=fopen(argv[2],"w"))==NULL) {
        printf("Unable to open %s.\n",argv[2]);
        exit(2);
    }
    printf("Here are the centers of gravity (XYZ):\n\n");
    fprintf(fp,"Here are the centers of gravity (XYZ):\n\n");

    Fbuffer=float_malloc(num_vox);
    if (num_bits_per_pixel==8)
        Bbuffer=unsigned_char_malloc(num_vox);
    else if (num_bits_per_pixel==16)
        Wbuffer=unsigned_short_malloc(num_vox);
    x_av=0.;
    y_av=0.;
    z_av=0.;
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
        /* FIND CENTER OF GRAVITY, CUMULATE INTO x_cum, y_cum, z_cum    */

        xx=0.;
        yy=0.;
        zz=0.;
        denom=0.;
        for (z=start-1;z<end;z++)
            for (y=0;y<ydim;y++)
                for (x=0;x<xdim;x++) {
                    xx+=(float)x*Fbuffer[((z*ydim)+y)*xdim+x];
                    yy+=(float)y*Fbuffer[((z*ydim)+y)*xdim+x];
                    zz+=(float)z*Fbuffer[((z*ydim)+y)*xdim+x];
                    denom+=Fbuffer[((z*ydim)+y)*xdim+x];
                }
        xx/=denom;
        yy/=denom;
        zz/=denom;

        printf("%g %g %g\n",xx,yy,zz);
        fprintf(fp,"%g %g %g\n",xx,yy,zz);

        x_av+=xx;
        y_av+=yy;
        z_av+=zz;
    }
    x_av/=(float)(argc-3);
    y_av/=(float)(argc-3);
    z_av/=(float)(argc-3);

    printf("\n\nAVERAGE CENTER OF GRAVITY:\n\n");
    fprintf(fp,"\nAVERAGE CENTER OF GRAVITY:\n\n");
    printf("%g %g %g\n",x_av,y_av,z_av);
    fprintf(fp,"%g %g %g\n",x_av,y_av,z_av);

    fclose(fp);

    /* ---------------------------------------------------------------- */
    /* DEALLOCATE MEMORY.                                               */

    free(Bbuffer);
    free(Wbuffer);
    free(Fbuffer);
    return 0;
}
