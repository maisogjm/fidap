#include <stdio.h>
#include <stdlib.h>
#include "analyze6.h"
unsigned char *unsigned_char_malloc(int num_items);
int read_unsigned_char_data(char *filename, unsigned char *array, int num_items);
int unsigned char_fwrite(char *filename, unsigned char *array, int num_items);

main (int argc, char *argv[])
{
    /* 09.08.93 : JM Maisog
        Converts a single-precision unsigned charing point flat file to
        unsigned 8-bit flat file.  RESCALES VALUES TO RANGE OF
        0 TO 255.
       06.28.96 : JMM
        Zooms a unsigned charing point image file by 2x in X and Y.
       07.19.96 : JMM
        Zooms an 8-bit image file by 2x in X and Y.
                                                                        */

    FILE *fHDR, *fIN, *fOUT;
    struct analyze_struct img_info;
    unsigned char *imgin;   /* Input voxel data.                        */
    unsigned char *imgout;  /* Output voxel data.                       */
    int input_num_vox;      /* Number of voxels in input image.         */
    int output_num_vox;     /* Number of voxels in output image.        */
    int vox_index;
    int read_status;
    int write_status;
    int xdim,ydim,zdim;     /* Dimensions of old image.                 */
    int Nxdim,Nydim,Nzdim;  /* Dimensions of new image.                 */
    int x,y,z,xx,yy,zz;     /* Indices into arrays.                     */


    /* **************************************************************** */
    /* VERIFY NUMBER OF COMMAND LINE ARGUMENTS.                         */
    if (argc != 4) {
        printf("zoomX2_8bit: need input .hdr file, input 8-bit .img file,\n");
        printf("             and output 8-bit .img file.\n");
        printf("             Output image will have dimensions:\n");
        printf("               New xdim = 2*(Old xdim)\n");
        printf("               New ydim = 2*(Old ydim)\n");
        printf("               New zdim = Old zdim\n");
        exit(1);
    }


    /* **************************************************************** */
    /* READ IN xdim, ydim & zdim FROM ANALYZE HEADER.                   */

    if ((fHDR = fopen(argv[1],"r")) == NULL) {
        printf("zoomX2_8bit: can't open %s\n",argv[1]);
        exit(1);
    }
    if(fread(&img_info, sizeof(struct analyze_struct), 1, fHDR) != 1){
        printf("zoomX2_8bit: error reading %s\n",argv[1]);
        img_info.hk.session_error = 2;
        exit(2);
    }
    fclose(fHDR);
    xdim = img_info.dime.dim[1];
    ydim = img_info.dime.dim[2];
    zdim = img_info.dime.dim[3];
    input_num_vox = xdim*ydim*zdim;
    printf("zoomX2_8bit: Input image dimensions are %d x %d x %d.\n",
                                                        xdim,ydim,zdim);
    printf("zoomX2_8bit: Number of voxels in input image = %d\n",input_num_vox);

    Nxdim = xdim*2;
    Nydim = ydim*2;
    Nzdim = zdim;
    output_num_vox = Nxdim*Nydim*Nzdim;
    printf("zoomX2_8bit: Output image dimensions are %d x %d x %d.\n",
                                                        Nxdim,Nydim,Nzdim);
    printf("zoomX2_8bit: Number of voxels in output image = %d\n",output_num_vox);

    /* **************************************************************** */
    /* ASSIGN MEMORY, READ IN FP DATA.                                  */
    imgin  = unsigned_char_malloc(input_num_vox);
    imgout = unsigned_char_malloc(output_num_vox);
    printf("Reading in data...\n");
    read_unsigned_char_data(argv[2],imgin,input_num_vox);
    printf("Now interpolating...\n");


    /* **************************************************************** */
    /* ZOOM X2.                                                         */

    for (z=0; z<zdim; z++)
        for (y=0; y<ydim; y++)
            for (x=0; x<xdim; x++) {
                imgout[((z*Nydim)+2*y)*(Nxdim)+x*2]=
                                            imgin[(z*ydim+y)*xdim+x];
                imgout[((z*Nydim)+2*y)*(Nxdim)+(x*2)+1]=
                                            imgin[(z*ydim+y)*xdim+x];

                imgout[((z*Nydim)+(2*y)+1)*(Nxdim)+x*2]=
                                            imgin[(z*ydim+y)*xdim+x];
                imgout[((z*Nydim)+(2*y)+1)*(Nxdim)+(x*2)+1]=
                                            imgin[(z*ydim+y)*xdim+x];
    }


    /* **************************************************************** */
    /* WRITE OUT ZOOMED IMAGE TO OUTPUT FILE                            */
    unsigned_char_fwrite(argv[3],imgout,output_num_vox);
}
