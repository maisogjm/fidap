#include <stdio.h>
#include <stdlib.h>
#include "analyze6.h"
float *float_malloc(int num_items);
int read_float_data(char *filename, float *array, int num_items);
int float_fwrite(char *filename, float *array, int num_items);

main (int argc, char *argv[])
{
    /* 09.08.93 : JM Maisog
        Converts a single-precision floating point flat file to
        unsigned 8-bit flat file.  RESCALES VALUES TO RANGE OF
        0 TO 255.
       06.28.96 : JMM
        Zooms a floating point image file by 2x in X and Y.
                                                                        */

    FILE *fHDR, *fIN, *fOUT;
    struct analyze_struct img_info;
    float *imgin;           /* Input voxel data.                        */
    float *imgout;          /* Output voxel data.                       */
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
        printf("zoomX2fpB: need input .hdr file, input floating point\n");
        printf("           .img file, and output floating point .img file.\n");
        printf("           Output image will have dimensions:\n");
        printf("             New xdim = 2*(Old xdim)\n");
        printf("             New ydim = 2*(Old ydim)\n");
        printf("             New zdim = Old zdim\n\n");
	printf("           Linear interpolation in plane will be done.\n");
        exit(1);
    }


    /* **************************************************************** */
    /* READ IN xdim, ydim & zdim FROM ANALYZE HEADER.                   */

    if ((fHDR = fopen(argv[1],"r")) == NULL) {
        printf("zoomX2fpB: can't open %s\n",argv[1]);
        exit(1);
    }
    if(fread(&img_info, sizeof(struct analyze_struct), 1, fHDR) != 1){
        printf("zoomX2fpB: error reading %s\n",argv[1]);
        img_info.hk.session_error = 2;
        exit(2);
    }
    fclose(fHDR);
    xdim = img_info.dime.dim[1];
    ydim = img_info.dime.dim[2];
    zdim = img_info.dime.dim[3];
    input_num_vox = xdim*ydim*zdim;
    printf("zoomX2fpB: Input image dimensions are %d x %d x %d.\n",
                                                        xdim,ydim,zdim);
    printf("zoomX2fpB: Number of voxels in input image = %d\n",input_num_vox);

    Nxdim = xdim*2;
    Nydim = ydim*2;
    Nzdim = zdim;
    output_num_vox = Nxdim*Nydim*Nzdim;
    printf("zoomX2fpB: Output image dimensions are %d x %d x %d.\n",
                                                        Nxdim,Nydim,Nzdim);
    printf("zoomX2fpB: Number of voxels in output image = %d\n",output_num_vox);

    /* **************************************************************** */
    /* ASSIGN MEMORY, READ IN FP DATA.                                  */
    imgin  = float_malloc(input_num_vox);
    imgout = float_malloc(output_num_vox);
    printf("Reading in data...\n");
    read_float_data(argv[2],imgin,input_num_vox);
    printf("Now interpolating...\n");


    /* **************************************************************** */
    /* ZOOM X2.                                                         */

        /* ************************************************************ */
        /* FIRST REPLICATE ALL ORIGINAL VOXELS.                         */
        for (z=0; z<zdim; z++)
            for (y=0; y<ydim; y++)
                for (x=0; x<xdim; x++)
                    imgout[((z*Nydim)+2*y)*(Nxdim)+x*2]=
                                                imgin[(z*ydim+y)*xdim+x];

        /* ************************************************************ */
        /* THEN, INTERPOLATE IN-PLANE VOXELS...                         */
        printf("zoomX2fpB: interpolating in-plane voxels...\n");
        for (z=0; z<Nzdim; z++) {

            /* ******************************************************** */
            /* ALONG ROWS...                                            */

            for (y=0; y<Nydim; y+=2) {
                for (x=1; x<=Nxdim-3; x+=2)
                    imgout[(z*Nydim+y)*Nxdim+x]=
                                (0.5)*  (imgout[(z*Nydim+y)*Nxdim+x-1]
                                        +imgout[(z*Nydim+y)*Nxdim+x+1]);
                imgout[(z*Nydim+y)*Nxdim+Nxdim-1]=0.;
            }

            /* ******************************************************** */
            /* ALONG COLUMNS...                                         */
            for (x=0; x<Nxdim; x+=2) {
                for (y=1; y<=Nydim-3; y+=2)
                    imgout[(z*Nydim+y)*Nxdim+x]=
                                (0.5)*  (imgout[(z*Nydim+y-1)*Nxdim+x]
                                        +imgout[(z*Nydim+y+1)*Nxdim+x]);
                imgout[(z*Nydim+Nydim-1)*Nxdim+x]=0.;
            }

            /* ******************************************************** */
            /* AND BETWEEN ROWS AND COLUMNS.                            */
            for (y=1; y<=Nydim-3; y+=2)
                for (x=1; x<=Nxdim-3; x+=2)
                    imgout[(z*Nydim+y)*Nxdim+x]=
                        (0.25)* (imgout[(z*Nydim+y-1)*Nxdim+x-1]
                                +imgout[(z*Nydim+y-1)*Nxdim+x+1]
                                +imgout[(z*Nydim+y+1)*Nxdim+x-1]
                                +imgout[(z*Nydim+y+1)*Nxdim+x+1]);


            /* ******************************************************** */
            /* ALSO REPLICATE EDGE ROWS AND COLUMNS.                    */
            for (y=0; y<Nydim; y++)
		imgout[(z*Nydim+y)*Nxdim+Nxdim-1]=imgout[(z*Nydim+y)*Nxdim-3];

            for (x=0; x<Nxdim; x++)
                imgout[(z*Nydim+Nydim-1)*Nxdim+x]=imgout[(z*Nydim+Nydim-3)*Nxdim+x];
        }


    /* **************************************************************** */
    /* WRITE OUT ZOOMED IMAGE TO OUTPUT FILE                            */
    float_fwrite(argv[3],imgout,output_num_vox);
}
