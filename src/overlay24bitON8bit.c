#include <stdio.h>
#include <stdlib.h>
#include "analyze6.h"

struct analyze_struct readhdr(char *argv[]);
unsigned char *unsigned_char_malloc(int num_items);
int read_unsigned_char_data(char *filename, unsigned char *array, int num_items);
int unsigned_char_fwrite(char *filename, unsigned char *array, int num_items);


/* ************************************************************************** */
main(int argc,char *argv[])
/*        Reads in an 8-bit and a 24-bit image volume, and overlays the 24-bit*/
/*        image volume onto the 8-bit image volume, producing 24-bit output.  */
{
        struct analyze_struct
                img_info;   /* Header information variable.                   */

        unsigned char *image8bit,*image24bit,*output24bit;
        FILE *fIN1, *fIN2, *fOUT;
        int i,num_vox;
        int xdim,ydim,zdim;
        int x,y,z,color;
        int vox_index8bit,vox_index24bit_red,vox_index24bit_green,vox_index24bit_blue;
        int read_status,write_status;


        /* **************************************************************** */
        /* CHECK COMMAND LINE ARGUMENTS.                                    */

        if (argc < 4) {
            printf("overlay24bitON8bit: need input header, input 8-bit image, input\n");
	    printf("                    24-bit image, and 24-bit output image filename.\n");
            return 1;
        }


        /* **************************************************************** */
        /* READ ANALYZE HEADER.                                             */

        img_info = readhdr(argv+1);

        xdim = img_info.dime.dim[1];
        ydim = img_info.dime.dim[2];
        zdim = img_info.dime.dim[3];
        num_vox = xdim*ydim*zdim;

        printf("overlay24bitON8bit: Number of voxels = %d\n",num_vox);


        /* **************************************************************** */
        /* ASSIGN MEMORY.                                                   */

        image8bit   = unsigned_char_malloc(num_vox);
        image24bit  = unsigned_char_malloc(num_vox*3);
        output24bit = unsigned_char_malloc(num_vox*3);


        /* **************************************************************** */
        /* READ IN DATA.                                                    */
        
        read_unsigned_char_data(argv[2], image8bit, num_vox);
        read_unsigned_char_data(argv[3], image24bit, num_vox*3);


        /* **************************************************************** */
        /* DO OVERLAY.                                                      */

        for (z=0; z<zdim; z++)
            for (y=0; y<ydim; y++)
                for (x=0; x<xdim; x++) {
                    vox_index8bit        = (z*ydim+y)*xdim+x;
                    vox_index24bit_red   = ((z*3+0)*ydim+y)*xdim+x;
                    vox_index24bit_green = ((z*3+1)*ydim+y)*xdim+x;
                    vox_index24bit_blue  = ((z*3+2)*ydim+y)*xdim+x;

                    if ((image24bit[vox_index24bit_red]==0)
                     && (image24bit[vox_index24bit_green]==0)
                     && (image24bit[vox_index24bit_blue]==0)) {
                        output24bit[vox_index24bit_red]=image8bit[vox_index8bit];
                        output24bit[vox_index24bit_green]=image8bit[vox_index8bit];
                        output24bit[vox_index24bit_blue]=image8bit[vox_index8bit];
                    }
                    else {
                        output24bit[vox_index24bit_red]=image24bit[vox_index24bit_red];
                        output24bit[vox_index24bit_green]=image24bit[vox_index24bit_green];
                        output24bit[vox_index24bit_blue]=image24bit[vox_index24bit_blue];
                    }
                }

        
        /* **************************************************************** */
        /* WRITE TO OUTPUT FILE.                                            */

        unsigned_char_fwrite(argv[4], output24bit, num_vox*3);

        free(image8bit);
        free(image24bit);
        free(output24bit);
}

