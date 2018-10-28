#include <stdio.h>
#include <stdlib.h>

unsigned short *unsigned_short_malloc(int num_items);
int read_unsigned_short_int_data(char *filename, unsigned short int *array, int num_items);
int *int_malloc(int num_items);
int read_float_data(char *filename, float *array, int num_items);
float *float_malloc(int num_items);

/* **************************************************************************** */
/* Returns mean value across all scans of a cluster, defined by
   cluster image input as first argument, and the cluster voxel value.          */
main(int argc, char *argv[])
{
    FILE *fIN;
    unsigned short int *cluster, cluster_val;
    int *vox_index_list, num_vox_list;
    int num_voxX2, num_vox;
    int image_index;
    float voxel_sum;
    float *imgin;
    int vox_index, vox_index2;


/* -------------------------------------------------------------------- */
/* Check input arguments.                                               */

    if (argc < 4) {
        printf("clustermeanFP: 16-bit cluster image, cluster value,\n");
        printf("               and names of floating point map file(s).\n");
        exit(1);
    }

/* -------------------------------------------------------------------- */
/*  Open cluster image file and count voxels.  Then allocate memory,    */
/*  and read in voxels.                                                 */

    if ((fIN=fopen(argv[1],"r")) == NULL) {
        printf("clustermeanFP: can't open %s\n",argv[1]);
        exit(2);
    }
    num_voxX2 = 0;
    while (fgetc(fIN)!=EOF)
        num_voxX2++;
    num_vox = num_voxX2/2;
    fclose(fIN);

    cluster = unsigned_short_malloc(num_vox);
    read_unsigned_short_int_data(argv[1], cluster, num_vox);


/* -------------------------------------------------------------------- */
/* Read in cluster voxel value.                                         */

    cluster_val = (unsigned short int) atoi(argv[2]);

/* -------------------------------------------------------------------- */
/* Determine voxel index values of voxels in cluster image which match  */
/* cluster_val.                                                         */

    num_vox_list = 0;
    for (vox_index=0; vox_index<num_vox; vox_index++)
        if (cluster[vox_index]==cluster_val)
            num_vox_list++;

    vox_index_list = int_malloc(num_vox_list);
    vox_index2=0;
    for (vox_index=0; vox_index<num_vox; vox_index++)
        if(cluster[vox_index]==cluster_val) {
            vox_index_list[vox_index2]=vox_index;
            vox_index2++;
        }

/* -------------------------------------------------------------------- */
/* Allocate memory for input floating point maps.                       */

    imgin=float_malloc(num_vox);
    
/* -------------------------------------------------------------------- */
/* Sequentially read in floating point data, calculate mean cluster     */
/* value across all images.                                             */

    voxel_sum=0.;
    for (image_index=3; image_index<argc; image_index++) {

        read_float_data(argv[image_index], imgin, num_vox);

	voxel_sum = 0.;
        for (vox_index2=0; vox_index2<num_vox_list; vox_index2++)
            voxel_sum+=imgin[vox_index_list[vox_index2]];

        printf("%g ",((float)voxel_sum)/((float)num_vox_list));

    }
    printf("\n");

    free(imgin);
    free(cluster);
    free(vox_index_list);

    exit(0);
}

