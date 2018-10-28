#include <stdio.h>
#include <stdlib.h>
#include "analyze6.h"

/* **************************************************************************** */
/* Returns mean value across all scans of a cluster, defined by
   cluster image input as first argument, and the cluster voxel value.          */
main(int argc, char *argv[])
{
    FILE *fp, *fIN;
    signed short *cluster, cluster_val;
    int *vox_index_list, num_vox_list;
    int num_voxX2, num_vox;
    int image_index;
    int voxel_sum;
    signed short *imgin;
    int vox_index, vox_index2, x, y, z;
    int xdim,ydim,zdim;
    int read_status;


/* -------------------------------------------------------------------- */
/* Check input arguments.                                               */

    if (argc < 4) {
        printf("clustermean: 16-bit cluster image, cluster value,\n");
        printf("             and names of 16-bit images file(s).\n");
        exit(1);
    }

/* -------------------------------------------------------------------- */
/*  Open cluster image file and count voxels.  Then allocate memory,    */
/*  and read in voxels.                                                 */

    if ((fIN=fopen(argv[1],"r")) == NULL) {
        printf("clustermean: can't open %s\n",argv[1]);
        exit(2);
    }
    num_voxX2 = 0;
    while (fgetc(fIN)!=EOF)
        num_voxX2++;
    num_vox = num_voxX2/2;

    cluster = (signed short int *) malloc(num_vox*sizeof(signed short int));
    if (cluster == NULL) {
        printf("clustermean: failed to malloc for cluster.\n");
        exit(3);
    }
    rewind(fIN);
    read_status=fread(cluster, sizeof(signed short int), num_vox, fIN);
    fclose(fIN);
    if (read_status != num_vox) {
        printf("%d voxels read in.  Should have been %d voxels.\n",
                                                read_status, num_vox);
        exit(4);
    }

/* -------------------------------------------------------------------- */
/* Read in cluster voxel value.                                         */

    cluster_val = (signed short) atoi(argv[2]);

/* -------------------------------------------------------------------- */
/* Determine voxel index values of voxels in cluster image which match  */
/* cluster_val.                                                         */

    num_vox_list = 0;
    for (vox_index=0; vox_index<num_vox; vox_index++)
        if (cluster[vox_index]==cluster_val)
            num_vox_list++;

    vox_index_list = (int *) malloc(sizeof(int)*num_vox_list);
    if (vox_index_list==NULL) {
        printf("clustermean: failed to malloc for vox_index_list.\n");
        exit(5);
    }
    vox_index2=0;
    for (vox_index=0; vox_index<num_vox; vox_index++)
        if(cluster[vox_index]==cluster_val) {
            vox_index_list[vox_index2]=vox_index;
            vox_index2++;
        }
/* -------------------------------------------------------------------- */
/* Allocate memory for input 16-bit images.                             */

    imgin=(signed short *)malloc(sizeof(signed short)*num_vox);
    if (imgin==NULL) {
        printf("clustermean: failed to malloc for imgin.\n");
        exit(6);
    }
    
/* -------------------------------------------------------------------- */
/* Sequentially read in 16-bit images, calculate mean cluster           */
/* value across all images.                                             */

    voxel_sum=0;
    for (image_index=3; image_index<argc; image_index++) {

        fp = fopen(*(argv+image_index),"r");
            read_status=fread(imgin, sizeof(signed short), num_vox, fp);
            if (read_status!=num_vox)
                printf("ERROR: read in only %d bytes from %s!\n",read_status,
                                                        *(argv+image_index));
        fclose(fp);

        for (vox_index2=0; vox_index2<num_vox_list; vox_index2++)
            voxel_sum+=imgin[vox_index_list[vox_index2]];

    }

    free(imgin);
    printf("%f\n",((float)voxel_sum)/((float)(num_vox_list*(argc-3))));

    exit(0);
}

