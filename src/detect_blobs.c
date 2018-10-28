#include <stdio.h>
#include <stdlib.h>
#include <math.h>

main (int argc, char *argv[])
{
    /* 04.19.94 : JM Maisog
	Takes 2 16-bit cluster images (output from blobs program)
	and finds which clusters between the two detect_blobs.
       10.27.94 : JMM
	Modified from overlap.c to detect what non-zero voxel-values
	are present in a 16-bit image, specifically meant for blobs
	output. */

    FILE *fIN;
    signed short int *imgin;
    int num_voxX2,max;
    int *table,tab_index;
    int num_vox;
    int vox_index;
    int read_status;

    /* **************************************************************** */
    /* VERIFY NUMBER OF COMMAND LINE ARGUMENTS.				*/
    if (argc != 2) {
	printf("detect_blobs: need 16-bit input file.\n");
	exit(1);
    }

    /* **************************************************************** */
    /* OPEN FILE FOR INPUT.						*/
    if ((fIN=fopen(argv[1],"r")) == NULL) {
	printf("detect_blobs: can't open %s\n",argv[1]);
	exit(2);
    }

    /* **************************************************************** */
    /* COUNT NUMBER OF BYTES IN INPUT FILE				*/
    num_voxX2 = 0;
    while (fgetc(fIN)!=EOF)
	num_voxX2++;

    num_vox = num_voxX2/2;
/*    printf("detect_blobs: Number of voxels = %d\n",num_vox); */

    /* **************************************************************** */
    /* ASSIGN MEMORY.							*/
    imgin = (signed short int *) malloc(num_vox*sizeof(signed short int));
    if (imgin == NULL) {
	printf("detect_blobs: failed to malloc for imgin.\n");
	exit(3);
    }

    /* **************************************************************** */
    /* REWIND INPUT FILE AND READ IN WORD DATA.				*/
    rewind(fIN);
    fIN=fopen(argv[1],"r");
    read_status=fread(imgin, sizeof(signed short int), num_vox, fIN);
    fclose(fIN);
    if (read_status != num_vox) {
        printf("%d voxels read in.  Should have been %d voxels.\n",
						read_status, num_vox);
        exit(5);
    }

    /* **************************************************************** */
    /* FIND MAXIMA OF CLUSTER IMAGE.					*/
    max=0;
    for (vox_index=0; vox_index<num_vox; vox_index++)
	max = (max<imgin[vox_index]) ? imgin[vox_index] : max;


    /* **************************************************************** */
    /* ALLOCATE MEMORY FOR TABLE OF CLUSTERS.				*/
    /* THEN INITIALIZE IT TO ZERO.					*/
    table = (int *) malloc(sizeof(int)*num_vox);
    if (table==NULL) {
	printf("detect_blobs: unable to malloc for table.\n");
	exit(8);
    }
    for (tab_index=0; tab_index<num_vox; tab_index++)
	table[tab_index]=0;

    /* **************************************************************** */
    /* COUNT VOXEL VALUES, THEN OUTPUT TO STANDARD OUT.			*/
    for (vox_index=0; vox_index<num_vox; vox_index++)
	table[imgin[vox_index]]=1;
    for (tab_index=1; tab_index<num_vox; tab_index++)
	if (table[tab_index]!=0)
	    printf("%d ",tab_index);
    printf("\n");
    free(imgin);
    free(table);
    exit(0);
}
