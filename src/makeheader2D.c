#include <stdio.h>
#include <stdlib.h>
#include "analyze6.h"

main(argc,argv)
/* makeheader4 : Program to make ANALYZE headers, based on information in the
   ANALYZE header, written by Jose' Ma. Maisog.
*/
int argc;
char **argv;
{
struct  analyze_struct hdr;
FILE *fp,*fopen();
int ret;

if ((argc != 15) && (argc != 16))
{
     printf("makeheader2D: Need filename, xdim, ydim, zdim, x pixdim, y pixdim, z pixdim,\n");
     printf("             volumes, bits/pixel, global max, global min, calibration max,\n");
     printf("             calibration min, AND datatype.\n");
     printf("             You can also add a scalefactor for SPM95 after the datatype.\n");
     printf("\n");
     printf("             CODE DATATYPE\n");
     printf("              0   Unknown\n");
     printf("              1   Binary (1 bit per voxel)\n");
     printf("              2   Unsigned character (8 bits per voxel)\n");
     printf("              4   Signed short (16 bits per voxel\n");
     printf("              8   Signed integer (32 bits per voxel\n");
     printf("             16   Floating point (32 bits per voxel)\n");
     printf("             32   Complex (64 bits per voxel)\n");
     printf("             64   Double (64 bits per voxel)\n");
     exit(1);
}
if ((fp=fopen(argv[1],"w"))==0) /*filename */
{
   printf ("unable to create: %s\n",argv[1]);
   exit(2);
}
memset(&hdr, 0, sizeof(struct analyze_struct));
hdr.hk.sizeof_hdr       = 348;
hdr.hk.data_type[0]     = '\0';
hdr.hk.db_name[0]       = '\0';
hdr.hk.extents          = 0;                 /* Now set to 0 to conform to SPM94 header.               */
hdr.hk.session_error    = 0;
hdr.hk.regular          = 'r';                   /* Indicates that all images and vols are same size   */
hdr.hk.hkey_un0         = '\0';

hdr.dime.dim[0]         = 2;                     /* Number of dimensions.  Here set to 4 instead of 3. */
hdr.dime.dim[1]         = atoi(argv[2]);         /* x dimension.                                       */
hdr.dime.dim[2]         = atoi(argv[3]);         /* y dimension.                                       */
hdr.dime.dim[3]         = atoi(argv[4]);         /* z dimension.                                       */
hdr.dime.dim[4]         = atoi(argv[8]);
hdr.dime.vox_units[0]   = '\0';
hdr.dime.cal_units[0]   = '\0';
hdr.dime.unused1        = 0;
hdr.dime.datatype       = atoi(argv[14]);
hdr.dime.bitpix         = atoi(argv[9]);
hdr.dime.dim_un0        = 0;
hdr.dime.pixdim[1]      = (double) atof(argv[5]); /* Pixel size in x dimension                          */
hdr.dime.pixdim[2]      = (double) atof(argv[6]); /* Pixel size in y dimension                          */
hdr.dime.pixdim[3]      = (double) atof(argv[7]); /* Pixel size in z dimension                          */
hdr.dime.pixdim[4]      = 0;                      /* Pixel size in time dimension.                      */
hdr.dime.pixdim[5]      = 0;
hdr.dime.pixdim[6]      = 0;
hdr.dime.pixdim[7]      = 0;
hdr.dime.pixdim[0]      = 0;
hdr.dime.vox_offset     = 0.;
if (argc==15)
    hdr.dime.funused1       = 1.;
else
    hdr.dime.funused1       = atof(argv[15]);
hdr.dime.funused2       = 0.;
hdr.dime.funused3       = 0.;
hdr.dime.cal_max        = (double) atof(argv[12]);
hdr.dime.cal_min        = (double) atof(argv[13]);
hdr.dime.compressed     = 0;
hdr.dime.verified       = 0.;
hdr.dime.glmax          = atoi(argv[10]);        /* Max pixel value for entire data base               */
hdr.dime.glmin          = atoi(argv[11]);        /* Min pixel value for entire data base               */

hdr.dh.descrip[0]       = '\0';
hdr.dh.aux_file[0]      = '\0';
hdr.dh.orient           = '0';                   /* Transverse unflipeed                               */
hdr.dh.originator[0]    = '\0';
hdr.dh.generated[0]     = '\0';
hdr.dh.scannum[0]       = '\0';
hdr.dh.patient_id[0]    = '\0';
hdr.dh.exp_date[0]      = '\0';
hdr.dh.exp_time[0]      = '\0';
hdr.dh.hist_un0[0]      = '\0';
hdr.dh.views            = 0;
hdr.dh.vols_added       = 0;
hdr.dh.start_field      = 0;
hdr.dh.field_skip       = 0;
hdr.dh.omax             = 0;
hdr.dh.omin             = 0;
hdr.dh.smax             = 0;
hdr.dh.smin             = 0;

if((ret=fwrite(&hdr,sizeof(struct analyze_struct),1,fp)) !=1)
{
    printf("error writing to %s %d \n",argv[1],ret);
    exit(3);
}
fclose(fp);
}
