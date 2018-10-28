#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "nrutil.h"
#include "analyze6.h"
/* #include "nr.h" */

#define MAXSTR 1024
struct analyze_struct readhdr(char *argv[]);
void fourn(float data[], unsigned long nn[], int ndim, int isign);
float *vector(long nl, long nh);
unsigned long *lvector(long nl, long nh);
void free_vector(float *v, long nl, long nh);
void free_lvector(unsigned long *v, long nl, long nh);

int main(int argc, char *argv[])
{
    /* 02.15.95 : Jose' Ma. Maisog, M.D.
       Reads in 16-bit scans, and convolves them with some kernel, e.g., a
       Gaussian smoothing kernel.  List of scans to be smoothed is read in
       from a text file, whose name is indicated on the command line.
       Convolution is done in the frequency domain with FFT's to make it fast.
       Uses C code from Numerical Recipes in C, 2nd ed. pp. 530-531.

       03.29.95 : JMM : Modified to enable program to be run on 8-bit, 16-bit,
       OR floating point images.                                                 
       
       05.22.95 : JMM : Added code to have program read in an 8-bit mask which
       will be applied to all images just before they're smoothed.
       
       08.31.95 : JMM : Added code to allow filenames to include directories.
       
       10.05.95 : JMM : Added method to correct smoothed image for edge artifact.
       
       09.23.96 : JMM : Modified naming of output files so that they would be
			placed in the current working directory, rather than
                        in the same directory as the input files.  This will give
                        the user more flexibility in the placement of the output.*/


                                    /* ***************************************** */
                                    /* FILE I/O.                                 */
    FILE *fp;                       /* General purpose file pointer.             */
    struct analyze_struct img_info; /* Header information variable.              */
    char **img_name;                /* Pointer to pointers to input image file   */
                                    /* names.  Indexing:                         */
                                    /*     img_name[scan][ASCII character]       */
    int tot_num_scans;              /* Total number of scans to be smoothed.     */
    int scan_index;                 /* Index for looping over scans.             */
    char dummy[MAXSTR];             /* Dummy pointer for counting number of      */
                                    /* entries in a text file.                   */
    char compress_command[MAXSTR];  /* Command to compress input scan file.      */
    char output_filename[MAXSTR];   /* Output file name.                         */
    int num_read_in;                /* Number of voxels read in.                 */
    int num_wrote_out;              /* Number of voxels written out.             */
    char compress_flag;             /* Indicates whether or not compression of   */
                                    /* input files is to be done after smoothing.*/
    int index,index2;               /* Indices into character strings.           */
    int length;                     /* Length of directory portion of pathname.  */
    char path[MAXSTR];              /* Directory portion of pathname.            */
    char file[MAXSTR];              /* Filename portion of pathname.             */

                                    /* ***************************************** */
                                    /* INPUT SCAN DATA.                          */
    int num_bits_per_pixel;         /* Bits per pixel in image.                  */
    unsigned char *Bbuffer;         /* Buffer to hold 8-bit image as read in     */
                                    /* from disk.                                */
    signed short int *Wbuffer;      /* Buffer to hold 16-bit image as read in    */
                                    /* from disk.                                */
    float *Fbuffer;                 /* Buffer to hold floating point image as    */
                                    /* read in from disk.                        */
    int x,y,z;                      /* 3D coordinates.                           */

    int xdim1, ydim1, zdim1;        /* Original dimensions of image data.        */
    int num_vox1;                   /* Original number of voxels.                */
    int vox_index1;                 /* Index for looping over voxels.            */

    int xdim2, ydim2, zdim2;        /* Padded dimensions, after doubling scan    */
                                    /* size in each dimension.  Extra padding    */
                                    /* with zeros, to avoid the "wrap-around"    */
                                    /* effect with circular convolution.         */
                                    /* xdim2=xdim2*2,ydim2=ydim2*2,zdim2=zdim2*2 */
    int Hxdim2, Hydim2, Hzdim2;     /* Half of xdim2, ydim2, and zdim2.          */
    int num_vox2;                   /* xdim2*ydim2*zdim2                         */
    int vox_index2;                 /* Index for looping over voxels.            */
    int num_to_z_pad;               /* Number of zero slices to pad with along   */
                                    /* the Z-axis.  At this time, in functional  */
                                    /* imaging the brain is never completely     */
                                    /* traversed along the Z-axis, whether       */
                                    /* acquisition was coronal, sagittal, or     */
                                    /* axial.  In smoothing, we don't want to    */
                                    /* smooth voxels in the first slice with     */
                                    /* voxels in the last slice, which is what   */
                                    /* will happen with circular convolution,    */
                                    /* which is what you get with FFT's, which   */
                                    /* we want to use to make convolution fast!  */
                                    /* To address this problem, we pad the image */
                                    /* volume with num_to_z_pad zero slices      */
                                    /* before the first slice.                   */
                                    /* Edge artifacts will be corrected for.     */

                                    /* ***************************************** */
                                    /* MASK.                                     */

    FILE *fMASK;                    /* File pointer to mask.                     */
    unsigned char *mask;            /* 8-bit mask to be applied to each volume   */
                                    /* before smoothing.                         */
    float *correction;              /* Image of correction factors to get rid of */
                                    /* edge artifact.                            */

                                    /* ***************************************** */
                                    /* SMOOTHING KERNEL.                         */
    float *kernel;                  /* Second floating point image, with which   */
                                    /* floating point copy of 16-bit image will  */
                                    /* be convolved.  This is the smoothing      */
                                    /* kernel.                                   */
    float sum;                      /* Sum of all voxels in gaussian smoothing   */
                                    /* kernel, should be near 1; if not,         */
                                    /* something's wrong.                        */
    int log2;                       /* Integer power of 2.                       */
    float pi;                       /* Pi.                                       */
    float dx,dy,dz;                 /* Distance from origin, in voxels, with     */
                                    /* wrap-around.  E.g., voxel at x=1 has dx   */
                                    /* of 1, and voxel at  x=xdim-1 also has dx  */
                                    /* of 1.                                     */
    float FWHMx, FWHMy, FWHMz;      /* FWHM's in X, Y, and Z.                    */
    float SDx, SDy, SDz;            /* SD's in X, Y, and Z.                      */
    float SDx_factor;               /* (1./(SDx*(float)sqrt(2.*pi)))             */
    float SDy_factor;               /* (1./(SDy*(float)sqrt(2.*pi)))             */
    float SDz_factor;               /* (1./(SDz*(float)sqrt(2.*pi)))             */
    float TWOxVARx_reciprocal;      /* (2.*SDx*SDx))                             */
    float TWOxVARy_reciprocal;      /* (2.*SDy*SDy))                             */
    float TWOxVARz_reciprocal;      /* (2.*SDz*SDz))                             */

                                    /* ***************************************** */
                                    /* FFT's.                                    */
    unsigned long *nn;              /* xdim2, ydim2, zdim2, for fourn.           */
    float *data1,*data2;            /* Internal variables used for FFTs and      */
                                    /* convolution.                              */
    float ad,bc,aMINUSb,cPLUSD,     /* For programming trick which may speed up  */
          aMINUSbTIMEScPLUSD;       /* computation of complex multiplication.    */
    float sign;                     /* Used to assign a sign to final value.     */


    /* ---------------------------------------------------------------- */
    /* CHECK INPUT ARGUMENTS.  Make sure correct number of command      */
    /* line arguments are present.                                      */

    if ((argc!=6) && (argc!=9)) {
        printf("convolve3D: need input header, input text file of scan filenames,\n");
        printf("            name of 8-bit mask, number of zero slices to pad with,\n");
        printf("            and EITHER:\n\n");
        printf("(1) name of input floating point file for smoothing filter,\n\n");
        printf("                            - OR -\n\n");
        printf("(2) FWHM of Gaussian smoothing filter in X, Y, and Z in voxel lengths,\n");
        printf("    and name of output floating point file for smoothing filter.\n");
        exit(1);
    }


#if 0
    /* ---------------------------------------------------------------- */
    /* QUERY USER WHETHER IMAGES ARE TO BE COMPRESSED AFTER SMOOTHING.  */

    printf("\nconvolve3D: do you want the original image files to be compressed\n");
    printf("            after smoothing <y/n>? ");
    scanf("%s",&compress_flag);
    while((compress_flag!='y') && (compress_flag!='n')) {
        printf("  Type 'y' or 'n': ");
        scanf("%s",&compress_flag);
    }
#endif
    compress_flag = 'n';

    /* ---------------------------------------------------------------- */
    /* READ IN HEADER INFO.                                             */
    /* In particular, get the 3D dimensions from the ANALYZE header.    */
    /* readhdr is a routine to read the ANALYZE header, and returns the */
    /* C structure stored in the ANALYZE header.                        */

    img_info = readhdr(argv+1);

    xdim1 = img_info.dime.dim[1];
    ydim1 = img_info.dime.dim[2];
    zdim1 = img_info.dime.dim[3];
    num_bits_per_pixel=img_info.dime.bitpix;
    num_vox1 = xdim1*ydim1*zdim1;


    /* ---------------------------------------------------------------- */
    /* CALCULATE PADDED DIMENSIONS OF INPUT IMAGE.  Need to ensure that */
    /* the scan  dimensions in 3D are all integer powers of 2.  If not, */
    /* scan will be padded to make them so.  This is done by padding    */
    /* the 3D image up so that each dimension is the smallest integer   */
    /* power of 2 greater than or equal to the original 3D image        */
    /* dimensions.  Note especially that zdim2 will be doubled, because */
    /* we need to pad with zero, as discussed above next to the         */
    /* declaration for the variable num_to_z_pad, and we also need a    */
    /* z-dimension which is an integer power of two.                    */

    log2 = (int)(log((double)xdim1)/log(2.));
    if ((int)pow(2.,log2)==xdim1)
        xdim2=xdim1;
    else
        xdim2=(int)pow(2,log2+1);
    printf("\n  Original scan x dimension = %d\n",xdim1);
    printf("  Next integer power of 2 >= %d is %d\n",xdim1,xdim2);
    printf("  After padding with zeros, internal x size is %d.\n",xdim2);
    Hxdim2=xdim2/2;

    log2 = (int)(log((double)ydim1)/log(2.));
    if ((int)pow(2.,log2)==ydim1)
        ydim2=ydim1;
    else
        ydim2=(int)pow(2,log2+1);
    printf("\n  Original scan y dimension = %d\n",ydim1);
    printf("  Next integer power of 2 >= %d is %d\n",ydim1,ydim2);
    printf("  After padding with zeros, internal y size is %d.\n",ydim2);
    Hydim2=ydim2/2;

    num_to_z_pad=atoi(argv[4]);
    log2 = (int)(log((double)(zdim1+num_to_z_pad))/log(2.));
    if ((int)pow(2.,log2)==zdim1+num_to_z_pad)
        zdim2=zdim1+num_to_z_pad;
    else
        zdim2=(int)pow(2,log2+1);
    printf("\n  Original scan z dimension = %d\n",zdim1);
    printf("  Will pad with %d zero slices before first slice.\n",num_to_z_pad);
    printf("  After padding, there will be %d slices.\n",zdim1+num_to_z_pad);
    printf("  Next integer power of 2 >= %d is %d\n",zdim1+num_to_z_pad,zdim2);
    printf("  After padding to next integer power of 2, internal z size is %d.\n",zdim2);
    Hzdim2=zdim2/2;

    num_vox2=xdim2*ydim2*zdim2;
    nn=lvector(1,3);
    nn[1]=(unsigned long)zdim2;
    nn[2]=(unsigned long)ydim2;
    nn[3]=(unsigned long)xdim2;


    /* ---------------------------------------------------------------- */
    /* COUNT NUMBER OF SCANS IN INPUT TEXT FILE.  THEN ALLOCATE MEMORY  */
    /* FOR img_name, AND THEN READ IN NAMES OF INPUT IMAGE FILES.       */
    /* List of scan filenames to be read in will be held in img_name[]. */

    if ((fp = fopen(argv[2],"r")) == NULL) {
        printf("convolve3D: can't open %s\n",argv[2]);
        exit(2);
    }
    tot_num_scans = 0;
    while(fscanf(fp,"%s",dummy)!=EOF)
        tot_num_scans++;

    img_name = (char **) malloc(tot_num_scans*sizeof(char *));
    if (img_name==NULL) {
        printf("convolve3D: error malloc'ing for img_name.\n");
        exit(3);
    }
    rewind(fp);
    for (scan_index=0; scan_index<tot_num_scans; scan_index++) {
        img_name[scan_index] = (char *) malloc(MAXSTR*sizeof(char));
        if (img_name[scan_index]==NULL) {
            printf("convolve3D: error malloc'ing for img_name[%d].\n",
                                                    scan_index);
            exit(4);
        }
        fscanf(fp,"%s\0",img_name[scan_index]);
    }
    fclose(fp);
    printf("\nconvolve3D: %d scan(s) to be smoothed:\n",tot_num_scans);
    for (scan_index=0; scan_index<tot_num_scans; scan_index++)
        printf("  %s\n",img_name[scan_index]);


    /* ---------------------------------------------------------------- */
    /* READ IN MASK.                                                    */

    printf("\nconvolve3D: reading in mask from %s...\n",argv[3]);
    mask = (unsigned char *)malloc(num_vox1*sizeof(unsigned char));
    if (mask==NULL) {
        printf("convolve3D: unable to malloc for mask.\n");
        exit(5);
    }
    if ((fMASK = fopen(argv[3],"r")) == NULL) {
        printf("convolve3D: can't open %s\n",argv[3]);
        exit(6);
    }
    num_read_in=fread(mask,sizeof(unsigned char),num_vox1,fMASK);
    fclose(fMASK);
    if(num_read_in!=num_vox1) {
        printf("convolve3D: %d voxels read in; should have been %d.\n",
                                        num_read_in,num_vox1);
        exit(7);
    }

    /* ---------------------------------------------------------------- */
    /* ALLOCATE MEMORY FOR FFT's AND FOR FLOATING POINT IMAGE ARRAY FOR */
    /* INTERNAL CALCULATIONS.  vector() is a Numerical Recipes in C     */
    /* routine which allocates memory to a one-offset pointer, not a    */
    /* zero-offset pointer.  num_vox2 is used because we're using the   */
    /* padded dimensions of the image.  The number of elements in the   */
    /* data vectors is num_vox2*2 because we need to specify both real  */
    /* and imaginary components.  Of course, for 16-bit raw fMR data,   */
    /* the imaginary component is zero for all voxels.                  */

    data1=vector(1,num_vox2*2);
    data2=vector(1,num_vox2*2);


    /* ---------------------------------------------------------------- */
    /* ALLOCATE MEMORY FOR SMOOTHING KERNEL, 2ND FLOATING POINT IMAGE.  */
    /* If a second floating point image file was given on the command   */
    /* line, read it in.  It will be used as the kernel with which the  */
    /* first input image will be convolved.  If instead three numbers   */
    /* were input, read  them in and calculate a gaussian smoothing     */
    /* kernel with those numbers as FWHM'S in X, Y, and Z.              */

    kernel      = (float *) malloc(sizeof(float)*num_vox2);
    if (kernel==NULL) {
        printf("convolve3D: unable to malloc for kernel.\n");
        exit(8);
    }
    if (argc==6) {
        printf("convolve3D: loading %s...\n",argv[5]);
        if ((fp = fopen(argv[5],"r")) == NULL) {
            printf ("convolve3D: can't open %s\n",argv[5]);
            exit(9);
        }
        num_read_in=fread(kernel,sizeof(float),num_vox2,fp);
        fclose(fp);
        if(num_read_in!=num_vox2) {
            printf("convolve3D: %d voxels read in; should have been %d.\n",
                                            num_read_in,num_vox2);
            exit(10);
        }
    }
    else {
        FWHMx=atof(argv[5]);
        FWHMy=atof(argv[6]);
        FWHMz=atof(argv[7]);
        SDx=FWHMx/(float)sqrt(8.*log(2.));
        SDy=FWHMy/(float)sqrt(8.*log(2.));
        SDz=FWHMz/(float)sqrt(8.*log(2.));
        pi = acos(0.)*2.;
        SDx_factor=(1./(SDx*(float)sqrt(2.*pi)));
        SDy_factor=(1./(SDy*(float)sqrt(2.*pi)));
        SDz_factor=(1./(SDz*(float)sqrt(2.*pi)));
        TWOxVARx_reciprocal=1./(2.*SDx*SDx);
        TWOxVARy_reciprocal=1./(2.*SDy*SDy);
        TWOxVARz_reciprocal=1./(2.*SDz*SDz);
        sum=0.;
        printf("\nconvolve3D: making gaussian smoothing kernel, FWHM=[%.3f, %.3f, %.3f]\n",
                                                                        FWHMx,FWHMy,FWHMz);
        for (z=0; z<zdim2; z++)
            for (y=0;y<ydim2; y++)
                for (x=0; x<xdim2; x++) {
                    vox_index2=(z*ydim2+y)*xdim2+x;

                    if (x<Hxdim2) dx=(float)x;
                    else          dx=(float)(xdim2-x);

                    if (y<Hydim2) dy=(float)y;
                    else          dy=(float)(ydim2-y);

                    if (z<Hzdim2) dz=(float)z;
                    else          dz=(float)(zdim2-z);

                    if ((FWHMx==0.) && (x!=0))
                        kernel[vox_index2]=0.;
                    else if ((FWHMy==0.) && (y!=0))
                        kernel[vox_index2]=0.;
                    else if ((FWHMz==0.) && (z!=0))
                        kernel[vox_index2]=0.;
                    else {
                        kernel[vox_index2]=1.;

                        if (FWHMx!=0.) 
                            kernel[vox_index2]*=
                                (SDx_factor
                                     *
                                exp(-dx*dx*TWOxVARx_reciprocal));

                        if (FWHMy!=0.) 
                            kernel[vox_index2]*=
                                (SDy_factor
                                     *
                                exp(-dy*dy*TWOxVARy_reciprocal));

                        if (FWHMz!=0.) 
                            kernel[vox_index2]*=
                                (SDz_factor
                                     *
                                exp(-dz*dz*TWOxVARz_reciprocal));
                    }
                    sum+=kernel[vox_index2];

/*                  Originally, the above complicated "if" structure
                    was the following relatively simple formula
                    for a 3D gaussian.  The "if" structure allows
                    for a FWHM of zero along one or more dimensions.
                      kernel[vox_index2]=
                        (1./(SDx*(float)sqrt(2.*pi)))
                                 *
                        exp(-dx*dx/(2.*SDx*SDx))
                                 *
                        (1./(SDy*(float)sqrt(2.*pi)))
                                 *
                        exp(-dy*dy/(2.*SDy*SDy))
                                 *
                        (1./(SDz*(float)sqrt(2.*pi)))
                                 *
                        exp(-dz*dz/(2.*SDz*SDz));
*/
                }

        printf("convolve3D: sum of all voxels in smoothing kernel = %f\n",sum);
        printf("            this should be nearly 1.0, otherwise something may be\n");
        printf("            wrong (e.g., FWHM too big in Z dimension).\n");


        /* NORMALIZE SUM OF SMOOTHING KERNEL VOXELS TO 1.0.             */
        /* So that smoothing doesn't cause a net brightening or dimming */
        /* of global intensity.                                         */
        printf("convolve3D: normalizing sum of all voxels in kernel to 1.0...\n");
        for (vox_index2=0; vox_index2<num_vox2; vox_index2++)
            kernel[vox_index2]/=sum;

        if ((fp = fopen(argv[8],"w")) == NULL) {
            printf("convolve3D: can't open %s\n",argv[8]);
            printf("convolve3D: will proceed with smoothing anyway.\n");
        }
        else {
            printf("convolve3D: writing smoothing kernel to floating point file %s.\n",
                                                                argv[8]);
            printf("convolve3D: dimensions of kernel will be [%d %d %d]\n",xdim2,ydim2,zdim2);
            num_wrote_out=fwrite(kernel,sizeof(float),num_vox2,fp);
            fclose(fp);
        }
    }


    /* ---------------------------------------------------------------- */
    /* CALCULATE FFT OF SMOOTHING KERNEL.                               */
    /* First, transfer data in kernel[] to data2[].                     */
    /* Then, send data2[] to routine fourn(), which returns the 3D FFT  */
    /* of data2[] in the same data space, data2[].                      */
    /* Now, note that the FFT as calculated by Numerical Recipes in C   */
    /* is different from the FFT as calculated by MATLAB.  They are     */
    /* complex conjugatest of each other.  This is because MATLAB fol-  */
    /* lows the engineering convention, where the sign of the exponent  */
    /* in the Fourier Transform is negative, whereas Numerical Recipes  */
    /* in C follows the physicists' convention, where the sign is       */
    /* positive.  As long as you stick to one convention, you're okay.  */
    /* So, I'll stick with the Numerical Recipes in C convention.       */

    for (vox_index2=0; vox_index2<num_vox2; vox_index2++) {
        data2[2*vox_index2+1]=kernel[vox_index2];
        data2[2*vox_index2+2]=0.;
    }
    free(kernel);

    if (argc==6)
        printf("convolve3D: calculating 3D FFT of %s...\n",argv[5]);
    else
        printf("convolve3D: calculating 3D FFT of gaussian smoothing kernel...\n");

    fourn(data2,nn,3,1);


    /* ---------------------------------------------------------------- */
    /* GENERATE THE MAP OF CORRECTION FACTORS.                          */

        /* ------------------------------------------------------------ */
        /* FIRST, ALLOCATE MEMORY FOR CORRECTION FACTOR IMAGE IN ZERO-  */
        /* PADDED IMAGE SPACE.                                          */

        printf("\nconvolve3D: making map of correction factors...\n");
        correction = (float *) malloc(sizeof(float)*num_vox2);
        if (correction==NULL) {
            printf("convolve3D: unable to malloc for correction image.\n");
            exit(11);
        }

        /* ------------------------------------------------------------ */
        /* THEN, ZERO OUT ZERO-PADDED IMAGE SPACE AND CORRECTION MAP,   */
        /* AND TRANSFER MASK INTO ZERO-PADDED IMAGE SPACE.              */

        for (vox_index2=0; vox_index2<num_vox2; vox_index2++) {
            data1[2*vox_index2+1]=0.;
            data1[2*vox_index2+2]=0.;
            correction[num_vox2]=0.;
        }
        for (z=0; z<zdim1; z++)
            for (y=0;y<ydim1; y++)
                for (x=0; x<xdim1; x++) {
                    vox_index1=(z*ydim1+y)*xdim1+x;
                    vox_index2=((z+num_to_z_pad)*ydim2+y)*xdim2+x;
                    data1[2*vox_index2+1]=(float)mask[vox_index1];
                }

        /* ------------------------------------------------------------ */
        /* THEN, CALCULATE FFT OF MASK.                                 */

        printf("\nconvolve3D: computing FFT of mask...\n");
        fourn(data1,nn,3,1);

        /* ------------------------------------------------------------ */
        /* THEN, CONVOLVE MASK WITH SMOOTHING KERNEL BY MULTIPLYING     */
        /* THEIR FFT's WITH EACH OTHER.  The use of the intermediate    */
        /* variables ad, bc, aMINUSb, cPLUSDcPLUSD, and                 */
        /* aMINUSbTIMEScPLUSd is just a trick to perform only three     */
        /* multiplications rather than four for every complex multi-    */
        /* plication performed.                                         */

        for (vox_index2=0; vox_index2<num_vox2; vox_index2++) {
            ad                 = data1[2*vox_index2+1]*data2[2*vox_index2+2];
            bc                 = data1[2*vox_index2+2]*data2[2*vox_index2+1];
            aMINUSb            = data1[2*vox_index2+1]-data1[2*vox_index2+2];
            cPLUSD             = data2[2*vox_index2+1]+data2[2*vox_index2+2];
            aMINUSbTIMEScPLUSD = aMINUSb*cPLUSD;

            /* REAL COMPONENT.                                          */
            data1[2*vox_index2+1]=aMINUSbTIMEScPLUSD-ad+bc;

            /* IMAGINARY COMPONENT.                                     */
            data1[2*vox_index2+2]=ad+bc;
        }

        /* ------------------------------------------------------------ */
        /* TAKE INVERSE FFT OF COMPLEX PRODUCT JUST CALCULATED, GIVING  */
        /* THE SMOOTHED MASK.  THIS IS THE MAP OF CORRECTION FACTORS.   */
        /* Need to divide by num_vox2 when performing inverse FT with   */
        /* the routine fourn.                                           */

        fourn(data1,nn,3,-1);

        for (vox_index2=0; vox_index2<num_vox2; vox_index2++) {
            data1[2*vox_index2+1]/=((float)(num_vox2));
            data1[2*vox_index2+2]/=((float)(num_vox2));
        }

        for (vox_index2=0; vox_index2<num_vox2; vox_index2++) {
            if (data1[2*vox_index2+1]<0) sign = -1.;
                else sign = 1.;
            correction[vox_index2]= sign *
                            (float)(sqrt((double)
                            (data1[2*vox_index2+1]*data1[2*vox_index2+1]
                                                    +
                             data1[2*vox_index2+2]*data1[2*vox_index2+2])
                                                    ));
        }

        /* ------------------------------------------------------------ */
        /* WRITE OUT CORRECTION MAP TO DISK.                            */
        /************** This section commented out since correction ************
	 ************** image is probably not of much immediate interest. ******

        printf("convolve3D: writing correction.img...\n");
        if ((fp = fopen("correction.img","w")) == NULL) {
            printf ("convolve3D: can't open correction.img.\n");
            exit(12);
        }
        num_wrote_out=fwrite(correction,sizeof(float),num_vox2,fp);
        fclose(fp);
        if(num_wrote_out!=num_vox2) {
            printf("convolve3D: %d voxels wrote out; should have been %d.\n",
                                            num_wrote_out,num_vox2);
            exit(13);
        }
        */

    /* ---------------------------------------------------------------- */
    /* ALLOCATE MEMORY FOR BUFFER.                                      */
    /* Make an internal buffer to hold input image data as read from    */
    /* disk.  May need to be converted to floating point for FFT.       */

    if (num_bits_per_pixel==8) {
        Bbuffer = (unsigned char *) malloc(sizeof(unsigned char)*num_vox1);
        if (Bbuffer==NULL) {
            printf("convolve3D: unable to malloc for Bbuffer.\n");
            exit(14);
        }
    }
    else if (num_bits_per_pixel==16) {
        Wbuffer = (signed short int *) malloc(sizeof(signed short int)*num_vox1);
        if (Wbuffer==NULL) {
            printf("convolve3D: unable to malloc for Wbuffer.\n");
            exit(15);
        }
    }
    else if (num_bits_per_pixel==32) {
        Fbuffer = (float *) malloc(sizeof(float)*num_vox1);
        if (Fbuffer==NULL) {
            printf("convolve3D: unable to malloc for Fbuffer.\n");
            exit(16);
        }
    }
    printf("num_bits_per_pixel = %d\n",num_bits_per_pixel);

    /* ---------------------------------------------------------------- */
    /* NOW LOOP OVER INPUT SCANS, SMOOTH EACH, AND WRITE OUT TO DISK.   */

    printf("Total number of scans = %d\n",tot_num_scans);

    for (scan_index=0; scan_index<tot_num_scans; scan_index++) {

	("scan # %d...\n",scan_index+1);

        /* ---------------------------------------------------------------- */
        /* LOAD NEXT SCAN INTO MEMORY, THEN COMPRESS IT.                    */

        if ((fp = fopen(img_name[scan_index],"r")) == NULL) {
            printf("convolve3D: can't open %s\n",img_name[scan_index]);
        }
        else {
            printf("\nconvolve3D: loading %s...\n",img_name[scan_index]);
            if (num_bits_per_pixel==8)
                num_read_in=fread(Bbuffer,sizeof(unsigned char),num_vox1,fp);
            else if (num_bits_per_pixel==16)
                num_read_in=fread(Wbuffer,sizeof(signed short int),num_vox1,fp);
            else if (num_bits_per_pixel==32)
                num_read_in=fread(Fbuffer,sizeof(float),num_vox1,fp);
            fclose(fp);
            if(num_read_in!=num_vox1) {
                printf("convolve3D: %d voxels read in; should have been %d.\n",
                                                num_read_in,num_vox1);
                exit(17);
            }
            if (compress_flag=='y') {
                sprintf(compress_command,"gzip %s",img_name[scan_index]);
                printf("convolve3D: %s...\n",compress_command);
                system(compress_command);
            }


            /* ---------------------------------------------------------------- */
            /* TRANSFER INPUT IMAGE DATA IN *buffer INTO *data1.                */
            /* Zero out data1 volume, then put contents of buffer in center of  */
            /* data1 volume.  This will in effect pad it with zeros.            */

            for (vox_index2=0; vox_index2<num_vox2; vox_index2++) {
                data1[2*vox_index2+1]=0.;
                data1[2*vox_index2+2]=0.;
            }
            printf("convolve3D: transfering %s in memory, with masking...\n",img_name[scan_index]);
            for (z=0; z<zdim1; z++)
                for (y=0;y<ydim1; y++)
                    for (x=0; x<xdim1; x++) {
                        vox_index1=(z*ydim1+y)*xdim1+x;
                        vox_index2=((z+num_to_z_pad)*ydim2+y)*xdim2+x;
                        if (num_bits_per_pixel==8)
                            data1[2*vox_index2+1]=(float)Bbuffer[vox_index1]*(float)mask[vox_index1];
                        else if (num_bits_per_pixel==16)
                            data1[2*vox_index2+1]=(float)Wbuffer[vox_index1]*(float)mask[vox_index1];
                        else if (num_bits_per_pixel==32)
                            data1[2*vox_index2+1]=Fbuffer[vox_index1]*(float)mask[vox_index1];
                    }


            /* ---------------------------------------------------------------- */
            /* FFT INPUT SCAN VOLUME.                                           */

            printf("convolve3D: 3D FFT'ing %s...\n",img_name[scan_index]);
            fourn(data1,nn,3,1);


            /* ---------------------------------------------------------------- */
            /* CONVOLVE THE TWO IMAGES IN SPECTRAL DOMAIN BY MULTIPLYING THE    */
            /* TWO FFT'S BY EACH OTHER.  Note indexing goes from 0 up to but    */
            /* not including num_vox2.  Complex multiplication.                 */
            /* Note: used to perform this the regular way, now commented out,   */
            /* but now using a programming trick to speed up computation of     */
            /* complex multiplication.  Can multiply two complex numbers using  */
            /* "only" 3 multiplications, 5 additions, and 2 negations.          */

            printf("convolve3D: smoothing %s...\n",img_name[scan_index]);
            for (vox_index2=0; vox_index2<num_vox2; vox_index2++) {

                ad                 = data1[2*vox_index2+1]*data2[2*vox_index2+2];
                bc                 = data1[2*vox_index2+2]*data2[2*vox_index2+1];
                aMINUSb            = data1[2*vox_index2+1]-data1[2*vox_index2+2];
                cPLUSD             = data2[2*vox_index2+1]+data2[2*vox_index2+2];
                aMINUSbTIMEScPLUSD = aMINUSb*cPLUSD;

                /* REAL COMPONENT.                                              */
                /* data3[2*vox_index2+1]=data1[2*vox_index2+1]*data2[2*vox_index2+1]
                                                           -
                                      data1[2*vox_index2+2]*data2[2*vox_index2+2]; */
                data1[2*vox_index2+1]=aMINUSbTIMEScPLUSD-ad+bc;

                /* IMAGINARY COMPONENT.                                         */
                /* data3[2*vox_index2+2]=data1[2*vox_index2+1]*data2[2*vox_index2+2]
                                                           +
                                      data1[2*vox_index2+2]*data2[2*vox_index2+1]; */
                data1[2*vox_index2+2]=ad+bc;
            }


            /* ---------------------------------------------------------------- */
            /* INVERSE FFT THE PRODUCT OF THE TWO FFT's.  THEN APPLY            */
            /* CORRECTIONS FOR EDGE ARTIFACTS, THEN PUT DATA BACK INTO buffer   */
            /* ARRAY.  Note use of xdim1, ydim1, and zdim1 rather than xdim2,   */
            /* ydim2, and zdim2 in the nested loops; this is to return an image */
            /* whose dimensions are the same as the original input image.  Note */
            /* division by num_vox2; this is necessary when using the fourn     */
            /* routine to do the inverse FFT.                                   */
            /*  The floor(sqrt(real^2 + imag^2)+0.5) construction is to round   */
            /* the magnitude sqrt(real^2 + imag^2) to the nearest integer.      */
            /* Ideally, the inverse FFT should have no imaginary component; it  */
            /* maybe just as well to simply set buffer[vox_index1] equal to     */
            /* data1[2*vox_index2+1], the real component.                       */

            fourn(data1,nn,3,-1);

            for (vox_index2=0; vox_index2<num_vox2; vox_index2++) {
                data1[2*vox_index2+1]/=((float)(num_vox2));
                data1[2*vox_index2+2]/=((float)(num_vox2));
            }

            for (z=0; z<zdim1; z++)
                for (y=0; y<ydim1; y++)
                    for (x=0; x<xdim1; x++) {
                        vox_index1=(z*ydim1+y)*xdim1+x;
                        vox_index2=((z+num_to_z_pad)*ydim2+y)*xdim2+x;
                        if (data1[2*vox_index2+1]<0) sign = -1.;
                            else sign = 1.;
                        if (mask[vox_index1]!=0) {
                            data1[2*vox_index2+1]/=correction[vox_index2];
                            data1[2*vox_index2+2]/=correction[vox_index2];
                        }
                        if (num_bits_per_pixel==8)
                            Bbuffer[vox_index1]=
                                        (unsigned char)
                                        floor((double)mask[vox_index1]*sqrt((double)
                                        (data1[2*vox_index2+1]*data1[2*vox_index2+1]
                                                                +
                                         data1[2*vox_index2+2]*data1[2*vox_index2+2])
                                                                )+0.5);
                        else if (num_bits_per_pixel==16)
                            Wbuffer[vox_index1]= sign *
                                        (signed short int)
                                        floor((double)mask[vox_index1]*sqrt((double)
                                        (data1[2*vox_index2+1]*data1[2*vox_index2+1]
                                                                +
                                         data1[2*vox_index2+2]*data1[2*vox_index2+2])
                                                                )+0.5);
                        else if (num_bits_per_pixel==32)
                            Fbuffer[vox_index1]= sign *
                                        (float)((double)mask[vox_index1]*sqrt((double)
                                        (data1[2*vox_index2+1]*data1[2*vox_index2+1]
                                                                +
                                         data1[2*vox_index2+2]*data1[2*vox_index2+2])
                                                                ));

                    }



            /* ---------------------------------------------------------------- */
            /* DETERMINE NAME OF OUTPUT FILE.  The output file will have the    */
            /* same name as the unsmoothed input file, with an "S" prepended to */
            /* the name.  E.g., if the name of the unsmoothed input file is     */
            /* fmri.img, the name of the smoothed output file will be Sfmri.img */

            length = strlen(img_name[scan_index]);
            index=length;
            while ((index>0) && (img_name[scan_index][index]!='/'))
                index--;
            if (index>0) {
                for (index2=0;index2<=index; index2++)
                    path[index2]=img_name[scan_index][index2];
                path[index2]='\0';
                for (index2=0;index2<length-index; index2++)
                    file[index2]=img_name[scan_index][index+index2+1];
                 sprintf(output_filename,"%sS%s",path,file);
/*                sprintf(output_filename,"S%s",file);*/
            }
            else
                sprintf(output_filename,"S%s",img_name[scan_index]);

            /* ---------------------------------------------------------------- */
            /* WRITE OUT OUTPUT TO DISK.                                        */
            /* Sfmri.img will be of the same dimensions as fmri.img, even       */
            /* internally it was padded with zeros.                             */

            printf("convolve3D: writing output to %s...\n",output_filename);
            if ((fp = fopen(output_filename,"wb")) == NULL) {
                printf("convolve3D: can't open %s\n",output_filename);
                exit(18);
            }
            if (num_bits_per_pixel==8)
                num_wrote_out=fwrite(Bbuffer,sizeof(unsigned char),num_vox1,fp);
            else if (num_bits_per_pixel==16)
                num_wrote_out=fwrite(Wbuffer,sizeof(signed short int),num_vox1,fp);
            else if (num_bits_per_pixel==32)
                num_wrote_out=fwrite(Fbuffer,sizeof(float),num_vox1,fp);
            fclose(fp);
            if(num_wrote_out!=num_vox1) {
                printf("convolve3D: %d voxels wrote out; should have been %d.\n",
                                                num_wrote_out,num_vox1);
                exit(19);
            }
        }
    }

    /* ---------------------------------------------------------------- */
    /* DEALLOCATE MEMORY.                                               */

    free(Bbuffer);
    free(Wbuffer);
    free(Fbuffer);
    free(img_name);
    free(mask);
/*    free(kernel); */
    free(correction);
    free_vector(data1,1,num_vox2*2);
    free_vector(data2,1,num_vox2*2);
    free_lvector(nn,1,3);

    return 0;
}
