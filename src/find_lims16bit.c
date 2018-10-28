void find_lims16bit(unsigned short int *mat, int xdim, int ydim, int zdim,
		int *xlo, int *xhi, int *ylo, int *yhi)

/* For use with fread method of reading in image data.			      */
/*	Takes as input a pointer mat to a matrix with dimensions xdim, ydim,  */
/*	and zdim.  Finds the maximum value in that matrix.  Finds the lowest  */
/*	and highest x and y coordinates	for all	entries	which are one third   */
/*	the value of the maximum or greater.				      */
/* 07.20.93 : JMM : Modified from find_lims2 to handle 16-bit data. */
{
    int x,y,z;
    int i;
    unsigned short int max=0;

    /* FIND MAX	VALUE.							      */
    for (i=0;i<=xdim*ydim*zdim-1;i++)
	max=(*(mat+i)>max) ? *(mat+i) : max;

    printf("xdim = %d, ydim = %d, zdim = %d\n",xdim,ydim,zdim);
    printf("maximum pixel value = %d\n",max);

    /* SET INITIAL VALUES OF xlo, xhi, ylo, yhi				      */
    *xlo=xdim-1;
    *xhi=0;
    *ylo=ydim-1;
    *yhi=0;
    
    /* FIND xlo, xhi, ylo, yhi						      */
    for (z=0;z<=zdim-1;z++)
	for (x=0;x<=xdim-1;x++)
	    for (y=0;y<=ydim-1;y++) {
		*xlo=(((int)*(mat+z*ydim*xdim+y*xdim+x)>=(int)max/3) && (x<*xlo)) ? x : *xlo;
		*xhi=(((int)*(mat+z*ydim*xdim+y*xdim+x)>=(int)max/3) && (x>*xhi)) ? x : *xhi;
		*ylo=(((int)*(mat+z*ydim*xdim+y*xdim+x)>=(int)max/3) && (y<*ylo)) ? y : *ylo;
		*yhi=(((int)*(mat+z*ydim*xdim+y*xdim+x)>=(int)max/3) && (y>*yhi)) ? y : *yhi;
	    }
}
