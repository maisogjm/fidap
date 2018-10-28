#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>

#define NMAX  50000
#define ALPHA 1.0
#define BETA  0.5
#define GAMMA 2.0
#define GET_PSUM for (j=1;j<=ndim;j++){for (i=1,sum=0.0;i<=mpts;i++)\
		      sum+=p[i][j];psum[j]=sum;}

double  amotry(double **p,double *y,double *psum,int ndim,
		double (*funk)(double *param),int ihi,int *nfunk,double fac);
double *vector(int nl, int nh);
void   nerror(char *error_text);
void   free_vector(double *v,int nl);
/* ************************************************************************** */
int amoebaD(double **p,double *y,int ndim,double ftol,double (*funk) (double *param),
								int *nfunk)
{
	int i,j,ilo,ihi,inhi,mpts=ndim+1;
	double ytry,ysave,sum,rtol,*psum;

	/* TIME VARIABLES, DELETE IF NECESSARY */
	char *s;
	time_t *tp;
	size_t smax = 10;
	tp       = (time_t *) malloc(4);
	s        = (char *)   malloc(26);
	printf(" TIME  NUM_F_EVALS  y[ihi]      y[ilo]     RTOL      BEST PARAMETER(1), PARAMETER(2), ...\n");


	psum=vector(1,ndim);
	*nfunk=0;
	GET_PSUM
	for(;;) {
	    ilo=1;
	    ihi=(y[1]>y[2]) ? (inhi=2,1) : (inhi=1,2);
	    for (i=1;i<=mpts;i++) {
		if (y[i]<y[ilo]) ilo=i;
		if (y[i]>y[ihi]) {
		    inhi=ihi;
		    ihi=i;
		} else if (y[i]>y[inhi])
		    if (i!=ihi)inhi=i;
	    }
	    rtol=2.0*fabs(y[ihi]-y[ilo])/(fabs(y[ihi])+fabs(y[ilo]));

	    /* TIME VARIABLES, DELETE IF NECESSARY */
	    time(tp);
	    strcpy(s, ctime(tp));
	    s[16]='\0';
	    printf(" %s %4d    %11.4f %11.4f %11.5f  ",s+11,*nfunk,y[ihi],y[ilo],rtol);
	    for (j = 1; j <= ndim; j++)
	        printf("%8.4f ",p[ilo][j]);
	    printf("\n");

	    if (rtol<ftol)break;
	    if (*nfunk>=NMAX) nerror("Too may iterations in AMOEBA.");
	    ytry=amotry(p,y,psum,ndim,funk,ihi,nfunk,-ALPHA);
	    if (ytry<=y[ilo]) {
		ytry=amotry(p,y,psum,ndim,funk,ihi,nfunk,GAMMA);
	    } else if (ytry>=y[inhi]) {
		ysave=y[ihi];
		ytry=amotry(p,y,psum,ndim,funk,ihi,nfunk,BETA);
		if (ytry>=ysave) {
		    for (i=1;i<=mpts;i++) {
			if (i!=ilo) {
			    for (j=1;j<=ndim;j++) {
				psum[j]=0.5*(p[i][j]+p[ilo][j]);
				p[i][j]=psum[j];
			    }
			    y[i]=(*funk)(psum);
			}
		    }
		    *nfunk+=ndim;
		    GET_PSUM
		}
	    }
	}
	free_vector(psum,1);
	return ilo;
}

/* ************************************************************************** */
double  amotry(double **p,double *y,double *psum,int ndim,
		double (*funk)(double *param),int ihi,int *nfunk,double fac)
{
	int j;
	double fac1,fac2,ytry,*ptry;

	ptry=vector(1,ndim);
	fac1=(1.0-fac)/ndim;
	fac2=fac1-fac;
	for (j=1;j<=ndim;j++) ptry[j]=psum[j]*fac1-p[ihi][j]*fac2;
	ytry=(*funk)(ptry);
	++(*nfunk);
	if (ytry<y[ihi]) {
	    y[ihi]=ytry;
	    for (j=1;j<=ndim;j++) {
		psum[j]+=ptry[j]-p[ihi][j];
		p[ihi][j]=ptry[j];
	    }
	}
	free_vector(ptry,1);
	return ytry;
}

/* ************************************************************************** */
double *vector(int nl, int nh)
{
	double *v;

	v=(double *) malloc((unsigned) (nh-nl+1)*sizeof(double));
	if (!v) nerror("Allocation failure in vector()");
	return v-nl;
}

/* ************************************************************************** */
void   nerror(char *error_text)
{
	fprintf(stderr,"Run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now exiting to UNIX...\n");
	exit(1);
}

/* ************************************************************************** */
void   free_vector(double *v,int nl)
{
	free((char*) (v+nl));
}

/* ************************************************************************** */
double **matrix(int nrl, int nrh, int ncl, int nch)
{
	int i;
	double **m;
	m=(double **) malloc((unsigned) (nrh-nrl+1)*sizeof(double*));
	if (!m) nerror("matrix: allocation 1 failure!\n");
	m-=nrl;
	for (i=nrl;i<=nrh;i++) {
	    m[i]=(double *) malloc((unsigned) (nch-ncl+1)*sizeof(double));
	    if (!m[i]) nerror("matrix: allocation 2 failure!\n");
	    m[i]-=ncl;
	}
	return m;
}
/* ************************************************************************** */
void free_matrix(double **m, int nrl, int nrh, int ncl, int nch)
{
	int i;
	for (i=nrh;i>=nrl;i--) free((char*) (m[i]+ncl));
	free((char*) (m+nrl));
}

