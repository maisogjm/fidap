#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "analyze6.h"

unsigned char *unsigned_char_malloc(int num_items);
int unsigned_char_fwrite(char *filename, unsigned char *array, int num_items);


/* ************************************************************************** */
main(int argc,char *argv[])
/*        Reads in an 8-bit and a 24-bit image volume, and overlays the 24-bit*/
/*        image volume onto the 8-bit image volume, producing 24-bit output.  */
{
        struct analyze_struct
                img_info;   /* Header information variable.                   */

        unsigned char *output24bit,background;
        FILE *fOUTPUT_IMAGE;
        int i,num_vox;
        int xdim,ydim;
        int x,y;
        int scale_type;           /* Flag indicating what sort of scaling   */
                                  /* is to be done.                         */
        float X[3],Y[3];          /* XY coordinates of the three vertices   */
                                  /* of the triangle.  These define three   */
                                  /* lines.                                 */

                                  /* LINEAR EQUATIONS.                      */
        float yint1,yint2,yint3;  /* Y-intercepts of the three lines.       */
        float slope1,slope2,slope3;/* Slopes of the three lines.            */

                                  /* VARIABLES USED TO DETERMINE WHETHER A  */
                                  /* POINT IS WITHIN THE TRIANGLE.          */
        float delta1,delta2,delta3,/* Difference between vertex and its     */
                                  /* projection on the line defined by the  */
                                  /* opposite side of the triangle.         */
                sign_pt1,         /* Difference between a point and its     */
                sign_pt2,         /* projection on the line defined by one  */
                sign_pt3;         /* of the sides of the triangle.          */

                                  /* DISTANCE CALCULATIONS.                 */
        float x4,y4;              /* Coordinates of the perpendicular       */
                                  /* projection of a point onto a line.     */
        float vdist1,vdist2,vdist3,/* Distances of the three vertices from  */
                                  /* the opposite side of the triangle.     */
              dist1,dist2,dist3,  /* Distances of a point from each side.   */
              sum_dists,          /* dist1+dist2+dist3                      */
              max_dist;           /* Maximum of dist1, dist2, and dist3.    */
        float m,b;                /* Slope and y-intercept used in finding  */
                                  /* projection (x4,y4) of a point onto a   */
                                  /* line defined by one of the sides of    */
                                  /* the triangle.                          */
        float denominator;        /* Denominator.                           */

                                  /* INDEXING INTO 24-BIT IMAGE VOLUME.     */
        int vox_index24bit_red,   /* Index position of the red component.   */
            vox_index24bit_green, /* Index position of the green component. */
            vox_index24bit_blue;  /* Index position of the blue component.  */

        int write_status;


        /* **************************************************************** */
        /* CHECK COMMAND LINE ARGUMENTS.                                    */

        if (argc != 12) {
            printf("make_triangle_palette: need xdim, ydim, x1, y1, x2, y2, x3, y3,\n");
            printf("                    background intensity value,\n");
            printf("                    24-bit output image filename,\n");
            printf("                    and scale type.\n\n");

            printf("(x1,y1), (x2,y2), and (x3,y3) are the coordinates of the\n");
            printf("vertices of the triangle.\n\n");
            printf("Scale types:\n");
            printf("  (1) mean(Vr, Vg, Vb)\n");
            printf("  (2) max(Vr, Vg, Vb)\n");
            printf("  (3) sqrt((Vr*Vr)+(Vg*Vg)+(Vb*Vb))\n");
            return 1;
        }

        /* **************************************************************** */
        /* READ IN COMMAND LINE ARGUMENTS.                                  */

        xdim = atoi(argv[1]);
        ydim = atoi(argv[2]);
        num_vox = xdim*ydim;
        printf("make_triangle_palette: Number of voxels = %d\n",num_vox);
        scale_type = atoi(argv[11]);

        X[0] = atof(argv[3]);
        Y[0] = atof(argv[4]);
        X[1] = atof(argv[5]);
        Y[1] = atof(argv[6]);
        X[2] = atof(argv[7]);
        Y[2] = atof(argv[8]);

        printf("Point #1 = (%g,%g)\n",X[0],Y[0]);
        printf("Point #2 = (%g,%g)\n",X[1],Y[1]);
        printf("Point #3 = (%g,%g)\n",X[2],Y[2]);
        background = (unsigned char)atoi(argv[9]);


        /* **************************************************************** */
        /* ASSIGN MEMORY.                                                   */

        printf("make_triangle_palette: allocating memory...\n");
        output24bit = unsigned_char_malloc(num_vox*3);


        /* **************************************************************** */
        /* CALCULATE SLOPES AND Y-INTERCEPTS OF THREE LINES OF TRIANGLE.    */
        /* ALSO CALCULATE delta1, delta2, delta3, WHICH WILL BE USED IN     */
        /* DETERMINING WHETHER A PARTICULAR POINT LIES WITHIN THE TRIANGLE  */
        /* OR NOT.                                                          */

        printf("make_triangle_palette: calculating slopes and intercepts...\n");
        slope1 = (Y[1]-Y[0])/(X[1]-X[0]);
        yint1  = Y[1]-(X[1]*slope1);
        delta1  = Y[2]-(X[2]*slope1+yint1);
        printf("Line 1 slope, intercept, and delta = %f, %f, and %f\n",slope1,yint1,delta1);

        slope2 = (Y[2]-Y[0])/(X[2]-X[0]);
        yint2  = Y[2]-(X[2]*slope2);
        delta2  = Y[1]-(X[1]*slope2+yint2);
        printf("Line 2 slope, intercept, and delta = %f, %f, and %f\n",slope2,yint2,delta2);

        slope3 = (Y[1]-Y[2])/(X[1]-X[2]);
        yint3  = Y[2]-(X[2]*slope3);
        delta3  = Y[0]-(X[0]*slope3+yint3);
        printf("Line 3 slope, intercept, and delta = %f, %f, and %f\n",slope3,yint3,delta3);


        /* **************************************************************** */
        /* CALCULATE DISTANCES OF EACH VERTEX FROM OPPOSITE SIDE.           */

        if (slope1!=0.) {
            m=-1./slope1;
            b=Y[2]-(X[2]*m);
            x4=(b-yint1)/(slope1-m);
            y4=x4*slope1+yint1;
            x4-=X[2];
            y4-=Y[2];
            vdist1=(float)sqrt((double)((x4*x4)+(y4*y4)));
        }
        else
            vdist1=Y[2]-Y[0];
        printf("vdist1 = %f\n",vdist1);

        if (slope2!=0.) {
            m=-1./slope2;
            b=Y[1]-(X[1]*m);
            x4=(b-yint2)/(slope2-m);
            y4=x4*slope2+yint2;
            x4-=X[1];
            y4-=Y[1];
            vdist2=(float)sqrt((double)((x4*x4)+(y4*y4)));
        }
        else
            vdist1=Y[1]-Y[0];
        printf("vdist2 = %f\n",vdist2);

        if (slope3!=0.) {
            m=-1./slope3;
            b=Y[0]-(X[0]*m);
            x4=(b-yint3)/(slope3-m);
            y4=x4*slope3+yint3;
            x4-=X[0];
            y4-=Y[0];
            vdist3=(float)sqrt((double)((x4*x4)+(y4*y4)));
        }
        else
            vdist3=Y[0]-Y[1];
        printf("vdist3 = %f\n",vdist3);


        /* **************************************************************** */
        /* MAKE PALETTE.                                                    */

        printf("make_triangle_palette: looping over voxels...\n");
        for (y=0; y<ydim; y++)
            for (x=0; x<xdim; x++) {

                vox_index24bit_red   = ((0*ydim)+y)*xdim+x;
                vox_index24bit_green = ((1*ydim)+y)*xdim+x;
                vox_index24bit_blue  = ((2*ydim)+y)*xdim+x;

                /* ******************************************************** */
                /* CHECK WHETHER CURRENT COORDINATES LIE WITHIN TRIANGLE.   */
                /* IF SO, CALCULATE DISTANCE FROM SIDES OF TRIANGLE, AND    */
                /* ASSIGN WEIGHTED INTENSITY BASED ON THIS DISTANCE.        */

                sign_pt1  = (float)y-((float)x*slope1+yint1);
                sign_pt2  = (float)y-((float)x*slope2+yint2);
                sign_pt3  = (float)y-((float)x*slope3+yint3);

                if ((sign_pt1*delta1>=0.)
                 && (sign_pt2*delta2>=0.)
                 && (sign_pt3*delta3>=0.)) {

                    /* **************************************************** */
                    /* DETERMINE DISTANCE FROM EACH SIDE OF TRIANGLE.       */

                    if (slope1!=0.) {
                        m=-1./slope1;
                        b=y-(x*m);
                        x4=(b-yint1)/(slope1-m);
                        y4=x4*slope1+yint1;
                        x4-=x;
                        y4-=y;
                        dist1=(float)sqrt((double)((x4*x4)+(y4*y4)));
                        dist1=vdist1-dist1;
                    }
                    else
                        dist1=y-Y[0];

                    if (slope2!=0.) {
                        m=-1./slope2;
                        b=y-(x*m);
                        x4=(b-yint2)/(slope2-m);
                        y4=x4*slope2+yint2;
                        x4-=x;
                        y4-=y;
                        dist2=(float)sqrt((double)((x4*x4)+(y4*y4)));
/*                        dist2=vdist2-dist2; */
                    }
                    else
                        dist2=y-Y[0];

                    if (slope3!=0.) {
                        m=-1./slope3;
                        b=y-(x*m);
                        x4=(b-yint3)/(slope3-m);
                        y4=x4*slope3+yint3;
                        x4-=x;
                        y4-=y;
                        dist3=(float)sqrt((double)((x4*x4)+(y4*y4)));
/*                        dist3=vdist3-dist3; */
                    }
                    else
                        dist3=y-Y[1];

                    if (scale_type == 1) 
                        denominator = dist1+dist2+dist3;
		    else if (scale_type == 2) {
                        denominator=dist1;
                        denominator=(dist2>denominator)?dist2:denominator;
                        denominator=(dist3>denominator)?dist3:denominator;
                    }
		    else if (scale_type == 3)
                        denominator=(float)sqrt((dist1*dist1)+(dist2*dist2)+(dist3*dist3));

                    output24bit[vox_index24bit_red]   = (unsigned char) floor((255.*dist1/denominator)+0.5);
                    output24bit[vox_index24bit_green] = (unsigned char) floor((255.*dist2/denominator)+0.5);
                    output24bit[vox_index24bit_blue]  = (unsigned char) floor((255.*dist3/denominator)+0.5);
                }

                /* ******************************************************** */
                /* IF NOT, GIVE THE BACKGROUND VALUE TO CURRENT COORDINATES.*/

                else {
                    output24bit[vox_index24bit_red]   = background;
                    output24bit[vox_index24bit_green] = background;
                    output24bit[vox_index24bit_blue]  = background;
                }
            }

        
        /* **************************************************************** */
        /* WRITE TO OUTPUT FILE.                                            */

        unsigned_char_fwrite(argv[10], output24bit, num_vox*3);

        free(output24bit);
}

