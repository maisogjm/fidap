#include <stdio.h>
#include <stdlib.h>

/* ************************************************************************** */
main(int argc, char *argv[])
{
        float Tx,Ty,Tz;         /* Talairach coordinates of interest.         */
        float Px,Py,Pz;         /* Where in image coordinates the Talairach   */
                                /* coordinates [Tz,Ty,Tz] map to.             */

        /* ****************************************************************** */
        /* CHECK COMMAND LINE ARGUMENTS*/
        if ((argc!=4) && (argc!=10)) {
            printf("tal2az: Need Talairach coordinates.\n");
            printf("\nYou can also include voxel dimensions and the ANALYZE coordinates\n");
            printf("of the origin (six extra arguments).\n");
            exit(1);
        }

        /* ****************************************************************** */
        /* DETERMINE TALAIRACH COORDINATES OF INTEREST.                       */

        Tx = atof(argv[1]);
        Ty = atof(argv[2]);
        Tz = atof(argv[3]);
        if (argc==4) {
            Px=(Tx/2.)+32.;
            Py=(Ty/2.)+52.;
            Pz=(Tz/4.)+7.;
        }
        if (argc==10) {
            Px=(Tx/atof(argv[4]))+atof(argv[7])-1;
            Py=(Ty/atof(argv[5]))+atof(argv[8])-1;
            Pz=(Tz/atof(argv[6]))+atof(argv[9])-1;
        }


        printf("%g %g %g\n",Px+1,Py+1,Pz+1);
}
