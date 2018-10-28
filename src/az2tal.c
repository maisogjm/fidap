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
            printf("az2tal: Need ANALYZE coordinates.\n");
            printf("\nYou can also include voxel dimensions and the ANALYZE coordinates\n");
            printf("of the origin (six extra arguments).\n");
            exit(1);
        }

        /* ****************************************************************** */
        /* DETERMINE TALAIRACH COORDINATES OF INTEREST.                       */

        Px = atof(argv[1]);
        Py = atof(argv[2]);
        Pz = atof(argv[3]);
        if (argc==4) {
            Tx=(Px-33.)*2.;
            Ty=(Py-53)*2.;
            Tz=(Pz-8.)*4.;
        }
        if (argc==10) {
            Tx=(Px-atof(argv[7]))*atof(argv[4]);
            Ty=(Py-atof(argv[8]))*atof(argv[5]);
            Tz=(Pz-atof(argv[9]))*atof(argv[6]);
        }

        printf("%g %g %g\n",Tx,Ty,Tz);
}
