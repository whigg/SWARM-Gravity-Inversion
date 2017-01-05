#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main (int argc, char *argv[])
{
    FILE *fpin, *fpout;
    int i, gps, flg;
    double acc[3], res[3];
    char line[400];

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
    if ( (fpin = fopen (argv[1],"r")) == NULL)
    {
        printf ("Cannot open acc file!\n");
        exit (0);
    }

    if ( (fpout = fopen (argv[2],"w")) == NULL)
    {
        printf ("Cannot open acc file!\n");
        exit (0);
    }


    i = 0;
    while (1)
    {
        if (fgets(line,400, fpin) ==NULL) break;
        sscanf (line, "%d%*s%lf%lf%lf%*f%*f%*f%lf%lf%lf%d",
            &gps, &acc[0], &acc[1], &acc[2], &res[0], &res[1], &res[2], &flg);
//        if (flg != 0)
//            continue;
//        if ((gps > 201787200 + 56200) && (gps < 201787200 + 56400))
//            continue;

        fprintf (fpout, "%10d %25.15e %25.15e %25.15e %25.15e %25.15e %25.15e %d \n",
            gps, acc[0], acc[1], acc[2], res[0], res[1], res[2], flg);

        i++;

    }


    fclose(fpin);
    fclose(fpout);

    return 0;

}
