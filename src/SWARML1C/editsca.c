#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main (int argc, char *argv[])
{
    FILE *fpin, *fpout;
    int n, i, gps2, gps1;
    double sca0[4] = {0}, sca2[4]={0}, sca1[4]={0}, diff, diffnew;
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


    gps1=0;gps2=0;
    i = 0;
    while (1)
    {
        if (fgets(line,400, fpin) ==NULL) break;
        sscanf (line, "%d%*s%*d%lf%lf%lf%lf",
            &gps2, &sca0[0], &sca0[1], &sca0[2], &sca0[3]);
        diff = fabs((sca0[0] - sca1[0])/(double)(gps2-gps1));
        if (diff > 0.001)
        {
            for (n = 0; n < 4; n++)
                sca2[n] = - sca0[n];
//            printf ("diff = %f\n", diff);
        }
        else 
        {
            for (n = 0; n < 4; n++)
                sca2[n] = sca0[n];
        }

        diffnew = fabs((sca2[0] - sca1[0])/(double)(gps2-gps1));
        fprintf (fpout, "%10d %20.15lf %20.15lf %20.15lf %20.15lf %f %20.15lf %20.15lf %20.15lf %20.15lf %f\n",
            gps2, sca2[0], sca2[1], sca2[2], sca2[3], diffnew, sca0[0], sca0[1], sca0[2], sca0[3], diff);

        gps1 = gps2;
        for (n = 0; n < 4; n++)
            sca1[n] = sca2[n];

        i++;

    }


    fclose(fpin);
    fclose(fpout);

    return 0;

}
