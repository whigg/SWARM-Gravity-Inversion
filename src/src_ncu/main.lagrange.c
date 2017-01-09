#include <stdio.h>
#include <stdlib.h>
#include <math.h>




double lagrange (double *y, int dim_y, int dim_x, double t, double *z)
{
    int i, j, k, m, dim, order = 2;
    double s;

    i = 0;
    while ((y[i * dim_x] < t) && (i < dim_y))
        i = i + 1;
    k = i - order;
    if (k < 0)
        k = 0;
    m = i + order - 1;
    if (m > dim_y - 1)
        m = dim_y - 1;

    for (dim = 0; dim < dim_x - 1; dim++)
    {
        z[dim] = 0;
    }

    for (i = k; i <= m; i++)
    {
        s = 1.0;
        for (j = k; j <= m; j++)
        {
            if (j != i)
            {
                s = s * (t - y[j * dim_x]) / (y[i * dim_x] - y[j * dim_x]);
            }
        }
        for (dim = 0; dim < dim_x - 1; dim++)
        {
            z[dim] = z[dim] + s * y[i * dim_x + dim + 1];
        }
    }
    return 0;
}




main (int argc, char *argv[])
{
    FILE *fpin, *fpout;
    int n, i, dim;
    double *sca_eph, gps0,gps1, gps2, gpsi, sca0[4] = {0}, sca2[4]={0}, sca1[4]={0}, diff;
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

    dim = 0;
    while (1)
    {
        if (fgets(line,400, fpin) ==NULL) break;
        dim ++;
    }
    rewind(fpin);

    sca_eph = (double *) calloc ( dim * 5, sizeof(double));

    gps1=0;gps2=0;
    i = 0;
    for (i=0;i<dim;i++)
    {
        if (fgets(line,400, fpin) ==NULL) break;
        sscanf (line, "%lf%*s%*d%lf%lf%lf%lf",
            &sca_eph[i * 5], &sca_eph[i * 5 + 1],&sca_eph[i * 5 + 2],&sca_eph[i * 5 + 3], &sca_eph[i * 5 + 4]);
        if (i == 0)
            gps0 = sca_eph[i * 5];

        gps2 = sca_eph[i * 5];
        for (n = 0; n < 4; n++)
            sca0[n] = sca_eph[i * 5 + 1 + n];

        diff = fabs((sca0[0] - sca1[0])/(gps2-gps1));
        if (diff > 0.01)
        {
            for (n = 0; n < 4; n++)
                sca_eph[i * 5 + 1 + n] = - sca0[n];
        }
        else 
        {
            for (n = 0; n < 4; n++)
                sca_eph[i * 5 + 1 + n] = sca0[n];
        }


        gps1 = gps2;
        for (n = 0; n < 4; n++)
            sca1[n] = sca_eph[i * 5 + 1 + n];

    }


    for (gpsi = gps0; gpsi < gps0 + 86400; gpsi = gpsi + 5)
    {
        lagrange (sca_eph, dim, 5, gpsi, sca2);


        fprintf (fpout, "%10d %20.15lf %20.15lf %20.15lf %20.15lf\n",
            (int)gpsi, sca2[0], sca2[1], sca2[2], sca2[3]);
    }

    fclose(fpin);
    fclose(fpout);


}
