/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/*

    Generate global geoid/gravity/... to see stripes from SHCs
    
    Version:    14 DEC 2011

    Copyright (c) 2010 SHANG Kun (shang.34@osu.edu) All Right Reserved 

    NOT FULLY VALIDATED -- SUBJECT TO CORRECTION

*/
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/

       
 
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "novas.h"
#include "grvts.h"


#define MAXLINE 300


/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
int main (int argc, char *argv[])
{
    FILE *fp_stdin, *fpout;	
    short int year, month, day, year0, month0, day0;
    int nmax, solvegrv, l, n, m, k, ind, zerotide = 0;
    double hour0, hour, t, t0, *coefa, zero=0;
    char line[MAXLINE], card[20], f_out[40], f_grv[40]={"\0"}, stdname[100];    
    
	if (argc == 1)
	{
		printf ("input file name: ");
		scanf ("%s", stdname);
	}
	else if (argc == 2)
	{
		sscanf (argv[1], "%s", stdname);
	}
	else
	{
        printf ("input argument error!\n");
        exit (0);
	}
    if ( (fp_stdin = fopen (stdname,"r")) == NULL)
    {
        printf ("Cannot open stdin file!\n");
        exit (0);
    }
//    printf ("\nreading input file \"%s\"...\n\n", stdname);

    while (feof(fp_stdin) == 0)
    {
        card[1] = '\0';
        fgets (line, MAXLINE, fp_stdin);
        sscanf (line, "%s", card);	
        
        if (strcmp (card,"INPUT") ==0)	
        {
            sscanf (line, "%*s%s", f_grv);	
        }
        if (strcmp (card, "NMAX") ==0)	
        {
            sscanf (line, "%*s%d", &nmax);
        }
        if (strcmp (card, "EPOCH") ==0)
        {
            sscanf (line, "%*s%hd%hd%hd%lf", &year, &month, &day, &hour);
        }
        if (strcmp (card, "EPOCH0") ==0)
        {
            sscanf (line, "%*s%hd%hd%hd%lf", &year0, &month0, &day0, &hour0);
        }
        if (strcmp (card, "OUTPUT") ==0)	
        {
            sscanf (line, "%*s%s", f_out);
        }
        if (strcmp (card, "ZEROTIDE") ==0)	
        {
            sscanf (line, "%*s%d", &zerotide);
        }
    }

    if ( (fpout = fopen (f_out,"w")) == NULL)
    {
        printf ("Cannot open fpout_a file!\n");
        exit (0);
    }

    t0 = julian_date (year0, month0, day0, hour0);

    t = julian_date (year, month, day, hour);


/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
    solvegrv = (nmax + 1) * (nmax + 1);

    coefa  = (double *) calloc (solvegrv, sizeof(double));

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/

//    opengrav (f_grv, coefa, GMA, nmax, 0);

//    openeigen (f_grv, coefa, nmax, (t-t0)/365.25);
    eig_open (f_grv, (t-t0)/365.25, nmax, nmax, coefa);
    unit0zero1 (coefa, nmax);

//    coefa[0] = 0; coefa[1] = 0; coefa[nmax + 1] = 0;
    
    if (zerotide == 1)
        coefa[2] = coefa[2] - 4.173576159e-9;
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/

    for (m = 0; m <= nmax; m ++)
    {
        l = nmax - m + 1;
        for (k = 0; k < l; k++)
        {
            n = k + m;
            if (n < 2) continue;
            if (m == 0)
            {
                fprintf(fpout, "%4d %4d %23.13e %23.13e \n",
                    n, m, coefa[k], zero);

            }
            else
            {
                ind = nmax + 1 + (2 * nmax - m + 2) * (m - 1);
                fprintf(fpout, "%4d %4d %23.13e %23.13e \n",
                    n, m, coefa[ind + n - m], coefa[ind + n - m + l]);
            }
        }
    }



            
    free(coefa);
    fclose(fpout);
    fclose(fp_stdin);


//    printf ("\nNormal end of ECHO!\n\npress any key to finish...\n");

    return 0;
    exit(0);

}
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
 



