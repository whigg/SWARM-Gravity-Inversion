
#include "coord.h"
#include "numrs.h"

#define MAXLINE 300


/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
int main (int argc, char *argv[])
{
    FILE *fp_stdin, *fpin_a, *fpin_b, *fpout_a, *fpout_b;	
    int year, month, day, days, mjds, mjde, gps_e, gps_i, i, n, 
        DT, NDATA, GPS_S;
    double tt, tjd[2], xi_a[3], vi_a[3], xi_b[3], vi_b[3], llh_a[3], llh_b[3],
        viv_a[3], viv_b[3], vix_a[3], vix_b[3], x_a[3], x_b[3], 
        v_a[3], v_b[3], ex_a[9], ev_a[9], ex_b[9], ev_b[9], tmp[9],
        exi_a[9], evi_a[9], exi_b[9], evi_b[9], eviv_a[9], evix_a[9], 
        eviv_b[9], evix_b[9], JD0;
    char line[MAXLINE], line_a[MAXLINE], line_b[MAXLINE], 
        card[20], f_eop[400], f_gnv1b_a[400], f_gnv1b_b[400], 
        f_gnv1c_a[400], f_gnv1c_b[400], stdname[100];    
    
    InfStruct info;


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

    while (feof(fp_stdin) == 0)
    {
        card[1] = '\0';
        fgets (line, MAXLINE, fp_stdin);
        sscanf (line, "%s", card);	

        if (strcmp (card, "EPOCH") ==0)
            sscanf (line, "%*s%d%d%d", &year, &month, &day);
        if (strcmp (card,"GNV1B") ==0)	
            sscanf (line, "%*s%s%s", f_gnv1b_a, f_gnv1b_b);	
        if (strcmp (card,"GNV1C") ==0)	
            sscanf (line, "%*s%s%s", f_gnv1c_a, f_gnv1c_b);	
        if (strcmp (card,"EOP") ==0)	
            sscanf (line, "%*s%s", f_eop);	
        if (strcmp (card, "DAYS") ==0)
            sscanf (line, "%*s%d", &days);
        if (strcmp (card, "DT") ==0)	
            sscanf (line, "%*s%d", &DT);

    }



/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
    if ( (fpin_a = fopen (f_gnv1b_a,"r")) == NULL)
    {
        printf ("Cannot open fpin_a file!\n");
        exit (0);
    }
    if ( (fpin_b = fopen (f_gnv1b_b,"r")) == NULL)
    {
        printf ("Cannot open fpin_b file!\n");
        exit (0);
    }

    if ( (fpout_a = fopen (f_gnv1c_a,"w")) == NULL)
    {
        printf ("Cannot open fpout_a file!\n");
        exit (0);
    }
    if ( (fpout_b = fopen (f_gnv1c_b,"w")) == NULL)
    {
        printf ("Cannot open fpout_b file!\n");
        exit (0);
    }


/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/

    JD0 = julian_date ((short)year,(short)month,(short)day,0);

    NDATA = days * 86400 / DT;


    mjds = (int)(JD0 - 2400000.5);
    mjde = mjds + 1;
    eop_open (f_eop, mjds, mjde);


    GPS_S = (int)((JD0 - T0) * 86400 + 0.5);
    gps_e = GPS_S + NDATA * DT;


    while (1)
    {
        if (fgets(line_a,300, fpin_a) ==NULL) break;
        if (fgets(line_b,300, fpin_b) ==NULL) break;

        sscanf (line_a, "%d", &gps_i);

        if ((gps_i - GPS_S) % DT == 0 
            && (gps_i - GPS_S) / DT < NDATA
            && (gps_i - GPS_S) / DT >= 0)
        {
            i = (gps_i - GPS_S) / DT;

            for (n = 0; n < 9; n ++)
            {
                ex_a[n] = 0;
                ev_a[n] = 0;
                ex_b[n] = 0;
                ev_b[n] = 0;
            }

            sscanf (line_a, "%*f%*s%*s%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf",
                &x_a[0], &x_a[1], &x_a[2], &ex_a[0], &ex_a[4], &ex_a[8],
                &v_a[0], &v_a[1], &v_a[2], &ev_a[0], &ev_a[4], &ev_a[8]);
            sscanf (line_b, "%*f%*s%*s%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf",
                &x_b[0], &x_b[1], &x_b[2], &ex_b[0], &ex_b[4], &ex_b[8],
                &v_b[0], &v_b[1], &v_b[2], &ev_b[0], &ev_b[4], &ev_b[8]);

            for (n = 0; n < 9; n ++)
            {
                ex_a[n] = ex_a[n] * ex_a[n];
                ev_a[n] = ev_a[n] * ev_a[n];
                ex_b[n] = ex_b[n] * ex_b[n];
                ev_b[n] = ev_b[n] * ev_b[n];
            }


            tt = gps_i - GPS_S + 19 + 32.184;

            tjd[0] = JD0;    tjd[1] = tt / 86400.0;

            getinfo (tjd, 2, &info);

            brmul(info.c_ei, x_a, 3,3,1, xi_a);
            brmul(info.c_ei, x_b, 3,3,1, xi_b);
            brmul(info.c_ei, v_a, 3,3,1, viv_a);
            brmul(info.c_ei, v_b, 3,3,1, viv_b);
            brmul(info.c_eidot, x_a, 3,3,1, vix_a);
            brmul(info.c_eidot, x_b, 3,3,1, vix_b);

            for (n = 0; n < 3; n ++)
            {
                vi_a[n] = viv_a[n] + vix_a[n];
                vi_b[n] = viv_b[n] + vix_b[n];
            }

            xyz2llr(x_a, llh_a);
            xyz2llr(x_b, llh_b);


            brmul(info.c_ei, ex_a, 3,3,3, tmp);
            brmul(tmp, info.c_ie,  3,3,3, exi_a);
            brmul(info.c_ei, ex_b, 3,3,3, tmp);
            brmul(tmp, info.c_ie,  3,3,3, exi_b);

            brmul(info.c_ei, ev_a, 3,3,3, tmp);
            brmul(tmp, info.c_ie,  3,3,3, eviv_a);
            brmul(info.c_ei, ev_b, 3,3,3, tmp);
            brmul(tmp, info.c_ie,  3,3,3, eviv_b);

            brmul(info.c_eidot, ex_a, 3,3,3, tmp);
            brmul(tmp, info.c_iedot,  3,3,3, evix_a);
            brmul(info.c_eidot, ex_b, 3,3,3, tmp);
            brmul(tmp, info.c_iedot,  3,3,3, evix_b);

            for (n = 0; n < 9; n ++)
            {
                evi_a[n] = eviv_a[n] + evix_a[n];
                evi_b[n] = eviv_b[n] + evix_b[n];
            }

            for (n = 0; n < 9; n ++)
            {
                exi_a[n] = sqrt(exi_a[n]);
                evi_a[n] = sqrt(evi_a[n]);
                exi_b[n] = sqrt(exi_b[n]);
                evi_b[n] = sqrt(evi_b[n]);
                ex_a[n] = sqrt(ex_a[n]);
                ev_a[n] = sqrt(ev_a[n]);
                ex_b[n] = sqrt(ex_b[n]);
                ev_b[n] = sqrt(ev_b[n]);
            }


            fprintf(fpout_a, 
                "%d A I %20.10f %20.10f %20.10f %20.14f %20.14f %20.14f %20.10f %20.14f %20.14f\t",
                gps_i, xi_a[0], xi_a[1], xi_a[2], vi_a[0], vi_a[1], vi_a[2], 
                llh_a[2], llh_a[0], llh_a[1]);
            fprintf(fpout_a, 
                "%10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e\n",
                exi_a[0], exi_a[4], exi_a[8], evi_a[0], evi_a[4], evi_a[8], 
                ex_a[0], ex_a[4], ex_a[8], ev_a[0], ev_a[4], ev_a[8]);
            fprintf(fpout_b, 
                "%d B I %20.10f %20.10f %20.10f %20.14f %20.14f %20.14f %20.10f %20.14f %20.14f\t",
                gps_i, xi_b[0], xi_b[1], xi_b[2], vi_b[0], vi_b[1], vi_b[2], 
                llh_b[2], llh_b[0], llh_b[1]);
            fprintf(fpout_b, 
                "%10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e\n",
                exi_b[0], exi_b[4], exi_b[8], evi_b[0], evi_b[4], evi_b[8], 
                ex_b[0], ex_b[4], ex_b[8], ev_b[0], ev_b[4], ev_b[8]);

            i++;
        }
    }

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/

  
    fclose(fpin_a);
    fclose(fpin_b);
    fclose(fpout_a);
    fclose(fpout_b);
    


    return 0;

}
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
 



