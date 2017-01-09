/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/*

    Copyright (c) 2012 Kun Shang (shang.34@osu.edu) All Right Reserved 

    NOT FULLY VALIDATED -- SUBJECT TO CORRECTION

*/
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/


/*#define MSDOS */
//#define LINUX



#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "novas.h"
#include "coord.h"
#include "grvts.h"
#include "numrs.h"

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
int main (int argc, char *argv[])
{
    FILE *fp_stdin, *fp_aotdat, *fp_aotbin;
    int i, n, year, month, day, step, step_aod, minor = 0, dim, dim_aod, 
        nmax, nmax_ots, nmax_ats1, nmax_ats2, nmax_aod, nmax_pts, mjds, mjde, 
        par, par_ots, par_aod, par_ats1, par_ats2, par_pts,
        ots_nper, ats1_nper, ats2_nper;
    char line[400], card[200], f_eop[200], f_adm[200], f_ots[200], f_opt[200],
        f_ats1[200],f_aod[200], f_ats2[200], f_out[200], stdname[100];
    double tt, tts, tte, jd0, jdtt[2], *coe_aot, *coe_ots, *coe_pts, 
        *coe_aod, *coe_ats1, *coe_ats2, *aot_eph, *aod_eph;

    InfStruct info;

    OTSperturb *ots_per=NULL, *ats1_per=NULL, *ats2_per=NULL;


/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
 
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
        fgets (line, 400, fp_stdin);
        sscanf (line, "%s", card);	
        if (strcmp (card, "EPOCH") ==0)
            sscanf (line, "%*s%d%d%d",  &year, &month, &day);
        if (strcmp (card, "STEP") ==0)
            sscanf (line, "%*s%d",      &step);
        if (strcmp (card, "NMAX") ==0)
            sscanf (line, "%*s%d",      &nmax);
        if (strcmp (card, "OTIDE") ==0)
            sscanf (line, "%*s%s%d",    f_ots, &nmax_ots);
        if (strcmp (card, "MINOR") ==0)
            sscanf (line, "%*s%d%s",    &minor, f_adm);
        if (strcmp (card, "ATIDE") ==0)
            sscanf (line, "%*s%s%d%s%d",f_ats1,&nmax_ats1,f_ats2,&nmax_ats2);
        if (strcmp (card, "AOD1B") ==0)
            sscanf (line, "%*s%s%d",    f_aod, &nmax_aod);
        if (strcmp (card, "PTIDE") ==0)
            sscanf (line, "%*s%s%d",    f_opt, &nmax_pts);
        if (strcmp (card, "EOP") ==0)
            sscanf (line, "%*s%s",      f_eop);
        if (strcmp (card, "OUTPUT") ==0)
            sscanf (line, "%*s%s",      f_out);
    }
    
    fclose(fp_stdin);
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/

    if (nmax_ots>nmax || nmax_aod>nmax || nmax_ats1>nmax || nmax_ats2>nmax)
    {
        printf ("nmax = %d is too small\n", nmax);
        exit (0);
    }

    if((fp_aotdat=fopen("aoteph.dat","w"))==NULL)
    {
        printf("Cannot write otide.dat!\n");
        exit(0);
    }
    if((fp_aotbin=fopen(f_out,"wb"))==NULL)
    {
        printf("Cannot write otide.dat!\n");
        exit(0);
    }


    jd0 = julian_date (year, month, day,0);

    tts = gps2tt (0);
    tte = tts + 86400.0;

    mjds = (int)(jd0 - 2400000.5);
    mjde = mjds + 1;
    eop_open (f_eop, mjds, mjde);


// the final matrix for output eph

//    nmax = nmax_ots;
    par = (nmax + 1) * (nmax + 1);
    dim = (int)(86400/step) + 1;
    coe_aot = (double *) calloc ( par, sizeof(double));
    aot_eph = (double *) calloc ((par + 1) * dim , sizeof(double));


// open ots file

    par_ots = (nmax_ots + 1) * (nmax_ots + 1);
    coe_ots = (double *) calloc ( par_ots, sizeof(double));
    ots_per = ts_open (f_ots, nmax_ots, &ots_nper);

    if (minor == 1)
        adm_open (f_adm);


// open aod file

    step_aod = 21600;
    dim_aod = (int)(86400/step_aod);
    par_aod = (nmax_aod + 1) * (nmax_aod + 1);
    coe_aod = (double *) calloc (par_aod, sizeof(double));
    aod_eph = (double *) calloc ( (par_aod + 1) * dim_aod , sizeof(double));

    aod_open (f_aod, nmax_aod, aod_eph);


// open ats file
//
    par_ats1 = (nmax_ats1 + 1) * (nmax_ats1 + 1);
    coe_ats1 = (double *) calloc (par_ats1, sizeof(double));
    ats1_per = ts_open (f_ats1, nmax_ats1, &ats1_nper);

    par_ats2 = (nmax_ats2 + 1) * (nmax_ats2 + 1);
    coe_ats2 = (double *) calloc (par_ats2, sizeof(double));
    ats2_per = ts_open (f_ats2, nmax_ats2, &ats2_nper);

// open pts file
//
    opt_open (f_opt, nmax_pts);
    par_pts = (nmax_pts + 1) * (nmax_pts + 1);
    coe_pts = (double *) calloc (par_pts, sizeof(double));



// remove S2 air tide from AOD1B data (AOD - ATS2)
//
    jdtt[0] = jd0;        
    tt = tts;
    for (i = 0; i < dim_aod; i++)
    {
        jdtt[1] = tt / 86400.0;   
        getinfo(jdtt, 2, &info);

        ts_read_minor(info.jdt, info.gmst, nmax_ats2, 0, ats2_nper, ats2_per, coe_ats2); 
        
        addsubcs (&aod_eph[i * (par_aod + 1) + 1], coe_ats2, nmax_aod, nmax_ats2, -1);

        tt = tt + step_aod;
    }
    

// AOT := OTS + ATS1 + ATS2 + (AOD - ATS2)
//
    jdtt[0] = jd0;        
    tt = tts;
    for (i = 0; i < dim; i++)
    {
        jdtt[1] = tt / 86400.0;   
        getinfo(jdtt, 2, &info);

        ts_read_minor(info.jdt, info.gmst, nmax_ots, minor, ots_nper, ots_per, coe_ots); 
        ts_read_minor(info.jdt, info.gmst, nmax_ats1, 0,  ats1_nper, ats1_per, coe_ats1); 
        ts_read_minor(info.jdt, info.gmst, nmax_ats2, 0,  ats2_nper, ats2_per, coe_ats2); 
        pt_read (info.mjdutc, info.xp, info.yp, nmax_pts, coe_pts); 


        lgr_order (aod_eph, 4, par_aod + 1 , tt, coe_aod, 1);

        memset(coe_aot,0, par * sizeof(double));

        addsubcs (coe_aot, coe_ots,  nmax, nmax_ots,  1);
        addsubcs (coe_aot, coe_ats1, nmax, nmax_ats1, 1);
        addsubcs (coe_aot, coe_ats2, nmax, nmax_ats2, 1);
        addsubcs (coe_aot, coe_aod,  nmax, nmax_aod,  1);
        addsubcs (coe_aot, coe_pts,  nmax, nmax_pts,  1);

        zero0zero1 (coe_aot, nmax);

        aot_eph[i * (par + 1)] = tt;
        for (n = 0; n < par; n++)
            aot_eph[i * (par + 1) + 1 + n] = coe_aot[n];
        tt = tt + step;
    }
    
        
// output
/*
    for (i = 0; i < dim; i++)
    {
        for (n = 0; n < par + 1; n++)
        {
            fprintf (fp_aotdat, "%e\t", 
                aot_eph[i * (par + 1) + n]);
        }
        fprintf (fp_aotdat, "\n");
    }
*/

    fwrite (aot_eph, sizeof(double) * (par + 1) * dim, 1, fp_aotbin);

    fclose(fp_aotdat);
    fclose(fp_aotbin);

    free (aot_eph);
    free (aod_eph);

    free (coe_aot);
    free (coe_ots);
    free (coe_aod);
    free (coe_ats1);
    free (coe_ats2);
    free (coe_pts);

    free (ots_per);
    free (ats1_per);
    free (ats2_per);

    eop_close();
    adm_close ();
    opt_close ();


    return 0;

}
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/

