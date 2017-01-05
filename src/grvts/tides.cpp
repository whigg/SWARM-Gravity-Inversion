
/*

------------------------------------------------------------------------

    Purpose: read and compute tide (ocean, atm., AOD tides, etc)
    Notes: 
    Programmer: Kun Shang @ 4.29.2014
    Functions:


        int opt_close ();

        int opt_open (char *f_pt, int nmax);


        int mtpole (double mjd, double xp, double yp, 
                double *m1, double *m2);

        int pt_read (double mjd, double xp, double yp, int nmax, 
                double *coef);

        int ts_close (OTSperturb *per);

        OTSperturb *ts_open (char *ts_name, int nmax, int *nper);
    
        int ts_read_minor(double jdt, double gmst, int nmax, int minor, 
                int nper, OTSperturb *per, double *coef);

        int adm_close ();

        int adm_open (char *adm_name);

        double arg2theta (double jdt, double gmst, int *n, 
                double *ang);

        double aod_open (char *file_aod, int nmax_aod, double *aod_eph);

        double stidecs_earth(InfStruct *info, double *stcs, 
                double c20pt, int anelastic, int freqdep, int poletide);

        double stfrqdep(double jdt, double gmst, 
                double *c20f, double *c21f, double *s21f, 
                double *c22f, double *s22f);


     Global variables:

        extern int NMAJ_OTS, NALL_OTS;
        extern double *OTADM;
        extern OTSconstit *ALL_OTS;



------------------------------------------------------------------------



*/




#ifndef _GRVTS_H_

#include "grvts.h"
    #define _GRVTS_H_
#endif




    int NMAJ_OTS = 18, NALL_OTS  = 256;
    double *OTADM, *OPTM1, *OPTM2;
    OTSconstit *ALL_OTS;


int opt_close ()
{
    free (OPTM1);
    free (OPTM2);
    return 0;
}


int opt_open (char *f_pt, int nmax)
{
    FILE *fp_grv;
    double cm1, sm1, cm2, sm2;
    int n, m, ic, is, l, ind;
    char string[200];

    if ((fp_grv = fopen (f_pt,"r")) == NULL)
    {
        printf ("Cannot open ocean pole tides file?\n");
        exit (0);
    }


    OPTM1 = (double *) calloc ( (nmax + 1) * (nmax + 1), sizeof(double));
    OPTM2 = (double *) calloc ( (nmax + 1) * (nmax + 1), sizeof(double));


    while (1)
    {
        if (fgets (string, 200, fp_grv) == NULL) break;
        sscanf (string, "%d%d%lf%lf%lf%lf", &n, &m, &cm1, &sm1, &cm2, &sm2); 
        if (n > nmax)
        { 
            continue;
        }
        else if (m == 0)
        {
            OPTM1[n] = cm1;
            OPTM2[n] = cm2;
        }
        else 
        {
            l = nmax - m + 1;
            ind = nmax + 1 + (2 * nmax - m + 2) * (m - 1);
            ic = ind + n - m;
            is = ind + n - m + l;
            OPTM1[ic] = cm1;
            OPTM1[is] = sm1;
            OPTM2[ic] = cm2;
            OPTM2[is] = sm2;
        }
    }

    fclose(fp_grv);
    return 0;
}




int mtpole (double mjd, double xp, double yp, 
        double *m1, double *m2)
{
    double t, xpm, ypm, mjd2000 = 51544.0, mjd2010 = 55197.0;

    t = (mjd - mjd2000)/365.25;

    if (mjd <= mjd2010)
    {
        xpm = 0.055974 + t * 1.8243e-3
                       + t * t * 0.18413e-3
                       + t * t * t * 0.007024e-3;
        ypm = 0.346346 + t * 1.7896e-3
                       - t * t * 0.10729e-3
                       - t * t * t * 0.000908e-3;
    }
    else 
    {
        xpm = 0.023513 + t * 7.6141e-3;
        ypm = 0.358891 - t * 0.6287e-3;
    }

    *m1 = xp - xpm;
    *m2 = - yp + ypm;

    return 0;
}







int pt_read (double mjd, double xp, double yp, int nmax, 
        double *coef)
{
    double m1, m2;
    int n, m, k, l, ind, ic, is;

    mtpole (mjd, xp, yp, &m1, &m2);
    

    for (m = 0; m <= nmax; m ++)
    {
        l = nmax - m + 1;
        for (k = 0; k < l; k++)
        {
            if (m==0)
            {
                n = k;
                coef[n] = OPTM1[n] * m1 + OPTM2[n] * m2;
            }
            else
            {
                ind = nmax + 1 + (2 * nmax - m + 2) * (m - 1);
                n = k + m;
                ic = ind + n - m;
                is = ind + n - m + l;
                coef[ic] = OPTM1[ic] * m1 + OPTM2[ic] * m2;
                coef[is] = OPTM1[is] * m1 + OPTM2[is] * m2;
            }
        }
    }

/*solid pole tide*/
    if (nmax >= 2)
    {
        n = 2; m = 1;
        l = nmax - m + 1;
        ind = nmax + 1 + (2 * nmax - m + 2) * (m - 1);

        ic = ind + n - m;
        is = ind + n - m + l;

        coef[ic] += -1.333e-9 * (m1 + 0.0115 * m2);
        coef[is] += -1.333e-9 * (m2 - 0.0115 * m1);
    }
/*solid pole tide*/

    return 0;
}






int ts_close (OTSperturb *per)
{
    free (per);
    return 0;
}


OTSperturb *ts_open (char *ts_name, int nmax, int *nper)
{
    FILE *fp_ot;
    int i, n;
    char string[200];

    OTSperturb *per;

    if ((fp_ot = fopen (ts_name,"r")) == NULL)
    {
        printf ("Cannot open tide file %s\n", ts_name);
        exit (0);
    }

    i = 0;
    while (1)
    {
        if (fgets (string, 200, fp_ot) == NULL) break;
        sscanf (string, "%*f%*s%d", &n);
        if (n > nmax) continue;
        i ++;
    }
    rewind(fp_ot);

    *nper = i;
//    printf ("NPER = %d\n", *nper);
    per = (OTSperturb *) calloc ( *nper, sizeof(OTSperturb));
    if ( per == NULL )    
    {
        printf("Malloc error!\n");
    }

    i = 0;
    while (1)
    {
        if (fgets (string, 200, fp_ot) == NULL) break;
        sscanf (string, "%*f%*s%d", &n);
        if (n > nmax) continue;
        sscanf (string, "%lf%s%d%d%lf%lf%lf%lf", 
            &per[i].ds, per[i].name, &per[i].n, &per[i].m, 
            &per[i].cp, &per[i].sp, &per[i].cm, &per[i].sm);    

        per[i].argn[0] = (int)(per[i].ds/100)%10;
        per[i].argn[1] = (int)(per[i].ds/10)%10 - 5;
        per[i].argn[2] = (int)(per[i].ds/1)%10 - 5;
        per[i].argn[3] = (int)(per[i].ds*10)%10 - 5;
        per[i].argn[4] = (int)(per[i].ds*100)%10 - 5;
        per[i].argn[5] = (int)(per[i].ds*1000)%10 - 5;
        per[i].fnm = 1e-11;

        i++;
    }

    fclose(fp_ot);

    return per;

}



int adm_close ()
{
    free (OTADM);
    free (ALL_OTS);
    return 0;
}



int adm_open (char *adm_name)
{
    FILE *fp_ad;    
    int i, j;

    double doodson[256] = {
    55.565,55.575,56.554,56.556,57.355,57.553,57.555,57.565,57.575,
    58.554,59.553,62.656,63.645,63.655,63.665,64.456,64.555,65.445,65.455,
    65.465,65.655,65.665,65.675,66.454,67.455,67.465,71.755,72.556,73.545,
    73.555,73.565,74.556,75.345,75.355,75.365,75.555,75.565,75.575,76.554,
    77.355,77.365,81.655,82.656,83.445,83.455,83.655,83.665,83.675,84.456,
    85.255,85.455,85.465,85.475,86.454,91.555,92.556,93.355,93.555,93.565,
    93.575,95.355,95.365,107.755,109.555,115.845,115.855,117.645,117.655,118.654,
    119.455,125.745,125.755,126.556,126.754,127.545,127.555,128.554,129.355,133.855,
    134.656,135.435,135.635,135.645,135.655,135.855,136.555,136.654,137.445,137.455,
    137.655,137.665,138.454,139.455,143.535,143.745,143.755,144.546,144.556,145.535,
    145.545,145.555,145.755,145.765,146.554,147.355,147.555,147.565,148.554,153.645,
    153.655,154.656,155.435,155.445,155.455,155.645,155.655,155.665,155.675,156.555,
    156.654,157.445,157.455,157.465,158.454,161.557,162.556,163.545,163.555,163.755,
    164.554,164.556,165.545,165.555,165.565,165.575,166.554,167.355,167.555,167.565,
    168.554,172.656,173.445,173.645,173.655,173.665,174.456,174.555,175.445,175.455,
    175.465,175.655,175.665,175.675,176.454,182.556,183.545,183.555,183.565,185.355,
    185.365,185.555,185.565,185.575,191.655,193.455,193.465,193.655,193.665,195.255,
    195.455,195.465,195.475,207.855,209.655,215.955,217.755,219.555,225.855,227.645,
    227.655,228.654,229.455,234.756,235.745,235.755,236.556,236.655,236.754,237.545,
    237.555,238.554,239.355,243.635,243.855,244.656,245.435,245.645,245.655,246.456,
    246.555,246.654,247.445,247.455,247.655,248.454,253.535,253.755,254.556,255.535,
    255.545,255.555,255.557,255.755,255.765,256.554,257.355,257.555,257.565,257.575,
    262.656,263.645,263.655,264.555,265.445,265.455,265.655,265.665,265.675,267.455,
    267.465,271.557,272.556,273.545,273.555,273.557,274.554,274.556,275.545,275.555,
    275.565,275.575,276.554,277.555,283.655,283.665,285.455,285.465,285.475,293.555,
    293.565,295.355,295.365,295.555,295.565,295.575,455.555};


    if ((fp_ad = fopen (adm_name,"r")) == NULL)
    {
        printf ("Cannot open admittance file?\n");
        exit (0);
    }
 

    ALL_OTS = (OTSconstit *) calloc ( NALL_OTS, sizeof(OTSconstit));
    for (i = 0; i < NALL_OTS; i ++)
    {
        ALL_OTS[i].ds = doodson[i];

        ALL_OTS[i].argn[0] = (int)(ALL_OTS[i].ds/100)%10;
        ALL_OTS[i].argn[1] = (int)(ALL_OTS[i].ds/10)%10 - 5;
        ALL_OTS[i].argn[2] = (int)(ALL_OTS[i].ds/1)%10 - 5;
        ALL_OTS[i].argn[3] = (int)(ALL_OTS[i].ds*10)%10 - 5;
        ALL_OTS[i].argn[4] = (int)(ALL_OTS[i].ds*100)%10 - 5;
        ALL_OTS[i].argn[5] = (int)(ALL_OTS[i].ds*1000)%10 - 5;
    }
 
   
    OTADM = (double *) calloc ( NALL_OTS * NMAJ_OTS, sizeof(double));
    for (j = 0; j < NMAJ_OTS; j++)
    {
        for (i = 0; i < NALL_OTS; i ++)
        {
            fscanf (fp_ad, "%lf", &OTADM[NALL_OTS * j + i]);
        }
    }

    fclose (fp_ad);    
        
    return 0; 
}





int ts_read_minor(double jdt, double gmst, int nmax, int minor, 
        int nper, OTSperturb *per, double *coef)
{
    double ang, cp, sp, cm, sm, fnm, ds, 
        cosang, sinang, fcos[18], fsin[18];
    int i,j, n,m, l, ind;



    for (i = 0; i < (nmax + 1) * (nmax + 1); i++)
    {
        coef[i] = 0;
    }

    if (minor == 1)
    {
        for (i = 0; i < NALL_OTS; i++)
        {
            arg2theta (jdt, gmst, ALL_OTS[i].argn, &ALL_OTS[i].ang);
            ALL_OTS[i].cosf = cos(ALL_OTS[i].ang);
            ALL_OTS[i].sinf = sin(ALL_OTS[i].ang);

        }

        for (j = 0; j < NMAJ_OTS; j++)
        {
            fcos[j] = 0;
            fsin[j] = 0;
            for (i = 0; i < NALL_OTS; i ++)
            {
                fcos[j] = fcos[j] + ALL_OTS[i].cosf * OTADM[NALL_OTS * j + i];
                fsin[j] = fsin[j] + ALL_OTS[i].sinf * OTADM[NALL_OTS * j + i];
            }
        }
    }

    

    for (i = 0; i < nper; i++)
    {
        if (per[i].n > nmax) 
//        if (per[i].n > nmax || per[i].n < 2) 
        {
            continue;
        }

        n = per[i].n;
        m = per[i].m;
        cp = per[i].cp;
        sp = per[i].sp;
        cm = per[i].cm;
        sm = per[i].sm;
        fnm = per[i].fnm;
        ds = per[i].ds;


        if (minor == 1)
        {
            if (ds == 55.565) j = 0;
            else if (ds == 55.575) j = 1;
            else if (ds == 56.554) j = 2;
            else if (ds == 57.555) j = 3;
            else if (ds == 65.455) j = 4;
            else if (ds == 75.555) j = 5;
            else if (ds == 85.455) j = 6;
            else if (ds == 93.555) j = 7;

            else if (ds == 135.655) j = 8;
            else if (ds == 145.555) j = 9;
            else if (ds == 163.555) j = 10;
            else if (ds == 165.555) j = 11;

            else if (ds == 235.755) j = 12;
            else if (ds == 245.655) j = 13;
            else if (ds == 255.555) j = 14;
            else if (ds == 273.555) j = 15;
            else if (ds == 275.555) j = 16;

            else if (ds == 455.555) j = 17;
            else 
            {
                printf ("error in otidecs: ds = %f not 18 major tides\n", ds);   
                exit (1);
            }
            cosang = fcos[j];
            sinang = fsin[j];
        }
        else if (minor == 0)
        {       
            arg2theta (jdt, gmst, per[i].argn, &ang);
            cosang = cos(ang);
            sinang = sin(ang);
        }
        else 
        {
            printf ("error in otidecs: minor = %d != 0 or 1\n", minor);   
            exit (1);
        }
    
      
        if (m == 0)
        {
            coef[n] = coef[n] + fnm * ((cp+cm) * cosang + (sp+sm) * sinang);
        }
        else
        {
            l = nmax - m + 1;
            ind = nmax + 1 + (2 * nmax - m + 2) * (m - 1);
            coef[ind + n - m] = coef[ind + n - m] + fnm * ((cp+cm) * cosang + (sp+sm) * sinang);
            coef[ind + n - m + l] = coef[ind + n - m + l] + fnm * ((sp-sm) * cosang - (cp-cm) * sinang);
        }
    }


    return 0;

}









double arg2theta (double jdt, double gmst, int *n, 
        double *ang)
{
    double t, a[5], b[6], theta, l, lp, F, D, Om;
    int i;
      
    t = (jdt - T0) / 36525.0;
    fund_args (t, a);
    l  = a[0];
    lp = a[1];
    F  = a[2];
    D  = a[3];
    Om = a[4];
//    Om = 0;

    b[1] = F + Om;
    b[0] = gmst + TWOPI / 2.0 - b[1];
    b[2] = b[1] - D;
    b[3] = b[1] - l;
//    b[4] = b[2];
    b[4] = - Om;
//    b[4] = 0;
    b[5] = b[1] - D - lp;

    theta = 0;
    for (i = 0; i < 6; i ++)
    {
        theta = theta + b[i] * n[i];
    }


    *ang = theta;

    return 0;
}







/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
double aod_open (char *file_aod, int nmax_aod, double *aod_eph)
{
    FILE *fp_aod;
    double c00, s00, c06, s06, c12, s12, c18, s18;
    int n,m, l, ind, par_aod;
    char string[400];

    if ((fp_aod = fopen (file_aod,"r")) == NULL)
    {
        printf ("Cannot open AOD file?\n");
        exit (0);
    }

    par_aod = (nmax_aod + 1) * (nmax_aod + 1);

    aod_eph[0 * (par_aod + 1)] = gps2tt(0);
    aod_eph[1 * (par_aod + 1)] = gps2tt(21600);
    aod_eph[2 * (par_aod + 1)] = gps2tt(21600 * 2);
    aod_eph[3 * (par_aod + 1)] = gps2tt(21600 * 3);

    while (1)
    {
        if (fgets (string, 400, fp_aod) == NULL) break;
        sscanf (string, "%d%d%lf%lf%lf%lf%lf%lf%lf%lf", 
            &n, &m, &c00, &s00, &c06, &s06, &c12, &s12, &c18, &s18);    

        if (n > nmax_aod )  continue;
        else if (m == 0)
        {
            aod_eph[0 * (par_aod + 1) + 1 + n] = c00;
            aod_eph[1 * (par_aod + 1) + 1 + n] = c06;
            aod_eph[2 * (par_aod + 1) + 1 + n] = c12;
            aod_eph[3 * (par_aod + 1) + 1 + n] = c18;
        }
        else 
        {
            l = nmax_aod - m + 1;
            ind = nmax_aod + 1 + (2 * nmax_aod - m + 2) * (m - 1);
            aod_eph[0 * (par_aod + 1) + 1 + ind + n - m] = c00;
            aod_eph[1 * (par_aod + 1) + 1 + ind + n - m] = c06;
            aod_eph[2 * (par_aod + 1) + 1 + ind + n - m] = c12;
            aod_eph[3 * (par_aod + 1) + 1 + ind + n - m] = c18;
            aod_eph[0 * (par_aod + 1) + 1 + ind + n - m + l] = s00;
            aod_eph[1 * (par_aod + 1) + 1 + ind + n - m + l] = s06;
            aod_eph[2 * (par_aod + 1) + 1 + ind + n - m + l] = s12;
            aod_eph[3 * (par_aod + 1) + 1 + ind + n - m + l] = s18;
        }
    }

    fclose(fp_aod);
    return 0;
}






double stidecs_earth(InfStruct *info, double *stcs, 
        double c20pt, int anelastic, int freqdep, int poletide)
//    nmax = 4;      
//    stcs = (double *) calloc ( (nmax + 1) * (nmax + 1), sizeof(double));
{
    short int moon = 9, earth = 2, sun = 10, n;

    double GMsun, gms2e, gmm2e, REk20, REk21, REk22, IMk21, IMk22,
        k20pa, k21pa, k22pa, k30, k31, k32, k33,
        ps[3], vs[3], pm[3], vm[3], tjd[2],
        pse[3], pme[3], llrs[3], llrm[3], pbar[4], t,
        p20m, p30m, p21m, p31m, p22m, p32m, p33m, 
        p20s, p30s, p21s, p31s, p22s, p32s, p33s,
        rerm, rers, c20, c30, c40, c21, s21, c22, s22, c31, s31, 
        c32, s32, c33, s33, c41, s41, c42, s42, 
        c20f, c21f, s21f, c22f, s22f, c21p, s21p, m1, m2, GM, radius;

    GM      = 398600.44180E+09;
    radius  = 6378136.6;

    GMsun   = 1.32712442076e20;
    gmm2e   = 0.0123000383;
    gms2e   = GMsun/GM;

    if (anelastic == 0)
    {
        REk20 =  0.29525;
        REk21 =  0.29470;
        REk22 =  0.29801;
        IMk21 =  0;
        IMk22 =  0;
        k20pa = -0.00087;
        k21pa = -0.00079;
        k22pa = -0.00057;
    }
    else
    {
        REk20 =  0.30190;
        REk21 =  0.29830;
        REk22 =  0.30102;
        IMk21 = -0.00144;
        IMk22 = -0.00130;
        k20pa = -0.00089;
        k21pa = -0.00080;
        k22pa = -0.00057;
    }

    k30   =  0.093;
    k31   =  0.093;
    k32   =  0.093;
    k33   =  0.094;






    tjd[0] = info->jd0;
    tjd[1] = info->tt/86400.0;

// Luni-solar ephemeris

    planet_ephemeris (tjd, sun, earth, ps, vs);
    planet_ephemeris (tjd, moon, earth, pm, vm);
    for (n = 0; n < 3; n++)
    {
        ps[n] = ps[n] * AU;
        pm[n] = pm[n] * AU;
    }
 
    
//    icrf2itrf(num, ps, pse);
//    icrf2itrf(num, pm, pme);

    brmul (info->c_ie, ps, 3, 3, 1, pse);   //inertial to fixed matrix gmat = rmat*tbt
    brmul (info->c_ie, pm, 3, 3, 1, pme);   //inertial to fixed matrix gmat = rmat*tbt


    xyz2llr(pse, llrs);
    xyz2llr(pme, llrm);

    t = sin(llrm[0] * DEG2RAD);
    lgdr(t, 3, 0, pbar); p20m = pbar[2]; p30m = pbar[3];
    lgdr(t, 3, 1, pbar); p21m = pbar[1]; p31m = pbar[2];
    lgdr(t, 3, 2, pbar); p22m = pbar[0]; p32m = pbar[1];
    lgdr(t, 3, 3, pbar); p33m = pbar[0];
    t = sin(llrs[0] * DEG2RAD);
    lgdr(t, 3, 0, pbar); p20s = pbar[2]; p30s = pbar[3];
    lgdr(t, 3, 1, pbar); p21s = pbar[1]; p31s = pbar[2];
    lgdr(t, 3, 2, pbar); p22s = pbar[0]; p32s = pbar[1];
    lgdr(t, 3, 3, pbar); p33s = pbar[0];


    rerm = radius / llrm[2];
    rers = radius / llrs[2];

// Frequency Independent Terms

// C20
    c20 = REk20/5.0 * ( gmm2e * pow(rerm, 3) * p20m 
                    + gms2e * pow(rers, 3) * p20s );
// C21/S21
    c21 = + REk21/5.0 * ( gmm2e * pow(rerm, 3) * p21m * cos(llrm[1] * DEG2RAD)
                    + gms2e * pow(rers, 3) * p21s * cos(llrs[1] * DEG2RAD) )
          + IMk21/5.0 * ( gmm2e * pow(rerm, 3) * p21m * sin(llrm[1] * DEG2RAD) 
                    + gms2e * pow(rers, 3) * p21s * sin(llrs[1] * DEG2RAD) );


    s21 = - IMk21/5.0 * ( gmm2e * pow(rerm, 3) * p21m * cos(llrm[1] * DEG2RAD)
                    + gms2e * pow(rers, 3) * p21s * cos(llrs[1] * DEG2RAD) )
          + REk21/5.0 * ( gmm2e * pow(rerm, 3) * p21m * sin(llrm[1] * DEG2RAD) 
                    + gms2e * pow(rers, 3) * p21s * sin(llrs[1] * DEG2RAD) );


// C22/S22
    c22 = + REk22/5.0 * ( gmm2e * pow(rerm, 3) * p22m * cos(llrm[1] * DEG2RAD * 2.0)
                    + gms2e * pow(rers, 3) * p22s * cos(llrs[1] * DEG2RAD * 2.0) )
          + IMk22/5.0 * ( gmm2e * pow(rerm, 3) * p22m * sin(llrm[1] * DEG2RAD * 2.0)
                    + gms2e * pow(rers, 3) * p22s * sin(llrs[1] * DEG2RAD * 2.0) );

    s22 = - IMk22/5.0 * ( gmm2e * pow(rerm, 3) * p22m * cos(llrm[1] * DEG2RAD * 2.0)
                    + gms2e * pow(rers, 3) * p22s * cos(llrs[1] * DEG2RAD * 2.0) )
          + REk22/5.0 * ( gmm2e * pow(rerm, 3) * p22m * sin(llrm[1] * DEG2RAD * 2.0)
                    + gms2e * pow(rers, 3) * p22s * sin(llrs[1] * DEG2RAD * 2.0) );


// C30
    c30 = k30/7.0 * ( gmm2e * pow(rerm, 4) * p30m 
                    + gms2e * pow(rers, 4) * p30s );
// C31/S31
    c31 = k31/7.0 * ( gmm2e * pow(rerm, 4) * p31m * cos(llrm[1] * DEG2RAD)
                    + gms2e * pow(rers, 4) * p31s * cos(llrs[1] * DEG2RAD) );
    s31 = k31/7.0 * ( gmm2e * pow(rerm, 4) * p31m * sin(llrm[1] * DEG2RAD)
                    + gms2e * pow(rers, 4) * p31s * sin(llrs[1] * DEG2RAD) );
// C32/S32
    c32 = k32/7.0 * ( gmm2e * pow(rerm, 4) * p32m * cos(llrm[1] * DEG2RAD * 2.0)
                    + gms2e * pow(rers, 4) * p32s * cos(llrs[1] * DEG2RAD * 2.0) );
    s32 = k32/7.0 * ( gmm2e * pow(rerm, 4) * p32m * sin(llrm[1] * DEG2RAD * 2.0)
                    + gms2e * pow(rers, 4) * p32s * sin(llrs[1] * DEG2RAD * 2.0) );
// C33/S33
    c33 = k33/7.0 * ( gmm2e * pow(rerm, 4) * p33m * cos(llrm[1] * DEG2RAD * 3.0)
                    + gms2e * pow(rers, 4) * p33s * cos(llrs[1] * DEG2RAD * 3.0) );
    s33 = k33/7.0 * ( gmm2e * pow(rerm, 4) * p33m * sin(llrm[1] * DEG2RAD * 3.0)
                    + gms2e * pow(rers, 4) * p33s * sin(llrs[1] * DEG2RAD * 3.0) );

// C40
    c40 = k20pa/5.0* ( gmm2e * pow(rerm, 3) * p20m 
                    + gms2e * pow(rers, 3) * p20s );
// C41/S41
    c41 = k21pa/5.0* ( gmm2e * pow(rerm, 3) * p21m * cos(llrm[1] * DEG2RAD)
                    + gms2e * pow(rers, 3) * p21s * cos(llrs[1] * DEG2RAD) );
    s41 = k21pa/5.0* ( gmm2e * pow(rerm, 3) * p21m * sin(llrm[1] * DEG2RAD)
                    + gms2e * pow(rers, 3) * p21s * sin(llrs[1] * DEG2RAD) );
// C42/S42
    c42 = k22pa/5.0* ( gmm2e * pow(rerm, 3) * p22m * cos(llrm[1] * DEG2RAD * 2.0)
                    + gms2e * pow(rers, 3) * p22s * cos(llrs[1] * DEG2RAD * 2.0) );
    s42 = k22pa/5.0* ( gmm2e * pow(rerm, 3) * p22m * sin(llrm[1] * DEG2RAD * 2.0)
                    + gms2e * pow(rers, 3) * p22s * sin(llrs[1] * DEG2RAD * 2.0) );

    
    
    
    stcs[0]  = 0; //c00; 
    stcs[1]  = 0; //c10; 
    stcs[2]  = c20; 
    stcs[3]  = c30; 
    stcs[4]  = c40; 

    stcs[5]  = 0; //c11; 
    stcs[6]  = c21; 
    stcs[7]  = c31; 
    stcs[8]  = c41; 
    stcs[9]  = 0; //s11; 
    stcs[10] = s21; 
    stcs[11] = s31; 
    stcs[12] = s41; 

    stcs[13] = c22; 
    stcs[14] = c32; 
    stcs[15] = c42; 
    stcs[16] = s22; 
    stcs[17] = s32; 
    stcs[18] = s42; 

    stcs[19] = c33; 
    stcs[20] = 0; //c43; 
    stcs[21] = s33; 
    stcs[22] = 0; //s43; 

    stcs[23] = 0; //c44; 
    stcs[24] = 0; //s44; 
  
// permanent tide
    stcs[2]  = stcs[2] - c20pt; 

    c20f = 0; c21f = 0; s21f = 0; c22f = 0; s22f = 0;
    c21p = 0; s21p = 0;

// Frequency Dependent Terms
    if (freqdep == 1)
        stfrqdep(info->jdt, info->gmst, &c20f, &c21f, &s21f, &c22f, &s22f);

// Solid Earth pole tide
    if (poletide == 1)
    {
        mtpole (info->mjd, info->xp, info->yp, &m1, &m2);
        c21p = -1.333e-9 * (m1 + 0.0115 * m2);
        s21p = -1.333e-9 * (m2 - 0.0115 * m1);
    }

    stcs[2]  += c20f; 
    stcs[6]  += c21f + c21p; 
    stcs[10] += s21f + s21p; 
    stcs[13] += c22f; 
    stcs[16] += s22f; 

    return 0;

}









double stfrqdep(double jdt, double gmst, 
        double *c20f, double *c21f, double *s21f, 
        double *c22f, double *s22f)
{
    static const double sets[71][8] = {
    { 0,5,5,5,6,5,  16.6e-12,  -6.7e-12 }, 
    { 0,5,5,5,7,5,  -0.1e-12,   0.1e-12 }, 
    { 0,5,6,5,5,4,  -1.2e-12,   0.8e-12 }, 
    { 0,5,7,5,5,5,  -5.5e-12,   4.3e-12 }, 
    { 0,5,7,5,6,5,   0.1e-12,  -0.1e-12 }, 
    { 0,5,8,5,5,4,  -0.3e-12,   0.2e-12 }, 
    { 0,6,3,6,5,5,  -0.3e-12,   0.7e-12 }, 
    { 0,6,5,4,4,5,   0.1e-12,  -0.2e-12 }, 
    { 0,6,5,4,5,5,  -1.2e-12,   3.7e-12 }, 
    { 0,6,5,4,6,5,   0.1e-12,  -0.2e-12 }, 
    { 0,6,5,6,5,5,   0.1e-12,  -0.2e-12 }, 
    { 0,7,3,5,5,5,   0.0e-12,   0.6e-12 }, 
    { 0,7,5,3,5,5,   0.0e-12,   0.3e-12 }, 
    { 0,7,5,5,5,5,   0.6e-12,   6.3e-12 }, 
    { 0,7,5,5,6,5,   0.2e-12,   2.6e-12 }, 
    { 0,7,5,5,7,5,   0.0e-12,   0.2e-12 }, 
    { 0,8,3,6,5,5,   0.1e-12,   0.2e-12 }, 
    { 0,8,5,4,5,5,   0.4e-12,   1.1e-12 }, 
    { 0,8,5,4,6,5,   0.2e-12,   0.5e-12 }, 
    { 0,9,3,5,5,5,   0.1e-12,   0.2e-12 }, 
    { 0,9,5,3,5,5,   0.1e-12,   0.1e-12 }, 
    { 1,2,5,7,5,5,  -0.1e-12,   0.0e-12 }, 
    { 1,2,7,5,5,5,  -0.1e-12,   0.0e-12 }, 
    { 1,3,5,6,4,5,  -0.1e-12,   0.0e-12 }, 
    { 1,3,5,6,5,5,  -0.7e-12,   0.1e-12 }, 
    { 1,3,7,4,5,5,  -0.1e-12,   0.0e-12 }, 
    { 1,4,5,5,4,5,  -1.3e-12,   0.1e-12 }, 
    { 1,4,5,5,5,5,  -6.8e-12,   0.6e-12 }, 
    { 1,4,7,5,5,5,   0.1e-12,   0.0e-12 }, 
    { 1,5,3,6,5,5,   0.1e-12,   0.0e-12 }, 
    { 1,5,5,4,4,5,   0.1e-12,   0.0e-12 }, 
    { 1,5,5,4,5,5,   0.4e-12,   0.0e-12 }, 
    { 1,5,5,6,5,5,   1.3e-12,  -0.1e-12 }, 
    { 1,5,5,6,6,5,   0.3e-12,   0.0e-12 }, 
    { 1,5,7,4,5,5,   0.3e-12,   0.0e-12 }, 
    { 1,5,7,4,6,5,   0.1e-12,   0.0e-12 }, 
    { 1,6,2,5,5,6,  -1.9e-12,   0.1e-12 }, 
    { 1,6,3,5,4,5,   0.5e-12,   0.0e-12 }, 
    { 1,6,3,5,5,5, -43.4e-12,   2.9e-12 }, 
    { 1,6,4,5,5,4,   0.6e-12,   0.0e-12 }, 
    { 1,6,4,5,5,6,   1.6e-12,  -0.1e-12 }, 
    { 1,6,5,3,4,5,   0.1e-12,   0.0e-12 }, 
    { 1,6,5,5,3,5,   0.1e-12,   0.0e-12 }, 
    { 1,6,5,5,4,5,  -8.8e-12,   0.5e-12 }, 
    { 1,6,5,5,5,5, 470.9e-12, -30.2e-12 }, 
    { 1,6,5,5,6,5,  68.1e-12,  -4.6e-12 }, 
    { 1,6,5,5,7,5,  -1.6e-12,   0.1e-12 }, 
    { 1,6,6,4,5,5,   0.1e-12,   0.0e-12 }, 
    { 1,6,6,5,4,4,  -0.1e-12,   0.0e-12 }, 
    { 1,6,6,5,5,4, -20.6e-12,  -0.3e-12 }, 
    { 1,6,6,5,5,6,   0.3e-12,   0.0e-12 }, 
    { 1,6,6,5,6,4,  -0.3e-12,   0.0e-12 }, 
    { 1,6,7,3,5,5,  -0.2e-12,   0.0e-12 }, 
    { 1,6,7,3,6,5,  -0.1e-12,   0.0e-12 }, 
    { 1,6,7,5,5,5,  -5.0e-12,   0.3e-12 }, 
    { 1,6,7,5,6,5,   0.2e-12,   0.0e-12 }, 
    { 1,6,8,5,5,4,  -0.2e-12,   0.0e-12 }, 
    { 1,7,3,6,5,5,  -0.5e-12,   0.0e-12 }, 
    { 1,7,3,6,6,5,  -0.1e-12,   0.0e-12 }, 
    { 1,7,5,4,4,5,   0.1e-12,   0.0e-12 }, 
    { 1,7,5,4,5,5,  -2.1e-12,   0.1e-12 }, 
    { 1,7,5,4,6,5,  -0.4e-12,   0.0e-12 }, 
    { 1,8,3,5,5,5,  -0.2e-12,   0.0e-12 }, 
    { 1,8,5,3,5,5,  -0.1e-12,   0.0e-12 }, 
    { 1,8,5,5,5,5,  -0.6e-12,   0.0e-12 }, 
    { 1,8,5,5,6,5,  -0.4e-12,   0.0e-12 }, 
    { 1,8,5,5,7,5,  -0.1e-12,   0.0e-12 }, 
    { 1,9,5,4,5,5,  -0.1e-12,   0.0e-12 }, 
    { 1,9,5,4,6,5,  -0.1e-12,   0.0e-12 }, 
    { 2,4,5,6,5,5,  -0.3e-12,   0.0e-12 }, 
    { 2,5,5,5,5,5,  -1.2e-12,   0.0e-12 }};

    double ang, c20 = 0, c21 = 0, s21 = 0, c22 = 0, s22 = 0;
    int i, nsets = 71, argn[6];

    for (i=0;i<nsets;i++)
    {
        argn[0] = (int)sets[i][0];
        argn[1] = (int)sets[i][1] - 5;
        argn[2] = (int)sets[i][2] - 5;
        argn[3] = (int)sets[i][3] - 5;
        argn[4] = (int)sets[i][4] - 5;
        argn[5] = (int)sets[i][5] - 5;
        
    
        arg2theta (jdt, gmst, argn, &ang);
        
        

//  C20 correction: Long period tidal constituent

        if(argn[0]==0) 
        {
           c20 = c20 + sets[i][6]*cos(ang) - sets[i][7]*sin(ang);
        }

//  C21/S21 correction: Diurnal period tidal constituent

        if(argn[0]==1) 
        {
            c21 = c21 + sets[i][6]*sin(ang) + sets[i][7]*cos(ang);
            s21 = s21 + sets[i][6]*cos(ang) - sets[i][7]*sin(ang);
        }

//  C22/S22 correction: Semi-diurnal period tidal constituent

        if(argn[0]==2) 
        {
            c22 = c22 + sets[i][6]*cos(ang);
            s22 = s22 - sets[i][6]*sin(ang);
        }
    }

    *c20f = c20;
    *c21f = c21;
    *s21f = s21;
    *c22f = c22;
    *s22f = s22;



    return 0;

}



