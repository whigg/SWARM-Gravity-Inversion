/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/*

  GRACE L1B to L1C & L2

  Version: 20 Dec 2011 
           fixed and stochastic constraint solution

  Copyright (c) 2011 Kun Shang (shang.34@osu.edu) All Right Reserved 

*/
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/

#ifndef __OMP_H__
    #include <omp.h>
#endif

//#include "l1b2l1c.h"
#include "grvts.h"
#include "coord.h"
#include "numrs.h"


//#define LFIX 4

int GPS_S, GPS_E, NDATA, DT;
double GMA[2];    

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
int main (int argc, char *argv[])
{
    FILE *fp_stdin, *fpin_l1c, *fpout_adj, *fpout_coef;  
    int i, n, nmax, mmax, slvls, l_old, ndat_max, factor_b1, factor_b2, 
        m, l, k, ind, num, *num_arc, *gps_arc, narc_max, *ind_m,
        LFIX = 0, ifix[4000], LSTK = 0, istk[4000], 
        nply, ncpr, nemp, ndt, narc, npd, *gpst, *gps_eph, 
        solvefor, solvegrv, solveemp;
    short int year, month, day;
    double *tp1, *csr, *gfz, *jpl, 
        *obs, GMR2, r2R, *coefcsr,  hour,
        *coefjpl, *coefgfz, *coefggm, JD0, t, Tprd, 
        *pt1, *pt2, *ai, *coefa, *coefe, *coefaam, factor_n, 
        *lambda, *lambdam, *ksiupd, *ksiupdm, *k0, *aksi, *l1cf, 
        *factor_a, *factor_b, factor_p, factor_r,
        *ksi, *ksim, *dksi, *dksiupd, r0[1], sigma0, pweight, *p0, *p0b0, 
        *knk0, *knk, *kksi, *id, *kfix, *kfixt,*kn, *nk, *nkknk,
        *atpa, *atpy, *cest, zero = 0.0, 
        omega, etpe, c, s, vfix[4000], vstk[4000],
        *llr1, *llr2, *l1c_eph, *l1ct, vb1, vb2;
    char line[500], card[20],
        f_aam[2][200]={"\0", "\0"}, 
         stdname[200], l1cfile[400];    
    time_t s0, s1, s2, s3, s4, s5, s6;
    long nl, solveforN, solvegrvN, solvefor2, ns, is, in;

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
// read input cards
    LFIX = 0;
    if (argc == 1 || argc == 2)
    {
        printf ("input file name: ");
        scanf ("%s%s", stdname, l1cfile);
    }
    else if (argc == 3)
    {
        sscanf (argv[1], "%s", stdname);
        sscanf (argv[2], "%s", l1cfile);
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
//    printf ("reading input file \"%s\"...\n", stdname);

    s0 = time(NULL);  
    while (feof(fp_stdin) == 0)
    {
        card[1] = '\0';
        if (fgets (line, 500, fp_stdin) == NULL) break;
        sscanf (line, "%s", card);  
        
        if (strcmp (card,"AAM") ==0)
        {
            sscanf (line, "%*s%s%s", f_aam[0], f_aam[1]);
        }
        if (strcmp (card, "DT") ==0)
        {
            sscanf (line, "%*s%d", &DT);
        }
        if (strcmp (card, "NMAX") ==0)  
        {
            sscanf (line, "%*s%d%d", &nmax, &mmax);
        }
        if (strcmp (card, "EMPR") ==0)  
        {
            sscanf (line, "%*s%d%d%d", &nply, &ncpr, &ndt);
        }
        
        if (strcmp (card, "SOLVE") ==0)  
        {
            sscanf (line, "%*s%d", &slvls);
        }


        if (strcmp (card, "FIX") ==0)   
        {
            sscanf (line, "%*s%d%d%lf%lf", &n, &m, &c, &s);
            if (m == 0) 
            {
                ifix[LFIX] = n;
                vfix[LFIX] = c;
                LFIX = LFIX + 1;
            }
            else 
            {
                l = nmax - m + 1;
                ind = nmax + 1 + (2 * nmax - m + 2) * (m - 1);

                ifix[LFIX] = ind + n - m;
                ifix[LFIX + 1] = ind + n - m + l;
                vfix[LFIX] = c;
                vfix[LFIX + 1] = s;
                LFIX = LFIX + 2;
            }
        }


        if (strcmp (card, "STK") ==0)   
        {
            sscanf (line, "%*s%d%d%lf%lf", &n, &m, &c, &s);
            if (m == 0) 
            {
                istk[LSTK] = n;
                vstk[LSTK] = c;
                LSTK = LSTK + 1;
            }
            else 
            {
                l = nmax - m + 1;
                ind = nmax + 1 + (2 * nmax - m + 2) * (m - 1);

                istk[LSTK] = ind + n - m;
                istk[LSTK + 1] = ind + n - m + l;
                vstk[LSTK] = c;
                vstk[LSTK + 1] = s;
                LSTK = LSTK + 2;
            }
        }

        if (strcmp (card, "FCT_P") ==0) 
        {
            sscanf (line, "%*s%lf", &factor_p);
        }
        if (strcmp (card, "FCT_B") ==0) 
        {
            sscanf (line, "%*s%d%d", &factor_b1, &factor_b2);
        }
        if (strcmp (card, "FCT_R") ==0) 
        {
            sscanf (line, "%*s%lf", &factor_r);
        }
        if (strcmp (card, "FCT_N") ==0)
        {
            sscanf (line, "%*s%lf", &factor_n);
        }
        if (strcmp (card, "PWEIGHT") ==0)   
        {
            sscanf (line, "%*s%lf", &pweight);
        }

    }


/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
// open output files 


    if ( (fpin_l1c = fopen (l1cfile,"r")) == NULL)
    {
        printf ("Cannot open fpout_a file!\n");
        exit (0);
    }

    if ( (fpout_coef = fopen ("COEF.txt","w")) == NULL)
    {
        printf ("Cannot open fpout_coef file!\n");
        exit (0);
    }

    if ( (fpout_adj = fopen ("adj.txt","w")) == NULL)
    {
        printf ("Cannot open fpout_adj file!\n");
        exit (0);
    }


    if (fgets(line,400, fpin_l1c) ==NULL) exit(0);
    sscanf (line, "%d", &GPS_S);
    NDATA = 1;
    while (1)
    {
        if (fgets(line,400, fpin_l1c) ==NULL) break;
        NDATA ++;
    }
    sscanf (line, "%d", &GPS_E);
    rewind(fpin_l1c);




//    NDATA = days * 86400 / DT;
    printf ("NDATA = %d\n", NDATA);

    l1c_eph = (double *) calloc (11 * NDATA, sizeof(double));
    gps_eph = (int *) calloc ( NDATA, sizeof(int));

    llr1 = (double *) calloc (3 * NDATA, sizeof(double));
    llr2 = (double *) calloc (3 * NDATA, sizeof(double));
    l1ct = (double *) calloc ( NDATA, sizeof(double));
    obs = (double *) calloc ( NDATA, sizeof(double));
    gpst = (int *) calloc ( NDATA, sizeof(int));
    tp1 = (double *) calloc ( NDATA, sizeof(double));
    csr = (double *) calloc ( NDATA, sizeof(double));
    jpl = (double *) calloc ( NDATA, sizeof(double));
    gfz = (double *) calloc ( NDATA, sizeof(double));

//    if (fgets(line,400, fpin_l1c) ==NULL) exit(0);
//    sscanf (line, "%d", &GPS_S);
//    rewind(fpin_l1c);

    printf("GPS_S = %d\n", GPS_S);
    JD0 = GPS_S / 86400.0 + T0;
    cal_date (JD0, &year, &month, &day, &hour);
    JD0 = julian_date (year,month,day,0);
    GPS_S = (int)((JD0 - T0) * 86400 + 0.5);
    printf("GPS_S = %d\tyear = %d\tmonth = %d\tday = %d\n", GPS_S, year, month, day);


    printf("GPS_E = %d\n", GPS_E);
    JD0 = GPS_E / 86400.0 + T0;
    cal_date (JD0, &year, &month, &day, &hour);
    day = day + 1;
    JD0 = julian_date (year,month,day,0);
    GPS_E = (int)((JD0 - T0) * 86400 + 0.5);
    printf("GPS_E = %d\tyear = %d\tmonth = %d\tday = %d\n", GPS_E, year, month, day);


    i = 0;
    while (1)
    {
        if (fgets(line,400, fpin_l1c) ==NULL) break;
        sscanf (line, "%d%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%*f%*f%*f%lf",
            &gps_eph[i], &l1c_eph[i * 11], &l1c_eph[i * 11 + 1], &l1c_eph[i * 11 + 2], 
            &l1c_eph[i * 11 + 3], &l1c_eph[i * 11 + 4], &l1c_eph[i * 11 + 5], &l1c_eph[i * 11 + 6],
            &l1c_eph[i * 11 + 7], &l1c_eph[i * 11 + 8], &l1c_eph[i * 11 + 9], &l1c_eph[i * 11 + 10]);
        i ++;
    }
    printf ("NDATA = %d\t i = %d\n", NDATA, i);

    solveemp = 0;
    nemp = nply + ncpr;

    if (nemp == 0)
    {
        n = 0;
        for (i = 0; i < NDATA; i ++)
        {
            if ((gps_eph[i] - GPS_S) % DT == 0)
            {
                gpst[n] = gps_eph[i];
                llr1[n * 3]     = l1c_eph[i * 11];
                llr1[n * 3 + 1] = l1c_eph[i * 11 + 1];
                llr1[n * 3 + 2] = l1c_eph[i * 11 + 2];
                llr2[n * 3]     = l1c_eph[i * 11 + 3];
                llr2[n * 3 + 1] = l1c_eph[i * 11 + 4];
                llr2[n * 3 + 2] = l1c_eph[i * 11 + 5];
                l1ct[n]         = l1c_eph[i * 11 + 6];
                csr[n]          = l1c_eph[i * 11 + 7];
                gfz[n]          = l1c_eph[i * 11 + 8];
                jpl[n]          = l1c_eph[i * 11 + 9];
                tp1[n]          = l1c_eph[i * 11 + 10];
                n++;
            }
        }
        NDATA = n;                    
        printf ("NDATA = %d where nemp == 0\n", NDATA);
    }


    if (nemp != 0)
    {
        npd = (int)(86400/ndt);
        narc_max = npd * (int)((GPS_E - GPS_S) / 86400.0);

        num_arc = (int *) calloc ( narc_max, sizeof(int));
        gps_arc = (int *) calloc ( narc_max, sizeof(int));
    
        ndat_max = (int)(ndt / DT);
        
        ind_m = (int *) calloc (ndat_max, sizeof(int));

        n = 0;
        m = 0;
        narc = 0;
        l_old = 0;
        for (i = 0; i < NDATA - 1; i ++)
        {
            if ((gps_eph[i] - GPS_S) % DT == 0)
            {
                ind_m[m] = i;
                m ++;
                l = (int)((gps_eph[i + 1] - GPS_S)/ndt);
                if (l < 0) l = 0;
                if (l > narc_max - 1) l = narc_max - 1;
                if (l == l_old)
                {
                    l_old = l;
                }
                else if (l > l_old)
                {
//                    JD0 = gps_eph[i] / 86400.0 + T0;
//                    cal_date (JD0, &year, &month, &day, &hour);
//                    printf ("l = %d\t l_old = %d\t year = %d month = %d day = %d hour = %f\t m = %d\n", 
//                            l, l_old, year, month, day, hour, m);
//                    l_old = l_old + 1;
                    if (m > ndat_max * 0.5 && m <= ndat_max)
                    {
                        num = m;
//                        m = 0;
//                        l = (int)((gps_eph[i] - GPS_S)/ndt);
//            if (l < 0) 
//            {   printf ("l = %d\n"); l = 0;}
//            if (l > narc_max - 1) 
//            {   printf ("l = %d\n"); l = narc_max - 1;}
                        num_arc[l_old] = narc;
                        gps_arc[l_old] = gps_eph[i];
                        narc ++;
                        for (k = 0; k < num; k ++)
                        {   
                            ind = ind_m[k];
                            gpst[n] = gps_eph[ind];
                            llr1[n * 3]     = l1c_eph[ind * 11];
                            llr1[n * 3 + 1] = l1c_eph[ind * 11 + 1];
                            llr1[n * 3 + 2] = l1c_eph[ind * 11 + 2];
                            llr2[n * 3]     = l1c_eph[ind * 11 + 3];
                            llr2[n * 3 + 1] = l1c_eph[ind * 11 + 4];
                            llr2[n * 3 + 2] = l1c_eph[ind * 11 + 5];
                            l1ct[n]         = l1c_eph[ind * 11 + 6];
                            csr[n]          = l1c_eph[ind * 11 + 7];
                            gfz[n]          = l1c_eph[ind * 11 + 8];
                            jpl[n]          = l1c_eph[ind * 11 + 9];
                            tp1[n]          = l1c_eph[ind * 11 + 10];
                            n++;
                        }
                    }
                    m = 0;
                    l_old = l_old + 1;
                }
            }
        }
        NDATA = n;                    
//        narc = narc - 1;
        solveemp = nemp * narc;
        printf ("NDATA = %d where nemp == %d\n", NDATA, nemp);
        printf ("narc = %d where narc_max == %d\n", narc, narc_max);
    }

    solvegrv = (nmax + 1) * (nmax + 1);
    solvefor = solvegrv + solveemp;


    coefaam  = (double *) calloc (solvegrv, sizeof(double));

    printf ("nmax = %d\n", nmax);
    printf ("solvefor = %d\n", solvefor);
    printf ("solvegrv = %d\n", solvegrv);

    solveforN = (long)(1.0 * solvefor * NDATA);
    solvegrvN = (long)(1.0 * solvegrv * NDATA);
    solvefor2 = (long)(1.0 * solvefor * solvefor);

    printf ("solveforN = %ld\t%fG\n", solveforN, solveforN /1024.0/1024.0/1024.0);
    printf ("solvegrvN = %ld\t%fG\n", solvegrvN, solvegrvN /1024.0/1024.0/1024.0);
    printf ("solvefor2 = %ld\t%fG\n", solvefor2, solvefor2 /1024.0/1024.0/1024.0);
    
    if (solveforN /1024.0/1024.0/1024.0 > 2) printf ("warning: solveforN > INT LIMIT!!!\n");
    if (solvegrvN /1024.0/1024.0/1024.0 > 2) printf ("warning: solvegrvN > INT LIMIT!!!\n");
    if (solvefor2 /1024.0/1024.0/1024.0 > 2) printf ("warning: solvefor2 > INT LIMIT!!!\n");

    GMA[0] = 398600.44150E+09;
    GMA[1] = 6378136.3;


//    opengrav (f_aam, coefaam, GMA, nmax, mmax, 0);
    grv_open (f_aam[0], f_aam[1], nmax, mmax, coefaam);
    zero0zero1 (coefaam, nmax);


//    printf ("succeed!\n");
//    fflush(stdout);

    atpy = (double *) calloc ( solvefor, sizeof(double));
    atpa = (double *) calloc ( solvefor2, sizeof(double));
    for (n = 0; n < solvefor; n++)
        atpy[n] = 0;
    for (nl = 0; nl < solvefor2; nl++)
        atpa[nl] = 0;

//    printf ("succeed!\n");
//    fflush(stdout);

    ai = (double *) calloc ( solveforN, sizeof(double));
    for (nl = 0; nl < solveforN; nl++)
        ai[nl] = 0;


/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
// open gravity field file


//    printf ("succeed!\n");
//    fflush(stdout);

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
// calibrate accelerometer data


    s1 = time(NULL);  
    printf("\n%5ld: seconds of reading L1C data\n", s1-s0);
    fflush(stdout);

#pragma omp parallel private(pt1, pt2, i, n, m, vb1, vb2, t, l, Tprd, in)
    {
        pt1  = (double *) calloc ( solvegrv, sizeof(double));
        pt2  = (double *) calloc ( solvegrv, sizeof(double));

        #pragma omp for
    for (i = 0; i < NDATA; i ++)
    {
        cs2gp_pt (&llr1[i * 3], coefaam, GMA[0], GMA[1], nmax, &vb1, pt1);
        cs2gp_pt (&llr2[i * 3], coefaam, GMA[0], GMA[1], nmax, &vb2, pt2);
//        obs[i]= vb2 - vb1;    // background model n=0~60
        obs[i]= l1ct[i];    // background model n=0~60
//        printf ("obs = %e\n", obs[i]);
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
//y=Ax+e; form A; add bias & trend estimation or not
/////////////////////////////////////////
//all the index should be long int
//
////////////////////////////////////////
        for(n = 0; n < solvegrv; n++)
        {
            in = (long)(1.0 * i * solvefor + n);
//            a1i[i * solvegrv + n] = pt2[n] - pt1[n];
            ai[in] = pt2[n] - pt1[n];
        }

        if (nemp != 0)
        {
            l = (int)((gpst[i] - GPS_S)/ndt);
            if (l < 0) l = 0;
            if (l > narc_max - 1) l = narc_max - 1;

            t = (gpst[i] - gps_arc[l]) / 86400.0;
            for (n = 0; n < nply; n++)
            {
                in = (long)(1.0 * i * solvefor + solvegrv + num_arc[l] * nemp + n);
                ai[in] = pow(t, n);
            }
            Tprd = tp1[i] / 86400.0;
            for (n = 0; n < ncpr; n = n + 2)
            {
                in = (long)(1.0 * i * solvefor + solvegrv + num_arc[l] * nemp + nply + n);
                ai[in] = sin (TWOPI / Tprd * t);
                in = (long)(1.0 * i * solvefor + solvegrv + num_arc[l] * nemp + nply + n + 1);
                ai[in] = cos (TWOPI / Tprd * t);
                Tprd = Tprd / 2.0;
            }
        }
        
    }
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/


        free(pt1);
        free(pt2);

    }

    s2 = time(NULL);  
    printf("\n%5ld: seconds of processing earth gravity and partial\n", s2-s1);
    fflush(stdout);




//    mt(ai, NDATA, solvefor, at);
//    brmul(at, ai, solvefor,NDATA,solvefor, atpa);
//    brmul(at,lhs,solvefor,NDATA,1,atpy);


/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
// HARD CODE!!!!!!  least-squre algorithm
// HARD CODE!!!!!!  least-squre algorithm
// HARD CODE!!!!!!  least-squre algorithm

    p0 = (double *) calloc ( solvefor, sizeof(double));
    p0b0 = (double *) calloc ( solvefor, sizeof(double));

    
///////////////////////////////////////////////////////////////
// start regularization
// 
    factor_a = (double *) calloc ( nmax + 1, sizeof(double));
    factor_b = (double *) calloc ( nmax + 1, sizeof(double));

//    pweight = 1/0.0001;
//    pweight = 1;


    GMR2 = 4e15;
    r2R = 6850/6378.136;

/////////////// parameter w.r.t degree n: factor_a
    for (n = 0; n <= nmax; n++)
    {
        factor_a[n] = pow(r2R, factor_p * (n+1));
    }


//////////////// parameter w.r.t. order m: factor_b
//    for (m = 0; m <= nmax; m++)
    for (m = factor_b1; m <= factor_b2; m++)
    {
//        factor_b[m] = m*m;
      factor_b[m] = pow (m, factor_n);
    }

//    for (m = factor_b1; m <= factor_b2; m++)
//    {
//        factor_b[m] = 1;
//    }

//    factor_b[14] = 1; factor_b[15] = 1; factor_b[16] = 1;  
//    factor_b[29] = 1; factor_b[30] = 1; factor_b[31] = 1;
//    factor_b[45] = 1; factor_b[46] = 1; factor_b[47] = 1;


    LSTK = 0;
    for (n = 0; n <= nmax; n++)
    {
        for (m = 0; m <= n; m++)
        {
            if (m == 0 && factor_b[m] !=0) 
            {
                istk[LSTK] = n;
                vstk[LSTK] = GMR2 * factor_a[n] * factor_b[m] * factor_r;
//                vstk[LSTK] = vstk[LSTK] * fnlgdr (0, n, m) * fnlgdr (0, n, m);
                LSTK = LSTK + 1;
            }
            else if (m != 0 && factor_b[m] !=0)
            {
                l = nmax - m + 1;
                ind = nmax + 1 + (2 * nmax - m + 2) * (m - 1);

                istk[LSTK] = ind + n - m;
                istk[LSTK + 1] = ind + n - m + l;
                vstk[LSTK] = GMR2 * factor_a[n] * factor_b[m] * factor_r;
                vstk[LSTK + 1] = GMR2 * factor_a[n] * factor_b[m] * factor_r;
//                vstk[LSTK] = vstk[LSTK] * fnlgdr (0, n, m) * fnlgdr (0, n, m);
//                vstk[LSTK + 1] = vstk[LSTK + 1] * fnlgdr (0, n, m) * fnlgdr (0, n, m);

//      if (m<=n-2)
//      {
//          vstk[LSTK] = 0;
//                  vstk[LSTK + 1] = 0;
//      }
 
                LSTK = LSTK + 2;
            }
        }
    }

//    printf("LSTK = %d\t from degree %d to %d\n", LSTK, factor_b1, factor_b2);
    for (n = 0; n<LSTK; n++)
    {
//        printf("istk = %d  vstk = %e\n", istk[n], vstk[n]);
    }




    for (n = 0; n < solvefor; n++)
    {
        p0[n] = 0;
    }

    for (n=0; n<LSTK; n++)
    {
        p0[istk[n]] = vstk[n];
    }


#pragma omp parallel for private(i,m,ns,is)
    for (n = 0; n < solvefor; n++)
    {
        p0b0[n] = p0[n] * coefaam[n];
        ns = (long)(1.0 * n * solvefor);
        for (i = 0; i < NDATA; i ++)
        {
            is = (long)(1.0 * i * solvefor);
            for (m = 0; m < solvefor; m++)
            {
                atpa[ns + m] +=  ai[is + n] * ai[is + m] * pweight;
            }
            atpy[n] +=  ai[is + n] * obs[i] * pweight;
        }
        atpa[ns + n] +=  p0[n];
        atpy[n] += p0b0[n]; 
    }

    s3 = time(NULL);  
    printf("\n%5ld: seconds of finish ATPA&ATPY\n", s3-s2);
    fflush(stdout);
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
//    for (n = 0; n < solvefor; n++)
//    {
//        atpa[n * solvefor + n] =  atpa[n * solvefor + n] + 1e-12;
//    }

    cest = (double *) calloc ( solvefor, sizeof(double));
    ksi = (double *) calloc ( solvefor, sizeof(double));

//    brinv (atpa,solvefor);  
/////////////////////////////////////////
// inversion may not work for n > 215
/////////////////////////////////////////
    if (slvls == 1)
    {
        printf ("Cholesky Decomposition is used\n");
        solvels_chol(atpa, solvefor, atpy, ksi, 0);
    }
    if (slvls == 2)
    {
        printf ("Gauss Elimination is used\n");
        solvegaus(atpa, solvefor, atpy, ksi);
    }

    s4 = time(NULL);   
    printf("\n%5ld: seconds of inverting N\n", s4-s3);
    fflush(stdout);

//    brmul (atpa,atpy,solvefor,solvefor,1,ksi);

/*
        if (slvls == 1)
            solvels_chol(BTrB, MSOL, BTry, dx, 0);
        if (slvls == 2)
            solvegaus(BTrB, MSOL, BTry, dx);
*/


////////////////////////////////////////////////
//estimate sigma0
//
/*
    brmul(atpy, ksi, 1, solvefor, 1, ctksi);
    brmul(obs, obs, 1, NDATA, 1, ytpy);
    brmul(coefaam, p0b0, 1, solvefor, 1, btpb);
    
    ytpy[0] = ytpy[0] * pweight;
    omega = ytpy[0] + btpb[0]- ctksi[0];
//    printf("omega = %30.20e\t%30.20e\t%30.20e\t%30.20e\n", ytpy[0], btpb[0], ctksi[0], omega);
*/    
    aksi = (double *) calloc ( NDATA, sizeof(double));
//    aks1i = (double *) calloc ( NDATA, sizeof(double));

//    brmul (ai,ksi,NDATA,solvefor,1,aksi);
    for (i = 0; i < NDATA; i++)
    {
        aksi[i] = 0;
        for (n = 0; n < solvefor; n++)
        {
            in = (long)(1.0 * i * solvefor + n);
            aksi[i] = aksi[i] + ai[in] * ksi[n];
        }
    }

//    brmul (a1i,ksi,NDATA,solvegrv,1,aks1i);
   
    etpe = 0;
    for (n=0; n<NDATA; n++)
    {
        etpe = etpe + (obs[n] - aksi[n]) * (obs[n] - aksi[n]);
    }

    etpe = etpe * pweight;
/*
    e0p0e0 = 0;
    for (n=0; n<solvefor; n++)
    {
        e0p0e0 = e0p0e0 + (ksi[n]-coefaam[n]) * (ksi[n]-coefaam[n]) * p0[n];
    }


    omega = etpe + e0p0e0;

    sigma0 = (etpe + e0p0e0 ) / (NDATA - solvefor + solvefor);
*/

    omega = etpe;

    sigma0 = (etpe) / (NDATA - solvefor);
    sigma0 = sqrt(sigma0);

//    printf("etpe = %30.20e\t%30.20e\t%30.20e\n", 
//        etpe, e0p0e0,  sigma0);
//    printf("omega = %30.20e\t%30.20e\t%30.20e\t%30.20e\t%30.20e\n", 
//        ytpy[0], btpb[0], ctksi[0], omega, sigma0);

//////////////////////////////////////////////////////////
//    printf("LFIX = %d\n", LFIX);

//  fixed constrant
//  *kifx, *kfixt, *kn, *knt, *nk, *kcest, *nkknk, *cupd,
    kfix = (double *) calloc ( solvefor * LFIX, sizeof(double));    
    kfixt = (double *) calloc ( solvefor * LFIX, sizeof(double));    
    kn = (double *) calloc ( solvefor * LFIX, sizeof(double));   
    nk = (double *) calloc ( solvefor * LFIX, sizeof(double));   
    nkknk = (double *) calloc ( solvefor * LFIX, sizeof(double));    
    ksiupd = (double *) calloc ( solvefor, sizeof(double));    
    knk = (double *) calloc ( LFIX*LFIX, sizeof(double));    
    knk0 = (double *) calloc ( LFIX*LFIX, sizeof(double));    
    id = (double *) calloc ( LFIX*LFIX, sizeof(double));    
    kksi = (double *) calloc ( LFIX, sizeof(double));    
    lambda = (double *) calloc ( LFIX, sizeof(double)); 




    
    k0 = (double *) calloc ( LFIX, sizeof(double));    



     
    for (n=0; n<LFIX; n++)
    {
        kfix[solvefor * n + ifix[n]] = 1.0;
        k0[n] = vfix[n];
    }
//    kfix[solvefor + 1] = 1.0;
//    kfix[solvefor * 2 + nmax + 1] = 1.0; 
//    kfix[solvefor * 3 + nmax + 1 + nmax] = 1.0;



    mt(kfix, LFIX, solvefor, kfixt);
    brmul (kfix,atpa,LFIX,solvefor,solvefor,kn);
    mt(kn, LFIX, solvefor, nk);
    brmul (kn,kfixt,LFIX,solvefor,LFIX,knk);
    for (n=0; n<LFIX*LFIX; n++)
    {
        knk0[n] = knk[n];
    }
//    brinv (knk,LFIX);   


//    solvels_chol(atpa, solvefor, atpy, ksi, 0);


//    brmul (atpa,kfixt,solvefor,solvefor,LFIX,nk);

    brmul (kfix,ksi,LFIX,solvefor,1,kksi);

    for (n=0; n<LFIX; n++)
    {
        kksi[n] = k0[n] - kksi[n];  //k0=[1 0 0 0]
    }
    
//    solvels_chol(knk, LFIX, kksi, lambda, 0);


    if (slvls == 1)
        solvels_chol(knk, LFIX, kksi, lambda, 0);
    if (slvls == 2)
        solvegaus(knk, LFIX, kksi, lambda);

//      solvegaus(atpa, solvefor, atpy, ksi);
//    brmul (knk,kksi, LFIX,LFIX,1,lambda);

    brmul(knk, knk0, LFIX, LFIX, LFIX, id);
//    for (n=0; n<LFIX; n++)
//    {
//        for (m=0; m<LFIX; m++)
//        {
//            printf("%10.5e ", id[n*LFIX + m]);
//        }
//        printf("\n");
//    }


    brmul (nk,lambda,solvefor,LFIX,1,ksiupd);
    
    for (n = 0; n < solvefor; n++)
    {
        cest[n] =  ksi[n] + ksiupd[n];
//        cest[n] =  ksi[n];
    }

//    brmul (ai,cest,NDATA,solvefor,1,aksi);

    for (i = 0; i < NDATA; i++)
    {
        aksi[i] = 0;
        for (n = 0; n < solvefor; n++)
        {
            in = (long)(1.0 * i * solvefor + n);
            aksi[i] = aksi[i] + ai[in] * cest[n];
        }
    }

    l1cf = (double *) calloc ( NDATA, sizeof(double));
//    brmul (a1i,cest,NDATA,solvegrv,1,aks1i);

    
    for (i = 0; i < NDATA; i++)
    {
        l1cf[i] = 0;
        for (n = 0; n < solvegrv; n++)
        {
            in = (long)(1.0 * i * solvefor + n);
            l1cf[i] = l1cf[i] + ai[in] * cest[n];
        }
    }

//    printf("cest done\n");


/////////////////////////////////////////////
//estimate sigma0
//
    dksiupd = (double *) calloc ( solvefor * solvefor, sizeof(double));
    dksi = (double *) calloc ( solvefor * solvefor, sizeof(double));
 
    brmul (nk,knk,solvefor,LFIX,LFIX,nkknk);
    brmul (nkknk,kn,solvefor,LFIX,solvefor,dksiupd);
    for (n = 0; n < solvefor; n++)
    {
        for (m = 0; m < solvefor; m++)
        {
            dksi[n * solvefor + m] = atpa[n * solvefor + m] - dksiupd[n * solvefor + m];
            if (dksi[n * solvefor + m] < 0) 
                dksi[n * solvefor + m] = 0;
        }
    }


    brmul (kksi,lambda,1,LFIX,1,r0);

//    sigma0 = (omega + r0[0] ) / (NDATA - solvefor + LFIX);
    sigma0 = sqrt ((omega + r0[0] ) / (NDATA - solvefor + LFIX));

//    printf("simga0 = %e\t omega = %e\t r0 = %e\t nml = %d\n", sigma0, omega, r0[0], NDATA - solvefor + LFIX );

////////////////////////////////////////////////////////
//
//
    s5 = time(NULL);   
    printf("\n%5ld: seconds of inverting CONSTRAINT\n", s5-s4);
    fflush(stdout);


/*
    dksi = (double *) calloc ( solvefor * solvefor, sizeof(double));

    for (n = 0; n < solvefor; n++)
    {
        cest[n] =  ksi[n];
    }
    for (n = 0; n < solvefor; n++)
    {
        for (m = 0; m < solvefor; m++)
        {
            dksi[n * solvefor + m] = atpa[n * solvefor + m];
        }
    }


*/
/*
    if (bias == 2)
    {
        efbias = cest[solvefor - 1];
        efbdot = cest[solvefor - 2];
    }
    else if (bias == 1)
    {
        efbias = cest[solvefor - 1];
        efbdot = zero;
    }
    else if (bias == 0)
    {
        efbias = zero;
        efbdot = zero;
    }
*/
    for (m = 0; m <= nmax; m ++)
    {
        l = nmax - m + 1;
        for (k = 0; k < l; k++)
        {
            if (m==0)
            {
                n = k + m;
                fprintf(fpout_coef, 
                    "%4d %4d %23.13e %23.13e %23.13e %23.13e\n", 
                    n, m, cest[k], zero,
                    sigma0 * sqrt(dksi[k * solvefor + k]), zero);
            }
            else
            {
                ind = nmax + 1 + (2 * nmax - m + 2) * (m - 1);
                n = k + m;
                fprintf(fpout_coef, 
                    "%4d %4d %23.13e %23.13e %23.13e %23.13e\n", 
                    n, m, 
                    cest[ind + n - m], 
                    cest[ind + n - m + l],
                    sigma0 * sqrt(dksi[(ind + n - m) * solvefor + ind + n - m]),
                    sigma0 * sqrt(dksi[(ind + n - m + l) * solvefor + ind + n - m + l]));
            }
        }
    }

//    for (n = solvegrv; n < solvefor; n ++)
//    {
//        printf ("%e\n", cest[n]);
//    }


    for (i = 0; i < NDATA; i++)
    {   
        fprintf (fpout_adj,
            "%d\t%20.10e\t%20.10e\t%20.10e\t%20.10e\t%20.10e\t%20.10e\t%20.10e\t%20.10e\n",  
            gpst[i], obs[i] - aksi[i], obs[i], aksi[i], l1cf[i],
            obs[i] - aksi[i] + l1cf[i], csr[i], gfz[i], jpl[i]);
    }

    s6 = time(NULL);   
    printf("\n%5ld: seconds of output data\n", s6-s4);
    fflush(stdout);

    exit(0);
//    fflush(0);

    free (gpst);
    free (ksi);
    free (ksim);
    free (lambda);
    free (lambdam);
    free (ksiupd);
    free (ksiupdm);
    free (dksi);
    free (dksiupd);

    free (ai);
    free (atpa);
    free (atpy);
    free (cest);
    free (coefa);
    free (coefe);
    free (coefaam);
    free (coefjpl);
    free (coefcsr);
    free (coefgfz);
    free (coefggm);
    free (obs);
    fclose(fpin_l1c);
    fclose(fpout_coef);
    fclose(fpout_adj);
   
    ephem_close();  /* remove this line for use with solsys version 2 */

    printf ("\nNormal end of ECHO!\n\npress any key to finish...\n");

    exit(0);

}
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
 


