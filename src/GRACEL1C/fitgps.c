/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/*

    Solar-system Model and Orbit Determination
    SMOD
  
    Version:    1008mars (20 Aug 2010)
    Version:    smod_1201 (1 Jan 2012)

    Copyright (c) 2010 SHANG Kun (shangkun@shao.ac.cn) All Right Reserved 
    Copyright (c) 2012 Kun Shang (shang.34@osu.edu) All Right Reserved 

    NOT FULLY VALIDATED -- SUBJECT TO CORRECTION

*/
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/


/*#define MSDOS */
#define LINUX



#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "eph_manager.h" 
#include "novas.h"
#include "smod.h"
#include "grvts.h"
#include "coord.h"
#include "numrs.h"

//#define orbformat  

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
int main (int argc, char *argv[])
{
    FILE *fp_stdin, *fp_simu, *fp_result, *fp_res, *fp_otbin, *fp_par;
    short int year, month, day, error, de_num;
    int iter, part, all, teout, accid = 0, n, i, j, gps, ndata, slvls, genrel, 
        maxiter, step_ot, step_at, mjds,mjde, indicator;
    double tt, tts, tte, sigma2, sigma, xsim[6], xobs[6], 
        start, finish, duration, xsm[6], xtm[6], nsp, nsv,  
        rmsold, rmsnew, rms2n, rmsv[6] = {0},
        *state, *bi, *bit, *btbi, *btyi, *BTry, *BTrB, *xbi,
        *p, *px, *dx, *x0, *xnew, *xold, weight[6], dyi[1],  
        jd_beg, jd_end, criteria, correl, step_orb, 
        gm[11]= {0};
    char lbl[1], line[400], card[20], f_eop[200], f_eph[200], f_ot[200],
        f_pt[200], f_aot[200], f_obs[200], 
        f_acc[200], f_sca[200], f_grv[2][200]={"\0", "\0"}, stdname[100];  

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
 

// 1AU = 1.5e11m = 0.15e12m = 0.15e15mm = 0.15e18micormeter

   
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
       	getch ();
        exit (0);
	}
    if ( (fp_stdin = fopen (stdname,"r")) == NULL)
    {
        printf ("Cannot open stdin file!\n");
       	getch();
        exit (0);
    }
//    printf ("\nreading input file \"%s\"...\n\n", stdname);


    MOBS = 0;  
    MSRP = 0;  
    MTK2 = 0;  
    MGCS = 0;  

    all = 1;
    teout = 0;
    while (feof(fp_stdin) == 0)
    {
        card[1] = '\0';
        fgets (line, 400, fp_stdin);
        sscanf (line, "%s", card);	

        if (strcmp (card, "planet") ==0)	
        {
            sscanf (line, "%*s%*s%d", &i);
//            i = i - 1;
            sscanf (line, "%*s%*s%*d%hd%lf", &PERB[i], &gm[i]);
        }

        if (strcmp (card, "indicator") ==0)
        {
            sscanf (line, "%*s%d", &indicator);
        }
        if (strcmp (card, "center") ==0)
        {
            sscanf (line, "%*s%d%lf%lf", &CT, &GMCT, &RCT);
//            CT = i - 1;
        }

        if (strcmp (card, "epoch") ==0)
        {
            sscanf (line, "%*s%hd%hd%hd", &year, &month, &day);
        }

        if (strcmp (card, "intstep") ==0)	
        {
            sscanf (line, "%*s%lf", &STEP_OR);
        }
        if (strcmp (card, "ation") ==0)	
        {
            sscanf (line, "%*s%lf%lf", &duration, &step_orb);
        }
        if (strcmp (card, "relativ") ==0)	
        {
   	        sscanf (line, "%*s%d", &genrel);
        }
        if (strcmp (card, "rmsconv") ==0)	
        {
   	        sscanf (line, "%*s%lf", &criteria);
        }
        if (strcmp (card, "itermax") ==0)	
        {
        sscanf (line, "%*s%d", &maxiter);
        }

        if (strcmp (card, "solvels") == 0)
        {
            sscanf (line, "%*s%d", &slvls);
        }

        if (strcmp (card, "gravity") ==0)   
        {
            sscanf (line, "%*s%d%d", &NMAX, &MMAX);
        }
        if (strcmp (card, "amratio") ==0)   
        {
            sscanf (line, "%*s%lf", &AMR);
        }
        if (strcmp (card, "nbody") ==0) 
        {
            sscanf (line, "%*s%d", &NBODY);
        }
        if (strcmp (card, "reltiv") ==0)    
        {
            sscanf (line, "%*s%d", &RELTIV);
        }
        if (strcmp (card, "permtide") ==0)
        {
            sscanf (line, "%*s%d%lf", &PERMT, &C20PERM);
        }
        if (strcmp (card, "bodytide") ==0)
        {
            sscanf (line, "%*s%d", &STIDE);
        }
        if (strcmp (card, "oceantide") ==0)
        {
            sscanf (line, "%*s%d%d%d", &OTIDE, &NOMAX, &step_ot);
        }
        if (strcmp (card, "atmostide") ==0)
        {
            sscanf (line, "%*s%d%d%d", &ATIDE, &NAMAX, &step_at);
        }

        if (strcmp (card, "opoletide") ==0)
        {
            sscanf (line, "%*s%d%d", &PTIDE, &NPMAX);
        }

        if (strcmp (card, "aodmax") ==0)
        {
            sscanf (line, "%*s%d%d", &AOD1B, &NDMAX);
        }

        if (strcmp (card,"ephfile") ==0)
        {
            sscanf (line, "%*s%s", f_eph);
        }
        if (strcmp (card,"aotfile") ==0)
        {
            sscanf (line, "%*s%s", f_aot);
        }
        if (strcmp (card,"eopfile") ==0)
        {
            sscanf (line, "%*s%s", f_eop);
        }
        if (strcmp (card,"grvfile") ==0)
        {
            sscanf (line, "%*s%s%s", f_grv[0], f_grv[1]);
        }
        if (strcmp (card,"oceanfile") ==0)
        {
            sscanf (line, "%*s%s", f_ot);
        }
        if (strcmp (card,"ptidefile") ==0)
        {
            sscanf (line, "%*s%s", f_pt);
        }

        if (strcmp (card,"obsfile") ==0)
        {
            sscanf (line, "%*s%s", f_obs);
        }

        if (strcmp (card,"accfile") ==0)
        {
            sscanf (line, "%*s%s%s", f_acc, lbl);
            if (strcmp (lbl,"A") ==0) accid = 1;
            if (strcmp (lbl,"B") ==0) accid = 2;
//            if (label == 'A') accid = 1;
//            if (label == 'B') accid = 2;
        }

        if (strcmp (card,"scafile") ==0)
        {
            sscanf (line, "%*s%s", f_sca);
        }

        if (strcmp (card, "acc_bias") ==0)
        {
            sscanf (line, "%*s%d%d", &MACC_BIAS, &MACC_DTBS);
        }
        if (strcmp (card, "acc_scal") ==0)
        {
            sscanf (line, "%*s%d%d", &MACC_SCAL, &MACC_DTSL);
        }


    }
    
    fclose(fp_stdin);
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/

    if ((error = ephem_open (f_eph, &jd_beg,&jd_end,&de_num)) != 0)
    {
        if (error == 1)
            printf ("JPL ephemeris file not found.\n");
        else
            printf ("Error reading JPL ephemeris file header.\n");
        getch();
        exit (error);
    }
    else
    {
//    printf ("JPL ephemeris DE%d open. Start JD = %10.2f  End JD = %10.2f\n",
//        de_num, jd_beg, jd_end);
//    printf ("\n");
    }


    start = clock();  
    if((fp_result=fopen("result.txt","w"))==NULL)
    {
        printf("Cannot write result.txt!\n");
      	getch();
        exit(0);
    }
    if((fp_simu=fopen(f_obs,"r"))==NULL)
    {
        printf("Cannot open obs data?\n");
      	getch();
        exit(0);
    }
    if((fp_res=fopen("forb.dat","w"))==NULL)
    {
        printf("Cannot write residua.dat!\n");
       	getch();
        exit(0);
    }

    if((fp_par=fopen("fpar.dat","w"))==NULL)
    {
        printf("Cannot write fpar.dat!\n");
        getch();
        exit(0);
    }


    JD0 = julian_date (year, month, day,0);

    AC_EPH  = (double *) calloc (86400 * 4, sizeof(double));
    SC_EPH  = (double *) calloc (17280 * 5, sizeof(double));

    readacc1 (f_acc);
    readsca1 (f_sca);

    cal_acc_1(accid);
//    cal2_acc();



/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/

    for (n = 0; n <= 10; n++)
        GMDE[n] = gm[n];


    if (fgets (line, 400, fp_simu) ==NULL) exit(0);
    sscanf(line, "%d%*s%*s%lf%lf%lf%lf%lf%lf", 
        &gps, &xsm[0], &xsm[1], &xsm[2], &xsm[3], &xsm[4], &xsm[5]);


//    utcs = gps2utc(gps);
//    utce = utcs + 86400.0;
//    tts = gps2tt(gps);
    tts = gps_GRACE2tt (JD0, gps);
    tte = tts + 86400.0;

    TT0 = tts;
//    mjds = (int)(JD0 - 2400000.5 + tts/86400.0) - 2;
//    NEOP = (int)(JD0 - 2400000.5 + tte/86400.0 - mjds) + 3;

//    EOPMT  = (double *) calloc (NEOP * 6, sizeof(double));
//    openeop (f_eop, mjds, NEOP, EOPMT);

    mjds = (int)(JD0 - 2400000.5);
    mjde = mjds + 1;
    eop_open (f_eop, mjds, mjde);


    if (NMAX >= 2)
    {
        COEF  = (double *) calloc ((NMAX + 1) * (NMAX + 1), sizeof(double));
        COEFG  = (double *) calloc ((NMAX + 1) * (NMAX + 1), sizeof(double));
//        opengrv (f_grv, COEFG, NMAX, MMAX);

        grv_open (f_grv[0], f_grv[1], NMAX, MMAX, COEFG);
        zero0zero1 (COEFG, NMAX);  

        NSMAX = 2;
        if (CT == 2 && STIDE > 1)
            NSMAX = 4;
        COEFS = (double *) calloc ( (NSMAX + 1) * (NSMAX + 1), sizeof(double));

//        DIM_OT = 17280;
        DIM_OT = (int)(86400/step_ot);
        OT_EPH  = (double *) calloc ( ( (NOMAX + 1) * (NOMAX + 1) + 1 ) * DIM_OT , sizeof(double));
        COEFO = (double *) calloc ( (NOMAX + 1) * (NOMAX + 1), sizeof(double));

        DIM_AT = (int)(86400/step_at);
        AT_EPH  = (double *) calloc ( ( (NAMAX + 1) * (NAMAX + 1) + 1 ) * DIM_AT , sizeof(double));
        COEFA = (double *) calloc ( (NAMAX + 1) * (NAMAX + 1), sizeof(double));
/*
        if (OTIDE == 1 || OTIDE == 2)
            openotcs_fes (f_ot);
        if (OTIDE == 3 || OTIDE == 4)
            openotcs_csr (f_ot);
        openoteph (tts, tte);
*/
        AOD_EPH  = (double *) calloc ( ( (NDMAX + 1) * (NDMAX + 1) + 1 ) * 5 , sizeof(double));
        COEFD  = (double *) calloc ( (NDMAX + 1) * (NDMAX + 1) , sizeof(double));
//        openaod (f_aod, NDMAX, NDMAX);



        OPTM1 = (double *) calloc ( (NPMAX + 1) * (NPMAX + 1), sizeof(double));
        OPTM2 = (double *) calloc ( (NPMAX + 1) * (NPMAX + 1), sizeof(double));
        COEFP = (double *) calloc ( (NPMAX + 1) * (NPMAX + 1), sizeof(double));
//        openopt (f_pt, NPMAX);

    }

    if((fp_otbin=fopen(f_aot,"rb"))==NULL)
    {
        printf("Cannot write oteph.bin!\n");
        getch();
        exit(0);
    }

/*
    for (i = 0; i < DIM_OT; i++)
    {
        for (n = 0; n < (NOMAX + 1) * (NOMAX + 1) + 1; n++)
        {
            fprintf (fp_otbin, "%e\t",   OT_EPH[i * ((NOMAX + 1) * (NOMAX + 1) + 1) + n]);
        }
        fprintf (fp_otbin, "\n");
    }
*/

    fread (OT_EPH, sizeof(double) * ( (NOMAX + 1) * (NOMAX + 1) + 1 ) * DIM_OT, 1, fp_otbin);

//    if((fp_atbin=fopen("ateph.bin","rb"))==NULL)
//    {
//        printf("Cannot write oteph.bin!\n");
//        getch();
//        exit(0);
//    }

//    fread (AT_EPH, sizeof(double) * ( (NAMAX + 1) * (NAMAX + 1) + 1 ) * DIM_AT, 1, fp_atbin);


//    exit(0);

//    lgr_order (OT_EPH, DIM_OT, (NOMAX + 1) * (NOMAX + 1) + 1 , tts, COEFO, 1);
//    lgr_order (OT_EPH, DIM_OT, (NOMAX + 1) * (NOMAX + 1) + 1 , tts + 60, COEFO, 1);
//    lgr_order (OT_EPH, DIM_OT, (NOMAX + 1) * (NOMAX + 1) + 1 , tts + 30, COEFO, 1);
//    lgr_order (AOD_EPH, 5, (NAMAX + 1) * (NAMAX + 1) + 1 , 0.125, COEFA, 1);

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
  
    getsolvefor();
  
    xbi = (double *) calloc (6 * MSOL, sizeof(double));
    bi = (double *) calloc ( MSOL, sizeof(double));
    bit = (double *) calloc ( MSOL, sizeof(double));
    btbi = (double *) calloc ( MSOL * MSOL, sizeof(double));
    btyi = (double *) calloc ( MSOL, sizeof(double));
    BTry = (double *) calloc ( MSOL, sizeof(double));
    BTrB = (double *) calloc ( MSOL * MSOL, sizeof(double));
    p = (double *) calloc ( MSOL * MSOL, sizeof(double));
    px = (double *) calloc ( MSOL, sizeof(double));
    dx = (double *) calloc ( MSOL, sizeof(double));
    x0 = (double *) calloc ( MSOL, sizeof(double));
    xnew = (double *) calloc ( MSOL, sizeof(double));
    xold = (double *) calloc ( MSOL, sizeof(double));
    state = (double *) calloc ( MSTA, sizeof(double));


    initsolvefor (xsm, x0);

//    printf ("%d\n", MACC);
//    printf ("%f\n", TT0);

    for (n = 0; n < MSOL; n ++)
        xnew[n] = x0[n];
    
//        updsolvefor (xnew);            
//        initsolvefor (xsm, xnew);

    rmsold = 2e10;
    rmsnew = 1e10;


    DIM_OR = (int)( (tte - tts) / STEP_OR) + 1;	
    OR_EPH = (double *) calloc ( DIM_OR * (MSTA + 1) , sizeof(double));	
        
    finish = clock();  
    duration = (double)(finish - start) / CLOCKS_PER_SEC; 
    fprintf (fp_result, "\nOD: Before estimation of orbit: \t%.4f seconds\n", 
        duration);

    for (iter = 0; iter < maxiter; iter ++)
    {

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/* initialize transition matrix and orbit trajectory, */
        rmsold = rmsnew;
        for (n = 0; n < MSOL; n ++)
        {
            xold[n] = xnew[n];
            dx[n] = 0;
        }
        for (n = 0; n < MSTA; n++)
            state[n] = 0;
        for (n = 0; n < 6; n++)
  	        state[n] = xold[n];
        for (n = 0; n < 6; n++)
            state[n * 6 + n + 6] = 1;
        
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/* integrate transition matrix and orbit trajectory*/
        
        for (i = 0; i < DIM_OR; i++)
        {
            tt = tts + STEP_OR * i;
//            lps = getlps (JD0 + utcs/86400.0);
//            tt = utc + (lps + 32.184);

            OR_EPH[i * (MSTA + 1)] = tt;
            for (n = 0; n < MSTA; n++)
                OR_EPH[i * (MSTA + 1) + 1 + n] = state[n]; 

            rkf78 (JD0, tt, STEP_OR, state, MSTA, fun_accel);

        }

        
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/* read observation and form BtB & Bty matrix*/

        for (n = 0; n < 6; n++)
            rmsv[n] = 0;
        rms2n = 0;
        sigma2 = 0;
        ndata = 0;
        for (n = 0; n < MSOL * MSOL; n++)
        {
            BTrB[n] =0;
        }
        for (n = 0; n < MSOL; n++)
        {
            BTry[n] = 0;
        }
        rewind(fp_simu);
        while (1)
        {
            if (fgets (line, 400, fp_simu) ==NULL) break;

            sscanf(line, "%d%*s%*s%lf%lf%lf%lf%lf%lf", 
                &gps, &xobs[0], &xobs[1], &xobs[2], 
                &xobs[3], &xobs[4], &xobs[5]);

//            tt = gps2tt(gps);
            tt = gps_GRACE2tt (JD0, gps);
            if (tt > tte) 
                continue;

//            weight = 1.0 / noise / noise; //range: m; velocity: m/s
            nsp = 3e-2;
            nsv = 3e-5;
            weight[0] = 1/nsp/nsp;              
            weight[1] = 1/nsp/nsp;              
            weight[2] = 1/nsp/nsp;              
            weight[3] = 1/nsv/nsv;              
            weight[4] = 1/nsv/nsv;              
            weight[5] = 1/nsv/nsv;              
         
            part = 1;
                    
            obs_mex (tt, xsim, part, xbi);

            for (i = 0; i < 6; i++)
            {
                dyi[0] = xobs[i]  - xsim[i];
                for (j = 0; j < MSOL; j ++)    
                    bi[j] = xbi[i * MSOL + j];

                rmsv[i] = rmsv[i] + dyi[0] * dyi[0];

                rms2n = rms2n + dyi[0] * dyi[0];
                sigma2 = sigma2 + dyi[0] * dyi[0] * weight[i];

                mt (bi, 1, MSOL, bit);
                brmul(bit,bi,MSOL,1,MSOL,btbi);
                brmul(bit,dyi,MSOL,1,1,btyi);
            
                for (n = 0; n < MSOL * MSOL; n++)
                {
                    BTrB[n] = BTrB[n] + btbi[n] * weight[i];
                }
                for (n = 0; n < MSOL; n++)
                {
                    BTry[n] = BTry[n] + btyi[n] * weight[i];
                }
                ndata++;
            }
        }
        rmsnew = sqrt (rms2n / ndata);
        sigma = sqrt (sigma2 / (ndata - MSOL));

        for (i = 0; i < 6; i++)
            rmsv[i] = sqrt (rmsv[i] / ndata);

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/* solve least-square */
        for (n = 0; n < MSOL; n++)
        {
            for (i = 0; i < MSOL; i++)
            {
                fprintf (fp_result, "%12.4e ", BTrB[n * MSOL + i]);
            }
            fprintf (fp_result, "%12.4e\n", BTry[n]);
        }
        fprintf (fp_result, "\n");


/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/

        if (slvls == 1)
            solvels_chol(BTrB, MSOL, BTry, dx, 0);
        if (slvls == 2)            
            solvegaus(BTrB, MSOL, BTry, dx);
   
        for (n = 0; n < MSOL; n++)
            xnew[n] = xold[n] + dx[n];
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/

        updsolvefor (xnew);            

        for (n = 0; n < MSOL; n++)
        {
            fprintf (fp_result, "%3d: ", n);
            for (i = 0; i <= n; i++)
            {
                correl = BTrB[i * MSOL + n] / 
                    sqrt (BTrB[n * MSOL + n] * BTrB[i * MSOL + i]);
                fprintf (fp_result, "%8.4f ", correl);
            }
            fprintf (fp_result, "\n");
        }
        fprintf (fp_result, "\n");
        
        for (n = 0; n < MSOL; n++)
        {
            for (i = 0; i <= n; i++)
            {
                fprintf (fp_result, "%12.4e ", BTrB[i*MSOL+n]);
            }
            fprintf (fp_result, "\n");

        }
        fprintf (fp_result, "\n");

        for (n = 0; n < MSOL; n++)
        {
            fprintf (fp_result, "%22.10f %22.10f %22.10f %22.10f %22.10f %22.10f\n", 
                x0[n], xold[n], xnew[n], x0[n] - xnew[n], dx[n], sqrt(BTrB[n*MSOL+n]));
        }
        finish = clock();  
        duration = (double)(finish - start) / CLOCKS_PER_SEC; 
        fprintf (fp_result, "\nOD: Finish estimation of orbit: \t%.4f seconds\n", 
            duration);
        fprintf (fp_result, "iter: %2d\t rms: %12.6f\t%8.2f%%\t sigma = %.4f\n\n\n", 
            iter, rmsnew, fabs((rmsold - rmsnew) / rmsold) * 100, sigma);
        fflush(fp_result);



        if ( fabs((rmsold - rmsnew) / rmsold) < criteria || iter == maxiter - 1)
        {
//            pt_orb (utcs, utce, 5, MSTA + 1);

            for (n = 0; n < 6; n++)
                xtm[n] = xnew[n];

            if (maxiter == 1)
            {
                for (n = 0; n < 6; n++)
                    xtm[n] = xold[n];
                if (MACC_BIAS != 0)
                    for (n = 0; n < MACC_PRBS; n ++)
                        MACC_EPBS[n] = 0;
                if (MACC_SCAL != 0)
                    for (n = 0; n < MACC_PRSL; n ++)
                        MACC_EPSL[n] = 0;
            }

            STEP_OR = 5;
            DIM_OR = (int)(86400.0 / STEP_OR);
            for (i = 0; i < DIM_OR; i++)
            {
                tt = tts + STEP_OR * i;
//                lps = getlps (JD0 + utcs/86400.0);
//                tt = utc + (lps + 32.184);

//                gps = tt - 32.184 - 19;
                fprintf (fp_res, "%d %s I %20.10f %20.10f %20.10f ",
                    (int) ((JD0 - T0) * 86400 + tt - 32.184 - 19), lbl, xtm[0], xtm[1], xtm[2]);
                fprintf (fp_res, "%20.14f %20.14f %20.14f %20.10f\n",
                    xtm[3], xtm[4], xtm[5], modvect(xtm));

                rkf78 (JD0, tt, STEP_OR, xtm, 6, fun_accel);

            }


            break;
        }

    }


    finish = clock();  
    duration = (double)(finish - start) / CLOCKS_PER_SEC; 
    fprintf (fp_result, "OD: Finish print orbit: \t%.4f seconds\n\n", 
            duration);
    
    fprintf (fp_result, "FITERR %10.2f %4d %2d %2d %12.8f %12.8f %12.8f %14.10f %14.10f %14.10f %12.8f %14.10f\n", 
        JD0, year, month, day,  rmsv[0], rmsv[1], rmsv[2], rmsv[3], rmsv[4], rmsv[5],
        modvect(&rmsv[0]), modvect(&rmsv[3]));
    fflush(fp_result);


    if (MACC_BIAS != 0)
    {
        for (i = 0; i < MACC_ARBS; i ++)
        {
            fprintf (fp_par, "BIAS  %d\t ", i * MACC_DTBS);
            for (n = 0; n < MACC_NOBS; n ++)
            {
                fprintf (fp_par, "%e\t ", MACC_EPBS[i * MACC_NOBS + n]);
            }
            fprintf (fp_par, "\n");
        }
    }


    if (MACC_SCAL != 0)
    {

        for (i = 0; i < MACC_ARSL; i ++)
        {
            fprintf (fp_par, "SCAL  %d\t ", i * MACC_DTSL);
            for (n = 0; n < MACC_NOSL; n ++)
            {
                fprintf (fp_par, "%e\t ", MACC_EPSL[i * MACC_NOSL + n]);
            }
            fprintf (fp_par, "\n");
        }

    }


//    for (n = 6; n < MSOL; n++)
//        printf ("%e\t %e\n",xnew[n], sqrt(BTrB[n*MSOL+n]));

//    exit(0);
    ephem_close();
    fclose(fp_result);
    fclose(fp_simu);
    fclose(fp_res);
    fclose(fp_par);

    free (AC_EPH);
    free (SC_EPH);
    free (OR_EPH); 
    free (EOPMT);
    free (COEF);
    free (COEFG);
    free (COEFS);
    free (OT_EPH);
    free (COEFO);
    free (AOD_EPH);
    free (COEFA);
    free (COEFD);


    free (state);
    free (xbi);
    free (bi); 
    free (bit); 
    free (btbi); 
    free (btyi); 
    free (BTry);
    free (BTrB); 
    free (p); 
    free (px); 
    free (dx); 
    free (x0); 
    free (xnew); 
    free (xold);
    exit(0);

}
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
