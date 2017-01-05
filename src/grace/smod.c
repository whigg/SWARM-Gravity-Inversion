
#ifndef _SMOD_H_
    #include "smod.h"
#endif
#ifndef _COORD_H_
    #include "coord.h"
#endif
#ifndef _GRVTS_H_
    #include "grvts.h"
#endif

#ifndef _l1b2l1c_H_
    #include "l1b2l1c.h"
#endif












void pt_orb (double ts_orb, double te_orb, double step_orb, int dim_eph)
{
    FILE *fp_fxyz, *fp_faei, *fp_frtn, *fp_fllh;
    int i, n;
    double tt, lps, utc, xtm[6], *eph, dist, velt, tp, rtn_p[3], rtn_v[3], 
        ele[6], llh[3], gps;
    double gm=398600.44150e+09;
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
    /*--print orbit--*/
    if((fp_fxyz=fopen("forb.xyz","w"))==NULL)
    {
        printf("Cannot write fort.xyz!\n");
        exit(0);
    }

    if((fp_faei=fopen("forb.aei","w"))==NULL)
    {
        printf("Cannot write fort.aei!\n");
        exit(0);
    }

    if((fp_frtn=fopen("forb.rtn","w"))==NULL)
    {
        printf("Cannot write fort.rtn!\n");
        exit(0);
    }
    if((fp_fllh=fopen("forb.llh","w"))==NULL)
    {
        printf("Cannot write fort.llh!\n");
        exit(0);
    }

        
    eph = (double *) calloc (dim_eph - 1, sizeof(double));

    i = 0;
    for (utc = ts_orb; utc <= te_orb; utc = ts_orb + step_orb * i)
    {
        lps = getlps (JD0 + ts_orb/86400.0);
        tt = utc + (lps + 32.184);

        lagrange (OR_EPH, DIM_OR, dim_eph, tt, eph);
        for (n = 0; n < 6; n++)
            xtm[n] = eph[n];


        dist = sqrt (xtm[0] * xtm[0] + xtm[1] * xtm[1] + xtm[2] * xtm[2]);
        velt = sqrt (xtm[3] * xtm[3] + xtm[4] * xtm[4] + xtm[5] * xtm[5]);

//        tp = xyz2aei(ele, &xtm[0], &xtm[3]);
        tp = xyz2aei(&xtm[0], &xtm[3], gm, ele);

        xyz2rtn(&xtm[0], &xtm[3], &xtm[0], rtn_p);
        xyz2rtn(&xtm[0], &xtm[3], &xtm[3], rtn_v);
    
        xyz2llr(xtm, llh);

        fprintf (fp_fxyz, "%14.4f  %14.6f  %26.14f  %26.14f  %26.14f  ",
            JD0, utc, xtm[0], xtm[1], xtm[2]);
        fprintf (fp_fxyz, "%24.16f  %24.16f  %24.16f  %16.4f  %14.6f \n", 
            xtm[3], xtm[4], xtm[5], dist, velt);

        fprintf (fp_faei, "%14.4f  %14.6f  %26.14f  %10.6f  %12.6f  ",
            JD0, utc, ele[0], ele[1], ele[2]);
        fprintf (fp_faei, "%12.6f  %12.6f  %12.4f  %12.4f \n", 
            ele[3], ele[4], ele[5], tp);

        fprintf (fp_frtn, "%14.4f  %14.6f  %16.4f  %16.4f  %16.4f  ",
            JD0, utc, rtn_p[0], rtn_p[1], rtn_p[2]);
        fprintf (fp_frtn, "%14.6f  %14.6f  %14.6f  %16.4f  %14.6f \n", 
            rtn_v[0], rtn_v[1], rtn_v[2], dist, velt);

        gps = tt - 32.184 - 19;
        fprintf (fp_fllh, "%d X I %25.15f %25.15f %25.15f  ",
           (int) ((JD0 - T0) * 86400 + gps), xtm[0], xtm[1], xtm[2]);
        fprintf (fp_fllh, "%25.18f %25.18f %25.18f ",
            xtm[3], xtm[4], xtm[5]);
        fprintf (fp_fllh, "%20.15f  %20.15f  %25.15f \n",
            llh[0], llh[1], llh[2]);

//        fprintf (fp_fllh, "%14.4f  %14.6f  %26.14f  %26.14f  %26.14f\n",
//            JD0, utc, llh[0], llh[1], llh[2] - RCT);
        i++;

    }
    fclose(fp_fxyz);
    fclose(fp_faei);
    fclose(fp_frtn);
    fclose(fp_fllh);
    free (eph);

    return;
}





    

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/*
 * geteop - interpolation of eop
 * mjd: double, input MJD
 * xp, yp, ut1_utc, dx, dy: output EOP
 
 * http://hpiers.obspm.fr/iers/eop/eopc04_05/eopc04_IAU2000.62-now
 * version: 20 Aug 2010
 */
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
void geteop (double mjd, double *xp, double *yp, 
               double *ut1_utc, double *dx, double *dy)
{
    int i, mjdi = 0;
    double x1 = 0, y1 = 0, dt1 = 0, dx1 = 0, dy1 = 0, 
           x2 = 0, y2 = 0, dt2 = 0, dx2 = 0, dy2 = 0; 

            
    for (i = 0; i < NEOP; i ++)
    {
        mjdi = (int)EOPMT[i * 6 + 0];

        if (mjdi == (int)mjd)
        {
            x1  = EOPMT[i * 6 + 1];
            y1  = EOPMT[i * 6 + 2];
            dt1 = EOPMT[i * 6 + 3];
            dx1 = EOPMT[i * 6 + 4];
            dy1 = EOPMT[i * 6 + 5];
            i++;
            x2  = EOPMT[i * 6 + 1];
            y2  = EOPMT[i * 6 + 2];
            dt2 = EOPMT[i * 6 + 3];
            dx2 = EOPMT[i * 6 + 4];
            dy2 = EOPMT[i * 6 + 5];

            break;
        }
    }    

    *xp      = x1 + (x2-x1) * (mjd - mjdi);
    *yp      = y1 + (y2-y1) * (mjd - mjdi);
    *ut1_utc = dt1 + (dt2-dt1) * (mjd - mjdi);
    *dx      = dx1 + (dx2-dx1) * (mjd - mjdi);
    *dy      = dy1 + (dy2-dy1) * (mjd - mjdi);

}








/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/*
* opengravfile ¨C 
* @param1: description of param1
* @param2: description of param2
*
* version: 20 Aug 2010
*/
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
double opengrv (char file_grv[2][200], double *coef, int nmax, int mmax)
{
    FILE *fp_grv;
    double c,s;
    int n,m, l, ind;
    char string[200], name[20];

    if ((fp_grv = fopen (file_grv[0],"r")) == NULL)
    {
        printf ("Cannot open gravity file?\n");
        exit (0);
    }


//    coef[0] = 1;    // include zero degree term
    coef[0] = 0;    // exclude zero degree term

    while (1)
    {
        if (fgets (string, 200, fp_grv) == NULL) break;
//        sscanf (string, "%d%d%lf%lf", &n, &m, &c, &s);    
        n = 0; m = 0;
        if (strlen(file_grv[1])==0)
        {
            sscanf (string, "%d%d%lf%lf", &n, &m, &c, &s);  
        }
        else 
        {
            sscanf (string, "%s", name);    
            if (strcmp (name,file_grv[1]) ==0)  
            {
                sscanf (string, "%*s%d%d%lf%lf", &n, &m, &c, &s);
//                printf ("n = %d m = %d c = %e s = %e\n", n, m, c, s);
            }
        }


//        if (n > nmax || n < 0)
        if (n > nmax || n < 2 || m > mmax)   // permanently exclude degree 1 @7/24/2012
            continue;
        else if (m == 0)
        {
            coef[n] = c;
        }
        else 
        {
            l = nmax - m + 1;
            ind = nmax + 1 + (2 * nmax - m + 2) * (m - 1);
            coef[ind + n - m] = c;
            coef[ind + n - m + l] = s;
//            coef[ind + n - m] = 0;
//            coef[ind + n - m + l] = 0;
        }
    }

//    printf ("coef[2] = %e\n", coef[2]);
    if (PERMT == 1)
    {
        coef[2]  = coef[2] + C20PERM; //tn32
//        coef[2]  = coef[2] - 4.1736e-9;  //tn36
    }

//    printf ("coef[2] = %e\n", coef[2]);
    fclose(fp_grv);
    return 0;
}












/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/*
* xyz2llh - xyz to latitude, longitude, height
* @param1: description of param1
* @param2: description of param2

* version: 20 Aug 2010
*/
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
void xyz2llh (double *vt, double *llh)
{
    double r, cosphi, phi, costhe, sinthe;
    r = sqrt (vt[0] * vt[0] + vt[1] * vt[1] + vt[2] * vt[2]);

    cosphi = vt[2] / r;
    phi = acos(cosphi) ;
    costhe = vt[0] / r / sin(phi);
    sinthe = vt[1] / r / sin(phi);
    llh[2] = r;
    llh[1] = chosephase(sinthe, costhe) * RAD2DEG;
    llh[0] = 90.0 - phi * RAD2DEG;
}





/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/*
* fun_pointmass - abandoned
* @param1: description of param1
* @param2: description of param2
*
* version: 20 Aug 2010
*/
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/

double earth_pointmass (double jd, double tdbs, double *x, double *f)
{
    double GM, radius, gmde, fnt[3], fgr[3], r, s2, rrd, a, b;
    int n, gamma;

    GM = GMCT;
    radius = RCT;

    gmde = GM * 86400.0 * 86400.0 / AU / AU / AU;
    gamma = 1;

    f[0] = x[3]; 
    f[1] = x[4]; 
    f[2] = x[5];

    r = sqrt (x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
    s2 = x[3] * x[3] + x[4] * x[4] + x[5] * x[5];
    rrd = x[0] * x[3] + x[1] * x[4] + x[2] * x[5];
	
    a = 2 * (1 + gamma) * gmde / r - gamma * s2;
    b = 2 * (1 + gamma) * rrd;
    for (n = 0; n < 3; n++)
        fgr[n] =  gmde / C_AUDAY / C_AUDAY / r / r / r 
        * ( a * x[n] + b * x[n+3] );

    fnt[0] = - gmde / (r*r*r) * x[0];
    fnt[1] = - gmde / (r*r*r) * x[1];
    fnt[2] = - gmde / (r*r*r) * x[2];

    for (n = 0; n < 3; n++)
    {
        f[3 + n] = fnt[n]; 
    }

	return 0;
}



double accel_point (double *tjd, double *x, double *fnt, double *fgr)
{
    double GM, radius, gmde, r, s2, rrd, a, b;
    int n, gamma;

    GM = GMCT;
    radius = RCT;

    gmde = GM * 86400.0 * 86400.0 / AU / AU / AU;
    gamma = 1;

    r = sqrt (x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
    s2 = x[3] * x[3] + x[4] * x[4] + x[5] * x[5];
    rrd = x[0] * x[3] + x[1] * x[4] + x[2] * x[5];
	
    a = 2 * (1 + gamma) * gmde / r - gamma * s2;
    b = 2 * (1 + gamma) * rrd;
    for (n = 0; n < 3; n++)
        fgr[n] =  gmde / C_AUDAY / C_AUDAY / r / r / r 
        * ( a * x[n] + b * x[n+3] );

    fnt[0] = - gmde / (r*r*r) * x[0];
    fnt[1] = - gmde / (r*r*r) * x[1];
    fnt[2] = - gmde / (r*r*r) * x[2];

	return 0;
}












double accel_pmiers (double *tjd, double *x, double *fnt, double *fgr)
{
    double GME, GMS, J, Jv[3], beta, gamma, r, v2, pv, pJ, a, b, p[3], v[3],
        pxv[3], vxJ[3], ps[3], vs[3], rs, vsxps[3], vsxpsxv[3],
        term1[3], term2[3], term3[3];
    int n;
    short int sun = 10;

    GME = GMCT; //m^3/s^2
    J = 9.8e8; //m^2/s
    gamma = 1;
    beta = 1;
    GMS = 1.32712442076e20; //m^3/s^2

    GME = GME * 86400.0 * 86400.0 / AU / AU / AU;
    GMS = GMS * 86400.0 * 86400.0 / AU / AU / AU;
    J = J * 86400.0 / AU / AU; 
    
    Jv[0] = 0; Jv[1] = 0; Jv[2] = J;
    p[0] = x[0]; p[1] = x[1]; p[2] = x[2];
    v[0] = x[3]; v[1] = x[4]; v[2] = x[5];

    r = modvect(p);
    v2 = dotvect(v, v);
    pv = dotvect(p, v);
    pJ = dotvect(p, Jv);
    crsvect(p, v, pxv);
    crsvect(v, Jv, vxJ);    

    planet_ephemeris (tjd, CT, sun, ps, vs);
    rs = modvect(ps);
    crsvect(vs, ps, vsxps);
    crsvect(vsxps, v, vsxpsxv);    

    a = 2 * (beta + gamma) * GME / r - gamma * v2;
    b = 2 * (1 + gamma) * pv;
    for (n = 0; n < 3; n++)
    {
        term1[n] = GME / C_AUDAY / C_AUDAY / r / r / r *
                ( a * p[n] + b * v[n] );
        term2[n] = GME / C_AUDAY / C_AUDAY / r / r / r * (1 + gamma) * 
                ( 3/r/r * pxv[n] * pJ + vxJ[n] );
        term3[n] = - GMS / C_AUDAY / C_AUDAY / rs / rs / rs * (1 + 2 * gamma) *
                vsxpsxv[n];

        fgr[n] = term1[n]
                + term2[n] + term3[n];
//        printf ("%15.12f\t%15.12f\t%15.12f\n", term1[n],term2[n],term2[n]);
    }


    fnt[0] = - GME / (r*r*r) * p[0];
    fnt[1] = - GME / (r*r*r) * p[1];
    fnt[2] = - GME / (r*r*r) * p[2];

    return 0;
}




double accel_ac_part (double *tjd, double *xic, double *acc)
{
    double tt, as[3], as0[3], qvec[4], c_is[9], c_si[9];
    int k, i, m, lbs = 0, lsl = 0;

    tt = tjd[1] * 86400.0;

//    lagrange (AC_EPH, DIM_AC, 4, tt, as);
//    lagrange (SC_EPH, DIM_SC, 5, tt, qvec);
//    lagrangelow (AC_EPH, DIM_AC, 4, tt, as0);
//    lagrangelow (SC_EPH, DIM_SC, 5, tt, qvec);

    lgr_order (AC_EPH, DIM_AC, 4, tt, as0, 4);
    lgr_order (SC_EPH, DIM_SC, 5, tt, qvec, 4);

    quat2mat_i2s (qvec, c_is);
    mt(c_is, 3, 3, c_si);

//    as[0] = as[0] + AX0 + AX1 * tjd[1];
//    as[1] = as[1] + AY0 + AY1 * tjd[1];
//    as[2] = as[2] + AZ0 + AZ1 * tjd[1];
    

    for (k = 0; k < 3; k ++)
    {
        as[k] = as0[k];
    }

    if (MACC_SCAL != 0)
    {
        lsl = (int)((tt - TT0)/MACC_DTSL);
        if (lsl < 0) lsl = 0;
        if (lsl > MACC_ARSL - 1) lsl = MACC_ARSL - 1;
        for (k = 0; k < 3; k ++)
        {
            as[k] = as0[k] * MACC_EPSL[lsl * MACC_NOSL + k];
        }
    }

    if (MACC_BIAS != 0)
    {
        lbs = (int)((tt - TT0)/MACC_DTBS);
        if (lbs < 0) lbs = 0;
        if (lbs > MACC_ARBS - 1) lbs = MACC_ARBS - 1;
        for (k = 0; k < 3; k ++)
        {
            for (i = 0; i < MACC_BIAS; i ++)
            {
                as[k] = as[k] + MACC_EPBS[lbs * MACC_NOBS + k + 3 * i] * pow (tjd[1], i);
            }
        }
    }
    
    brmul(c_si, as, 3,3,1, acc);


    if (MACC_BIAS != 0)
    {
        for (k = 0; k < 3 * MACC_PRBS; k ++)
        {
            DADBS[k] = 0;
        }

        for (k = 0; k < 3; k ++)
        {
            for (i = 0; i < MACC_BIAS; i ++)
            {
                for (m = 0; m < 3; m ++)
                {
                    DADBS[k * MACC_PRBS + lbs * MACC_NOBS + i * 3 + m] = c_si[k * 3 + m] * pow (tjd[1], i);
                }
            }
        }
    }

    if (MACC_SCAL != 0)
    {
        for (k = 0; k < 3 * MACC_PRSL; k ++)
        {
            DADSL[k] = 0;
        }

        for (k = 0; k < 3; k ++)
        {
            for (m = 0; m < 3; m ++)
            {
                DADSL[k * MACC_PRSL + lsl * MACC_NOSL + m] = c_si[k * 3 + m] * as0[m];
            }
        }
    }






    return 0;

}



/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
void readsca1 (char *infile)
{
    FILE *fpSCA;
    int i, gps;
    char line[200];


/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
    if ( (fpSCA = fopen (infile,"r")) == NULL)
    {
        printf ("Cannot open sca file!\n");
        exit (0);
    }

    i = 0;
    while (1)
    {
        if (fgets(line,200, fpSCA) ==NULL) break;
//        sscanf (line, "%d%*s%*d%lf%lf%lf%lf",
        sscanf (line, "%d%lf%lf%lf%lf",
            &gps, &SC_EPH[i * 5 + 1], &SC_EPH[i * 5 + 2], 
            &SC_EPH[i * 5 + 3], &SC_EPH[i * 5 + 4]);
        SC_EPH[i * 5] = gps - (JD0 - T0) * 86400.0 + 19 + 32.184;
        i++;
    }


    fclose(fpSCA);

    DIM_SC = i;

    return;
}



void readacc1 (char *infile)
{
    FILE *fpACC;
    int n_acca, i, gps;
    char line[400];

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
    if ( (fpACC = fopen (infile,"r")) == NULL)
    {
        printf ("Cannot open acc file!\n");
        exit (0);
    }


    i = 0;
    while (1)
    {
        if (fgets(line,400, fpACC) ==NULL) break;
        n_acca ++;
//        sscanf (line, "%d%*s%lf%lf%lf",
        sscanf (line, "%d%lf%lf%lf",
            &gps, &AC_EPH[i * 4 + 1], &AC_EPH[i * 4 + 2], &AC_EPH[i * 4 + 3]);
        AC_EPH[i * 4] = gps - (JD0 - T0) * 86400.0 + 19 + 32.184;
        i++;

    }


    fclose(fpACC);

    DIM_AC = i;

    return;

}















/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
double cal_acc_1(int accid)
{

    double scale_x_A, scale_y_A, scale_z_A, scale_x_B, scale_y_B, scale_z_B,
        bias_x_A, bias_y_A, bias_z_A, bias_x_B, bias_y_B, bias_z_B, 
        mjd, mjd0, mjdm, gps, gpss, mjds;
    int i;

    scale_x_A = 0.9595;
    scale_y_A = 0.9797;
    scale_z_A = 0.9485;
    scale_x_B = 0.9465;
    scale_y_B = 0.9842;
    scale_z_B = 0.9303;

    gpss = AC_EPH[0] + (JD0 - T0) * 86400.0 - 19 - 32.184;
    mjds = gpss / 86400.0 + T0 - 2400000.5;

    if (mjds > 55562) // Jan. 1, 2011
    {
        scale_x_A = scale_x_A * 0.98;
        scale_z_A = scale_z_A * 0.98;
        scale_x_B = scale_x_B * 0.98;
        scale_z_B = scale_z_B * 0.98;
    }



    mjdm = 52705; //March 7, 2003

if (accid == 1)
{

    for (i = 0; i < DIM_AC; i++)
    {
        gps = AC_EPH[i * 4] + (JD0 - T0) * 86400.0 - 19 - 32.184;
        mjd = gps / 86400.0 + T0 - 2400000.5;
        if ( mjd < mjdm)
        {
            mjd0 = 52532;
            bias_x_A = - 1.106 
                       + 2.233e-4 * ( mjd - mjd0 ) 
                       + 2.5e-7   * ( mjd - mjd0 ) * ( mjd - mjd0 );
            bias_y_A =   27.042 
                       + 4.46e-3  * ( mjd - mjd0 ) 
                       + 1.1e-6   * ( mjd - mjd0 ) * ( mjd - mjd0 );                       
            bias_z_A = - 0.5486 
                       - 1.139e-6 * ( mjd - mjd0 ) 
                       + 1.7e-7   * ( mjd - mjd0 ) * ( mjd - mjd0 );
        }
        else
        {
            mjd0 = 53736;
            bias_x_A = - 1.2095 
                       - 4.128e-5 * ( mjd - mjd0 ) 
                       + 9.7e-9   * ( mjd - mjd0 ) * ( mjd - mjd0 );
            bias_y_A =   29.3370 
                       + 6.515e-4 * ( mjd - mjd0 ) 
                       - 3.9e-7   * ( mjd - mjd0 ) * ( mjd - mjd0 );                       
            bias_z_A = - 0.5606 
                       - 2.352e-6 * ( mjd - mjd0 ) 
                       + 3.8e-9   * ( mjd - mjd0 ) * ( mjd - mjd0 );
        }

        AC_EPH[i * 4 + 1] = bias_x_A * 1e-6 + scale_x_A * AC_EPH[i * 4 + 1];
        AC_EPH[i * 4 + 2] = bias_y_A * 1e-6 + scale_y_A * AC_EPH[i * 4 + 2];
        AC_EPH[i * 4 + 3] = bias_z_A * 1e-6 + scale_z_A * AC_EPH[i * 4 + 3];
/*        
        if (i == 0)
        {
            x0 = AC_EPH[1];
            y0 = AC_EPH[2];
            z0 = AC_EPH[3];
        }
        AC_EPH[i * 4 + 1] = AC_EPH[i * 4 + 1] - x0;
        AC_EPH[i * 4 + 2] = AC_EPH[i * 4 + 2] - y0;
        AC_EPH[i * 4 + 3] = AC_EPH[i * 4 + 3] - z0;
*/

    }

}

if (accid == 2)
{    

    for (i = 0; i < DIM_AC; i++)
    {
        gps = AC_EPH[i * 4] + (JD0 - T0) * 86400.0 - 19 - 32.184;
        mjd = gps / 86400.0 + T0 - 2400000.5;
        if ( mjd < mjdm)
        {
            mjd0 = 52532;
            bias_x_B = - 0.5647 
                       - 7.788e-5 * ( mjd - mjd0 ) 
                       + 2.4E-7   * ( mjd - mjd0 ) * ( mjd - mjd0 );
            bias_y_B =   7.5101 
                       + 7.495E-3 * ( mjd - mjd0 ) 
                       - 9.6E-6   * ( mjd - mjd0 ) * ( mjd - mjd0 );                       
            bias_z_B = - 0.8602 
                       + 1.399E-4 * ( mjd - mjd0 ) 
                       + 2.5E-7   * ( mjd - mjd0 ) * ( mjd - mjd0 );
        }
        else
        {
            mjd0 = 53736;
            bias_x_B = - 0.6049 
                       - 1.982E-5 * ( mjd - mjd0 ) 
                       + 3.5E-9   * ( mjd - mjd0 ) * ( mjd - mjd0 );
            bias_y_B =   10.6860 
                       + 1.159E-3 * ( mjd - mjd0 ) 
                       - 4.3E-7   * ( mjd - mjd0 ) * ( mjd - mjd0 );                       
            bias_z_B = - 0.7901 
                       + 4.783E-5 * ( mjd - mjd0 ) 
                       - 6.5E-9   * ( mjd - mjd0 ) * ( mjd - mjd0 );
        }


        AC_EPH[i * 4 + 1] = bias_x_B * 1e-6 + scale_x_B * AC_EPH[i * 4 + 1];
        AC_EPH[i * 4 + 2] = bias_y_B * 1e-6 + scale_y_B * AC_EPH[i * 4 + 2];
        AC_EPH[i * 4 + 3] = bias_z_B * 1e-6 + scale_z_B * AC_EPH[i * 4 + 3];
/*
        if (i == 0)
        {
            x0 = AC_EPH[1];
            y0 = AC_EPH[2];
            z0 = AC_EPH[3];
        }
        AC_EPH[i * 4 + 1] = AC_EPH[i * 4 + 1] - x0;
        AC_EPH[i * 4 + 2] = AC_EPH[i * 4 + 2] - y0;
        AC_EPH[i * 4 + 3] = AC_EPH[i * 4 + 3] - z0;
*/
    }
}
   
    return 0;

}
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/

























/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/*
* fun_fullstate -transition matrix(36), orbit(6), sensitivity matrix(6*DYNPAR)
* @param1: description of param1
* @param2: description of param2
*
* version: 20 Aug 2010
*/
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
void fun_accel (int dim, double jd, double tt, double *state, double *fstate)
{
    int n, i,k, part;
    double tjd[2], xic[6], dfdx[36], dxdx0[36],
        dadr1[9],
        dadr2[9],
        dadr3[9], dadsrpb[3], dadsrpt[3],
        dadr4[9], dadk2[3], ac[3], 
        acc[3], dadr[9], 
        fxic[6], fdxdx0[36], ap[3], an[3], ar[3], ag[3],
        *dadp, *dxdp, *dfdpp, *dfdp, *fdxdp;

    tjd[0] = jd;    tjd[1] = tt / 86400.0;
//    tjd[0] = jd;    tjd[1] = tt;

    if (dim < 6)
    {
        printf ("error: fun_accel input dim < 6!\n");
        exit (0);
    }
    else if (dim == 6)
        part = 0;
    else if (dim > 6)
    {
        part = 1;    
    }
    for (n = 0; n < 6; n++)
    {
        xic[n] = state[n];
    }   

/* acc, partial to xyz: dadr, partial to parameters dadp*/
//    accel_ntrel (tjd, xic, part, acc1, dadr1, dadp1);
//    accel_nonsp (tjd, xic, part, acc2, dadr2, dadp2);
//    accel_radpr (tjd, xic, part, acc3, dadr3, dadp3);
/*todo: air drag acc & partial to vxvyvz dadv*/

    accel_pm_part (tjd, xic, ap, part, dadr1);
    accel_nb_part (tjd, xic, an, part, dadr2);
    accel_sr_part (tjd, xic, ar, part, dadr3, dadsrpb, dadsrpt);
    accel_gt_part (tjd, xic, ag, part, dadr4, dadk2);

//    printf ("%e\t%e\t%e\n", ac[0], ac[1], ac[2]);

    accel_ac_part (tjd, xic, ac);

//    printf ("%e\t%e\t%e\n", ac[0], ac[1], ac[2]);
//    printf ("%e\t%e\t%e\n", an[0], an[1], an[2]);


    for (n = 0; n <= 2; n++)
    {
//        ac[n] = 0;
        acc[n] = ap[n] + an[n] + ar[n] + ag[n] + ac[n];
    }

//    acc[0] = acc[0] + AX0 + AX1 * tt;
//    acc[1] = acc[1] + AY0 + AY1 * tt;
//    acc[2] = acc[2] + AZ0 + AZ1 * tt;

    fxic[0] = xic[3];
    fxic[1] = xic[4];
    fxic[2] = xic[5];
    fxic[3] = acc[0];
    fxic[4] = acc[1];
    fxic[5] = acc[2];

    for (n = 0; n < 6; n++)
    {
        fstate[n] = fxic[n];
    }

    if (part == 0)
    {
        return;
    }

    
    for (n = 0; n < 36; n++)
    {
        dxdx0[n] = state[n + 6];
    }

    for (n = 0; n <= 8; n++)
    {
        dadr[n] = dadr1[n] + dadr2[n] + dadr3[n] + dadr4[n];
//        dadr[n] = dadr1[n];
    }


    for (n = 0; n < 36; n++)
    {
        dfdx[n] = 0;
    }
    dfdx[3]  = 1; 
    dfdx[10] = 1; 
    dfdx[17] = 1;
    for (n = 0; n < 3; n++)
    {
        dfdx[n + 18] = dadr[n];
        dfdx[n + 24] = dadr[n + 3];
        dfdx[n + 30] = dadr[n + 6];
    }
    brmul(dfdx, dxdx0, 6, 6, 6, fdxdx0);


    for (n = 0; n < 36; n++)
    {
        fstate[n + 6] = fdxdx0[n];
    }


    if (MDYN == 0)
        return;
///////////////////////stop here for GRACE w/o acc parameters//////////////////

    dadp = (double *) calloc ( 3 * MDYN, sizeof(double));
    dxdp = (double *) calloc ( 6 * MDYN, sizeof(double));
    dfdpp = (double *) calloc ( 6 * MDYN, sizeof(double));
    dfdp = (double *) calloc ( 6 * MDYN, sizeof(double));
    fdxdp = (double *) calloc ( 6 * MDYN, sizeof(double));


    for (n = 0; n < 6 * MDYN; n++)
    {
        dxdp[n] = state[n + 42];
    }




    i = 0;
    if (MSRP > 0)
    {     
        for (n = 0; n < 3; n++)
        {
            dadp[n * MDYN + i] = dadsrpb[n];
        }
        i++;
    }
    if (MSRP > 1)
    {     
        for (n = 0; n < 3; n++)
        {
            dadp[n * MDYN + i] = dadsrpt[n];
        }
        i++;
    }
    if (MTK2 > 0)
    {     
        for (n = 0; n < 3; n++)
        {
            dadp[n * MDYN + i] = dadk2[n];
        }
        i++;
    }
    if (MGCS > 0)
    {     
        for (k = 0; k < MGCS; k ++)
        {
            for (n = 0; n < 3; n++)
            {
                dadp[n * MDYN + i] = CSinfo[k].dadcs[n];
            }
            i++;
        }
    }


    if (MACC_BIAS != 0)
    {
        for (k = 0; k < MACC_PRBS; k ++)
        {
            for (n = 0; n < 3; n++)
            {
                dadp[n * MDYN + i] = DADBS[n * MACC_PRBS + k];
            }
            i++;
        }
    }

    if (MACC_SCAL != 0)
    {
        for (k = 0; k < MACC_PRSL; k ++)
        {
            for (n = 0; n < 3; n++)
            {
                dadp[n * MDYN + i] = DADSL[n * MACC_PRSL + k];
            }
            i++;
        }
    }






    brmul(dfdx, dxdp, 6, 6, MDYN, dfdpp);
    for (n = 0; n < 3 * MDYN; n++)
    {
        dfdp[n] = 0;
        dfdp[n + 3 * MDYN] = dadp[n];
    }
    for (n = 0; n < 6 * MDYN; n++)
    {
        fdxdp[n] = dfdpp[n] + dfdp[n];
    }
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/


    for (n = 0; n < 6 * MDYN; n++)
    {
        fstate[n + 42]= fdxdp[n];
    }

    free (dadp);
    free (dxdp);
    free (dfdpp);
    free (dfdp);
    free (fdxdp);

    return;
}










double accel_pm_part (double *tjd, double *x, double *acc, int part, double *dadr)
{
    double GME, GMS, J, Jv[3], beta, gamma, r, v2, pv, pJ, a, b, p[3], v[3],
        pxv[3], vxJ[3], ps[3], vs[3], xsc[6], rs, vsxps[3], vsxpsxv[3],
        unit[9], ppt[9], r5, r3, fgr[3], term1[3], term2[3], term3[3];
    int n;
    short int sun = 10;

    GME = GMCT; //m^3/s^2

//    GME = GME * 86400.0 * 86400.0 / AU / AU / AU;
    p[0] = x[0]; p[1] = x[1]; p[2] = x[2];
    v[0] = x[3]; v[1] = x[4]; v[2] = x[5];

    r = modvect(p);

    acc[0] = - GME / (r*r*r) * p[0];
    acc[1] = - GME / (r*r*r) * p[1];
    acc[2] = - GME / (r*r*r) * p[2];

    
    if (part == 1)
    {
        unit[0] = 1; unit[1] = 0; unit[2] = 0;
        unit[3] = 0; unit[4] = 1; unit[5] = 0;
        unit[6] = 0; unit[7] = 0; unit[8] = 1;
        r5 = pow (r, 5);
        r3 = pow (r, 3);
        brmul (p, p, 3,1,3, ppt);
        for (n = 0; n <= 8; n++)
        {
            dadr[n] = 3 * GME * ppt[n] / r5
            - GME * unit[n] / r3;
        }
    }



    if (RELTIV == 0)
        return 0;

    J = 9.8e8; //m^2/s
    gamma = 1;
    beta = 1;
    GMS = 1.32712442076e20; //m^3/s^2
//    GMS = GMS * 86400.0 * 86400.0 / AU / AU / AU;
//    J = J * 86400.0 / AU / AU; 
    
    Jv[0] = 0; Jv[1] = 0; Jv[2] = J;

    v2 = dotvect(v, v);
    pv = dotvect(p, v);
    pJ = dotvect(p, Jv);
    crsvect(p, v, pxv);
    crsvect(v, Jv, vxJ);    

//    planet_ephemeris (tjd, CT, sun, ps, vs);
    get_ephemeris (tjd, CT, sun, xsc);
    for (n = 0; n < 3; n++)
    {
        ps[n] = xsc[n];
        vs[n] = xsc[n + 3];
    }

    rs = modvect(ps);
    crsvect(vs, ps, vsxps);
    crsvect(vsxps, v, vsxpsxv);    

    a = 2 * (beta + gamma) * GME / r - gamma * v2;
    b = 2 * (1 + gamma) * pv;
    for (n = 0; n < 3; n++)
    {
        term1[n] = GME / C / C / r / r / r *
                ( a * p[n] + b * v[n] );
        term2[n] = GME / C / C / r / r / r * (1 + gamma) * 
                ( 3/r/r * pxv[n] * pJ + vxJ[n] );
        term3[n] = - GMS / C / C / rs / rs / rs * (1 + 2 * gamma) *
                vsxpsxv[n];

        fgr[n] = term1[n]
                + term2[n] + term3[n];
//        printf ("%15.12f\t%15.12f\t%15.12f\n", term1[n],term2[n],term2[n]);
    }


    acc[0] = acc[0] + fgr[0];
    acc[1] = acc[1] + fgr[1];
    acc[2] = acc[2] + fgr[2];

    return 0;




}






int get_ephemeris (double tjd[2], int to, int from, double *x)
{
    double jd0 = 2451545.00000000, tdbj2000, fromTtoS[6] = {0}, pos[3], vel[3];
    int n;
    short int center, target;


    if (from <= 12 && to <= 12)
    {
        center = (short int)from;
        target = (short int)to;
        planet_ephemeris (tjd, target, center, pos, vel);
        x[0] = pos[0] * AU;
        x[1] = pos[1] * AU;
        x[2] = pos[2] * AU;
        x[3] = vel[0] * AU / 86400.0;
        x[4] = vel[1] * AU / 86400.0;
        x[5] = vel[2] * AU / 86400.0;
    }
    else if (from == I_TITAN || to == I_TITAN) //titan
    {

        tdbj2000 = ((tjd[0] - jd0) + tjd[1]) * 86400.0;

//        spkezr_c ("SATURN", tdbj2000, "J2000", "NONE", "TITAN", fromTtoS, &lt);
/*
Procedure

   void spkezr_c ( ConstSpiceChar     *targ,
                   SpiceDouble         et,
                   ConstSpiceChar     *ref,
                   ConstSpiceChar     *abcorr,
                   ConstSpiceChar     *obs,
                   SpiceDouble         starg[6],
                   SpiceDouble        *lt        )
 
 
   Return the state (position and velocity) of a target body 
   relative to an observing body, optionally corrected for light 
   time (planetary aberration) and stellar aberration. 
*/

        if (from == I_TITAN) center = (short int) to;
        else center = (short int) from;
        planet_ephemeris (tjd, center, 5, pos, vel);

        for (n = 0; n < 3; n++)
        {    
            x[n] = pos[n] * AU + fromTtoS[n] * 1000.0;
            x[n + 3] = vel[n] * AU / 86400.0 + fromTtoS[n + 3] * 1000.0;
        }
        if (to == I_TITAN)
        {
            for (n = 0; n < 6; n++)
            {
                x[n] = - x[n];
            }
        }    
    }
    else 
    {
        printf ("error in get_ephemeris: from = %d\t to = %d\n", from, to);
    }
    return 0;
}









double accel_nb_part (double *tjd, double *xic, double *acc, int part, double *dadr)
{
    int n;
    short int ssbary = 11;
    double xcb[6], xib[6], acb[3], aib[3], dadrc[9], dadri[9];

    if (NBODY == 0)
    {
        if (part == 1)
        {
            for (n = 0; n <= 8; n++)
            {
                dadr[n] = 0;
            }
        }

        for (n = 0; n <= 2; n++)
        {
            acc[n] = 0;
        }

        return 0;
    }


//    planet_ephemeris (tjd, CT, ssbary, &xcb[0], &xcb[3]);
    get_ephemeris (tjd, CT, ssbary, xcb);
    f_bcrs (tjd, xcb, CT, acb, part, dadrc);
    for (n = 0; n <= 5; n++)
    {
        xib[n] = xic[n] + xcb[n];
    }
    f_bcrs (tjd, xib, CT, aib, part, dadri);

    for (n = 0; n <= 2; n++)
    {
        acc[n] = aib[n] - acb[n];
    }

    if (part == 1)
    {
        for (n = 0; n <= 8; n++)
        {
            dadr[n] = dadri[n];
        }
    }

    return 0;
}











/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
double f_bcrs (double *tjd, double *xi, int exclude, 
                   double *acc, int part, double *dadr)
{
    double fnt[3], fgr[3], xj[11][6], xij[11][6], rij[11], xjk[6], rjk, 
        xddj[3], sumil, sumjk, sdi2, sdj2, rdirdj, rrrdr2, rjirdd, 
        rij5, rij3, xijt[9], gra, grb, beta, gamma, unit[9];
    short int ssbary, l, k, j, n, flag_gr;
    
    ssbary = 11;
    gamma = 1.0;
    beta = 1.0;
    unit[0] = 1; unit[1] = 0; unit[2] = 0;
    unit[3] = 0; unit[4] = 1; unit[5] = 0;
    unit[6] = 0; unit[7] = 0; unit[8] = 1;

    for (j = 0; j <= 10; j++)
    {
//        planet_ephemeris (jd, j, ssbary, &xj[j][0], &xj[j][3]);
        get_ephemeris (tjd, j, ssbary, xj[j]);
        for (n = 0; n < 6; n++)
        {
            xij[j][n] = xi[n] - xj[j][n];
        }
        rij[j] = sqrt (xij[j][0] * xij[j][0] 
            + xij[j][1] * xij[j][1] + xij[j][2] * xij[j][2]);
    }
    
    flag_gr = 0;
    for (n = 0; n < 3; n ++)
        fnt[n] = 0;
    for (j = 0; j <= 10; j++)
    {
        if (PERB[j] == 2)
            flag_gr = 1;
        if (PERB[j] == 0)
            continue;
        if (j == exclude)
            continue;
        for (n = 0; n < 3; n++)
            fnt[n] = fnt[n] 
            - GMDE[j] / (rij[j] * rij[j] * rij[j]) * xij[j][n];
    }

    if (part == 1)
    {
        for (n = 0; n <= 8; n++)
        {
            dadr[n] = 0;
        }
        for (j = 0; j <= 10; j++)
        {
            if (j == exclude)
                continue;
            rij5 = pow (rij[j], 5);
            rij3 = pow (rij[j], 3);
            brmul (xij[j], xij[j], 3,1,3, xijt);
            for (n = 0; n <= 8; n++)
            {
                dadr[n] = dadr[n] + 3 * GMDE[j] * xijt[n] / rij5
                - GMDE[j] * unit[n] / rij3;
            }
        }
    }

    if (flag_gr == 0)
    {
        for (n = 0; n < 3; n++)
            acc[n] =  fnt[n];
        return 0;
    }

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/

    sdi2 = xi[3] * xi[3] + xi[4] * xi[4] + xi[5] * xi[5];
    sumil = 0;
    for (l = 0; l < 11; l ++)
    {
        if ( l == exclude)
            continue;
        if (PERB[l] != 2)
            continue;
        sumil = sumil + GMDE[l] / rij[l];
    }

    for (n = 0; n < 3; n ++)
        fgr[n] = 0;
    for (j = 0; j < 11; j ++)
    {
        if (PERB[j] != 2)
            continue;
        if (j == exclude)
            continue;
        sumjk = 0;
        for (n = 0; n < 3; n ++)
            xddj[n] = 0;
        for (k = 0; k < 11; k ++)
        {
            if (k == j)	
                continue;	//k!=j
            if (PERB[k] != 2)
                continue;
            for (n = 0; n < 3; n++)
                xjk[n] = xj[j][n] - xj[k][n];
            rjk = sqrt (xjk[0] * xjk[0] + xjk[1] * xjk[1] + xjk[2] * xjk[2]);
            sumjk = sumjk + GMDE[k] / rjk;
            for (n = 0; n < 3; n ++)
                xddj[n] = xddj[n] - GMDE[k] / (rjk * rjk * rjk) * xjk[n];
        }
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
        sdj2 = xj[j][3] * xj[j][3] + xj[j][4] * xj[j][4] 
            + xj[j][5] * xj[j][5];
        rdirdj = xi[3] * xj[j][3] + xi[4] * xj[j][4] + xi[5] * xj[j][5];
        rrrdr2 = pow( ( xij[j][0] * xj[j][3] + xij[j][1] * xj[j][4] 
            + xij[j][2] * xj[j][5]) / rij[j], 2);
        rjirdd = - ( xij[j][0] * xddj[0] + xij[j][1] * xddj[1] 
            + xij[j][2] * xddj[2]);
        
        gra = - 2 * (beta + gamma) * sumil - (2 * beta -1) * sumjk 
            + gamma * sdi2 + (1 + gamma) * sdj2
            - 2 * (1 + gamma) * rdirdj - 1.5 * rrrdr2 + 0.5 * rjirdd;

        grb = xij[j][0] * ((2+2*gamma) * xi[3] - (1+2*gamma) * xj[j][3])
            + xij[j][1] * ((2+2*gamma) * xi[4] - (1+2*gamma) * xj[j][4])
            + xij[j][2] * ((2+2*gamma) * xi[5] - (1+2*gamma) * xj[j][5]);

        for (n = 0; n < 3; n ++)
        {
            fgr[n] = fgr[n] 
                + GMDE[j] / (rij[j] * rij[j] * rij[j]) 
                * ( - xij[j][n]) * gra / C / C
                + GMDE[j] / (rij[j] * rij[j] * rij[j]) 
                * xij[j][n + 3] * grb / C / C
                + GMDE[j] / rij[j] * (3 + 4 * gamma) * 0.5   
                * xddj[n] / C / C;
        }
    }

    for (n = 0; n < 3; n++)
        acc[n] =  fgr[n] + fnt[n];
    return 1;
}



















double accel_sr_part (double *tjd, double *xic, double *acc, int part,
        double *dadr, double *dadsrpb, double *dadsrpt)
{
    double j, c1, rsp, usp[3], xis[6], xsc[6], f, 
        xist[9], unit[9], rsp3;
    short int n, sun;


    if (AMR == 0)
    {
        if (part == 1)
        {
            for (n = 0; n <= 8; n++)
            {
                dadr[n] = 0;
            }
            for (n = 0; n <= 3; n++)
            {
                dadsrpb[n] = 0;
                dadsrpt[n] = 0;
            }
        }

        for (n = 0; n <= 2; n++)
        {
            acc[n] = 0;
        }

        return 0;
    }


    sun = 10;
    unit[0] = 1; unit[1] = 0; unit[2] = 0;
    unit[3] = 0; unit[4] = 1; unit[5] = 0;
    unit[6] = 0; unit[7] = 0; unit[8] = 1;

//    planet_ephemeris (tjd, sun, CT, &xsc[0], &xsc[3]);    
    get_ephemeris (tjd, sun, CT, xsc);
    for (n = 0; n <= 5; n++)
    {
        xis[n] = xic[n] - xsc[n];
    }
    rsp = sqrt (xis[0] * xis[0] + xis[1] * xis[1] + xis[2] * xis[2]);
    usp[0] = xis[0] / rsp;
    usp[1] = xis[1] / rsp;
    usp[2] = xis[2] / rsp;

    j  = 1352.5;   //kg/s3
//    j  = 1359.4;   //kg/s3
//    m  = SATMASS;     //kg
//    ap = SATAREA;        //m2
    c1 = j / C * AU * AU;   //kg/s2/m*au*au
    f  = c1 * AMR  / rsp / rsp;    
//    f  = c1 * ap / m  / rsp / rsp;    
//kg/s2/m*au*au * m2 / kg / au  / au = m/s2
//    f = f / AU * 86400.0 * 86400.0;


//    printf ("SRPB = %f\t SRPT = %f\n", SRPB, SRPT);
    for (n = 0; n < 3; n++)
    {
       acc[n] = f * usp[n] * (1 + SRPB + SRPT * tjd[1]);
    }

    if (part == 0)
        return 1;

    rsp3 = rsp * rsp * rsp;
    brmul (xis, xis, 3,1,3, xist);
    for (n = 0; n <= 8; n++)
    {
        dadr[n] = - f * (3 * xist[n] / rsp3 - unit[n] / rsp) * 
            (1 + SRPB + SRPT * tjd[1]);
    }

    for (n = 0; n <= 2; n++)
    {
        dadsrpb[n] = f * usp[n];      
        dadsrpt[n] = f * usp[n] * tjd[1];   
    }

    return 0;
}




















double accel_gt_part (double *tjd, double *xic, double *acc, int part,
        double *dadr, double *dadk2)
{
    int n,k, blst[12], nb = 0;
    double GM, radius, pi[3], pe[3], llr[3], ae[3], dadk2e[3], 
        dadre[9], dadres[9], dadrei[9];
    InfStruct info;

    GM = GMCT;
    radius = RCT;

    for (n = 0; n <= 2; n++)
    {
        acc[n] = 0;
    }

    if (NMAX < 2)
    {
        if (part == 1)
        {
            for (n = 0; n <= 8; n++)
            {
                dadr[n] = 0;
            }
            for (n = 0; n <= 3; n++)
            {
                dadk2[n] = 0;
            }
        }

        return 0;
    }

    for (n = 0; n < 3; n++)
    {
        pi[n] = xic[n];
    }

    getinfo(tjd, 2, &info);

    brmul(info.c_ie, pi, 3, 3, 1, pe);  
    xyz2llh(pe, llr);

    for (n = 0; n < (NMAX + 1) * (NMAX + 1); n++)
    {   
        COEF[n] = COEFG[n];
    }


    if (STIDE != 0)
    {


        if (STIDE >= 3 && CT == 2)
            stidecs_earth(&info, COEFS, C20PERM, 1,1,0);
//            stidecs_Anelastic(&info, 1, COEFS);
        else if (STIDE == 2 && CT == 2)
            stidecs_earth(&info, COEFS, C20PERM, 0,0,0);
//            stidecs(tjd, info.c_ie, 1, COEFS);
        else if (CT == 2)
        {
            blst[0] = 10; blst[1] = 9; nb = 2;
            stidecs_k2 (&info, K2, COEFS, blst, nb);
        }
        else if (CT == 9)
        {
            blst[0] = 10; blst[1] = 2; nb = 2;
            stidecs_k2 (&info, K2, COEFS, blst, nb);
        }
        else if (CT == 20)
        {
            blst[0] = 10; blst[1] = 5; nb = 2;
            stidecs_k2 (&info, K2, COEFS, blst, nb);
        }
        else 
        {
            blst[0] = 10; nb = 1;
            stidecs_k2 (&info, K2, COEFS, blst, nb);
        }
        
//        sumcs (COEF, COEFS, NMAX, NSMAX);
        addsubcs (COEF, COEFS, NMAX, NSMAX, 1);
    
    }


    if (OTIDE != 0) 
    {
//        otidecs_csr(&info, NOMAX, COEFO);
        lgr_order (OT_EPH, DIM_OT, (NOMAX + 1) * (NOMAX + 1) + 1 , info.tt, COEFO, 4);
//        sumcs (COEF, COEFO, NMAX, NOMAX);
        addsubcs (COEF, COEFO, NMAX, NOMAX, 1);
    }

/*
    if (ATIDE != 0) 
    {
        lgr_order (AT_EPH, DIM_AT, (NAMAX + 1) * (NAMAX + 1) + 1 , info.tt, COEFA, 4);
        sumcs (COEF, COEFA, NMAX, NAMAX);
    }


    if (AOD1B != 0) 
    {
        t = info.gps / 86400.0;
        lgr_order (AOD_EPH, DIM_AO, (NDMAX + 1) * (NDMAX + 1) + 1 , t, COEFD, 1);
        sumcs (COEF, COEFD, NMAX, NDMAX);
    }


    if (PTIDE != 0)
    {
        poletide (&info, NPMAX, COEFP);
        sumcs (COEF, COEFP, NMAX, NPMAX);
    }
*/
    cs2ada (llr, COEF,GM, radius, NMAX, ae, part, dadre, 0);
    brmul(info.c_ei, ae, 3, 3, 1, acc);


    if (part == 0)
        return 1;


    brmul(dadre, info.c_ie, 3, 3, 3, dadrei);			
    brmul(info.c_ei, dadrei, 3, 3, 3, dadr);	

    if (MTK2 == 1)
    {
        stidecs_k2 (&info, 1, COEFS, blst, nb);
        cs2ada (llr, COEFS,GM, radius, NSMAX, dadk2e, 0, dadres, 0);
        brmul(info.c_ei, dadk2e, 3, 3, 1, dadk2);
    }

    if (MGCS > 0)
    {
        for (k = 0; k < MGCS; k ++)
        {
            brmul(info.c_ei, CSinfo[k].dadcse, 3, 3, 1, CSinfo[k].dadcs);
        }    
    }
    return 0;

}




//    nmax = 4;      
//    stcs = (double *) calloc ( (nmax + 1) * (nmax + 1), sizeof(double));
double stidecs_k2(InfStruct *info, double k2, double *stcs, int *body, int nbody)
{
    int sun, i;

    double xs[6], gms2e, tjd[2], 
        pse[3], llrs[3], pbar[4], t,
        p20s, p30s, p21s, p31s, p22s, p32s, p33s,
        rers, c20, c21, s21, c22, s22;

    tjd[0] = info->jd0;
    tjd[1] = info->tt/86400.0;


    c20 = 0; c21 = 0; s21 = 0; c22 = 0; s22 = 0;
    for (i = 0; i < nbody; i++)
    {
        sun = body[i]; 
        if (sun > 12)
            continue;
    
        get_ephemeris (tjd, sun, CT, xs);
//        planet_ephemeris (tjd, sun, earth, ps, vs);

        brmul (info->c_ie, xs, 3, 3, 1, pse);  
    
        xyz2llh(pse, llrs);

        t = sin(llrs[0] * DEG2RAD);
        lgdr(t, 3, 0, pbar); p20s = pbar[2]; p30s = pbar[3];
        lgdr(t, 3, 1, pbar); p21s = pbar[1]; p31s = pbar[2];
        lgdr(t, 3, 2, pbar); p22s = pbar[0]; p32s = pbar[1];
        lgdr(t, 3, 3, pbar); p33s = pbar[0];


        gms2e  =  GMDE[sun]/GMCT;
        rers = RCT / llrs[2];

        c20 += k2/5.0 * ( gms2e * pow(rers, 3) * p20s );
        c21 += k2/5.0 * ( gms2e * pow(rers, 3) * p21s * cos(llrs[1] * DEG2RAD) );
        s21 += k2/5.0 * ( gms2e * pow(rers, 3) * p21s * sin(llrs[1] * DEG2RAD) );
        c22 += k2/5.0 * ( gms2e * pow(rers, 3) * p22s * cos(llrs[1] * DEG2RAD * 2.0) );
        s22 += k2/5.0 * ( gms2e * pow(rers, 3) * p22s * sin(llrs[1] * DEG2RAD * 2.0) );
    }
    
    stcs[0]  = 0; //c00; 
    stcs[1]  = 0; //c10; 
    stcs[2]  = c20; 
    stcs[3]  = 0; 
    stcs[4]  = c21; 

    stcs[5]  = 0; //s11; 
    stcs[6]  = s21; 
    stcs[7]  = c22; 
    stcs[8]  = s22; 
    
    return 0;

}






















double accel_slrad (double *tjd, double *xic, double *acc)
{
    double j, c1, rsp, usp[3], xis[6], xsc[6], f, 
        unit[9];
    short int n, sun;

    sun = 10;
    unit[0] = 1; unit[1] = 0; unit[2] = 0;
    unit[3] = 0; unit[4] = 1; unit[5] = 0;
    unit[6] = 0; unit[7] = 0; unit[8] = 1;

    planet_ephemeris (tjd, sun, CT, &xsc[0], &xsc[3]);
    for (n = 0; n <= 5; n++)
    {
        xis[n] = xic[n] - xsc[n];
    }
    rsp = sqrt (xis[0] * xis[0] + xis[1] * xis[1] + xis[2] * xis[2]);
    usp[0] = xis[0] / rsp;
    usp[1] = xis[1] / rsp;
    usp[2] = xis[2] / rsp;

    j  = 1352.5;   //kg/s3
//    j  = 1359.4;   //kg/s3
//    m  = SATMASS;     //kg
//    ap = SATAREA;        //m2
    c1 = j / C * 1 * 1;   //kg/s2/m*au*au
    f  = c1 * AMR  / rsp / rsp;    
//    f  = c1 * ap / m  / rsp / rsp;    
//kg/s2/m*au*au * m2 / kg / au  / au = m/s2
    f = f / AU * 86400.0 * 86400.0;

    acc[0] = f * usp[0];
    acc[1] = f * usp[1];
    acc[2] = f * usp[2]; 


    return 0;
}





double accel_nbody (double *tjd, double *xic, double *fnt, double *fgr)
{
    int n;
    short int ssbary = 11;

    double xcb[6], xib[6], fnti[3], fntb[3], fgri[3], fgrb[3];

    planet_ephemeris (tjd, CT, ssbary, &xcb[0], &xcb[3]);
    force_bcrs (tjd, xcb, CT, fntb, fgrb);
    for (n = 0; n <= 5; n++)
    {
        xib[n] = xic[n] + xcb[n];
    }
    force_bcrs (tjd, xib, CT, fnti, fgri);
//    force_bcrs (tjd, xib, 99, fnti, fgri);
    for (n = 0; n <= 2; n++)
    {
        fnt[n] = fnti[n] - fntb[n];
        fgr[n] = fgri[n] - fgrb[n];

    }

    return 0;
}


/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
double force_bcrs (double *jd, double *xi, short int exclude, 
                   double *fnt, double *fgr)
{
    double xj[11][6], xij[11][6], rij[11], xjk[6], rjk, 
        xddj[3], sumil, sumjk, sdi2, sdj2, rdirdj, rrrdr2, rjirdd, gm[11], GMDE[11],
        gra, grb, beta, gamma, unit[9],gm2de;
    short int ssbary, l, k, j, n;
    
    gm[0] =   2.203208082807623e+13;
    gm[1] =      3.248586038641429e+14;
    gm[2] =     398600.44150E+09;
    gm[3] =     4.28283719012840e+13;
    gm[4] =      1.267127698227696e+17;
    gm[5] =     3.794062664949063e+16;
    gm[6] =      5.794549096929744e+15;
    gm[7] =     6.836534169987595e+15;
    gm[8] =    9.816009029289940e+11;
    gm[9] =      4.902801056E+12;
    gm[10] =      1.32712442076e20;

    gm2de = 86400.0 * 86400.0 / AU / AU / AU;
    for (n = 0; n <= 10; n++)
        GMDE[n] = gm[n] * gm2de;

    ssbary = 11;
//    ssbary = 10;
    gamma = 1.0;
    beta = 1.0;
    unit[0] = 1; unit[1] = 0; unit[2] = 0;
    unit[3] = 0; unit[4] = 1; unit[5] = 0;
    unit[6] = 0; unit[7] = 0; unit[8] = 1;

    for (j = 0; j <= 10; j++)
    {
        planet_ephemeris (jd, j, ssbary, &xj[j][0], &xj[j][3]);
        for (n = 0; n < 6; n++)
        {
            xij[j][n] = xi[n] - xj[j][n];
        }
        rij[j] = sqrt (xij[j][0] * xij[j][0] 
            + xij[j][1] * xij[j][1] + xij[j][2] * xij[j][2]);
    }
    
    for (n = 0; n < 3; n ++)
        fnt[n] = 0;
    for (j = 0; j <= 10; j++)
    {
        if (j == exclude)
            continue;
        for (n = 0; n < 3; n++)
            fnt[n] = fnt[n] 
            - GMDE[j] / (rij[j] * rij[j] * rij[j]) * xij[j][n];
    }


/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/

    sdi2 = xi[3] * xi[3] + xi[4] * xi[4] + xi[5] * xi[5];
    sumil = 0;
    for (l = 0; l < 11; l ++)
    {
        if ( l == exclude)
            continue;
        sumil = sumil + GMDE[l] / rij[l];
    }

    for (n = 0; n < 3; n ++)
        fgr[n] = 0;
    for (j = 0; j < 11; j ++)
    {
        if (j == exclude)
            continue;
        sumjk = 0;
        for (n = 0; n < 3; n ++)
            xddj[n] = 0;
        for (k = 0; k < 11; k ++)
        {
            if (k == j)	
                continue;	//k!=j
            for (n = 0; n < 3; n++)
                xjk[n] = xj[j][n] - xj[k][n];
            rjk = sqrt (xjk[0] * xjk[0] + xjk[1] * xjk[1] + xjk[2] * xjk[2]);
            sumjk = sumjk + GMDE[k] / rjk;
            for (n = 0; n < 3; n ++)
                xddj[n] = xddj[n] - GMDE[k] / (rjk * rjk * rjk) * xjk[n];
        }
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
        sdj2 = xj[j][3] * xj[j][3] + xj[j][4] * xj[j][4] 
            + xj[j][5] * xj[j][5];
        rdirdj = xi[3] * xj[j][3] + xi[4] * xj[j][4] + xi[5] * xj[j][5];
        rrrdr2 = pow( ( xij[j][0] * xj[j][3] + xij[j][1] * xj[j][4] 
            + xij[j][2] * xj[j][5]) / rij[j], 2);
        rjirdd = - ( xij[j][0] * xddj[0] + xij[j][1] * xddj[1] 
            + xij[j][2] * xddj[2]);
        
        gra = - 2 * (beta + gamma) * sumil - (2 * beta -1) * sumjk 
            + gamma * sdi2 + (1 + gamma) * sdj2
            - 2 * (1 + gamma) * rdirdj - 1.5 * rrrdr2 + 0.5 * rjirdd;

        grb = xij[j][0] * ((2+2*gamma) * xi[3] - (1+2*gamma) * xj[j][3])
            + xij[j][1] * ((2+2*gamma) * xi[4] - (1+2*gamma) * xj[j][4])
            + xij[j][2] * ((2+2*gamma) * xi[5] - (1+2*gamma) * xj[j][5]);

        for (n = 0; n < 3; n ++)
        {
            fgr[n] = fgr[n] 
                + GMDE[j] / (rij[j] * rij[j] * rij[j]) 
                * ( - xij[j][n]) * gra / C_AUDAY / C_AUDAY
                + GMDE[j] / (rij[j] * rij[j] * rij[j]) 
                * xij[j][n + 3] * grb / C_AUDAY / C_AUDAY
                + GMDE[j] / rij[j] * (3 + 4 * gamma) * 0.5   
                * xddj[n] / C_AUDAY / C_AUDAY;
        }
    }
    return 1;
}















/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/*
* lagrange interpolation order = 6, 2*order points
* @param1: description of param1
* @param2: description of param2
* todo
    order = input parameter    
* version: 20 Aug 2010
*/
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
double lagrange (double *y, int dim_y, int dim_x, double t, double *z)
{ 
    int i, j, k, m, dim, order = 8;
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
















double obs_alt (double jd, double utc, double *obs, int part, double *bmat)
{
    int n, lps, i;
    double r, h, ref, *eph, *dxdp = NULL, *dodpo = NULL, *dodpd = NULL, 
        xc2[6], dxdx0[36], dodx[6], 
        dodx0[6], tt, dx[3];

    ref = RCT; //should be topography height

    lps = getlps (JD0 + utc/86400.0);
    tt = utc + (lps + 32.184);

//    tjd[0] = JD0;    tjd[1] = tt / 86400.0;
//    get_ephemeris (tjd, 2, CT, xsc);
 
    if (part == 0)
    {
        lagrange (OR_EPH, DIM_OR, 7, tt, xc2);
    }
    if (part == 1)
    {
        eph = (double *) calloc (42 + 6 * MDYN, sizeof(double));
        lagrange (OR_EPH, DIM_OR, 42 + 6 * MDYN + 1, tt, eph);
        for (n = 0; n < 6; n++)
            xc2[n] = eph[n];
        for (n = 0; n < 36; n++)
            dxdx0[n] = eph[n + 6];
        if (MDYN > 0)
        {
            dxdp = (double *) calloc (6 * MDYN, sizeof(double));
            dodpd = (double *) calloc (MDYN, sizeof(double));
            for (n = 0; n < 6 * MDYN; n++)
                dxdp[n] = eph[n + 42];
        }
        if (MOBS > 0)
            dodpo = (double *) calloc (MOBS, sizeof(double));
        free (eph);
    }

    h = modvect(xc2) - RCT;

    for (n = 0; n < 3; n++)
        dx[n] = xc2[n];
    r = modvect(dx);

    *obs = r - ref + BASB + BAST * utc;

    if (part == 0)
        return h;

    for (n = 0; n < 3; n++)
    {
        dodx[n] = (xc2[n])/r;
        dodx[n + 3] = 0;
    }

    brmul (dodx, dxdx0, 1, 6, 6, dodx0);
    for (n = 0; n < 6; n++)
        bmat[n] = dodx0[n];
    
    if (MEST == 0) 
        return h;

    i = 0;
    if (MOBS > 0)
    {     
        dodpo[i] = 1;
        i++;
    }
    if (MOBS > 1)
    {
        dodpo[i] = utc;
        i++;
    }

    if (MDYN > 0)
    {
        brmul (dodx, dxdp, 1, 6, MDYN, dodpd);
    }
    if (MOBS > 0)
    for (n = 0; n < MOBS; n++)
        bmat[6 + n] = dodpo[n];
    if (MDYN > 0)
    for (n = 0; n < MDYN; n++)
        bmat[6 + MOBS + n] = dodpd[n];
    
    if (MDYN > 0)
    {
        free (dxdp);
        free (dodpd);
    }
    if (MOBS > 0)
        free (dodpo);
    return h;


}



double obs_mex (double tt, double *xsim, int part, double *xbmat)
{
    int n, i;
    double h, *eph, *dodpo = NULL, *dodpd = NULL, dodx0[36], xc2[6];

//    gps0 = (JD0 - T0) * 86400.0;
//    tt = gps - gps0 + 19 + 32.184;
//    lps = getlps (JD0 + tt / 86400.0);
//    utc = gps - gps0 - (lps - 19);

//    lps = getlps (JD0 + utc/86400.0);
//    tt = utc + (lps + 32.184);

    if (part == 0)
    {
        lagrange (OR_EPH, DIM_OR, 7, tt, xc2);
    }
    if (part == 1)
    {
        eph = (double *) calloc (42 + 6 * MDYN, sizeof(double));
        lagrange (OR_EPH, DIM_OR, 42 + 6 * MDYN + 1, tt, eph);
        for (n = 0; n < 6; n++)
            xsim[n] = eph[n];
        for (n = 0; n < 36; n++)
//            dxdx0[n] = eph[n + 6];
            dodx0[n] = eph[n + 6];
        if (MDYN > 0)
        {
//            dodp = (double *) calloc (6 * MDYN, sizeof(double));
            dodpd = (double *) calloc (6 * MDYN, sizeof(double));
            for (n = 0; n < 6 * MDYN; n++)
                dodpd[n] = eph[n + 42];
        }
        if (MOBS > 0)
            dodpo = (double *) calloc (6 * MOBS, sizeof(double));
        free (eph);
    }

    h = modvect(xsim) - RCT;
        
//    for (n = 0; n < 3; n++)
//        dv[n] = xc2[n + 3] - xsc[n + 3];
//    v = modvect(dv);

//    *obs = v + BASB + BAST * utc;

    if (part == 0)
        return h;

//    for (n = 0; n < 3; n++)
//    {
//        dodx[n] = 0;
//        dodx[n + 3] = (xc2[n + 3] - xsc[n + 3])/v;
//    }

//    brmul (dodx, dxdx0, 1, 6, 6, dodx0);

    for (i = 0; i < 6; i++)
    {
        for (n = 0; n < 6; n++)
            xbmat[i * MSOL + n] = dodx0[i * 6 + n];
    }
    
    if (MEST == 0) 
        return h;

//    i = 0;
//    if (MOBS > 0)
//    {     
//        dodpo[i] = 1;
//        i++;
//    }
//    if (MOBS > 1)
//    {
//        dodpo[i] = utc;
//        i++;
//    }

//    if (MDYN > 0)
//    {
//        brmul (dodx, dxdp, 1, 6, MDYN, dodpd);
//    }
    if (MOBS > 0)
    {
        for (i = 0; i < 6; i++)
        {
            for (n = 0; n < MOBS; n++)
                xbmat[i * MSOL + 6 + n] = dodpo[i * MOBS + n];
        }
    }

    if (MDYN > 0)
    {
        for (i = 0; i < 6; i++)
        {
            for (n = 0; n < MDYN; n++)
                xbmat[i * MSOL + 6 + MOBS + n] = dodpd[i * MDYN + n];
        }

    }
    if (MDYN > 0)
    {
//        free (dxdp);
        free (dodpd);
    }
    if (MOBS > 0)
        free (dodpo);
    return h;


}










double obs_vel (double jd, double utc, double *obs, int part, double *bmat)
{
    int n, lps, i;
    double v, h, *eph, *dxdp = NULL, *dodpo = NULL, *dodpd = NULL, 
        xc2[6], dxdx0[36], dodx[6], 
        dodx0[6], tt, tjd[2], xsc[6], dv[3];

    lps = getlps (JD0 + utc/86400.0);
    tt = utc + (lps + 32.184);

    tjd[0] = JD0;    tjd[1] = tt / 86400.0;
    get_ephemeris (tjd, 2, CT, xsc);
 
    if (part == 0)
    {
        lagrange (OR_EPH, DIM_OR, 7, tt, xc2);
    }
    if (part == 1)
    {
        eph = (double *) calloc (42 + 6 * MDYN, sizeof(double));
        lagrange (OR_EPH, DIM_OR, 42 + 6 * MDYN + 1, tt, eph);
        for (n = 0; n < 6; n++)
            xc2[n] = eph[n];
        for (n = 0; n < 36; n++)
            dxdx0[n] = eph[n + 6];
        if (MDYN > 0)
        {
            dxdp = (double *) calloc (6 * MDYN, sizeof(double));
            dodpd = (double *) calloc (MDYN, sizeof(double));
            for (n = 0; n < 6 * MDYN; n++)
                dxdp[n] = eph[n + 42];
        }
        if (MOBS > 0)
            dodpo = (double *) calloc (MOBS, sizeof(double));
        free (eph);
    }

    h = modvect(xc2) - RCT;
        
    for (n = 0; n < 3; n++)
        dv[n] = xc2[n + 3] - xsc[n + 3];
    v = modvect(dv);

    *obs = v + BASB + BAST * utc;

    if (part == 0)
        return h;

    for (n = 0; n < 3; n++)
    {
        dodx[n] = 0;
        dodx[n + 3] = (xc2[n + 3] - xsc[n + 3])/v;
    }

    brmul (dodx, dxdx0, 1, 6, 6, dodx0);
    for (n = 0; n < 6; n++)
        bmat[n] = dodx0[n];
    
    if (MEST == 0) 
        return h;

    i = 0;
    if (MOBS > 0)
    {     
        dodpo[i] = 1;
        i++;
    }
    if (MOBS > 1)
    {
        dodpo[i] = utc;
        i++;
    }

    if (MDYN > 0)
    {
        brmul (dodx, dxdp, 1, 6, MDYN, dodpd);
    }
    if (MOBS > 0)
    for (n = 0; n < MOBS; n++)
        bmat[6 + n] = dodpo[n];
    if (MDYN > 0)
    for (n = 0; n < MDYN; n++)
        bmat[6 + MOBS + n] = dodpd[n];
    
    if (MDYN > 0)
    {
        free (dxdp);
        free (dodpd);
    }
    if (MOBS > 0)
        free (dodpo);
    return h;


}









double obs_dsn (double jd, double utc, double *obs, int part, double *bmat)
{
    int n, lps, i;
    double r, h, *eph, *dxdp = NULL, *dodpo = NULL, *dodpd = NULL, 
        xc2[6], dxdx0[36], dodx[6], 
        dodx0[6], tt, tjd[2], xsc[6], dx[3];

    lps = getlps (JD0 + utc/86400.0);
    tt = utc + (lps + 32.184);

    tjd[0] = JD0;    tjd[1] = tt / 86400.0;
    get_ephemeris (tjd, 2, CT, xsc);
 
    if (part == 0)
    {
        lagrange (OR_EPH, DIM_OR, 7, tt, xc2);
    }
    if (part == 1)
    {
        eph = (double *) calloc (42 + 6 * MDYN, sizeof(double));
        lagrange (OR_EPH, DIM_OR, 42 + 6 * MDYN + 1, tt, eph);
        for (n = 0; n < 6; n++)
            xc2[n] = eph[n];
        for (n = 0; n < 36; n++)
            dxdx0[n] = eph[n + 6];
        if (MDYN > 0)
        {
            dxdp = (double *) calloc (6 * MDYN, sizeof(double));
            dodpd = (double *) calloc (MDYN, sizeof(double));
            for (n = 0; n < 6 * MDYN; n++)
                dxdp[n] = eph[n + 42];
        }
        if (MOBS > 0)
            dodpo = (double *) calloc (MOBS, sizeof(double));
        free (eph);
    }

    h = modvect(xc2) - RCT;
        
    for (n = 0; n < 3; n++)
        dx[n] = xc2[n] - xsc[n];
    r = modvect(dx);

    *obs = r + BASB + BAST * utc;

    if (part == 0)
        return h;

    for (n = 0; n < 3; n++)
    {
        dodx[n] = (xc2[n] - xsc[n])/r;
        dodx[n + 3] = 0;
    }

    brmul (dodx, dxdx0, 1, 6, 6, dodx0);
    for (n = 0; n < 6; n++)
        bmat[n] = dodx0[n];
    
    if (MEST == 0) 
        return h;

    i = 0;
    if (MOBS > 0)
    {     
        dodpo[i] = 1;
        i++;
    }
    if (MOBS > 1)
    {
        dodpo[i] = utc;
        i++;
    }

    if (MDYN > 0)
    {
        brmul (dodx, dxdp, 1, 6, MDYN, dodpd);
    }
    if (MOBS > 0)
    for (n = 0; n < MOBS; n++)
        bmat[6 + n] = dodpo[n];
    if (MDYN > 0)
    for (n = 0; n < MDYN; n++)
        bmat[6 + MOBS + n] = dodpd[n];
    
    if (MDYN > 0)
    {
        free (dxdp);
        free (dodpd);
    }
    if (MOBS > 0)
        free (dodpo);
    return h;


}







void getsolvefor ()
{
    int k;
       
//    MOBS = 0;   //2; 
//    MSRP = 0;   //2; 
//    MTK2 = 0;   //1; 
//    MGCS = 1;   //6;

//    MACC_DT = 3600;

    MACC_PRBS = 0;
    MACC_PRSL = 0;

    if (MACC_BIAS != 0)
    {
        MACC_ARBS = (int)(86400 / MACC_DTBS);
        MACC_NOBS = 3 * MACC_BIAS;
        MACC_PRBS = MACC_NOBS * MACC_ARBS;
        MACC_EPBS = (double *) calloc (MACC_PRBS, sizeof(double));
        DADBS = (double *) calloc (3 * MACC_PRBS, sizeof(double));
    }


    if (MACC_SCAL != 0)
    {
        MACC_ARSL = (int)(86400 / MACC_DTSL);
        MACC_NOSL = 3 * MACC_SCAL;
        MACC_PRSL = MACC_NOSL * MACC_ARSL;
        MACC_EPSL = (double *) calloc (MACC_PRSL, sizeof(double));
        DADSL = (double *) calloc (3 * MACC_PRSL, sizeof(double));
        for (k = 0; k < MACC_PRSL; k ++)
        {
            MACC_EPSL[k] = 1;
        }
    }

    MACC = MACC_PRBS + MACC_PRSL;
    MDYN = MSRP + MTK2 + MGCS + MACC;                // dim of sensitivity matrix
    MSOL = MOBS + MDYN + 6;                       // dim of regress matrix
    MEST = MOBS + MDYN;  //= MOBS + MDYN - 6

    MSTA = 42 + MDYN * 6;

/*
    if (MGCS > 0)    
    {
        CSinfo = (CSStruct *) calloc ( MGCS, sizeof(CSStruct));

        CSinfo[0].n = 2; CSinfo[0].m = 0; CSinfo[0].cs = 0;
//        CSinfo[1].n = 3; CSinfo[1].m = 0; CSinfo[1].cs = 0; 
    }
*/
    return;

}








void initsolvefor (double *xsm, double *x)
{
    int i, k, n, m, ind, l, ic, is, label;

    for (k = 0; k < 6; k ++)
    {
        x[k] = xsm[k];
    }
    i = 6;
    if (MOBS > 0)
    {
        x[i] = BASB;
        i++;
    }
    if (MOBS > 1)
    {
        x[i] = BAST;
        i++;
    }
    if (MSRP > 0)
    {
        x[i] = SRPB;
        i++;
    }
    if (MSRP > 1)
    {
        x[i] = SRPT ;
        i++;
    }
    if (MTK2 > 0)
    {
        x[i] = K2;
        i++;
    }

    if (MGCS > 0)
    for (k = 0; k < MGCS; k ++)
    {
        n = CSinfo[k].n; m = CSinfo[k].m; label = CSinfo[k].cs;

        if (m == 0)
        {
            x[i] = COEFG[n] + CSinfo[k].initv;
            COEFG[n] = x[i];
        }
        else
        {
            l = NMAX - m + 1;
            ind = NMAX + 1 + (2 * NMAX - m + 2) * (m - 1);
            ic = ind + n - m;
            is = ind + n - m + l;
            if (label == 1)
            {
                x[i] = COEFG[ic] + CSinfo[k].initv;
                COEFG[ic] = x[i];
            }

            if (label == -1)
            {
                x[i] = COEFG[is] + CSinfo[k].initv;
                COEFG[is] = x[i];
            }
        }
        i++;
    }



    if (MACC_BIAS != 0)
    {
        for (k = 0; k < MACC_PRBS; k ++)
        {
            x[i] = MACC_EPBS[k];
            i++;
        }
    }

    if (MACC_SCAL != 0)
    {
        for (k = 0; k < MACC_PRSL; k ++)
        {
            x[i] = MACC_EPSL[k];
            i++;
        }
    }


    return;

}






void updsolvefor (double *x)
{
    int i, k, n, m, ind, l, ic, is, label;

    i = 6;
    if (MOBS > 0)
    {
        BASB = x[i];
        i++;
    }
    if (MOBS > 1)
    {
        BAST = x[i];
        i++;
    }
    if (MSRP > 0)
    {
        SRPB = x[i];
        i++;
    }
    if (MSRP > 1)
    {
        SRPT = x[i];
        i++;
    }
    if (MTK2 > 0)
    {
        K2 = x[i];
        i++;
    }

    if (MGCS > 0)
    for (k = 0; k < MGCS; k ++)
    {
        n = CSinfo[k].n; m = CSinfo[k].m; label = CSinfo[k].cs;

        if (m == 0)
        {
            COEFG[n] = x[i];
        }
        else
        {
            l = NMAX - m + 1;
            ind = NMAX + 1 + (2 * NMAX - m + 2) * (m - 1);
            ic = ind + n - m;
            is = ind + n - m + l;
            if (label == 1)
            {
                COEFG[ic] = x[i];
            }

            if (label == -1)
            {
                COEFG[is] = x[i];
            }
        }
        i++;
    }


    if (MACC_BIAS != 0)
    {
        for (k = 0; k < MACC_PRBS; k ++)
        {
            MACC_EPBS[k]  = x[i];
            i++;
        }
    }

    if (MACC_SCAL != 0)
    {
        for (k = 0; k < MACC_PRSL; k ++)
        {
            MACC_EPSL[k] = x[i];
            i++;
        }
    }




    return;



}




/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/*
* simula_phase - simulate total phase count observable
* @param1: description of param1
* @param2: description of param2
* todo: 
        1 one-way doppler deltat accumlated error
* version: 20 Aug 2010
*/
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
double simula_phase (double utc3, double utc0, double *station3, 
                     short int uplink, double *station1, short int genrel,
                     double *calculable, double *azimuth, double *elevation, 
                     short int part, double *bmat)
{
    double txice[7], txics[7], *bmats, *bmate, deltat;
    int n;
    
    real128 lts2[3], lte2[3];

    bmats = (double *) calloc ( SLOVEFOR, sizeof(double));
    bmate = (double *) calloc ( SLOVEFOR, sizeof(double));


    ltsolution (utc0, station3, uplink, station1, genrel, lts2, 
        azimuth, elevation, part, bmats, txics);
    ltsolution (utc3, station3, uplink, station1, genrel, lte2, 
        azimuth, elevation, part, bmate, txice);
  
    *calculable = (lte2[2] - lts2[2]) * C;

    if (uplink == 0)    //one-way doppler deltat, time correction
    {
        delta_tdb (txice, txics, &deltat);
        *calculable = *calculable + deltat * (txice[0] - txics[0]) * C;
    }

    if (part == 1)
    {
        for (n = 0; n < 6 + DYNPAR; n++)
        {
            bmat[n] = bmate[n] - bmats[n];
        }
        bmat[8] = 1;
        bmat[9] = utc3;
    }
    *calculable = *calculable + BIAS + DBIA * utc3;
    
    free (bmats); 
    free (bmate);

    return 0;
}


/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/*
* simula_dople - simulate doppler observable
* @param1: description of param1
* @param2: description of param2
* todo: 
        1 one-way doppler deltat accumlated error
* version: 20 Aug 2010
*/
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
double simula_dople (double utc3, double tc, double *station3, 
                     short int uplink, double *station1, short int genrel,
                     double *calculable, double *azimuth, double *elevation, 
                     short int part, double *bmat)
{
    double txice[7], txics[7], *bmats, *bmate, deltat;
    double dop_old, dop_new;
    int n;
    real128 lts2[3], lte2[3], dlt;

    bmats = (double *) calloc ( SLOVEFOR, sizeof(double));
    bmate = (double *) calloc ( SLOVEFOR, sizeof(double));


    ltsolution (utc3 + tc / 2, station3, uplink, station1, genrel, lte2, 
        azimuth, elevation, part, bmate, txice);
    ltsolution (utc3 - tc / 2, station3, uplink, station1, genrel, lts2, 
        azimuth, elevation, part, bmats, txics);

    dop_old = (double) ((lte2[1] - lts2[1]) / (real128)tc);

    dlt = lte2[2] - lts2[2];
    dop_new = (double) (dlt / (real128)tc * (real128)C);
    *calculable = dop_new;

    if (uplink == 0)     //one-way doppler deltat, proper time correction
    {
        delta_tdb (txice, txics, &deltat);
        *calculable = *calculable + deltat * C;
    }

    if (part == 1)
    {
        for (n = 0; n < 6 + DYNPAR; n++)
        {
            bmat[n] = (bmate[n] - bmats[n]) / tc;
        }
        bmat[8] = 1;
        bmat[9] = utc3;
    }
    *calculable = *calculable + BIAS + DBIA * utc3;
    free (bmats); 
    free (bmate);
    return 0;
}



/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/*
* delta_tdb - one-way doppler deltat, proper time correction
* txice: [0]:  satellite TDB time(s), [1]~[7]satellite coordinates(AU, day)
* @param2: description of param2
*
* version: 20 Aug 2010
*/
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/

double delta_tdb (double *txice, double *txics, double *deltat)
{
    double ie, is, t, ied, isd, tjde[2], tjds[2], xcbe[6], xcbs[6],
        xibe[6], xibs[6];
    short int ssbary, n;
    
    ssbary = 11;
    t = txice[0] - txics[0];
    tjde[0] = JD0;
    tjde[1] = txice[0] / 86400.0;
    tjds[0] = JD0;
    tjds[1] = txics[0] / 86400.0;

    planet_ephemeris (tjde, CENTER, ssbary, &xcbe[0], &xcbe[3]);
    for (n = 0; n < 6; n++)
        xibe[n] = txice[n + 1] + xcbe[n];
    planet_ephemeris (tjds, CENTER, ssbary, &xcbs[0], &xcbs[3]);
    for (n = 0; n < 6; n++)
        xibs[n] = txics[n + 1] + xcbe[n];

    delta_iid (tjde, xibe, &ie, &ied);
    delta_iid (tjds, xibs, &is, &isd);

    *deltat = 1.0 / 2.0 * (ie + is) * t - 1.0 / 12.0 * (ied - isd) * t * t;
    return 0;
}




/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/*
* delta_tdb -  one-way light-time proper time correction
* txice: [0]:  satellite TDB time(s), [1]~[7]satellite coordinates(AU, day)
* @param2: description of param2
* todo: 
        1 Uobl, time correction due to non-spherical potential 
* version: 20 Aug 2010
*/
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/

double delta_iid (double *jd, double *xi, double *ii, double *id)
{
    double uu, vv, L, ud, vd,
    rdd[3], xj[11][6], xij[11][6], rij[11], rijd[11]; 
    short int ssbary, j, n;
    
    ssbary = 11;
    L = 1.550520e-8;

    uu = 0;
    ud = 0;
    for (n = 0; n < 3; n ++)
        rdd[n] = 0;
    for (j = 0; j <= 10; j++)
    {
        if (PERB[j] == 0)
            continue;
        planet_ephemeris (jd, j, ssbary, &xj[j][0], &xj[j][3]);
        for (n = 0; n < 6; n++)
        {
            xij[j][n] = xi[n] - xj[j][n];
        }
        rij[j] = sqrt (xij[j][0] * xij[j][0] 
            + xij[j][1] * xij[j][1] + xij[j][2] * xij[j][2]);
        rijd[j] = (xij[j][0] * xij[j][3] + xij[j][1] * xij[j][4] 
            + xij[j][2] * xij[j][5]) / rij[j];
        uu = uu + GMDE[j] / rij[j];
        ud = ud - GMDE[j] / rij[j] / rij[j] * rijd[j];
        for (n = 0; n < 3; n++)
            rdd[n] = rdd[n] 
            - GMDE[j] / (rij[j] * rij[j] * rij[j]) * xij[j][n];
    }

    vv = xi[3] * xi[3] + xi[4] * xi[4] + xi[5] * xi[5];
    vd = 2.0 * (xi[3] * rdd[0] + xi[4] * rdd[1] + xi[5] * rdd[2]);

    *ii = (uu + vv / 2.0) / C_AUDAY / C_AUDAY - L;
    *id = (ud + vd / 2.0) / C_AUDAY / C_AUDAY / 86400.0;

    return 0;
}


/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/*
* simula_range - simulate doppler observable
* @param1: description of param1
* @param2: description of param2
* todo: 

* version: 20 Aug 2010
*/
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
double simula_range (double utc3, double *station3, short int uplink, 
                     double *station1, short int genrel, 
                     double *calculable, double *azimuth, double *elevation, 
                     short int part, double *bmat)
{
    double txic[7];
    real128 lt[3];

    ltsolution (utc3, station3, uplink, station1, genrel, lt, 
        azimuth, elevation, part, bmat, txic);

    if (uplink == 0)    //one-way range: only work for near-earth satellite
    {
        *calculable = (double) lt[0];
    }
    if (uplink == 1)    //two/three-way range: work for deep space
    {
        *calculable = (double) (lt[2] * (real128)C);
    }

    if (part == 1)
    {
        bmat[8] = 1;    //DBIAS
        bmat[9] = utc3;
    }
    *calculable = *calculable + BIAS + DBIA * utc3;
    return 0;
}






/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/*

* ltsolution - light time solution
* @param:
    utc3        : unit: day; 
    station3    : receive station;
    uplink      : no uplink == 0; yes uplink == 1 
    station1    : transmit station;		
    genrel      : no general relativity correction == 0; yes == 1
    calculable  : ;
    azimuth     : ;
    elevation   : ;
    partial     : no partial == 0 ; yes partial == 1 
    bmat        : (partial == 0: satellite coordinates(6), partial == 1: partial)
    txic        : satellite coordinates(t2)

* 
* version: 20 Aug 2010

*/
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/

double ltsolution (double utc_3, double *station3, short int uplink,  
                   double *station1, short int genrel, real128 *lt, 
                   double *azimuth, double *elevation, short int part,
                   double *bmat, double *txic)
{
    double re3fi[3], re3[3], re1fi[3], re1[3], vec[3], xe3[6], re3n[3], 
        xe1[6], re1n[3], xc2[6], 
        secdiff, 
        ra, dec, zd, az, secdiff3, secdiff1 = 0,
        utc_1 = 0, ut1_3, ut1_1, tt_3, tt_1, tdb_3, tdb_2, tdb_1, 
        ut1_utc, xp, yp, xp3, yp3, dx, dy, delta_t3, t, elong, u, v, 
        dxdx0[36], dxdp[6 * DYNPAR], dodx[6], dodp[DYNPAR], 
        dodx0[6], dodpp[DYNPAR], eph[42 + 6 * DYNPAR], 
        te[9], llh3[3], llh1[3];
    real128 tao231, tao232, tao121 = 0, tao122, taoerr, r23,  r12 = 0, xb3[6], xb2[6], xb1[6];
    int n, flag, dim_par;
	
    taoerr   = 1.0e-12L; //1nanosec;
//    taoerr   = 1.0e-8L; //1nanosec;
    dim_par = 42 + 6 * DYNPAR + 1;

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/*--light time iteration 2 -> 3 --*/

    re3fi[0] = station3[0];
    re3fi[1] = station3[1];
    re3fi[2] = station3[2];
    xyz2llh(re3fi, llh3);

    elong = llh3[1] * DEG2RAD;
    u = sqrt (re3fi[0] * re3fi[0] + re3fi[1] * re3fi[1]) / 1000.0;
    v = re3fi[2] / 1000.0;

/*--time scales transformation --*/
    geteop (utc_3, &xp3, &yp3, &ut1_utc, &dx, &dy);	
    delta_t3 = 32.184 + LEAPSECS - ut1_utc;
    ut1_3 = utc_3 + ut1_utc;
    tt_3 = utc_3 + (LEAPSECS + 32.184); 	
    secdiff3 = iauDtdb (JD0, tt_3 / 86400.0, ut1_3 / 86400.0, elong, u, v);
    tdb_3 = tt_3 + secdiff3;

/*--station coordinate interpolation--*/    
    lagrange (TE_EPH, DIM_TE, 10, utc_3, te);        
    brmul (te, re3fi, 3, 3, 1, re3);
    lagrange (TE_EPH, DIM_TE, 10, utc_3 + 1.0, te);        
    brmul (te, re3fi, 3, 3, 1, re3n);
    for (n = 0; n < 3; n++)
    {
        xe3[n] = re3[n] / AU;
        xe3[n + 3] = (re3n[n] - re3[n]) / AU * 86400;
    }

/*--satellite coordinate interpolation--*/    
    if (part == 0)
    {
        lagrange (OR_EPH, DIM_OR, 7, tdb_3, xc2);
    }
    if (part == 1)
    {
        lagrange (OR_EPH, DIM_OR, dim_par, tdb_3, eph);
        for (n = 0; n < 6; n++)
            xc2[n] = eph[n + 36];
    }

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/*--iteration--*/    
    r23 = lt_form (tdb_3, tdb_3, xe3, xc2, genrel, xb3, xb2);
    tao231 = (real128)r23 * (real128)AU_SEC;
    tdb_2 = tdb_3 - (double) tao231;

    flag = -1;
    do
    {
        flag++;
        tao232 = tao231;
        if (part == 0)
        {
            lagrange (OR_EPH, DIM_OR, 7, tdb_2, xc2);
        }
        if (part == 1)
        {
            lagrange (OR_EPH, DIM_OR, dim_par, tdb_2, eph);
            for (n = 0; n < 6; n++)
                xc2[n] = eph[n + 36];
        }
        
        r23 = lt_form (tdb_3, tdb_2, xe3, xc2, genrel, xb3, xb2);
        tao231 = (real128)r23 * (real128)AU_SEC;
        tdb_2 = tdb_3 - (double) tao231;
    }while (fabsl (tao232-tao231) > taoerr);



/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/*--light time iteration 1 -> 2--*/    


    if (uplink == 1)
    {
        re1fi[0] = station1[0];
        re1fi[1] = station1[1];
        re1fi[2] = station1[2];
        xyz2llh(re1fi, llh1);

        elong = llh1[1] * DEG2RAD;
        u = sqrt (re1fi[0] * re1fi[0] + re1fi[1] * re1fi[1]) / 1000.0;
        v = re1fi[2] / 1000.0;

    /*--time scales transformation --*/
        tdb_1 = tdb_2; //unit: s
        tdb2tt (JD0 + tdb_1 / 86400.0, &t, &secdiff);
        tt_1 = tdb_1 - secdiff;
        utc_1 = tt_1 - (LEAPSECS + 32.184); 	
        geteop (utc_1, &xp, &yp, &ut1_utc, &dx, &dy);	
        ut1_1 = utc_1 + ut1_utc;	        
        secdiff1 = iauDtdb (JD0, tdb_1 / 86400.0, ut1_1 / 86400.0, 
            elong, u, v);
        tt_1 = tdb_1 - secdiff1;
        utc_1 = tt_1 - (LEAPSECS + 32.184); 	

    /*--station coordinate interpolation--*/    
        lagrange (TE_EPH, DIM_TE, 10, utc_1, te);        
        brmul (te, re1fi, 3, 3, 1, re1);
        lagrange (TE_EPH, DIM_TE, 10, utc_1 + 1.0, te);        
        brmul (te, re1fi, 3, 3, 1, re1n);
        for (n = 0; n < 3; n++)
        {
            xe1[n] = re1[n] / AU;
            xe1[n + 3] = (re1n[n] - re1[n]) / AU * 86400;
        }

        r12 = lt_form (tdb_1, tdb_2, xe1, xc2, genrel, xb1, xb2);
        tao121 = (real128)r12 * (real128)AU_SEC;
        tdb_1 = tdb_2 - (double) tao121;
     
    /*--iteration--*/    
        flag = -1;
        do
        {
            flag++;
            tao122 = tao121;

        /*--time scales transformation --*/
            tdb2tt (JD0 + tdb_1 / 86400.0, &t,&secdiff);
            tt_1 = tdb_1 - secdiff;
            utc_1 = tt_1 - (LEAPSECS + 32.184); 
            geteop (utc_1, &xp, &yp, &ut1_utc, &dx, &dy);	

            ut1_1 = utc_1 + ut1_utc;
            secdiff1 = iauDtdb (JD0, tdb_1 / 86400.0, ut1_1 / 86400.0, 
                elong, u, v);
            tt_1 = tdb_1 - secdiff1;
            utc_1 = tt_1 - (LEAPSECS + 32.184); 	

        /*--station coordinate interpolation--*/    
            lagrange (TE_EPH, DIM_TE, 10, utc_1, te);        
            brmul (te, re1fi, 3, 3, 1, re1);
            lagrange (TE_EPH, DIM_TE, 10, utc_1 + 1.0, te);        
            brmul (te, re1fi, 3, 3, 1, re1n);
            for (n = 0; n < 3; n++)
            {
                xe1[n] = re1[n] / AU;
                xe1[n + 3] = (re1n[n] - re1[n]) / AU * 86400;
            }
            
            r12 = lt_form (tdb_1, tdb_2, xe1, xc2, genrel, xb1, xb2);
            tao121 = (real128)r12 * (real128)AU_SEC;
            tdb_1 = tdb_2 - (double) tao121;
        }while (fabsl (tao122-tao121) > taoerr);
    }

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/*-- partial time: tdb2 --*/    
    if (part == 0)
    {
        lagrange (OR_EPH, DIM_OR, 7, tdb_2, bmat);

    }
    if (part == 1)
    {
        lagrange (OR_EPH, DIM_OR, dim_par, tdb_2, eph);

        for (n = 0; n < 36; n++)
            dxdx0[n] = eph[n];
        for (n = 0; n < 6 * DYNPAR; n++)
            dxdp[n] = eph[n + 42];
        lt_part (xb3, xb2, xb1, uplink, dodx, dodp);
        brmul (dodx, dxdx0, 1, 6, 6, dodx0);
        brmul (dodx, dxdp, 1, 6, DYNPAR, dodpp);
        for (n = 0; n < DYNPAR; n++)
            dodp[n] = dodpp[n] + dodp[n];

        for (n = 0; n < 3; n++)
        {
            bmat[n] = dodx0[n];
            bmat[n + 3] = dodx0[n + 3] * 86400.0;
        }
        bmat[6] = dodp[0] * AU; // l/c: au, 
        bmat[7] = dodp[1] * AU * 86400.0; // l/(c*d-1): au*d, 
    }

    txic[0] = tdb_2;
    for (n = 1; n < 7; n++)	
    {
        txic[n] = xc2[n - 1];
    }

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/*--calculate light time observable --*/
    if (uplink == 0)
    {
        lt[0] = (real128)r23 * (real128)AU;			
        lt[1] = ((real128)utc_3 - (real128)tdb_2) * (real128)C;
        lt[2] = (real128)tao231 - (real128)secdiff3 - (real128)(LEAPSECS + 32.184);
    }

    if (uplink == 1)
    {
        lt[0] = ((real128)r12 + (real128)r23) * (real128)AU;			
        lt[1] = ((real128)utc_3 - (real128)utc_1) * (real128)C;			
        lt[2] = (real128)tao231 + (real128)tao121 + (real128)secdiff1 - (real128)secdiff3;			
//        lt[2] = tao231 + tao121 + secdiff1 - secdiff3;			
    }

    for (n = 0; n < 3; n++)
        vec[n] = xb2[n] - xb3[n];
    vector2radec (vec, &ra,&dec);

    azelev (ut1_3 / 86400.0 + JD0, delta_t3, ACCURACY,
               xp3, yp3, llh3, ra, dec, &zd, &az);

    *azimuth = az;
    *elevation = 90.0 - zd;

    return 0;
}



/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/*
* lt_part - partial of light time 
* @param1: description of param1
* @param2: description of param2
*
* version: 20 Aug 2010
*/
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
double lt_part (real128 *xb3, real128 *xb2, real128 *xb1, int uplink, 
                 double *dodx, double *dodp)
{
    double r23, r12, p23, p12, pt2, pt1, rp12, rpt1, dt2dx[3], dt1dx[3];
    int n;
    
    r23 = sqrt ((xb2[0] - xb3[0]) * (xb2[0] - xb3[0]) 
        + (xb2[1] - xb3[1]) * (xb2[1] - xb3[1])
        + (xb2[2] - xb3[2]) * (xb2[2] - xb3[2]));

    p23 = ((xb3[0] - xb2[0]) * xb2[3] + (xb3[1] - xb2[1]) * xb2[4]
        + (xb3[2] - xb2[2]) * xb2[5]) / r23;

    pt2 = (1 - p23 / C_AUDAY);

    dt2dx[0] = (xb3[0] - xb2[0]) / r23 / pt2;
    dt2dx[1] = (xb3[1] - xb2[1]) / r23 / pt2;
    dt2dx[2] = (xb3[2] - xb2[2]) / r23 / pt2;

    if (uplink == 0)
    {
        dodx[0] = - dt2dx[0];
        dodx[1] = - dt2dx[1];
        dodx[2] = - dt2dx[2];
    }
    if (uplink == 1)
    {
        r12 = sqrt ((xb2[0] - xb1[0]) * (xb2[0] - xb1[0]) 
            + (xb2[1] - xb1[1]) * (xb2[1] - xb1[1])
            + (xb2[2] - xb1[2]) * (xb2[2] - xb1[2]));

        p12 = ((xb2[0] - xb1[0]) * xb1[3] + (xb2[1] - xb1[1]) * xb1[4]
            + (xb2[2] - xb1[2]) * xb1[5]) / r12;
        
        pt1 = (1 - p12 / C_AUDAY);

        rp12 = ((xb2[0] - xb1[0]) * (xb2[3] - xb1[3]) 
              + (xb2[1] - xb1[1]) * (xb2[4] - xb1[4]) 
              + (xb2[2] - xb1[2]) * (xb2[5] - xb1[5])) / r12;

        rpt1 = (1 - (rp12 + p12) / C_AUDAY);


        for (n = 0; n < 3; n++)
        {
            dt1dx[n] = (dt2dx[n] * rpt1 - (xb2[n] - xb1[n]) / r12) / pt1;
        }

        for (n = 0; n < 3; n++)
        {
            dodx[n] = - dt1dx[n];
        }

    }

    dodx[3] = 0; 
    dodx[4] = 0; 
    dodx[5] = 0;
    
    dodp[0] = 0; //
    dodp[1] = 0; //
    
    return 0;
}



/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/*
* lt_form - calculate the light time equation
* @param:
    tdb3, re3[3]    : station time, coordinates (unit:AU)
    tdb2, rp2[3]    : satellite time, coordinates (AU)
    genrel          : 
    *rs3            : output: station coordinates to SSB (AU)
    *rs2            : output: satellite coordinates to SSB (AU)
    return          : light time solution (AU)
* version: 20 Aug 2010
*/
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
real128 lt_form (double tdb3, double tdb2, double *re3, double *rp2, 
                int genrel, real128 *rs3, real128 *rs2)
{
    double gamma, tjd2[2], tjd3[2], re[3], ve[3], 
        rp[3], vp[3], xe[6], xp[6];
    real128 rlight, rgen, ri, rj, r12, rse[3], rsp[3], r23[3], rlt[3];
    short int earth = 2, sun = 10, j;
    int n;
    gamma = 1;
    tjd2[0] = JD0;
    tjd3[0] = JD0;
    tjd2[1] = tdb2 / 86400.0;
    tjd3[1] = tdb3 / 86400.0;

    planet_ephemeris (tjd3, earth, sun, &xe[0], &xe[3]);
    for (n = 0; n < 6; n++)
        rs3[n] = (real128)re3[n] + (real128)xe[n];
    ri = sqrtl (rs3[0] * rs3[0] + rs3[1] * rs3[1] + rs3[2] * rs3[2]);

    planet_ephemeris (tjd2, CENTER, sun, &xp[0], &xp[3]);
    for (n = 0; n < 6; n++)
        rs2[n] = (real128)rp2[n] + (real128)xp[n];
    rj = sqrtl (rs2[0] * rs2[0] + rs2[1] * rs2[1] + rs2[2] * rs2[2]);


    for (n = 0; n < 3; n++)
        r23[n] = (real128)xe[n] - (real128)xp[n];
    for (n = 0; n < 3; n++)
        rlt[n] = (real128)r23[n] + (real128)re3[n] - (real128)rp2[n];


//    rlight = sqrtl (((real128)rs3[0] - (real128)rs2[0]) * ((real128)rs3[0] - (real128)rs2[0]) 
//        + ((real128)rs3[1] - (real128)rs2[1]) * ((real128)rs3[1] - (real128)rs2[1])
//        + ((real128)rs3[2] - (real128)rs2[2]) * ((real128)rs3[2] - (real128)rs2[2]));

     rlight = sqrtl( rlt[0] * rlt[0] +  rlt[1] * rlt[1] +  rlt[2] * rlt[2]);
     rlight = sqrtl ((rs3[0] - rs2[0]) * (rs3[0] - rs2[0]) 
        + (rs3[1] - rs2[1]) * (rs3[1] - rs2[1])
        + (rs3[2] - rs2[2]) * (rs3[2] - rs2[2]));


    if (genrel == 1)
    {
        rgen = (1L + (real128)gamma) * (real128)GMDE[10] / (real128)C_AUDAY / (real128)C_AUDAY 
            * logl ((ri + rj + rlight 
            + (1L + (real128)gamma) * (real128)GMDE[10] / (real128)C_AUDAY / (real128)C_AUDAY) 
            / (ri + rj - rlight 
            + (1L + (real128)gamma) * (real128)GMDE[10] / (real128)C_AUDAY / (real128)C_AUDAY));
        rlight = rlight + rgen;
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
        for (j = 0; j <= 9; j++)
        {
//			if (j ==9) continue;		
            planet_ephemeris (tjd3, earth, j, re, ve);
            for (n = 0; n < 3; n++)
                rse[n] = (real128)re3[n] + (real128)re[n];
            ri = sqrtl (rse[0] * rse[0] + rse[1] * rse[1] + rse[2] * rse[2]);

            planet_ephemeris (tjd2, CENTER, j, rp, vp);
            for (n = 0; n < 3; n++)
                rsp[n] = (real128)rp2[n] + (real128)rp[n];
            rj = sqrtl (rsp[0] * rsp[0] + rsp[1] * rsp[1] + rsp[2] * rsp[2]);
	
            r12 = sqrtl((rse[0] - rsp[0]) * (rse[0] - rsp[0]) 
                + (rse[1] - rsp[1]) * (rse[1] - rsp[1])
                + (rse[2] - rsp[2]) * (rse[2] - rsp[2]));
            rgen = (1L + (real128)gamma) * (real128)GMDE[j] / (real128)C_AUDAY / (real128)C_AUDAY 
                * logl ((ri + rj + r12 ) / (ri + rj - r12));
            rlight = rlight + rgen;
        }
    }
    return rlight;
}





void azelev (double jd_ut1, double delta_t, short int accuracy,
              double x, double y, double *llh, double ra,
              double dec, double *zd, double *az)

{

   double sinlat, coslat, sinlon, coslon, sindc, cosdc, sinra, cosra,
      uze[3], une[3], uwe[3], uz[3], un[3], uw[3], p[3], pz, pn, pw,
      proj;

/*
   Preliminaries.
*/

   sinlat = sin (llh[0] * DEG2RAD);
   coslat = cos (llh[0] * DEG2RAD);
   sinlon = sin (llh[1] * DEG2RAD);
   coslon = cos (llh[1] * DEG2RAD);
   sindc = sin (dec * DEG2RAD);
   cosdc = cos (dec * DEG2RAD);
   sinra = sin (ra * 15.0 * DEG2RAD);
   cosra = cos (ra * 15.0 * DEG2RAD);

/*
   Set up orthonormal basis vectors in local Earth-fixed system.

   Define vector toward local zenith in Earth-fixed system (z axis).
*/
   uze[0] = coslat * coslon;
   uze[1] = coslat * sinlon;
   uze[2] = sinlat;

/*
   Define vector toward local north in Earth-fixed system (x axis).
*/

   une[0] = -sinlat * coslon;
   une[1] = -sinlat * sinlon;
   une[2] = coslat;

/*
   Define vector toward local west in Earth-fixed system (y axis).
*/

   uwe[0] = sinlon;
   uwe[1] = -coslon;
   uwe[2] = 0.0;

/*
   Obtain vectors in celestial system.

   Rotate Earth-fixed orthonormal basis vectors to celestial system
   (wrt equator and equinox of date).
*/

   ter2cel (jd_ut1,0.0,delta_t,1,accuracy,1,x,y,uze, uz);
   ter2cel (jd_ut1,0.0,delta_t,1,accuracy,1,x,y,une, un);
   ter2cel (jd_ut1,0.0,delta_t,1,accuracy,1,x,y,uwe, uw);

/*
   Define unit vector 'p' toward object in celestial system
   (wrt equator and equinox of date).
*/

   p[0] = cosdc * cosra;
   p[1] = cosdc * sinra;
   p[2] = sindc;

/*
   Compute coordinates of object wrt orthonormal basis.

   Compute components of 'p' - projections of 'p' onto rotated
   Earth-fixed basis vectors.
*/

   pz = p[0] * uz[0] + p[1] * uz[1] + p[2] * uz[2];
   pn = p[0] * un[0] + p[1] * un[1] + p[2] * un[2];
   pw = p[0] * uw[0] + p[1] * uw[1] + p[2] * uw[2];

/*
   Compute azimuth and zenith distance.
*/

   proj = sqrt (pn * pn + pw * pw);

   if (proj > 0.0)
      *az = -atan2 (pw, pn) * RAD2DEG;

   if (*az < 0.0)
      *az += 360.0;

   if (*az >= 360.0)
      *az -= 360.0;

   *zd = atan2 (proj, pz) * RAD2DEG;
}





/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/*
* fun_pointmass - abandoned
* @param1: description of param1
* @param2: description of param2
*
* version: 20 Aug 2010
*/
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/

double fun_pointmass (double tdbs, double *x, double *f)
{
    double fnt[3], fgr[3], r, s2, rrd, a, b;
    int n, gamma;

    gamma = 1;

    f[0] = x[3]; 
    f[1] = x[4]; 
    f[2] = x[5];

    r = sqrt (x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
    s2 = x[3] * x[3] + x[4] * x[4] + x[5] * x[5];
    rrd = x[0] * x[3] + x[1] * x[4] + x[2] * x[5];
	
    a = 2 * (1 + gamma) * GMDE[CENTER] / r - gamma * s2;
    b = 2 * (1 + gamma) * rrd;
    for (n = 0; n < 3; n++)
        fgr[n] =  GMDE[CENTER] / C_AUDAY / C_AUDAY / r / r / r 
        * ( a * x[n] + b * x[n+3] );

    fnt[0] = - GMDE[CENTER] / (r*r*r) * x[0];
    fnt[1] = - GMDE[CENTER] / (r*r*r) * x[1];
    fnt[2] = - GMDE[CENTER] / (r*r*r) * x[2];

    for (n = 0; n < 3; n++)
    {
        f[3 + n] = fnt[n] + fgr[n]; 
    }

	return 0;
}



/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/*
* fun_fullaccel - 
* @param1: description of param1
* @param2: description of param2
*
* version: 20 Aug 2010
*/
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
double fun_fullaccel (double tdbs, double *xic, double *fxic)
{
    int n;
    short int part = 0;
    double tjd[2], acc1[3], acc2[3], acc3[3], acc[3], dum1[1], dum2[1];

    tjd[0] = JD0;
    tjd[1] = tdbs;

    accel_ntrel (tjd, xic, part, acc1, dum1, dum2);
    accel_nonsp (tjd, xic, part, acc2, dum1, dum2);
    accel_radpr (tjd, xic, part, acc3, dum1, dum2);

    for (n = 0; n <= 2; n++)
    {
        acc[n] = acc1[n] + acc2[n] + acc3[n];
    }

    fxic[0] = xic[3];
    fxic[1] = xic[4];
    fxic[2] = xic[5];
    fxic[3] = acc[0];
    fxic[4] = acc[1];
    fxic[5] = acc[2];

    return 0;
}


/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/*
* fun_fullstate -transition matrix(36), orbit(6), sensitivity matrix(6*DYNPAR)
* @param1: description of param1
* @param2: description of param2
*
* version: 20 Aug 2010
*/
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
double fun_fullstate (double tdbs, double *state, double *fstate)
{
    int n;
    short int part = 1;
    double tjd[2], xic[6], dfdx[36], dxdx0[36], dfdp[6 * DYNPAR], 
        dfdpp[6 * DYNPAR], dxdp[6 * DYNPAR],
        acc1[3], dadr1[9], dadp1[3 * DYNPAR],
        acc2[3], dadr2[9], dadp2[3 * DYNPAR],
        acc3[3], dadr3[9], dadp3[3 * DYNPAR],
        acc[3], dadr[9], dadp[3 * DYNPAR],
        fxic[6], fdxdx0[36], fdxdp[6 * DYNPAR];

    tjd[0] = JD0;
    tjd[1] = tdbs;

    for (n = 0; n < 36; n++)
    {
        dxdx0[n] = state[n];
    }
    for (n = 0; n < 6; n++)
    {
        xic[n] = state[n + 36];
    }
    for (n = 0; n < 6 * DYNPAR; n++)
    {
        dxdp[n] = state[n + 42];
    }

/* acc, partial to xyz: dadr, partial to parameters dadp*/
    accel_ntrel (tjd, xic, part, acc1, dadr1, dadp1);
    accel_nonsp (tjd, xic, part, acc2, dadr2, dadp2);
    accel_radpr (tjd, xic, part, acc3, dadr3, dadp3);
/*todo: air drag acc & partial to vxvyvz dadv*/

    for (n = 0; n <= 2; n++)
    {
        acc[n] = acc1[n] + acc2[n] + acc3[n];
    }
    for (n = 0; n <= 8; n++)
    {
        dadr[n] = dadr1[n] + dadr2[n] + dadr3[n];
    }
    for (n = 0; n < 3 * DYNPAR; n++)
    {
        dadp[n] = dadp1[n] + dadp2[n] + dadp3[n];
    }

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/

    for (n = 0; n < 36; n++)
    {
        dfdx[n] = 0;
    }
    dfdx[3]  = 1; 
    dfdx[10] = 1; 
    dfdx[17] = 1;
    for (n = 0; n < 3; n++)
    {
        dfdx[n + 18] = dadr[n];
        dfdx[n + 24] = dadr[n + 3];
        dfdx[n + 30] = dadr[n + 6];
    }
    brmul(dfdx, dxdx0, 6, 6, 6, fdxdx0);

    fxic[0] = xic[3];
    fxic[1] = xic[4];
    fxic[2] = xic[5];
    fxic[3] = acc[0];
    fxic[4] = acc[1];
    fxic[5] = acc[2];

    brmul(dfdx, dxdp, 6, 6, DYNPAR, dfdpp);
    for (n = 0; n < 3 * DYNPAR; n++)
    {
        dfdp[n] = 0;
        dfdp[n + 3 * DYNPAR] = dadp[n];
    }
    for (n = 0; n < 6 * DYNPAR; n++)
    {
        fdxdp[n] = dfdpp[n] + dfdp[n];
    }

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/

    for (n = 0; n < 36; n++)
    {
        fstate[n] = fdxdx0[n];
    }
    for (n = 0; n < 6; n++)
    {
        fstate[n + 36] = fxic[n];
    }
    for (n = 0; n < 6 * DYNPAR; n++)
    {
        fstate[n + 42]= fdxdp[n];
    }

    return 0;
}




/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/*
* accel_ntrel - Newtonian + Relativistic acceleration
* @param1: description of param1
* @param2: description of param2
*
* version: 20 Aug 2010
*/
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
double accel_ntrel (double *tjd, double *xic, short int part, 
                    double *acc, double *dadr, double *dadp)
{
    int n;
    short int ssbary = 11;
    double xcb[6], acb[3], xib[6], aib[3], dadr1[9], dum[9], 
        dadp0[3*DYNPAR], dadp1[3*DYNPAR];

    planet_ephemeris (tjd, CENTER, ssbary, &xcb[0], &xcb[3]);
    accel_bcrs (tjd, xcb, part, CENTER, acb, dum, dadp0);
    for (n = 0; n <= 5; n++)
    {
        xib[n] = xic[n] + xcb[n];
    }
    accel_bcrs (tjd, xib, part, 99, aib, dadr1, dadp1);
    for (n = 0; n <= 2; n++)
    {
        acc[n] = aib[n] - acb[n];
    }

    if (part == 1)
    {
        for (n = 0; n <= 8; n++)
        {
            dadr[n] = dadr1[n];
        }
        for (n = 0; n <= 3 * DYNPAR - 1; n++)
        {
            dadp[n] = dadp1[n] - dadp0[n];
        }
    }
    return 0;
}






/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/*
* accel_bcrs - Newtonian + Relativistic acceleration
* @param1: description of param1
* @param2: description of param2
* todo: 
        1 partial to parameters
* version: 20 Aug 2010
*/
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
double accel_bcrs (double *jd, double *xi, short int part, short int exclude, 
                   double *acc, double *dadr, double *dadp)
{
    double fnt[3], fgr[3], xj[11][6], xij[11][6], rij[11], xjk[6], rjk, 
        xddj[3], sumil, sumjk, sdi2, sdj2, rdirdj, rrrdr2, rjirdd, 
        rij5, rij3, xijt[9], gra, grb, beta, gamma, unit[9];
    short int ssbary, l, k, j, n, flag_gr;
    
    ssbary = 11;
    gamma = 1.0;
    beta = 1.0;
    unit[0] = 1; unit[1] = 0; unit[2] = 0;
    unit[3] = 0; unit[4] = 1; unit[5] = 0;
    unit[6] = 0; unit[7] = 0; unit[8] = 1;

    for (j = 0; j <= 10; j++)
    {
        planet_ephemeris (jd, j, ssbary, &xj[j][0], &xj[j][3]);
        for (n = 0; n < 6; n++)
        {
            xij[j][n] = xi[n] - xj[j][n];
        }
        rij[j] = sqrt (xij[j][0] * xij[j][0] 
            + xij[j][1] * xij[j][1] + xij[j][2] * xij[j][2]);
    }
    
    flag_gr = 0;
    for (n = 0; n < 3; n ++)
        fnt[n] = 0;
    for (j = 0; j <= 10; j++)
    {
        if (PERB[j] == 2)
            flag_gr = 1;
        if (PERB[j] == 0)
            continue;
        if (j == exclude)
            continue;
        for (n = 0; n < 3; n++)
            fnt[n] = fnt[n] 
            - GMDE[j] / (rij[j] * rij[j] * rij[j]) * xij[j][n];
    }

    if (part == 1)
    {
        for (n = 0; n <= 3 * DYNPAR - 1; n++)
        {
            dadp[n] = 0;
        }
        for (n = 0; n <= 8; n++)
        {
            dadr[n] = 0;
        }
        for (j = 0; j <= 10; j++)
        {
            if (j == exclude)
                continue;
            rij5 = pow (rij[j], 5);
            rij3 = pow (rij[j], 3);
            brmul (xij[j], xij[j], 3,1,3, xijt);
            for (n = 0; n <= 8; n++)
            {
                dadr[n] = dadr[n] + 3 * GMDE[j] * xijt[n] / rij5
                - GMDE[j] * unit[n] / rij3;
            }
        }
    }

    if (flag_gr == 0)
    {
        for (n = 0; n < 3; n++)
            acc[n] =  fnt[n];
        return 0;
    }

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/

    sdi2 = xi[3] * xi[3] + xi[4] * xi[4] + xi[5] * xi[5];
    sumil = 0;
    for (l = 0; l < 11; l ++)
    {
        if ( l == exclude)
            continue;
        if (PERB[l] != 2)
            continue;
        sumil = sumil + GMDE[l] / rij[l];
    }

    for (n = 0; n < 3; n ++)
        fgr[n] = 0;
    for (j = 0; j < 11; j ++)
    {
        if (PERB[j] != 2)
            continue;
        if (j == exclude)
            continue;
        sumjk = 0;
        for (n = 0; n < 3; n ++)
            xddj[n] = 0;
        for (k = 0; k < 11; k ++)
        {
            if (k == j)	
                continue;	//k!=j
            if (PERB[k] != 2)
                continue;
            for (n = 0; n < 3; n++)
                xjk[n] = xj[j][n] - xj[k][n];
            rjk = sqrt (xjk[0] * xjk[0] + xjk[1] * xjk[1] + xjk[2] * xjk[2]);
            sumjk = sumjk + GMDE[k] / rjk;
            for (n = 0; n < 3; n ++)
                xddj[n] = xddj[n] - GMDE[k] / (rjk * rjk * rjk) * xjk[n];
        }
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
        sdj2 = xj[j][3] * xj[j][3] + xj[j][4] * xj[j][4] 
            + xj[j][5] * xj[j][5];
        rdirdj = xi[3] * xj[j][3] + xi[4] * xj[j][4] + xi[5] * xj[j][5];
        rrrdr2 = pow( ( xij[j][0] * xj[j][3] + xij[j][1] * xj[j][4] 
            + xij[j][2] * xj[j][5]) / rij[j], 2);
        rjirdd = - ( xij[j][0] * xddj[0] + xij[j][1] * xddj[1] 
            + xij[j][2] * xddj[2]);
        
        gra = - 2 * (beta + gamma) * sumil - (2 * beta -1) * sumjk 
            + gamma * sdi2 + (1 + gamma) * sdj2
            - 2 * (1 + gamma) * rdirdj - 1.5 * rrrdr2 + 0.5 * rjirdd;

        grb = xij[j][0] * ((2+2*gamma) * xi[3] - (1+2*gamma) * xj[j][3])
            + xij[j][1] * ((2+2*gamma) * xi[4] - (1+2*gamma) * xj[j][4])
            + xij[j][2] * ((2+2*gamma) * xi[5] - (1+2*gamma) * xj[j][5]);

        for (n = 0; n < 3; n ++)
        {
            fgr[n] = fgr[n] 
                + GMDE[j] / (rij[j] * rij[j] * rij[j]) 
                * ( - xij[j][n]) * gra / C_AUDAY / C_AUDAY
                + GMDE[j] / (rij[j] * rij[j] * rij[j]) 
                * xij[j][n + 3] * grb / C_AUDAY / C_AUDAY
                + GMDE[j] / rij[j] * (3 + 4 * gamma) * 0.5   
                * xddj[n] / C_AUDAY / C_AUDAY;
        }
    }

    for (n = 0; n < 3; n++)
        acc[n] =  fgr[n] + fnt[n];
    return 1;
}




/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/*
* accel_radpr - solar radiation press & partial to srp coefficients 
* @param1: description of param1
* @param2: description of param2
* todo: 
        1 earth shadow 
* version: 20 Aug 2010
*/
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
double accel_radpr (double *tjd, double *xic, short int part, 
                    double *acc, double *dadr, double *dadp)
{
    double j, c1, ap, m, rsp, usp[3], xis[6], xsc[6], f, 
        xist[9], unit[9], rsp3;
    short int n, sun;

    sun = 10;
    unit[0] = 1; unit[1] = 0; unit[2] = 0;
    unit[3] = 0; unit[4] = 1; unit[5] = 0;
    unit[6] = 0; unit[7] = 0; unit[8] = 1;

    planet_ephemeris (tjd, sun, CENTER, &xsc[0], &xsc[3]);
    for (n = 0; n <= 5; n++)
    {
        xis[n] = xic[n] - xsc[n];
    }
    rsp = sqrt (xis[0] * xis[0] + xis[1] * xis[1] + xis[2] * xis[2]);
    usp[0] = xis[0] / rsp;
    usp[1] = xis[1] / rsp;
    usp[2] = xis[2] / rsp;

    j  = 1352.5;   //kg/s3
    m  = SATMASS;     //kg
    ap = SATAREA;        //m2
    c1 = j / C * 1 * 1;   //kg/s2/m*au*au
    f  = c1 * ap / m  / rsp / rsp;    
//kg/s2/m*au*au * m2 / kg / au  / au = m/s2
    f = f / AU * 86400.0 * 86400.0;

//    acc[0] = f * usp[0] * (1 + CONS + DCON * tjd[1]);
//    acc[1] = f * usp[1] * (1 + CONS + DCON * tjd[1]);
//    acc[2] = f * usp[2] * (1 + CONS + DCON * tjd[1]);
    acc[0] = f * usp[0];
    acc[1] = f * usp[1];
    acc[2] = f * usp[2]; 

    if (part == 0)
        return 1;

    rsp3 = rsp * rsp * rsp;
    brmul (xis, xis, 3,1,3, xist);
    for (n = 0; n <= 8; n++)
        dadr[n] = - f * (3 * xist[n] / rsp3 - unit[n] / rsp) ;
//        * (1 + CONS + DCON * tjd[1]);

//    for (n = 0; n <= 2; n++)
//    {
//        dadp[n * DYNPAR] = f * usp[n];      
//        dadp[n * DYNPAR + 1] = dadp[n * DYNPAR] * tjd[1];   
//        dadp[n * DYNPAR] = 0;      
//        dadp[n * DYNPAR + 1] = 0;   

//    }

        for (n = 0; n <= 3 * DYNPAR - 1; n++)
        {
            dadp[n] = 0;
        }


    return 0;
}






//    nmax = 4;      
//    stcs = (double *) calloc ( (nmax + 1) * (nmax + 1), sizeof(double));
double stidecs_old(double *tjd, double gma1, double k2, 
               double *c20, double *c21, double *s21, double *c22, double *s22)
{

    double gms2e;                   
    double gmm2e;
    

//    short int moon = 9, earth = 2, sun = 10;
    short int moon = 2, sun = 10;

    double ps[3], vs[3], pm[3], vm[3],
        pse[3], pme[3], llrs[3], llrm[3], pbar[4], t,
        p20m, p30m, p21m, p31m, p22m, p32m, p33m, 
        p20s, p30s, p21s, p31s, p22s, p32s, p33s,
        rerm, rers, tb[9], tbt[9];


// Luni-solar ephemeris

    planet_ephemeris (tjd, sun, CENTER, ps, vs);
    planet_ephemeris (tjd, moon, CENTER, pm, vm);

    iau_pns (tjd, tb, CENTER);	
    mt (tb, 3, 3, tbt);			
    brmul (tbt,ps,3,3,1,pse);	
    brmul (tbt,pm,3,3,1,pme);	

    
    xyz2llh(pse, llrs);
    xyz2llh(pme, llrm);

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


    gms2e  =  GMDE[sun]/GMDE[CENTER];                   
//    gmm2e  =  GMDE[moon]/GMDE[CENTER]; 
    gmm2e = 0;

    rerm = gma1 / llrm[2];
    rers = gma1 / llrs[2];

// Frequency Independent Terms

// C20
    *c20 = k2/5.0 * ( gmm2e * pow(rerm, 3) * p20m 
                    + gms2e * pow(rers, 3) * p20s );
// C21/S21
    *c21 = k2/5.0 * ( gmm2e * pow(rerm, 3) * p21m * cos(llrm[1] * DEG2RAD)
                    + gms2e * pow(rers, 3) * p21s * cos(llrs[1] * DEG2RAD) );
    *s21 = k2/5.0 * ( gmm2e * pow(rerm, 3) * p21m * sin(llrm[1] * DEG2RAD) 
                    + gms2e * pow(rers, 3) * p21s * sin(llrs[1] * DEG2RAD) );
// C22/S22
    *c22 = k2/5.0 * ( gmm2e * pow(rerm, 3) * p22m * cos(llrm[1] * DEG2RAD * 2.0)
                    + gms2e * pow(rers, 3) * p22s * cos(llrs[1] * DEG2RAD * 2.0) );
    *s22 = k2/5.0 * ( gmm2e * pow(rerm, 3) * p22m * sin(llrm[1] * DEG2RAD * 2.0)
                    + gms2e * pow(rers, 3) * p22s * sin(llrs[1] * DEG2RAD * 2.0) );
  
    
    return 0;

}




/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/*
* accel_nonsp - non-spherical force 
* @param1: description of param1
* @param2: description of param2
* version: 20 Aug 2010
*/
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
double accel_nonsp (double *tjd, double *xic, short int part, 
                    double *acc, double *dadr, double *dadp)
{
    int n, m, i;
    static int flag = 0;
//    static double *cn0, *cnm, *snm, gma[2];

    double xfc[3] = {0}, tb[9] = {0}, tbt[9] ={0}, rmat[9] = {0}, 
        gmat[9] = {0}, gmatt[9] = {0}, r, rxy, sinf, cosf, sinl, cosl, 
        lamta, *pn, *pnm, *pnp, *pnmp, *pnpp, *pnmpp,
        frj[3] = {0}, frcs[3] = {0}, fr[3] = {0}, *cosml, *sinml, *aprn,
        peprp[27], pepr[27], prtpx[9], prtpy[9], prtpz[9], pgtpx[9], 
        pgtpy[9], pgtpz[9],	pgx[3], pgy[3], pgz[3], part1[9], prjpx[3], 
        prjpy[3], prjpz[3], prcspx[3], prcspy[3], prcspz[3], prpr[9], 
        gtpr[9], part2[9], prtpxx[9], prtpyy[9], prtpzz[9];

    double dfd2r[3], dfd2[3], dfdkr[3], dfdk[3], k2, c20, c21, s21, c22, s22,
        unit, nup, ndown;

    if (part == 1)
    {
        for (n = 0; n <= 3 * DYNPAR - 1; n++)
        {
            dadp[n] = 0;
        }
        for (n = 0; n <= 8; n++)
        {
            dadr[n] = 0;
        }
    }

    if (GRAVDEGREE < 2)
    {
        for (n = 0; n <= 2; n++)
        {
            acc[n] = 0;
        }
        return 1;
    }

    if (flag != 9)			//
    {
        cn0  = (double *) calloc (GRAVDEGREE, sizeof(double));
        cnm  = (double *) calloc (GRAVDEGREE * GRAVDEGREE, sizeof(double));
        snm  = (double *) calloc (GRAVDEGREE * GRAVDEGREE, sizeof(double));
        opengravfile (cn0, cnm, snm, gma);
        flag = 9;
    }


    k2 = CONS;
    stidecs_old(tjd, gma[1], k2, &c20, &c21, &s21, &c22, &s22);

//    cn0[1] = (-8.745054708184200e-04 + c20 + DCON * 1.0e-16 * tjd[1]) * sqrt(5);
    n = 2;
    cn0[n-1] = (j2 + c20 + DCON * 1.0e-8) * sqrt(2*n+1);
    n = 3;
    cn0[n-1] = (j3 + BIAS * 1.0e-8) * sqrt(2*n+1);
    n = 4;
    cn0[n-1] = (j4 + DBIA * 1.0e-8) * sqrt(2*n+1);


    n = 2; m = 1;
    {
            nup = 1.0;
            ndown = 1.0;
            for (i = 1; i <= n - m; i++)
                nup = nup * i;
            for (i = 1; i <= n + m; i++)
                ndown = ndown * i;
            unit = sqrt (2 * (2*n+1.0)*nup/ndown);
            cnm[(n-1)*GRAVDEGREE + (m-1)] = (jc21 + c21)* unit;
            snm[(n-1)*GRAVDEGREE + (m-1)] = (js21 + s21)* unit;
    }

    n = 2; m = 2;
    {
            nup = 1.0;
            ndown = 1.0;
            for (i = 1; i <= n - m; i++)
                nup = nup * i;
            for (i = 1; i <= n + m; i++)
                ndown = ndown * i;
            unit = sqrt (2 * (2*n+1.0)*nup/ndown);
            cnm[(n-1)*GRAVDEGREE + (m-1)] = (jc22 + c22)* unit;
            snm[(n-1)*GRAVDEGREE + (m-1)] = (js22  + s22)* unit;
    }

    //    cn0[2] = -1.188691064601560e-05 * sqrt(7) * (1 + DCON);

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/*****************************rotation matrix ********************************/
    iau_pns (tjd, tb, CENTER);	//tb
    mt (tb, 3, 3, tbt);			//tbt
    brmul (tbt,xic,3,3,1,xfc);	//rb
	
    r = sqrt (xfc[0] * xfc[0] + xfc[1] * xfc[1] + xfc[2] * xfc[2]);	
    //define up-east-north system 
    rxy  = sqrt (xfc[0] * xfc[0] + xfc[1] * xfc[1]);
    sinf = xfc[2] / r;
    cosf = rxy / r;
    sinl = xfc[1] / rxy;
    cosl = xfc[0] / rxy;

    rmat[0] = cosf * cosl;		//from fixed to up-east-north system: rmat
    rmat[1] = cosf * sinl;
    rmat[2] = sinf;
    rmat[3] = -sinl;
    rmat[4] = cosl;
    rmat[5] = 0;
    rmat[6] = -sinf * cosl;
    rmat[7] = -sinf * sinl;
    rmat[8] = cosf;

    brmul (rmat,tbt,3,3,3,gmat);	//inertial to fixed matrix gmat = rmat*tbt
    mt (gmat, 3, 3, gmatt);		//fixed to inertial matrix gmatt

    lamta = chosephase (sinl, cosl); //rad

    cosml = (double *) calloc ( GRAVDEGREE, sizeof(double)); //cos(m*lamta)
    sinml = (double *) calloc ( GRAVDEGREE, sizeof(double)); //sin(m*lamta)
    aprn = (double *) calloc ( GRAVDEGREE, sizeof(double));  //sin(m*lamta)

    for (m = 1; m <= GRAVDEGREE; m++)
    {
        cosml[m-1] = cos(m*lamta);
        sinml[m-1] = sin(m*lamta);
    }
	
    for (n = 1; n <= GRAVDEGREE; n++)
    {
        aprn[n-1] = pow (gma[1] / r, n);
    }

/******************************************************************/

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/****************Legendre Polynomial**********************************/
    pn    = (double *) calloc ( GRAVDEGREE, sizeof(double));				
    //Pn
    pnp   = (double *) calloc ( GRAVDEGREE, sizeof(double));				
    //Pn'
    pnm   = (double *) calloc ( GRAVDEGREE * GRAVDEGREE, sizeof(double));	
    //secf*Pmn
    pnmp  = (double *) calloc ( GRAVDEGREE * GRAVDEGREE, sizeof(double));	
    //cosf*Pmn'
    pnpp  = (double *) calloc ( GRAVDEGREE, sizeof(double));				
    //Pn''
    pnmpp = (double *) calloc ( GRAVDEGREE * GRAVDEGREE, sizeof(double));	
    //cos2fPmn''

    pn[0]   = sinf; 
    pnp[0]  = 1; 
    pnpp[0] = 0;
    pn[1]   = 3.0/2.0*sinf*sinf - 1.0/2.0;
    pnp[1]  = sinf + 2 * sinf;
    pnpp[1] = 3;
    
    for (n = 3; n <= GRAVDEGREE; n++)
    {
        pn[n-1] = (2 * n - 1.0) / n * sinf * pn[n-2] 
            - (n - 1.0) / n * pn[n-3]; //tmd!!!
        pnp[n-1] = sinf * pnp[n-2] + n * pn[n-2];
        pnpp[n-1] = sinf * pnpp[n-2] + (n+1) * pnp[n-2];
    }

    pnm[0] = 1; //secfP11 = 1
    for (n = 2; n <= GRAVDEGREE; n++)
    {
        pnm[(n-1) * GRAVDEGREE + n - 1] 
            = (2 * n - 1.0) * cosf * pnm[(n - 2) * GRAVDEGREE + (n - 2)];
    }
	
    pnm[GRAVDEGREE] = (2 * 2.0 - 1.0) / (2 - 1.0) * sinf * pnm[0];	
    //secfP21 = pnm[GRAVDEGREE]
    for (n = 3; n <= GRAVDEGREE; n++)
    {
        for (m = 1; m < n; m++)
        {
            pnm[(n-1) * GRAVDEGREE + (m-1)] = (2 * n - 1.0) / (n-m) 
                * sinf * pnm[(n-2) * GRAVDEGREE + (m-1)]
                - (n + m - 1.0) / (n - m) 
                * pnm[(n - 3) * GRAVDEGREE + (m - 1)];
//			printf ("%d\t%d\t%f\n", n, m, pnm[(n-1)*n2 + (m-1)]);
        }
    }

    pnmp[0] = -sinf * pnm[0];		//cosfP11'
    for (n = 2; n <= GRAVDEGREE; n++)
    {
        for (m = 1; m <= n; m++)
        {
            pnmp[(n - 1) * GRAVDEGREE + (m - 1)] = 
                - n * sinf * pnm[(n - 1) * GRAVDEGREE + (m - 1)] 
                + (n + m) * pnm[(n - 2) * GRAVDEGREE + (m - 1)];
        }
    }
    
    pnmpp[0] = sinf * pnmp[0] / cosf - pnm[0] * cosf;	
    //cos2fP11''
    pnmpp[GRAVDEGREE] = sinf * pnmp[GRAVDEGREE] / cosf 
        - pnm[GRAVDEGREE] * cosf - 3 * sinf * pnm[GRAVDEGREE + 1];	
    //cos2fP21'' = pnmpp[GRAVDEGREE]
    pnmpp[GRAVDEGREE + 1] = - 2 * pnm[GRAVDEGREE + 1] * cosf;	
    //cos2fP22'' = pnmpp[GRAVDEGREE+1]
    for (n = 3; n <= GRAVDEGREE; n++)
    {
        for (m = 1; m <= n; m++)
        {
            if (m == 1)
            {
                pnmpp[(n - 1) * GRAVDEGREE + (m - 1)] = 
                    + sinf * pnmp[(n - 1) * GRAVDEGREE + (m - 1)] / cosf 
                    - pnm[(n - 1) * GRAVDEGREE + (m - 1)] * cosf
                    - 3 * sinf * pnm[(n - 1) * GRAVDEGREE + (m - 1) + 1] 
                    + pnm[(n - 1) * GRAVDEGREE + (m - 1) + 2] * cosf;
            }
            else
            {
                pnmpp[(n-1)*GRAVDEGREE + (m-1)] = - (n - 2) * sinf 
                    * pnmp[(n - 1) * GRAVDEGREE + (m-1)] / cosf
                    + (n + m) * pnmp[(n - 2) * GRAVDEGREE + (m-1)] / cosf 
                    - n * pnm[(n - 1) * GRAVDEGREE + (m - 1)] * cosf;
            }
        }
    }
/*****************Legendre Polynomial**********************************/
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/

/**********************************************************/
    for (n = 1; n <= GRAVDEGREE; n++)
    {
        frj[0] = frj[0] + (-cn0[n - 1]) * aprn[n - 1] * (n + 1) * pn[n - 1];
        frj[1] = frj[1] + 0;
        frj[2] = frj[2] + (-cn0[n - 1]) * aprn[n - 1] * (-cosf) * pnp[n - 1];
    }
    for (n = 1; n <= GRAVDEGREE; n++)
    {
        for (m = 1; m <= n; m++)
        {
            if ( n == GRAVDEGREE && m > GRAVORDER)
            {
//				printf ("%d\t%d\n",n,m);
                break;
            }
            frcs[0] = frcs[0] + aprn[n - 1] 
                * ( - (n + 1)) * pnm[(n - 1) * GRAVDEGREE + (m - 1)] * cosf
                * (cnm[(n - 1) * GRAVDEGREE + (m - 1)] * cosml[m - 1] 
                + snm[(n - 1) * GRAVDEGREE + (m - 1)] * sinml[m - 1]);

            frcs[1] = frcs[1] + aprn[n - 1] 
                * m * pnm[(n - 1) * GRAVDEGREE + (m - 1)] 
                * (-cnm[(n - 1) * GRAVDEGREE + (m - 1)] * sinml[m - 1] 
                + snm[(n - 1) * GRAVDEGREE + (m - 1)] * cosml[m - 1]);
            frcs[2] = frcs[2] + aprn[n-1] 
                * pnmp[(n - 1) * GRAVDEGREE + (m - 1)] 
                * (cnm[(n - 1) * GRAVDEGREE + (m - 1)] * cosml[m - 1] 
                + snm[(n - 1) * GRAVDEGREE + (m - 1)] * sinml[m - 1]);
        }
    }

    for (n = 0; n < 3; n++)
    {
        fr[n] = (frj[n] + frcs[n]) * gma[0] / r / r;
    }

    brmul(gmatt,fr,3,3,1,acc);  //from fixed acc to inertial acc

    if (part == 0)
    {
        free (pn);
        free (pnp);
        free (pnm);
        free (pnmp);
        free (pnpp);
        free (pnmpp);
        free (cosml);
        free (sinml);
        free (aprn);
        return 1;
    }


/*************************************************************/

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/************************partial part1***********************************/

    peprp[0]  = 0;	
    peprp[1]  = - sinl / r;					
    peprp[2]  = - sinf * cosl / r;	//(573)
    peprp[3]  = 0;	
    peprp[4]  = cosl / r;					
    peprp[5]  = - sinf * sinl / r;	//(574)		//11.10change
    peprp[6]  = 0;	
    peprp[7]  = 0;							
    peprp[8]  = cosf / r;			//(575)
    peprp[9]  = 0;	
    peprp[10] = - cosl / r / cosf;			
    peprp[11] = 0;					//(576)
    peprp[12] = 0;	
    peprp[13] = - sinl / r / cosf;			
    peprp[14] = 0;					//(577)
    peprp[15] = 0;	
    peprp[16] = 0;							
    peprp[17] = 0;					//(578)
    peprp[18] = 0;	
    peprp[19] = sinf * sinl / r / cosf;		
    peprp[20] = - cosf * cosl / r;	//(579)
    peprp[21] = 0;	
    peprp[22] = - sinf * cosl / r / cosf;	
    peprp[23] = - cosf * sinl / r;	//(580)
    peprp[24] = 0;	
    peprp[25] = 0;							
    peprp[26] = - sinf / r;			//(581)

    brmul(peprp, gmat, 9, 3, 3, pepr);		//571

    for (n = 0; n < 9; n++)					//570
    {
        prtpx[n] = pepr[n*3];
        prtpy[n] = pepr[n*3+1];
        prtpz[n] = pepr[n*3+2];
    }

    mt (prtpx, 3, 3, prtpxx);			//OH MY GOD!!!11.11
    mt (prtpy, 3, 3, prtpyy);			//OH MY GOD!!!11.11
    mt (prtpz, 3, 3, prtpzz);			//OH MY GOD!!!11.11


    brmul(tb, prtpxx, 3, 3, 3, pgtpx);		//568
    brmul(tb, prtpyy, 3, 3, 3, pgtpy);  
    brmul(tb, prtpzz, 3, 3, 3, pgtpz);  

    brmul(pgtpx, fr, 3, 3, 1, pgx);			//558 first term	//11.10 change
    brmul(pgtpy, fr, 3, 3, 1, pgy);  
    brmul(pgtpz, fr, 3, 3, 1, pgz);  

    for (n = 0; n < 3; n++)	
    {
        part1[n*3] = pgx[n]; 
        part1[n*3+1] = pgy[n]; 
        part1[n*3+2] = pgz[n];
    }

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/************************partial part2*************************************/

    for (n = 0; n < 3; n++)
    {
        prjpx[n]  = 0;	
        prjpy[n]  = 0;	
        prjpz[n]  = 0;
        prcspx[n] = 0;	
        prcspy[n] = 0;	
        prcspz[n] = 0;
    }

    for (n = 1; n <= GRAVDEGREE; n++)
    {
        prjpx[0] = prjpx[0] - (-cn0[n - 1]) * aprn[n - 1] 
            * (n + 1) * pn[n - 1] * (n + 2);        //561
        prjpx[2] = prjpx[2] - (-cn0[n - 1]) * aprn[n - 1] 
            * (-cosf) * pnp[n - 1] * (n + 2);       //561
        prjpz[0] = prjpz[0] + (-cn0[n - 1]) * aprn[n - 1] 
            * (n + 1) * cosf * pnp[n - 1];      //563
		prjpz[2] = prjpz[2]	+ (-cn0[n - 1]) * aprn[n - 1] 
            * ( sinf * pnp[n - 1] - cosf * cosf * pnpp[n - 1] );//563
    }
	
    for (n = 1; n <= GRAVDEGREE; n++)
    {
        for (m = 1; m <= n; m++)
        {
            if ( n == GRAVDEGREE && m > GRAVORDER)
            {
//				printf ("%d\t%d\n",n,m);
                break;
            }
//from 564 to 566
            prcspx[0] = prcspx[0] - aprn[n - 1] 
                * ( - (n + 1)) * pnm[(n - 1) * GRAVDEGREE + (m - 1)] * cosf
                * (cnm[(n - 1) * GRAVDEGREE + (m - 1)] * cosml[m - 1] 
                + snm[(n - 1) * GRAVDEGREE + (m - 1)] * sinml[m - 1]) 
                * (n + 2);			
            prcspx[1] = prcspx[1] - aprn[n - 1] 
                * m * pnm[(n - 1) * GRAVDEGREE + (m - 1)] 
                * ( - cnm[(n - 1) * GRAVDEGREE + (m - 1)] * sinml[m - 1] 
                + snm[(n - 1) * GRAVDEGREE + (m - 1)] * cosml[m - 1]) 
                * (n + 2);
            prcspx[2] = prcspx[2] - aprn[n - 1] 
                * pnmp[(n - 1) * GRAVDEGREE + (m - 1)] 
                * (cnm[(n - 1) * GRAVDEGREE + (m - 1)] * cosml[m - 1] 
                + snm[(n - 1) * GRAVDEGREE + (m - 1)] * sinml[m - 1])
                * (n + 2);

            prcspy[0] = prcspy[0] + m * aprn[n - 1] 
                * (n + 1) * pnm[(n - 1) * GRAVDEGREE + (m - 1)]
                * (cnm[(n - 1) * GRAVDEGREE + (m - 1)] * sinml[m - 1] 
                - snm[(n - 1) * GRAVDEGREE + (m - 1)] * cosml[m - 1]);
            prcspy[1] = prcspy[1] + m * aprn[n - 1] 
                * m * pnm[(n - 1) * GRAVDEGREE + (m - 1)] / cosf
                * ( - cnm[(n - 1) * GRAVDEGREE + (m - 1)] * cosml[m - 1] 
                - snm[(n - 1) * GRAVDEGREE + (m - 1)] * sinml[m - 1]);				
            prcspy[2] = prcspy[2] + m * aprn[n - 1] 
                * pnmp[(n - 1) * GRAVDEGREE + (m - 1)] / cosf
                * ( - cnm[(n - 1) * GRAVDEGREE + (m - 1)] * sinml[m - 1] 
                + snm[(n - 1) * GRAVDEGREE + (m - 1)] * cosml[m - 1]);

            prcspz[0] = prcspz[0] + aprn[n - 1] 
                * ( - (n + 1)) * pnmp[(n - 1) * GRAVDEGREE + (m - 1)]
                * (cnm[(n - 1) * GRAVDEGREE + (m - 1)] * cosml[m - 1] 
                + snm[(n - 1) * GRAVDEGREE + (m - 1)] * sinml[m - 1]);
            prcspz[1] = prcspz[1] + aprn[n - 1] 
                * m * (sinf * pnm[(n - 1) * GRAVDEGREE + (m - 1)] 
                + pnmp[(n - 1) * GRAVDEGREE + (m - 1)]) / cosf
                * ( - cnm[(n - 1) * GRAVDEGREE + (m - 1)] * sinml[m - 1] 
                + snm[(n - 1) * GRAVDEGREE + (m - 1)] * cosml[m - 1]);			
            prcspz[2] = prcspz[2] + aprn[n - 1] 
                * ( pnmpp[(n - 1) * GRAVDEGREE + (m - 1)] 
                - pnmp[(n - 1) * GRAVDEGREE + (m - 1)] * sinf / cosf)
                * ( cnm[(n - 1) * GRAVDEGREE + (m - 1)] * cosml[m - 1] 
                + snm[(n - 1) * GRAVDEGREE + (m - 1)] * sinml[m - 1]);
        }
    }

    for (n = 0; n < 3; n++)	
    {
        prpr[n*3]   = (prjpx[n] + prcspx[n]) * gma[0] / r / r / r;
        prpr[n*3+1] = (prjpy[n] + prcspy[n]) * gma[0] / r / r / r;
        prpr[n*3+2] = (prjpz[n] + prcspz[n]) * gma[0] / r / r / r;
    }

    brmul(prpr, gmat, 3, 3, 3, gtpr);			
    brmul(gmatt, gtpr, 3, 3, 3, part2);	

    for (n = 0; n <= 8; n++)
    {
        dadr[n] = part1[n] + part2[n]; 
    }

/*****************************************************************************/




    n = 2;
    dfd2r[0] =  (- 1.0e-8 * sqrt(2*n+1)) * aprn[n - 1] * (n + 1) * pn[n - 1];
    dfd2r[1] =  0;
    dfd2r[2] =  (- 1.0e-8 * sqrt(2*n+1)) * aprn[n - 1] * (-cosf) * pnp[n - 1];

    for (n = 0; n < 3; n++)
    {
        dfd2r[n] = dfd2r[n] * gma[0] / r / r;
    }
  
    brmul(gmatt,dfd2r,3,3,1,dfd2);  


    for (n = 0; n <= 2; n++)
    {
//        dadp[n * DYNPAR] = dfd2[n];      
//        dadp[n * DYNPAR + 1] = dfd3[n];   
        dadp[n * DYNPAR + 1] = dfd2[n];   
    }

    n = 3;
    dfd2r[0] =  (- 1.0e-8 * sqrt(2*n+1)) * aprn[n - 1] * (n + 1) * pn[n - 1];
    dfd2r[1] =  0;
    dfd2r[2] =  (- 1.0e-8 * sqrt(2*n+1)) * aprn[n - 1] * (-cosf) * pnp[n - 1];

    for (n = 0; n < 3; n++)
    {
        dfd2r[n] = dfd2r[n] * gma[0] / r / r;
    }
  
    brmul(gmatt,dfd2r,3,3,1,dfd2);  


    for (n = 0; n <= 2; n++)
    {
//        dadp[n * DYNPAR] = dfd2[n];      
//        dadp[n * DYNPAR + 1] = dfd3[n];   
        dadp[n * DYNPAR + 2] = dfd2[n];   
    }

    n = 4;
    dfd2r[0] =  (- 1.0e-8 * sqrt(2*n+1)) * aprn[n - 1] * (n + 1) * pn[n - 1];
    dfd2r[1] =  0;
    dfd2r[2] =  (- 1.0e-8 * sqrt(2*n+1)) * aprn[n - 1] * (-cosf) * pnp[n - 1];

    for (n = 0; n < 3; n++)
    {
        dfd2r[n] = dfd2r[n] * gma[0] / r / r;
    }
  
    brmul(gmatt,dfd2r,3,3,1,dfd2);  


    for (n = 0; n <= 2; n++)
    {
//        dadp[n * DYNPAR] = dfd2[n];      
//        dadp[n * DYNPAR + 1] = dfd3[n];   
        dadp[n * DYNPAR + 3] = dfd2[n];   
    }


    k2 = 1;
    stidecs_old(tjd, gma[1], k2, &c20, &c21, &s21, &c22, &s22);

    n = 2;
    cn0[n-1] = (c20) * sqrt(2*n+1);

    dfdkr[0] =  (-cn0[n - 1]) * aprn[n - 1] * (n + 1) * pn[n - 1];
    dfdkr[1] =  0;
    dfdkr[2] =  (-cn0[n - 1]) * aprn[n - 1] * (-cosf) * pnp[n - 1];


    n = 2; m = 1;
    {
            nup = 1.0;
            ndown = 1.0;
            for (i = 1; i <= n - m; i++)
                nup = nup * i;
            for (i = 1; i <= n + m; i++)
                ndown = ndown * i;
            unit = sqrt (2 * (2*n+1.0)*nup/ndown);
            cnm[(n-1)*GRAVDEGREE + (m-1)] = (  c21)* unit;
            snm[(n-1)*GRAVDEGREE + (m-1)] = (  s21)* unit;

            dfdkr[0] = dfdkr[0] + aprn[n - 1] 
                * ( - (n + 1)) * pnm[(n - 1) * GRAVDEGREE + (m - 1)] * cosf
                * (cnm[(n - 1) * GRAVDEGREE + (m - 1)] * cosml[m - 1] 
                + snm[(n - 1) * GRAVDEGREE + (m - 1)] * sinml[m - 1]);

            dfdkr[1] = dfdkr[1] + aprn[n - 1] 
                * m * pnm[(n - 1) * GRAVDEGREE + (m - 1)] 
                * (-cnm[(n - 1) * GRAVDEGREE + (m - 1)] * sinml[m - 1] 
                + snm[(n - 1) * GRAVDEGREE + (m - 1)] * cosml[m - 1]);
            dfdkr[2] = dfdkr[2] + aprn[n-1] 
                * pnmp[(n - 1) * GRAVDEGREE + (m - 1)] 
                * (cnm[(n - 1) * GRAVDEGREE + (m - 1)] * cosml[m - 1] 
                + snm[(n - 1) * GRAVDEGREE + (m - 1)] * sinml[m - 1]);

    }

    n = 2; m = 2;
    {
            nup = 1.0;
            ndown = 1.0;
            for (i = 1; i <= n - m; i++)
                nup = nup * i;
            for (i = 1; i <= n + m; i++)
                ndown = ndown * i;
            unit = sqrt (2 * (2*n+1.0)*nup/ndown);
            cnm[(n-1)*GRAVDEGREE + (m-1)] = (   c22)* unit;
            snm[(n-1)*GRAVDEGREE + (m-1)] = (   s22)* unit;

            dfdkr[0] = dfdkr[0] + aprn[n - 1] 
                * ( - (n + 1)) * pnm[(n - 1) * GRAVDEGREE + (m - 1)] * cosf
                * (cnm[(n - 1) * GRAVDEGREE + (m - 1)] * cosml[m - 1] 
                + snm[(n - 1) * GRAVDEGREE + (m - 1)] * sinml[m - 1]);

            dfdkr[1] = dfdkr[1] + aprn[n - 1] 
                * m * pnm[(n - 1) * GRAVDEGREE + (m - 1)] 
                * (-cnm[(n - 1) * GRAVDEGREE + (m - 1)] * sinml[m - 1] 
                + snm[(n - 1) * GRAVDEGREE + (m - 1)] * cosml[m - 1]);
            dfdkr[2] = dfdkr[2] + aprn[n-1] 
                * pnmp[(n - 1) * GRAVDEGREE + (m - 1)] 
                * (cnm[(n - 1) * GRAVDEGREE + (m - 1)] * cosml[m - 1] 
                + snm[(n - 1) * GRAVDEGREE + (m - 1)] * sinml[m - 1]);

    }


    for (n = 0; n < 3; n++)
    {
        dfdkr[n] = dfdkr[n] * gma[0] / r / r;
    }
  
    brmul(gmatt,dfdkr,3,3,1,dfdk);  


    for (n = 0; n <= 2; n++)
    {
        dadp[n * DYNPAR] = dfdk[n];      
//        dadp[n * DYNPAR + 1] = dfd3[n];   
//        dadp[n * DYNPAR + 1] = dfd2[n];   
    }



/*****************************************************************************/

    free (pn);
    free (pnp);
    free (pnm);
    free (pnmp);
    free (pnpp);
    free (pnmpp);
    free (cosml);
    free (sinml);
    free (aprn);
    return 0;
}




/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/*
* opengravfile - open gravity field
* @param1: description of param1
* @param2: description of param2
*
* version: 20 Aug 2010
*/
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
double opengravfile (double *cn0, double *cnm, double *snm, double *gma)
{
    FILE *fp_gra;
    double value, c,s, nup, ndown, unit;
    int n,m, i;
    char string[100], name[20];

    if ((fp_gra = fopen (FILE_GRV,"r")) == NULL)
    {
        printf ("Cannot open gravity file?\n");
        exit (0);
    }

    while (feof (fp_gra) == 0)
    {
        fgets (string, 100, fp_gra);
        sscanf (string, "%s%lf", name, &value);	
        if (strcmp (name,"Gm") ==0)	
        {
            gma[0] = value / AU / AU / AU * 86400.0 * 86400.0;
        }
        if (strcmp (name,"RefDistance") ==0)
        {
            gma[1] = value / AU;
        }
        if (strcmp (name,"BEGIN") ==0)	
            break;
    }
   
    while (feof(fp_gra) == 0)
    {
        n = 999;
        m = 999;
        fgets (string, 100, fp_gra);
        sscanf (string, "%d%d%lf%lf", &n, &m, &c, &s);	
        if (n > GRAVDEGREE)
            continue;
        else if (m == 0)
        {
            unit = sqrt (2 * n + 1.0);
            cn0[n-1] = c * unit;
            if (n == 2) j2 = c;
            if (n == 3) j3 = c;
            if (n == 4) j4 = c;

        }
        else 
        {
            if (n == 2 && m == 1) {jc21 = c; js21 = s; }
            if (n == 2 && m == 2) {jc22 = c; js22 = s; }

            nup = 1.0;
            ndown = 1.0;
            for (i = 1; i <= n - m; i++)
                nup = nup * i;
            for (i = 1; i <= n + m; i++)
                ndown = ndown * i;
            unit = sqrt (2 * (2*n+1.0)*nup/ndown);
            cnm[(n-1)*GRAVDEGREE + (m-1)] = c * unit;
            snm[(n-1)*GRAVDEGREE + (m-1)] = s * unit;
        }
    }

    fclose(fp_gra);
    return 0;
}



