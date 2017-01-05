
/*
------------------------------------------------------------------------

    Purpose: transformation between time and space coordinates
    Notes: 
    Programmer: Kun Shang @ 4.29.2014
    Functions:


        double gps_GRACE2utc (double jd0_day, double gps_GRACE);

        double gps_GRACE2tt (double jd0_day, double gps_GRACE);

        double gps2utc (double jd0_day, double gps_sec);

        double gps2tt (double gps_sec);

        int eop_open (char *eop_name, int mjd_beg, int mjd_end);

        int eop_read (double mjdt, double *xp, double *yp, 
            double *ut1_utc, double *dx, double *dy);
    
        int eop_close ();

        double lgr_order (double *y, int dim_y, int dim_x, 
            double t, double *z, int order);

        double getinfo(double *jdtt, int center, InfStruct *info);

        void cip2icrf (double *jdtt, double dx, double dy, 
                double *vt, double *vc);

        void tod2icrf (double *jdtt, double *vt, double *vc);

        void icrf2tod (double *jdtt, double *vc, double *vt);

        int getlps (double jdutc);

        double iau_pns (double *jd, double *te, int cent);

        double iau_s (double *jd, double *tes, int cent);
        
        void iau_pn (double *jd, double *tes, int cent);

        void in2pa(double *jd, double *te);

        short int mbf2cel (double *jd_tdb, double *te);

        void in2me (double *jd, double *te, short int derivation);

        void me2pa (double *te);


     Global variables:

        extern double *EOPEPH;
        extern int NEOP;



------------------------------------------------------------------------
*/

#ifndef _COORD_H_

#include "coord.h"
    #define _COORD_H_
#endif



    double *EOPEPH;
    int NEOP;




/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/

double gps2utc (double jd0_day, double gps_sec)
{
    double utc, tt;
    int lps;

    tt = gps_sec + 19 + 32.184;
    lps = getlps (jd0_day + tt / 86400.0);
    utc = gps_sec + 19 - lps;

    return utc;
}


double gps_GRACE2utc (double jd0_day, double gps_GRACE)
{
    double gps0, utc, tt;
    int lps;

    gps0 = (jd0_day - T0) * 86400.0;
    tt = gps_GRACE - gps0 + 19 + 32.184;
    lps = getlps (jd0_day + tt / 86400.0);
    utc = gps_GRACE - gps0 - (lps - 19);

    return utc;
}




double gps2tt (double gps_sec)
{
    double tt;

    tt = gps_sec + 19 + 32.184;

    return tt;
}


double gps_GRACE2tt (double jd0_day, double gps_GRACE)
{
    double gps0, tt;

    gps0 = (jd0_day - T0) * 86400.0;
    tt = gps_GRACE - gps0 + 19 + 32.184;

    return tt;
}




/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/*
    from mjd_beg - 1 to mjd_end + 1
    e.g. mjd_beg = 3, mjd_end = 4 
      => begein with 2, end with 5 (included)
      => NEOP = 4 (2,3,4,5)
*/
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
int eop_open (char *eop_name, int mjd_beg, int mjd_end)
{
    FILE *fp_eop;
    int i, mjdi;
    char string[160];

    if ((fp_eop = fopen (eop_name,"r")) == NULL)
    {
        printf ("Cannot open eop file?\n");
        exit (0);
    }

    NEOP = mjd_end - mjd_beg + 3;

    EOPEPH  = (double *) calloc (NEOP * 6, sizeof(double));

    while (feof(fp_eop) == 0)
    {
        fgets (string, 160, fp_eop);
        sscanf (string, "%*d%*d%*d%d", &mjdi);
        if (mjdi == mjd_beg - 2)
        {
            for (i = 0; i < NEOP; i ++)
            {
                fgets (string, 160, fp_eop);
                sscanf (string, "%*d%*d%*d%lf%lf%lf%lf%*f%lf%lf",
                    &EOPEPH[i * 6 + 0], &EOPEPH[i * 6 + 1], &EOPEPH[i * 6 + 2],
                    &EOPEPH[i * 6 + 3], &EOPEPH[i * 6 + 4], &EOPEPH[i * 6 + 5]);
            }
            break;
        }
    }

    fclose (fp_eop);
    return 0;


}



int eop_read (double mjdt, double *xp, double *yp, 
        double *ut1_utc, double *dx, double *dy)
{
    double eop[5];
    int mjds, mjde, n;

    mjds = EOPEPH[0];
    mjde = EOPEPH[6 * (NEOP - 1)];

    if (mjdt <= mjds)
    {
        printf ("error in eop_read: mjdt(%f) <= mjds(%d)\n", mjdt, mjds);
        exit (0);
    }

    if (mjdt >= mjde)
    {
        printf ("error in eop_read: mjdt(%f) >= mjds(%d)\n", mjdt, mjde);
        exit (0);
    }


    lgr_order (EOPEPH, NEOP, 6, mjdt, eop, 1);

    for (n = 0; n < 5; n++)
    {
        if (eop[n] > 2.0)
        {
            printf ("error in eop_read: eop[%d] = %f > 2.0\n", n, eop[n]);
            exit(0);
        }
    }

    *xp      = eop[0];
    *yp      = eop[1];
    *ut1_utc = eop[2];
    *dx      = eop[3];
    *dy      = eop[4];

    return 0;

}






int eop_close ()
{
    free (EOPEPH);
    return 0;
}







double getinfo(double *jdtt, int center, InfStruct *info)
{
    int n;
    double gmsth, angvel, ux[3] = {1,0,0}, uy[3] = {0,1,0}, uz[3] = {0,0,1}, 
        tx[3], ty[3], tz[3], xp, yp, ut1_utc, dx, dy, 
        t, st_dot, theta_dot, we[3];


    info->jd0 = jdtt[0];
    info->tt = jdtt[1] * 86400.0;

    info->jdtt[0] = jdtt[0];    
    info->jdtt[1] = jdtt[1];    

    info->jdt = info->jdtt[0] + info->jdtt[1];
//    info->jdt = info->jd0 + info->tt / 86400.0;


    info->mjd = info->jdt - 2400000.5;

    info->leaps = getlps (info->jd0);

    info->utc = info->tt - info->leaps - 32.184;
        
    info->gps = info->tt - 32.184 - 19;


    info->mjdutc = info->jd0 - 2400000.5 + info->utc/86400.0;
    eop_read (info->mjdutc, &xp, &yp, &ut1_utc, &dx, &dy);  


    info->xp      = xp;
    info->yp      = yp;
    info->ut1_utc = ut1_utc;
    info->dx      = dx;
    info->dy      = dy;


    info->deltat = 32.184 + info->leaps - info->ut1_utc;
    info->ut1 = info->utc + info->ut1_utc;

    sidereal_time (info->jd0, info->ut1/86400.0, info->deltat,0,1,1, &gmsth);
//    sidereal_time (info->jd0, info->tt/86400.0, info->deltat,0,1,1, &gmsth);
//    sidereal_time (info->jd0, info->tt/86400.0, 0,0,1,1, &gmsth);

    info->gmst = gmsth / 24 * 360.0 * DEG2RAD;

    if (center == 2)
    {
        cel_pole (info->jdt, 2, info->dx * 1e3, info->dy * 1e3);
    
        cel2ter (info->jd0, info->ut1 / 86400.0, info->deltat, 1, 1, 0,
            info->xp, info->yp, ux, tx);
        cel2ter (info->jd0, info->ut1 / 86400.0, info->deltat, 1, 1, 0,
            info->xp, info->yp, uy, ty);
        cel2ter (info->jd0, info->ut1 / 86400.0, info->deltat, 1, 1, 0,
            info->xp, info->yp, uz, tz);
    
        for (n = 0; n < 3; n++)
        {
            info->c_ie[n*3] = tx[n];
            info->c_ie[n*3+1] = ty[n];
            info->c_ie[n*3+2] = tz[n];
        }
        
        mt(info->c_ie, 3, 3, info->c_ei);
    
    
        ter2cel_dot (info->jd0, info->ut1 / 86400.0, info->deltat, 1, 1, 0,
            info->xp, info->yp, ux, tx);
        ter2cel_dot (info->jd0, info->ut1 / 86400.0, info->deltat, 1, 1, 0,
            info->xp, info->yp, uy, ty);
        ter2cel_dot (info->jd0, info->ut1 / 86400.0, info->deltat, 1, 1, 0,
            info->xp, info->yp, uz, tz);
    
        for (n = 0; n < 3; n++)
        {
            info->c_eidot[n*3] = tx[n];
            info->c_eidot[n*3+1] = ty[n];
            info->c_eidot[n*3+2] = tz[n];
        }
    
        mt(info->c_eidot, 3, 3, info->c_iedot);
    
    
        t = (info->jdt - T0) / 36525.0;
    
        st_dot =
                   (((( - 5 * 0.0000000368   * t
                        - 4 * 0.000029956  ) * t
                        - 4 * 0.00000044   ) * t
                        + 2 * 1.3915817    ) * t
                        +  4612.156534     );
    
    
        theta_dot = 1.00273781191135448 * TWOPI / 86400.0 
                + st_dot / 206264.806 / 36525.0 / 86400.0;
    
//        we[0] = 0; we[1] = 0; we[2] = ANGVEL;
        we[0] = 0; we[1] = 0; we[2] = theta_dot;
    
        cip2icrf (info->jdtt, info->dx, info->dy, we, info->wi);
    }
    else
    {
        angvel = iau_pns(info->jdtt, info->c_ei, center);
        mt(info->c_ei, 3, 3, info->c_ie);     
        we[0] = 0; we[1] = 0; we[2] = angvel;
        tod2icrf (info->jdtt, we, info->wi);
    }


    return 0;
}


void cip2icrf (double *jdtt, double dx, double dy, 
        double *vt, double *vc)
{
    double jd_tdb, v3[3], v4[3];

    jd_tdb = jdtt[0] + jdtt[1];
    cel_pole (jd_tdb, 2, dx * 1e3, dy * 1e3);
    nutation (jd_tdb,-1,1,vt, v3);
    precession (jd_tdb,v3,T0, v4);
    frame_tie (v4,-1, vc);

}




void tod2icrf (double *jdtt, double *vt, double *vc)
{
    double tepn[9];

    iau_pn (jdtt, tepn, 3);        //IAU inertial to J2000 inertial
    brmul (tepn,vt,3,3,1,vc);

}

void icrf2tod (double *jdtt, double *vc, double *vt)
{
    double tepn[9], te[9];

    iau_pn (jdtt, tepn, 3);        //IAU inertial to J2000 inertial
    mt(tepn, 3, 3, te);
    brmul (te,vc,3,3,1,vt);
}










/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/*
 *
 * getlps - get the leap seconds value for input JD
 * 
 * jdutc: double, Julian Day of UTC 
 * return: short int, leap seconds
 *
 * version: Mar 2013
 * 
 */
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
int getlps (double jdutc)
{
/*
 *
1972 JAN  1 =JD 2441317.5  TAI-UTC=  10.0       S + (MJD - 41317.) X 0.0      S
1972 JUL  1 =JD 2441499.5  TAI-UTC=  11.0       S + (MJD - 41317.) X 0.0      S
1973 JAN  1 =JD 2441683.5  TAI-UTC=  12.0       S + (MJD - 41317.) X 0.0      S
1974 JAN  1 =JD 2442048.5  TAI-UTC=  13.0       S + (MJD - 41317.) X 0.0      S
1975 JAN  1 =JD 2442413.5  TAI-UTC=  14.0       S + (MJD - 41317.) X 0.0      S
1976 JAN  1 =JD 2442778.5  TAI-UTC=  15.0       S + (MJD - 41317.) X 0.0      S
1977 JAN  1 =JD 2443144.5  TAI-UTC=  16.0       S + (MJD - 41317.) X 0.0      S
1978 JAN  1 =JD 2443509.5  TAI-UTC=  17.0       S + (MJD - 41317.) X 0.0      S
1979 JAN  1 =JD 2443874.5  TAI-UTC=  18.0       S + (MJD - 41317.) X 0.0      S
1980 JAN  1 =JD 2444239.5  TAI-UTC=  19.0       S + (MJD - 41317.) X 0.0      S
1981 JUL  1 =JD 2444786.5  TAI-UTC=  20.0       S + (MJD - 41317.) X 0.0      S
1982 JUL  1 =JD 2445151.5  TAI-UTC=  21.0       S + (MJD - 41317.) X 0.0      S
1983 JUL  1 =JD 2445516.5  TAI-UTC=  22.0       S + (MJD - 41317.) X 0.0      S
1985 JUL  1 =JD 2446247.5  TAI-UTC=  23.0       S + (MJD - 41317.) X 0.0      S
1988 JAN  1 =JD 2447161.5  TAI-UTC=  24.0       S + (MJD - 41317.) X 0.0      S
1990 JAN  1 =JD 2447892.5  TAI-UTC=  25.0       S + (MJD - 41317.) X 0.0      S
1991 JAN  1 =JD 2448257.5  TAI-UTC=  26.0       S + (MJD - 41317.) X 0.0      S
1992 JUL  1 =JD 2448804.5  TAI-UTC=  27.0       S + (MJD - 41317.) X 0.0      S
1993 JUL  1 =JD 2449169.5  TAI-UTC=  28.0       S + (MJD - 41317.) X 0.0      S
1994 JUL  1 =JD 2449534.5  TAI-UTC=  29.0       S + (MJD - 41317.) X 0.0      S
1996 JAN  1 =JD 2450083.5  TAI-UTC=  30.0       S + (MJD - 41317.) X 0.0      S
1997 JUL  1 =JD 2450630.5  TAI-UTC=  31.0       S + (MJD - 41317.) X 0.0      S
1999 JAN  1 =JD 2451179.5  TAI-UTC=  32.0       S + (MJD - 41317.) X 0.0      S
2006 JAN  1 =JD 2453736.5  TAI-UTC=  33.0       S + (MJD - 41317.) X 0.0      S
2009 JAN  1 =JD 2454832.5  TAI-UTC=  34.0       S + (MJD - 41317.) X 0.0      S
2012 JUL  1 =JD 2456109.5  TAI-UTC=  35.0       S + (MJD - 41317.) X 0.0      S
 *
*/

    short int lps;
    double jd = jdutc;

    if (jd >= 2456109.5)
        lps = 35;
    else if (jd >= 2454832.5)
        lps = 34;
    else if (jd >= 2453736.5)
        lps = 33;
    else if (jd >= 2451179.5)
        lps = 32;
    else if (jd >= 2450630.5)
        lps = 31;
    else if (jd >= 2450083.5)
        lps = 30;
    else if (jd >= 2449534.5)
        lps = 29;
    else if (jd >= 2449169.5)
        lps = 28;
    else if (jd > 2448804.5)
        lps = 27;
    else if (jd >= 2448257.5)
        lps = 26;
    else if (jd >= 2447892.5)
        lps = 25;
    else if (jd >= 2447161.5)
        lps = 24;
    else if (jd >= 2446247.5)
        lps = 23;
    else if (jd >= 2445516.5)
        lps = 22;
    else if (jd >= 2445151.5)
        lps = 21;
    else if (jd >= 2444786.5)
        lps = 20;
    else 
    {
        printf ("No leapsecond configured before 1981 JUL  1 =JD 2444786.5\n");
        exit (0);
    }
    return lps;
  

}

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/





/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/*
* iau_pns - planet fixed to J2000 inertial (for gravity field)
        Report of the IAU/IAGWorking Group on cartographic
        coordinates and rotational elements: 2006
* @param1: description of param1
* @param2: description of param2
* todo: 
* version: 20 Aug 2010
*/
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
double iau_pns (double *jd, double *te, int cent) 
{
    double tes[9] = {0}, tepn[9] ={0}, angvel;

    if (cent == 9)
    {
        mbf2cel (jd, te);
        angvel = 13.17635815 * DEG2RAD / 86400;
        in2pa (jd, tes);  mt (tes, 3, 3, te);  // need JPLEPH to uncomment
    }
    else
    {
        cent = cent +1;
        angvel = iau_s (jd, tes, cent);          //IAU fixed to IAU inertial
        iau_pn (jd, tepn, cent);        //IAU inertial to J2000 inertial
        brmul (tepn,tes,3,3,3,te);
    }
    return angvel;
}


/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/*
* iau_s - from IAU fixed to IAU inertial, true-of-date equator and equinox
* @param1: description of param1
* @param2: description of param2
*
* version: 20 Aug 2010
*/
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
double iau_s (double *jd, double *tes, int cent)
{
    double d, str = 0, cosst, sinst, rate = 0;
    
    d = jd[0] - 2451545.0;
    d = d + jd[1];

    switch (cent)       //sun0, mercury1, ..., pluto9 
    {
        case 0 : str = 84.176 + 14.1844000 * d; rate = 14.1844000; break;
        case 1 : str = 329.548 + 6.1385025 * d; rate =  6.1385025; break;
        case 2 : str =  160.20 - 1.4813688 * d; rate = -1.4813688; break;
        case 3 : str = 190.147 + 360.9856235 * d; rate = 360.9856235; break;
        case 4 : str = 176.63 + 350.89198226 * d; rate = 350.89198226; break;    
        case 5 : str = 284.95 + 870.5366420 * d; rate = 870.5366420; break;
        case 6 : str =  38.90 + 810.7939024 * d; rate = 810.7939024; break;
        case 7 : str = 203.81 - 501.1600928 * d; rate = - 501.1600928; break;
        case 8 : str = 253.18 + 536.3128492 * d - 
                 0.48 * sin ((357.85 + 52.316 * d / 36525.0 ) * DEG2RAD);
                rate = 536.3128492;
                break;
        case 9 : str = 237.305 - 56.3625225 * d; rate = - 56.3625225;  break;
        case I_TITAN : str = 186.5855 + 22.5769768 * d; rate = 22.5769768; break;
    }

    cosst = cos (str * DEG2RAD);
    sinst = sin (str * DEG2RAD);
    
    tes[0] = cosst;
    tes[1] = -sinst;
    tes[2] = 0;
    tes[3] = sinst;
    tes[4] = cosst;
    tes[5] = 0;
    tes[6] = 0;
    tes[7] = 0;
    tes[8] = 1;

    rate = rate * DEG2RAD / 86400.0;

    return rate;
}


/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/*
* iau_pn - from IAU inertial (for all planets) to J2000 inertial
* @param1: description of param1
* @param2: description of param2
*
* version: 20 Aug 2010
*/
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
void iau_pn (double *jd, double *tes, int cent)
{   
    double ra0 = 0, dec0 = 0, jcent, cr, sr, cd, sd;
    
    jcent = jd[0] - 2451545.0;
    jcent = (jcent + jd[1]) / 36525.0;
    
    switch (cent)       //sun0, mercury1, ..., pluto9 
    {
    case 0 : 
        ra0  = 286.13;
        dec0 = 63.87; 
        break;      
    case 1 : 
        ra0  = 281.01 - 0.033 * jcent;
        dec0 = 61.45 - 0.005 * jcent; 
        break;          
    case 2 : 
        ra0  = 272.76;
        dec0 = 67.16; 
        break;  
    case 3 : 
        ra0  = 0.00 - 0.641 * jcent;
        dec0 = 90.0 - 0.557 * jcent; 
        break;  
    case 4 : 
        ra0  = 317.68143 - 0.1061 * jcent;
        dec0 = 52.88650 - 0.0609 * jcent; 
        break;  
    case 5 : 
        ra0  = 268.05 - 0.009 * jcent;
        dec0 = 64.49 + 0.003 * jcent; 
        break;  
    case 6 : 
        ra0  = 40.589 - 0.036 * jcent;
        dec0 = 83.537 - 0.004 * jcent; 
        break;
    case 7 : 
        ra0  = 257.311;
        dec0 = -15.175; 
        break;
    case 8 : 
        ra0  = 299.36 + 0.70 * sin ((357.85 + 52.316 * jcent) * DEG2RAD);
        dec0 = 43.46 - 0.51 * cos ((357.85 + 52.316 * jcent) * DEG2RAD); 
        break;
    case 9 : 
        ra0 = 313.02;
        dec0 = 9.09; 
        break;
    case I_TITAN : 
        ra0 = 39.4827;
        dec0 = 83.4279;        
        break;
    }

    cr = cos (ra0 * DEG2RAD);
    sr = sin (ra0 * DEG2RAD);
    cd = cos (dec0 * DEG2RAD);
    sd = sin (dec0 * DEG2RAD);

    tes[0] = -sr;
    tes[1] = -cr * sd;
    tes[2] = cr * cd;
    tes[3] = cr;
    tes[4] = -sr * sd;
    tes[5] = sr * cd;
    tes[6] = 0;
    tes[7] = cd;
    tes[8] = sd;
    return;
}






/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/*
* in2pa -   from inertial to moon fixed (PA)
* @param1: description of param1
* @param2: description of param2
*
* version: 20 Aug 2010
*/
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
void in2pa(double *jd, double *te)
{
    double lib[6] = {0}, tb1[9], tb2[9], tb3[9], tb32[9];
    long int target, center;

    target = 15;
    center = 0;

    dpleph_(jd, &target, &center, lib);
//    DPLEPH(jd, &target, &center, lib);

    rotmatz (lib[0], tb1, 0);
    rotmatx (lib[1], tb2, 0);
    rotmatz (lib[2], tb3, 0);

    brmul(tb3, tb2, 3, 3, 3, tb32);
    brmul(tb32, tb1, 3, 3, 3, te);
}



/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/*
* mbf2cel - simulate doppler observable
* @param1: description of param1
* @param2: description of param2
*
* version: 20 Aug 2010
*/
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
short int mbf2cel (double *jd_tdb, double *te)
/*
------------------------------------------------------------------------

   PURPOSE:
      This function rotates a vector from the moon body-fixed system 
      to the celestial system.  

   REFERENCES:
      P. Kenneth Seidelmann et. al. (2007). Report of the IAU/IAGWorking 
      Group on cartographic coordinates and rotational elements: 2006

   INPUT
   ARGUMENTS:
      jd_tdb[2] (double)
         TDB Julian date.
         High-order part (jd_tdb[0]) & Low-order part (jd_tdb[0]).
      method (short int)
         Selection for method
            = 0 ... IAU report formulae
            = 1 ... NASA/JPL DE/LE ephemeris 
      ref_sys (short int)
         Reference system in which moon body-fixed system is given
            = 0 ... Mean Earth/polar axis (ME) system
            = 1 ... Principal Axis (PA) system
      derivation (short int)
         Seclection derivation of parameters
            = 0 ... No derivation, vecc is normal
            = 1 ... fisrt parameter derivation, vecc is derivation
            = 2 ... second 
            = 3 ... third
      vecm[3] (double)
         Position vector referred to moon body-fixed system

   OUTPUT
   ARGUMENTS:
      vecc[3] (double)
         Position vector referred to ICRF axes (celestial system)

   RETURNED
   VALUE:
      =  0  ... everything is ok.
      =  1  ... invalid value of 'ref_sys'
      =  2  ... invalid value of 'method'

   GLOBALS
   USED:

   FUNCTIONS
   CALLED:

   VER./DATE/
   PROGRAMMER:
      V1.0/03-10/ (SHAO).

   NOTES:

------------------------------------------------------------------------
*/
{
    short int error = 0;
    double tb1[9], tb2[9], tbt1[9], tbt2[9];

/*
IAU report formulae
*/
    me2pa(tb1);
    mt(tb1, 3, 3, tbt1);

    in2me(jd_tdb, tb2, 0);  
    mt(tb2, 3, 3, tbt2);
    brmul(tbt2,tbt1,3,3,3,te);  //ME2ICRF

    return error;
}


/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/*
* in2me - from inertial to moon fixed (ME)
* @param1: description of param1
* @param2: description of param2
*
* version: 20 Aug 2010
*/
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
void in2me (double *jd, double *te, short int derivation)
{
    double ra, dec, w, lib[3], d, T, tb1[9], tb2[9], tb3[9], tb32[9], 
        E1, E2, E3, E4, E5, E6, E7, E8, E9, E10, E11, E12, E13;
    
    d = jd[0] - 2451545.0 + jd[1];
    T = d / 36525.0;

    E1  = 125.045 - 0.0529921 * d;  
    E2  = 250.089 - 0.1059842 * d;  
    E3  = 260.008 + 13.0120009 * d;
    E4  = 176.625 + 13.3407154 * d; 
    E5  = 357.529 + 0.9856003 * d;  
    E6  = 311.589 + 26.4057084 * d;
    E7  = 134.963 + 13.0649930 * d; 
    E8  = 276.617 + 0.3287146 * d;  
    E9  = 34.226 + 1.7484877 * d;
    E10 = 15.134 - 0.1589763 * d;   
    E11 = 119.743 + 0.0036096 * d;  
    E12 = 239.961 + 0.1643573 * d;
    E13 = 25.053 + 12.9590088 * d;

    ra  = 269.9949 + 0.0031 * T - 3.8787 * sin (E1 * DEG2RAD) 
        - 0.1204 * sin (E2 * DEG2RAD) + 0.0700 * sin (E3 * DEG2RAD) 
        - 0.0172 * sin (E4 * DEG2RAD) + 0.0072 * sin (E6 * DEG2RAD) 
        - 0.0052 * sin (E10 * DEG2RAD) + 0.0043 * sin (E13 * DEG2RAD);
    dec = 66.5392 + 0.0130 * T + 1.5419 * cos (E1 * DEG2RAD) 
        + 0.0239 * cos (E2 * DEG2RAD) - 0.0278 * cos (E3 * DEG2RAD) 
        + 0.0068 * cos (E4 * DEG2RAD) - 0.0029 * cos (E6 * DEG2RAD)
        + 0.0009 * cos (E7 * DEG2RAD) + 0.0008 * cos (E10 * DEG2RAD) 
        - 0.0009 * cos (E13 * DEG2RAD);
    w   = 38.3213 + 13.17635815 * d - 1.4e-12 * d * d 
        + 3.5610 * sin (E1 * DEG2RAD) + 0.1208 * sin (E2 * DEG2RAD) 
        - 0.0642 * sin (E3 * DEG2RAD) + 0.0158 * sin (E4 * DEG2RAD)
        + 0.0252 * sin (E5 * DEG2RAD) - 0.0066 * sin (E6 * DEG2RAD) 
        - 0.0047 * sin (E7 * DEG2RAD) - 0.0046 * sin (E8 * DEG2RAD) 
        + 0.0028 * sin (E9 * DEG2RAD) + 0.0052 * sin (E10 * DEG2RAD)
        + 0.0040 * sin (E11 * DEG2RAD) + 0.0019 * sin (E12 * DEG2RAD) 
        - 0.0044 * sin (E13 * DEG2RAD);

    lib[0] = (90.0 + ra) * DEG2RAD;
    lib[1] = (90.0 - dec) * DEG2RAD;
    lib[2] = w * DEG2RAD;
    
    rotmatz (lib[0], tb1, 0);
    rotmatx (lib[1], tb2, 0);
    rotmatz (lib[2], tb3, 0);

    if (derivation == 1)
        rotmatz (lib[0], tb1, 1);
    if (derivation == 2)
        rotmatx (lib[1], tb2, 1);
    if (derivation == 3)
        rotmatz (lib[2], tb3, 1);
    
    brmul(tb3, tb2, 3, 3, 3, tb32);
    brmul(tb32, tb1, 3, 3, 3, te);
}


/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/*
* me2pa - simulate doppler observable
* @param1: description of param1
* @param2: description of param2
*
* version: 20 Aug 2010
*/
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
void me2pa (double *te)
{
    double tb1[9], tb2[9], tb3[9], tb32[9];

    rotmatx ( 0.1462 * ASEC2RAD, tb1, 0);
    rotmaty (79.0768 * ASEC2RAD, tb2, 0);
    rotmatz (63.8986 * ASEC2RAD, tb3, 0);
    brmul(tb3, tb2, 3, 3, 3, tb32);
    brmul(tb32, tb1, 3, 3, 3, te);
}






