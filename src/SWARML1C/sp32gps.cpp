#include <fstream>
#include <iostream>
#include <sstream>
#include <string>


#include "coord.h"
#include "numrs.h"


using namespace std;


//#define MAXLINE 300

//short int ACCURACY = 1;

//int LAG_ORD = 8;

//int NEOP;
//double *EOPMT;

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
int main (int argc, char *argv[])
{
    int year, month, day, hour, min, days, mjds;
    int gps_e, gps_i, dim, i, n, gpsweek, vel;
    double jdi, sec, xe[3], ve[3], xve[6], xi[3], vi[3], utc, tt, 
        tjd[2], xi_a[3], vi_a[3], xi_b[3], vi_b[3], llh_a[3], llh_b[3], 
        viv[3], vix[3], x_a[3], x_b[3], v_a[3], v_b[3], 
        *gnv_a, *gnv_b, ex_a[9], ev_a[9], ex_b[9], ev_b[9], tmp[9], *orb_eph,
        exi_a[9], evi_a[9], exi_b[9], evi_b[9], eviv_a[9], evix_a[9], eviv_b[9], evix_b[9];
//    char line[MAXLINE];  
    string stdname, card, fsp3, fgpsi, fgpse, feop;  
    time_t s0, s1, s2, s3;
    int DT, NDATA, GPS_S;

    InfStruct info;

    if (argc < 2) {
        cout << "Usage: " << argv[0] << " <input file name>" << endl;
        return EXIT_FAILURE;
    }

    ifstream if_std, if_sp3;
    ofstream of_gpsi, of_gpse; 
    if_std.open(argv[1]); 
//    ostringstream line;
//    stringstream sline; 
    string line;
    while(!if_std.eof()) {
        getline(if_std, line);
        cout << line;
        istringstream iline(line);
        iline >> card;
        if (card == "INPUT") {iline >> fsp3; if_sp3.open(fsp3.c_str()); }
        if (card == "OUTPUT") {iline >> fgpsi >> fgpse; of_gpsi.open(fgpsi.c_str()); of_gpse.open(fgpse.c_str()); }
        if (card == "EOP") {iline >> feop;}
        if (card == "EPOCH") {iline >> year >> month >> day; }
    }
    

    JD0 = julian_date ((short)year,(short)month,(short)day,0);

//    NDATA = days * 86400 / DT;

    mjds = (int)(JD0 - 2400000.5);
    mjde = mjds + 1;
    eop_open (feop, mjds - 1, mjde + 1);





    GPS_S = (int)((JD0 - T0) * 86400 + 0.5);
    gps_e = GPS_S + days * 86400;

//    printf ("year = %d\nJD0 = %f\njdgps0 = %f\nGPS_S = %d",
//        year, JD0, jdgps0, GPS_S);


    if (vel == 0) dim = 4;
    if (vel == 1) dim = 7;

    orb_eph = (double *) calloc ( NDATA * dim, sizeof(double));


    for (i = 0; i < 22; i++)
        fgets(line,300, fpin);

    i = 0;
    while (1)
    {
        if (fgets(line,300, fpin) ==NULL) break;

        sscanf (line, "%*s%d%d%d%d%d%lf", 
            &year, &month, &day, &hour, &min, &sec);
    
        jdi = julian_date ((short)year,(short)month,(short)day,0);

        orb_eph[i*dim] = (int)((jdi - T0) * 86400 + 0.5) + hour * 3600 + min * 60 + sec;

        if (fgets(line,300, fpin) ==NULL) break;
        sscanf (line, "%*s%lf%lf%lf", 
            &orb_eph[i * dim + 1],&orb_eph[i * dim + 2],&orb_eph[i * dim + 3]);
        
        if (vel == 1)
        {
            if (fgets(line,300, fpin) ==NULL) break;
            sscanf (line, "%*s%lf%lf%lf", 
                &orb_eph[i * dim + 4],&orb_eph[i * dim + 5],&orb_eph[i * dim + 6]);
        }
//        printf ("%lf\t%lf\t%lf\t%lf\n",orb_eph[i*4], orb_eph[i * 4 + 1],orb_eph[i * 4 + 2],orb_eph[i * 4 + 3]);
        i++;
    }

    printf ("i = %d\t NDATA = %d\n", i, NDATA);

    for (gps_i = GPS_S; gps_i < gps_e; gps_i = gps_i + 5)
    {
        if (vel == 0)
        {
            lag_deri (orb_eph, i, dim, gps_i, xe, ve);
            for (n = 0; n < 3; n ++)
            {
                xe[n] = xe[n] * 1000;
                ve[n] = ve[n] * 1000;
            }
        }
        if (vel == 1)
        {
            lagrange (orb_eph, i, dim, gps_i, xve);
            for (n = 0; n < 3; n ++)
            {
                xe[n] = xve[n] * 1000;
                ve[n] = xve[n + 3] * 0.1;
            }
        }

        tt = gps_i - GPS_S + 19 + 32.184;
        tjd[0] = JD0;    tjd[1] = tt / 86400.0;
        getinfo (tjd, &info);

        brmul(info.c_ei, xe, 3,3,1, xi);

        brmul(info.c_ei, ve, 3,3,1, viv);
        brmul(info.c_eidot, xe, 3,3,1, vix);

        for (n = 0; n < 3; n ++)
            vi[n] = viv[n] + vix[n];

        fprintf (fpout1, "%10d X E %20.15lf %20.15lf %20.15lf %20.15lf %20.15lf %20.15lf\n",
            (int)gps_i, xe[0], xe[1], xe[2], ve[0], ve[1], ve[2]);
        fprintf (fpout2, "%10d X I %20.15lf %20.15lf %20.15lf %20.15lf %20.15lf %20.15lf\n",
            (int)gps_i, xi[0], xi[1], xi[2], vi[0], vi[1], vi[2]);
    }



/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/

  
    fclose(fpin);
    fclose(fpout1);
    fclose(fpout2);
    

//    printf ("\nNormal end of ECHO!\n\npress any key to finish...\n");

    exit(0);

}
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
 



short int getlps (double jd)
{

    short int lps;


    if (jd > 2456109.5)
        lps = 35;
    else if (jd > 2454832.5)
        lps = 34;
    else if (jd > 2453736.5)
        lps = 33;
    else if (jd > 2451179.5)
        lps = 32;
    else if (jd > 2450630.5)
        lps = 31;
    else if (jd > 2450083.5)
        lps = 30;
    else if (jd > 2449534.5)
        lps = 29;
    else if (jd > 2449169.5)
        lps = 28;
    else if (jd > 2448804.5)
        lps = 27;
    else if (jd > 2448257.5)
        lps = 26;
    else if (jd > 2447892.5)
        lps = 25;
    else if (jd > 2447161.5)
        lps = 24;
    else if (jd > 2446247.5)
        lps = 23;
    else if (jd > 2445516.5)
        lps = 22;
    else if (jd > 2445151.5)
        lps = 21;
    else if (jd > 2444786.5)
        lps = 20;
    else
    {
        printf ("No leapsecond configured before 1981 JUL  1 =JD 2444786.5\n");
        exit (0);
    }
    return lps;


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
    int i, mjdi;
    double x1, y1, dt1, dx1, dy1, x2, y2, dt2, dx2, dy2; 

            
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





double getinfo(double *tjd, InfStruct *info)
{
    int n;
    double gmsth, ux[3] = {1,0,0}, uy[3] = {0,1,0}, uz[3] = {0,0,1}, 
        tx[3], ty[3], tz[3], xp, yp, ut1_utc, dx, dy;


    info->jd0 = tjd[0];
    info->tt = tjd[1] * 86400.0;

    info->jdt = info->jd0 + info->tt / 86400.0;
        
    info->leaps = getlps (info->jdt);

    info->utc = info->tt - info->leaps - 32.184;

    info->mjd = info->jd0 - 2400000.5 + info->utc/86400.0;
    geteop (info->mjd, &xp, &yp, &ut1_utc, &dx, &dy);	

    info->xp      = xp;
    info->yp      = yp;
    info->ut1_utc = ut1_utc;
    info->dx      = dx;
    info->dy      = dy;


    info->deltat = 32.184 + info->leaps - info->ut1_utc;
    info->ut1 = info->utc + info->ut1_utc;

    sidereal_time (info->jd0, info->ut1/86400.0, info->deltat,0,1,1, &gmsth);

    info->gmst = gmsth / 24 * 360.0 * DEG2RAD;


    cel_pole (info->jdt, 2, info->dx * 1e3, info->dy * 1e3);

    ter2cel (info->jd0, info->ut1 / 86400.0, info->deltat, 1, 1, 0,
        info->xp, info->yp, ux, tx);
    ter2cel (info->jd0, info->ut1 / 86400.0, info->deltat, 1, 1, 0,
        info->xp, info->yp, uy, ty);
    ter2cel (info->jd0, info->ut1 / 86400.0, info->deltat, 1, 1, 0,
        info->xp, info->yp, uz, tz);

    for (n = 0; n < 3; n++)
    {
        info->c_ei[n*3] = tx[n];
        info->c_ei[n*3+1] = ty[n];
        info->c_ei[n*3+2] = tz[n];
    }
    
    mt(info->c_ei, 3, 3, info->c_ie);





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




    return 0;
}




void mt (double *a, int m, int n, double *b)
{
    int i, j;
    for (i = 0; i <= m - 1; i++)
    {
        for (j = 0; j <= n - 1; j++)
            b[j * m + i] = a[i * n + j];
    }
    return;
}


/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/*
* brmul - 
* @param1: description of param1
* @param2: description of param2
*
* version: 20 Aug 2010
*/
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
void brmul (double *a, double *b, int m,int n, int k,double *c)
{ 
    int i, j, l, u;
    for (i = 0; i <= m - 1; i++)
    {
        for (j = 0; j <= k - 1; j++)
        {
            u = i * k + j; 
            c[u] = 0.0;
            for (l = 0; l <= n - 1; l++)
                c[u] = c[u] + a[i * n + l] * b[l * k + j];
        }
    }
    return;
}





















void openeop (char file_eop[200], int mjds, int num, double *eopmat)
{
    FILE *fp_eop;
    int i, mjdi;
    char string[160];
   
    if ((fp_eop = fopen (file_eop,"r")) == NULL)
    {
        printf ("Cannot open eop file?\n");
        exit (0);
    }
        
//    for (i = 0; i < 13;i++)
//        fgets (string, 160, fp_eop);
    while (feof(fp_eop) == 0)
    {
        fgets (string, 160, fp_eop);
        sscanf (string, "%*d%*d%*d%d", &mjdi);
        if (mjdi == mjds - 1)
        {
            for (i = 0; i < num; i ++)
            {
                fgets (string, 160, fp_eop);
                sscanf (string, "%*d%*d%*d%lf%lf%lf%lf%*lf%lf%lf",
                    &eopmat[i * 6 + 0], &eopmat[i * 6 + 1], &eopmat[i * 6 + 2], 
                    &eopmat[i * 6 + 3], &eopmat[i * 6 + 4], &eopmat[i * 6 + 5]);
//                printf("mjd = %f\n", eopmat[i * 6 + 0]);
            }
            break;
        }
    }    

    fclose (fp_eop);
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
    int i, j, k, m, dim, order = LAG_ORD;
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




double lag_deri (double *y, int dim_y, int dim_x, double t, double *z, double *dz)
{
    int l,i, j, k, m, ms, me, n, dim, order = LAG_ORD;
    double lm, dlm;

    i = 0;
    while ((y[i * dim_x] < t) && (i < dim_y))
        i = i + 1;
    ms = i - order;
    if (ms < 0)
        ms = 0;
    me = i + order - 1;
    if (me > dim_y - 1)
        me = dim_y - 1;

    for (dim = 0; dim < dim_x - 1; dim++)
    {
        z[dim] = 0;
        dz[dim] = 0;
    }

    for (m = ms; m <= me; m ++)
    {
        lm = 1.0;
        for (j = ms; j <= me; j ++)
        {
            if (j != m)
            {
                lm = lm * (t - y[j * dim_x]) / (y[m * dim_x] - y[j * dim_x]);
            }
        }
        for (dim = 0; dim < dim_x - 1; dim++)
        {
            z[dim] = z[dim] + lm * y[m * dim_x + dim + 1];
        }
    }



    for (m = ms; m <= me; m ++)
    {
        dlm = 0;
        for (l = ms; l <= me; l++)
        {
            if (l != m)
            {
                lm = 1.0;
                for (j = ms; j <= me; j ++)
                {
                    if ((j != m) && (j != l))
                    {
                        lm = lm * (t - y[j * dim_x]) / (y[m * dim_x] - y[j * dim_x]);
                    }
                }
                dlm = dlm + lm / (y[m * dim_x] - y[l * dim_x]);
            }
        }
        for (dim = 0; dim < dim_x - 1; dim++)
        {
            dz[dim] = dz[dim] + dlm * y[m * dim_x + dim + 1];
        }
    }





    return 0;

}







/*
 *
 *

 *
 *

http://www.phys.ufl.edu/~coldwell/wsteve/FDERS/lagrange.for

      SUBROUTINE UELAG(NL,X,MBEG,NSKIP,ALAG,XDAT)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION ALAG(*),XDAT(*)
      DO M=1,NL
        IF(M+MBEG.EQ.NSKIP)THEN
          ALAG(M)=0
        ELSE
          ALAG(M)=1
          DO J=1,NL
            IF(J.NE.M.AND.J+MBEG.NE.NSKIP)THEN
              ALAG(M)=ALAG(M)*(X-XDAT(J+MBEG))/(XDAT(M+MBEG)-
     2        XDAT(J+MBEG))
            ENDIF
          ENDDO
        ENDIF
      ENDDO
      RETURN
      END
      SUBROUTINE UDLAG(NL,X,MBEG,NSKIP,DLAG,XDAT)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION DLAG(*),XDAT(*)
      DO M=1,NL
        DLAG(M)=0
        IF(M+MBEG.NE.NSKIP)THEN
          DO L=1,NL
            IF(L.NE.M.AND.L+MBEG.NE.NSKIP)THEN
              PROD=1
              DO J=1,NL
                IF(J.NE.M.AND.J+MBEG.NE.NSKIP.AND.J.NE.L)THEN
                  PROD=PROD*(X-XDAT(J+MBEG))/(XDAT(M+MBEG)-
     2             XDAT(J+MBEG))
                ENDIF
              ENDDO
              DLAG(M)=DLAG(M)+PROD/(XDAT(M+MBEG)-XDAT(L+MBEG))
            ENDIF
          ENDDO
        ENDIF
      ENDDO
      RETURN
      END
      SUBROUTINE UDDLAG(NL,X,MBEG,NSKIP,DDLAG,XDAT)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION DDLAG(*),XDAT(*)
      DO M=1,NL
        DDLAG(M)=0
        IF(M+MBEG.NE.NSKIP)THEN
          DO L=1,NL
            IF(L.NE.M.AND.L+MBEG.NE.NSKIP)THEN
              SUM=0
              DO N=1,NL
                IF(N.NE.L.AND.N.NE.M.AND.N+MBEG.NE.NSKIP)THEN
                  PROD=1
                  DO J=1,NL
                    IF(J.NE.M.AND.J+MBEG.NE.NSKIP.AND.J.NE.L.AND.
     2               J.NE.N)THEN
                      PROD=PROD*(X-XDAT(J+MBEG))/(XDAT(M+MBEG)-
     2                XDAT(J+MBEG))
                    ENDIF
                  ENDDO
                  SUM=SUM+PROD/(XDAT(M+MBEG)-XDAT(N+MBEG))
                ENDIF
              ENDDO
              DDLAG(M)=DDLAG(M)+SUM/(XDAT(M+MBEG)-XDAT(L+MBEG))
            ENDIF
          ENDDO
        ENDIF
      ENDDO
      RETURN
      END



*/



double itrf2icrf(double jd, double utc, double *vt, double *vc)
{
    double xp = 0, yp = 0, ut1_utc = 0, dx = 0, dy = 0, delta_t, ut1, tt;

    geteop (utc, &xp, &yp, &ut1_utc, &dx, &dy);	

    delta_t = 32.184 + LEAPSECS - ut1_utc;
    ut1 = utc + ut1_utc;
    tt = utc + (LEAPSECS + 32.184); 	

    cel_pole (jd + tt / 86400.0, 2, dx * 1e3, dy * 1e3);

    ter2cel (jd, ut1 / 86400.0, delta_t, 1, ACCURACY, 0,
        xp, yp, vt, vc);

    return 0;


}




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


double chosephase (double sinvalue, double cosvalue)
{
    double sv = sinvalue, cv = cosvalue;
    if (sv >= 0 && cv >= 0) 
        return (asin (sv));
    if (sv > 0 && cv < 0) 
        return (acos (cv));
    if (sv < 0 && cv < 0) 
        return ( - asin (sv) + TWOPI / 2.0);
    else 
        return (asin (sv) + TWOPI);
}


/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
double itrf2icrf_dot(double jd, double utc, double *vt, double *rt, double *vc)
{
    double xp = 0, yp = 0, ut1_utc = 0, dx = 0, dy = 0, delta_t, ut1, tt,
        vc_1[3], vc_2[3];

    geteop (utc, &xp, &yp, &ut1_utc, &dx, &dy);	

    delta_t = 32.184 + LEAPSECS - ut1_utc;
    ut1 = utc + ut1_utc;
    tt = utc + (LEAPSECS + 32.184); 	

    cel_pole (jd + tt / 86400.0, 2, dx * 1e3, dy * 1e3);

    ter2cel (jd, ut1 / 86400.0, delta_t, 1, ACCURACY, 0,
        xp, yp, vt, vc_1);

    ter2cel_dot (jd, ut1 / 86400.0, delta_t, 1, ACCURACY, 0,
        xp, yp, rt, vc_2);

        
    vc[0] = vc_1[0] + vc_2[0];
    vc[1] = vc_1[1] + vc_2[1];
    vc[2] = vc_1[2] + vc_2[2];

    return 0;

}


/********sidereal_time */

short int sidereal_time_dot (double jd_high, double jd_low,
                         double delta_t,short int gst_type,
                         short int method, short int accuracy,

                         double *gst, double *gst_dot)
/*
------------------------------------------------------------------------

   PURPOSE:
      Computes the Greenwich apparent sidereal time, at Julian date
      'jd_high' + 'jd_low'.

   REFERENCES:
      Kaplan, G. (2005), US Naval Observatory Circular 179.

   INPUT
   ARGUMENTS:
      jd_high (double)
         High-order part of UT1 Julian date.
      jd_low (double)
         Low-order part of UT1 Julian date.
      delta_t (double)
         Difference TT-UT1 at 'jd_high'+'jd_low', in seconds
         of time.
      gst_type (short int)
         = 0 ... compute Greenwich mean sidereal time
         = 1 ... compute Greenwich apparent sidereal time
      method (short int)
         Selection for method
            = 0 ... CIO-based method
            = 1 ... equinox-based method
      accuracy (short int)
         Selection for accuracy
            = 0 ... full accuracy
            = 1 ... reduced accuracy

   OUTPUT
   ARGUMENTS:
      *gst (double)
         Greenwich apparent sidereal time, in hours.

   RETURNED
   VALUE:
      (short int)
         = 0         ... everything OK
         = 1         ... invalid value of 'accuracy'
         = 2         ... invalid value of 'method'
         > 10, < 30  ... 10 + error from function 'cio_rai'

   GLOBALS
   USED:
      T0, RAD2DEG

   FUNCTIONS
   CALLED:
      era                novas.c
      tdb2tt             novas.c
      e_tilt             novas.c
      cio_location       novas.c
      cio_basis          novas.c
      nutation           novas.c
      precession         novas.c
      frame_tie          novas.c
      fabs               math.h
      fmod               math.h
      atan2              math.h

   VER./DATE/
   PROGRAMMER:
      V1.0/06-92/TKB (USNO/NRL Optical Interfer.) Translate Fortran.
      V1.1/08-93/WTH (USNO/AA) Update to C programing standards.
      V1.2/03-98/JAB (USNO/AA) Expand documentation.
      V1.3/08-98/JAB (USNO/AA) Match flow of the Fortran counterpart.
      V2.0/09-03/JAB (USNO/AA) Incorporate the 2nd-reference changes.
      V2.1/08-04/JAB (USNO/AA) Incorporate the 1st-reference changes.
      V2.2/12-05/WKP (USNO/AA) Updated error handling.
      V2.3/01-06/WKP (USNO/AA) Changed 'mode' to 'method' and 'accuracy'.
      V2.4/04-06/JAB (USNO/AA) Use precession-in-RA terms in mean
                               sidereal time from third reference.
      V2.5/07-06/JAB (USNO/AA) Implement 'cio_location' construct.
      V2.6/06-08/WKP (USNO/AA) Changed value of direction argument in
                               call to 'nutation' from 1 to -1 for
                               consistency.

   NOTES:
      1. The Julian date may be split at any point, but for highest
      precision, set 'jd_high' to be the integral part of the Julian
      date, and set 'jd_low' to be the fractional part.
      2. This function is the C version of NOVAS Fortran routine
      'sidtim'.

------------------------------------------------------------------------
*/
{
   short int error = 0;
   short int ref_sys;

   static double ee;
   static double jd_last = -99.0;
   double unitx[3] = {1.0, 0.0, 0.0};
   double jd_ut, jd_tt, jd_tdb, tt_temp, t, theta, a, b, c, d,
      ra_cio, x[3], y[3], z[3], w1[3], w2[3], eq[3], ha_eq, st, st_dot,
      secdiff, eqeq;

/*
   Invalid value of 'accuracy'.
*/

   if ((accuracy < 0) || (accuracy > 1))
      return (error = 1);

/*
   Time argument for precession and nutation components of sidereal
   time is TDB.  First approximation is TDB = TT, then refine.
*/

   jd_ut = jd_high + jd_low;
   jd_tt = jd_ut + (delta_t / 86400.0);
   jd_tdb = jd_tt;
   tdb2tt (jd_tdb, &tt_temp,&secdiff);
   jd_tdb = jd_tt + (secdiff / 86400.0);

   t = (jd_tdb - T0) / 36525.0;

/*
   Compute the Earth Rotation Angle.  Time argument is UT1.
*/

   theta = era (jd_high, jd_low);

/*
   Compute the equation of the equinoxes if needed, depending upon the
   input values of 'gst_type' and 'method'.  If not needed, set to zero.
*/

   if (((gst_type == 0) && (method == 0)) ||       /* GMST; CIO-TIO */
       ((gst_type == 1) && (method == 1)))         /* GAST; equinox */
   {
      if (fabs (jd_tdb - jd_last) > 1.0e-8)
      {
         e_tilt (jd_tdb,accuracy, &a,&b,&ee,&c,&d);
         jd_last = jd_tdb;
      }
      eqeq = ee * 15.0;
   }
    else
   {
      eqeq = 0.0;
   }

/*
   Compute Greenwich sidereal time depending upon input values of
   'method' and 'gst_type'.
*/

   switch (method)
   {
      case (0):

/*
   Use 'CIO-TIO-theta' method.  See Circular 179, Section 6.5.4.
*/

/*
   Obtain the basis vectors, in the GCRS, of the celestial intermediate
   system.
*/

         if ((error = cio_location (jd_tdb,accuracy, &ra_cio,
            &ref_sys)) != 0)
         {
            *gst = 99.0;
            return (error += 10);
         }

         cio_basis (jd_tdb,ra_cio,ref_sys,accuracy, x,y,z);

/*
   Compute the direction of the true equinox in the GCRS.
*/

         nutation (jd_tdb,-1,accuracy,unitx, w1);
         precession (jd_tdb,w1,T0, w2);
         frame_tie (w2,-1, eq);

/*
   Compute the hour angle of the equinox wrt the TIO meridian
   (near Greenwich, but passes through the CIP and TIO).
*/

         ha_eq = theta - atan2 ((eq[0] * y[0] + eq[1] * y[1] +
            eq[2] * y[2]), (eq[0] * x[0] + eq[1] * x[1] +
            eq[2] * x[2])) * RAD2DEG;

/*
   For mean sidereal time, subtract the equation of the equinoxes.
*/

         ha_eq -= (eqeq / 240.0);

         ha_eq = fmod (ha_eq, 360.0) / 15.0;
         if (ha_eq < 0.0)
            ha_eq += 24.0;
         *gst = ha_eq;
         break;

      case (1):

/*
   Use equinox method.  See Circular 179, Section 2.6.2.
*/

/*
   Precession-in-RA terms in mean sidereal time taken from third
   reference, eq. (42), with coefficients in arcseconds.
*/

         st = eqeq + 0.014506 +
               (((( -    0.0000000368   * t
                    -    0.000029956  ) * t
                    -    0.00000044   ) * t
                    +    1.3915817    ) * t
                    + 4612.156534     ) * t;

/*
   Form the Greenwich sidereal time.
*/


         *gst = fmod ((st / 3600.0 + theta), 360.0) / 15.0;


//         st_dot = (((( - 5 * 0.0000000368 * t - 4 * 0.000029956 ) * t - 3 * 0.00000044 ) * t + 2 * 1.3915817 ) * t + 4612.156534 );

         st_dot = 
               (((( - 5 * 0.0000000368   * t
                    - 4 * 0.000029956  ) * t
                    - 4 * 0.00000044   ) * t
                    + 2 * 1.3915817    ) * t
                    +  4612.156534     );

         
         
         *gst_dot = 1.00273781191135448 * TWOPI / 86400.0 + st_dot / 206264.806 / 36525.0 / 86400.0;


         if (*gst < 0.0)
            *gst += 24.0;
         break;

/*
   Invalid value of 'method'.
*/

      default:
         *gst = 99.0;
         error = 2;
         break;
   }

   return (error);
}



/********ter2cel */

short int ter2cel_dot (double jd_ut_high, double jd_ut_low, double delta_t,
                   short int method, short int accuracy, short int option,
                   double x, double y, double *vect,

                   double *vecc)
/*
------------------------------------------------------------------------


------------------------------------------------------------------------
*/
{
   short int error = 0;
   short int rs, j;

   double jd_ut1, jd_tt, dummy, secdiff, jd_tdb, gast, gast_dot, r_cio, theta,
   v1[3], v2[3], v3[3], v4[3], xx[3], yy[3], zz[3];

/*
   Invalid value of 'accuracy'.
*/

   if ((accuracy < 0) || (accuracy > 1))
      return (error = 1);

/*
   Compute the TT Julian date corresponding to the input UT1 Julian
   date.
*/

   jd_ut1 = jd_ut_high + jd_ut_low;
   jd_tt = jd_ut1 + (delta_t / 86400.0);

/*
   Compute the TDB Julian date corresponding to the input UT1 Julian
   date.
*/

   jd_tdb = jd_tt;
   tdb2tt (jd_tdb, &dummy,&secdiff);
   jd_tdb = jd_tt + secdiff / 86400.0;

   switch (method)
   {
      case (0):

/*
   'CIO-TIO-THETA' method.

   See second reference, eq. (3) and (4).

   Apply polar motion, transforming the vector to the terrestrial
   intermediate system.
*/

         if ((x == 0.0) && (y == 0.0))
         {
            v1[0] = vect[0];
            v1[1] = vect[1];
            v1[2] = vect[2];
         }
          else
            wobble (jd_tdb,x,y,vect, v1);

/*
   Obtain the basis vectors, in the GCRS, of the celestial intermediate
   system.
*/

      if ((error = cio_location (jd_tdb,accuracy, &r_cio,&rs)) != 0)
         return (error += 10);

      if ((error = cio_basis (jd_tdb,r_cio,rs,accuracy, xx,yy,zz)) != 0)
         return (error += 20);

/*
   Compute and apply the Earth rotation angle, 'theta', transforming the
   vector to the celestial intermediate system.
*/

         theta = era (jd_ut_high,jd_ut_low);
         spin (-theta,v1, v2);

/*
   Transform the vector from the celestial intermediate system to the
   GCRS.
*/

         vecc[0] = xx[0] * v2[0] + yy[0] * v2[1] + zz[0] * v2[2];
         vecc[1] = xx[1] * v2[0] + yy[1] * v2[1] + zz[1] * v2[2];
         vecc[2] = xx[2] * v2[0] + yy[2] * v2[1] + zz[2] * v2[2];
         break;

      case (1):

/*
   Equinox mode.

   Apply polar motion.
*/

         if ((x == 0.0) && (y == 0.0))
         {
            for (j = 0; j < 3; j++)
            {
               v1[j] = vect[j];
            }
         }
          else
            wobble (jd_tdb,x,y,vect, v1);

/*
   Apply Earth rotation.
*/

         sidereal_time_dot (jd_ut_high,jd_ut_low,delta_t,1,1,accuracy, &gast, &gast_dot);
         spin_dot (-gast * 15.0, -gast_dot, v1, v2);

/*
   'option' = 1 skips remaining transformations.
*/

         if (option == 1)
         {
            vecc[0] = v2[0];
            vecc[1] = v2[1];
            vecc[2] = v2[2];
         }
          else
         {

/*
   Apply precession, nutation, and frame tie.
*/

            nutation (jd_tdb,-1,accuracy,v2, v3);
            precession (jd_tdb,v3,T0, v4);
            frame_tie (v4,-1, vecc);
         }
         break;

/*
   Invalid value of 'method'.
*/

      default:
         error = 2;
         break;
   }

   return (error);
}



/********spin */

void spin_dot (double angle, double angle_dot, double *pos1,

           double *pos2)
/*
------------------------------------------------------------------------

   PURPOSE:
      This function transforms a vector from one coordinate system
      to another with same origin and axes rotated about the z-axis.


------------------------------------------------------------------------
*/
{
   static double ang_last = -999.0;
   static double xx, yx, zx, xy, yy, zy, xz, yz, zz;
   double angr, cosang, sinang;

   if (fabs (angle - ang_last) >= 1.0e-12)
   {
      angr = angle * DEG2RAD;
      cosang = cos (angr);
      sinang = sin (angr);

/*
   Rotation matrix follows.
*/

//      xx =  cosang;
//      yx =  sinang;
//      zx =  0.0;
//      xy =  -sinang;
//      yy =  cosang;
//      zy =  0.0;
//      xz =  0.0;
//      yz =  0.0;
//      zz =  1.0;

      xx =  -sinang;
      yx =  cosang;
      zx =  0.0;
      xy =  -cosang;
      yy =  -sinang;
      zy =  0.0;
      xz =  0.0;
      yz =  0.0;
      zz =  0.0;

      ang_last = angle;
   }

/*
   Perform rotation.
*/

   pos2[0] = (xx * pos1[0] + yx * pos1[1] + zx * pos1[2]) * angle_dot;
   pos2[1] = (xy * pos1[0] + yy * pos1[1] + zy * pos1[2]) * angle_dot;
   pos2[2] = (xz * pos1[0] + yz * pos1[1] + zz * pos1[2]) * angle_dot;

   return;
}




void xyz2rtn(double *x, double *v, double *xyz, double *rtn)
{ 
    double scal_x, scal_v, vr[3], vn[3], vt[3];

    scal_x=sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
    scal_v=sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);

// c...unit vector in R direction
    vr[0]=x[0]/scal_x;
    vr[1]=x[1]/scal_x;
    vr[2]=x[2]/scal_x;
 
// c...unit direction in N direction
    vn[0]=(vr[1]*v[2]-vr[2]*v[1])/scal_v;
    vn[1]=(vr[2]*v[0]-vr[0]*v[2])/scal_v;
    vn[2]=(vr[0]*v[1]-vr[1]*v[0])/scal_v;
        
// c...unit direction in T direction
    vt[0]=(vn[1]*vr[2]-vn[2]*vr[1]);
    vt[1]=(vn[2]*vr[0]-vn[0]*vr[2]);
    vt[2]=(vn[0]*vr[1]-vn[1]*vr[0]);


// drtn(i,4)=dsqrt(dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2])
    rtn[0]=xyz[0]*vr[0]+xyz[1]*vr[1]+xyz[2]*vr[2];
    rtn[1]=xyz[0]*vt[0]+xyz[1]*vt[1]+xyz[2]*vt[2];
    rtn[2]=xyz[0]*vn[0]+xyz[1]*vn[1]+xyz[2]*vn[2];

    return;
}
