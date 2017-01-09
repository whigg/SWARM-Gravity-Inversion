
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/

#ifndef _SL1B_H_
    #define _SL1B_H_


    #ifndef __STDIO_H__
        #include <stdio.h>
    #endif

    #ifndef __MATH_H__
        #include <math.h>
    #endif

    #ifndef __STRING_H__
        #include <string.h>
    #endif

    #ifndef __STDLIB_H__
        #include <stdlib.h>
    #endif

    #ifndef __CTYPE_H__
        #include <ctype.h>
    #endif

    #ifndef __TIME_H__
        #include <time.h>
    #endif

    #ifndef _NOVAS_H_
        #include "novas.h"
    #endif

    #ifndef _NOVASCON_H_
        #include "novascon.h"
    #endif

    #ifndef _NUTATION_H_
        #include "nutation.h"
    #endif



    #ifndef __EPHMAN_H_
        #include "eph_manager.h"
    #endif

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/

    short int ACCURACY, LEAPSECS, INTERP;
    double JD0;
    char FILE_EOP[80];
//    int NDATA, DT, GPS_S;

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
    
    



    typedef struct InfStruct
    {
        int leaps;
        double jd0;
        double jdt;
        double mjd;
        double gmst;
        double deltat;
        double utc;
        double ut1;
        double tt;
        double xp;
        double yp;
        double ut1_utc;
        double dx;
        double dy;
        double c_ie[9];
        double c_ei[9];
        double c_iedot[9];
        double c_eidot[9];
    }InfStruct;





    void mt (double *a, int m, int n, double *b);
    void brmul (double *a, double *b, int m,int n, int k,double *c);


    double getinfo(double *tjd, InfStruct *info);

    short int getlps (double jd);

    void xyz2llh (double *vt, double *llh);

    double lagrange (double *y, int dim_y, int dim_x, double t, double *z);
    double lag_deri (double *y, int dim_y, int dim_x, double t, double *z, double *dz);
    void openeop (char file_eop[200], int mjds, int num, double *eopmat);

    void geteop (double utcsec, double *xp, double *yp, 
                double *ut1_utc, double *dx, double *dy);

    double itrf2icrf (double jd, double utc, double *vt, double *vc);

    double chosephase (double sinvalue, double cosvalue);

    double potential (double tjd, double *xic, double *xfc, 
                    double *ptt, double *rp);

    void xyz2rtn(double *x, double *v, double *xyz, double *rtn);

    short int sidereal_time_dot (double jd_high, double jd_low,
                             double delta_t,short int gst_type,
                             short int method, short int accuracy,
                             double *gst, double *gst_dot);

    short int ter2cel_dot (double jd_ut_high, double jd_ut_low, 
        double delta_t, short int method, short int accuracy, 
        short int option, double x, double y, double *vect, double *vecc);

    void spin_dot (double angle, double angle_dot, double *pos1,
               double *pos2);

    double itrf2icrf_dot(double jd, double utc, double *vt,
        double *rt, double *vc);


#endif


/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/

