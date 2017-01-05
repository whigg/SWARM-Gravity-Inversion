/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/

#ifndef _COORD_H_
    #define _COORD_H_


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

    #ifndef _EPHMAN_H_
        #include "eph_manager.h"
    #endif

//    #ifndef _GRVTS_H_
//        #include "grvts.h"
//    #endif

    #ifndef _NUMRS_H_
        #include "numrs.h"
    #endif


/* Reference epoch (J2000), Julian Date */
    #define DJ00 (2451545.0)

/* Days per Julian millennium */
    #define DJM (365250.0)

/* Degrees to radians */
    #define DD2R (1.745329251994329576923691e-2)

/* 2Pi */
    #define D2PI (6.283185307179586476925287)



        #define I_TITAN 99

        extern double *EOPEPH;
        extern int NEOP;


        typedef struct
        {
            int leaps;
            double jd0;
            double jdt;
            double jdtt[2];
            double mjd;
            double mjdutc;
            double gps;
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
            double wi[3];
        }InfStruct;
    
    



/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/*tedot.c*/
        short int ter2cel_dot (double jd_ut_high, double jd_ut_low, 
                    double delta_t, short int method, 
                    short int accuracy, short int option,
                    double x, double y, double *vect,
                    double *vecc);

        short int sidereal_time_dot (double jd_high, double jd_low,
                    double delta_t,short int gst_type,
                    short int method, short int accuracy,
                    double *gst, double *gst_dot);

        void spin_dot (double angle, double angle_dot, double *pos1,
                    double *pos2);


/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/*xyz.c*/
        int aei2xyz (double ele[6], double gm, 
                    double pos[3], double vel[3]);

        double kepler (double M,double e);
        
        double xyz2aei (double pos[3], double vel[3], double gm,
                        double ele[6]);

        void xyz2llr (double *vt, double *llr);
            
        void llr2xyz (double *llr, double *vt);

        double chosephase (double sinvalue, double cosvalue);

        void xyz2rtn (double *x, double *v, double *xyz, double *rtn);


/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/*xform.c*/

        double gps_GRACE2utc (double jd0_day, double gps_GRACE);

        double gps_GRACE2tt (double jd0_day, double gps_GRACE);

        double gps2utc (double jd0_day, double gps_sec);

        double gps2tt (double gps_sec);

        int eop_open (char *eop_name, int mjd_beg, int mjd_end);

        int eop_read (double mjdt, double *xp, double *yp, 
            double *ut1_utc, double *dx, double *dy);
    
        int eop_close ();

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

        void dpleph_ (double *tjdtdb, long int *targe, long int *center,
                double *posvel);
    
 
/*dtdb.c*/
        double iauDtdb(double date1, double date2,
               double ut, double elong, double u, double v);
   
    
#endif
    

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/

