
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/

#ifndef _SMOD_H_
    #define _SMOD_H_

    #ifndef REAL128
       #define REAL128
       typedef long double real128;
    #endif

    #ifndef __STDIO_H__
        #include <stdio.h>
    #endif

    #ifdef MSDOS
        #ifndef __CONIO_H__
            #include <conio.h>
        #endif
    #endif

    #ifdef LINUX
        void getch(void) {}
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

    #ifndef __EPHMAN_H_
        #include "eph_manager.h"
    #endif

#ifndef _COORD_H_
    #include "coord.h"
#endif
#ifndef _GRVTS_H_
    #include "grvts.h"
#endif


/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/




    #define DYNPAR 4
    
    short int CENTER, GRAVDEGREE, GRAVORDER,
        ACCURACY, LEAPSECS, PERB[11];
    int DIM_TE, DIM_OR, DIM_AC, DIM_SC, DIM_OT, DIM_AT;
    double JD0, STEP_OR, STEP_TE, *TE_EPH, *OR_EPH, *AC_EPH, *SC_EPH, GMDE[11], 
        SRPJ, ASRP, CONS, DCON, BIAS, DBIA, SATMASS, SATAREA;
    char FILE_GRV[80], FILE_EOP[80];
    double *cn0, *cnm, *snm, gma[2];
    double j2, j3, j4, jc21, js21, jc22, js22;


    int PERMT, STIDE, OTIDE,ATIDE, PTIDE, AOD1B, NBODY, RELTIV, NFES, NMAX, MMAX, NSMAX, NOMAX, NEOP, NAMAX, NPMAX, NDMAX,
        CT;
    double C20PERM, AMR, *COEF, *COEFG, *COEFS, *COEFO, *COEFP, *COEFA, *COEFD, *AOD_EPH, *OT_EPH, *AT_EPH, *EOPMT, *OPTM1, *OPTM2, GMCT, RCT;

    int IBASB, IBAST, ISRPB, ISRPT, IK2;
    
    double SRPB, SRPT, BASB, BAST, K2, AX0, AY0, AZ0, AX1, AY1, AZ1;
    
    int MEST, SLOVEFOR, MSTA, MOBS, MSRP, MGCS, MACC, MTK2, MDYN, MSOL;
    int MACC_DTBS, MACC_DTSL, MACC_NOBS, MACC_NOSL, MACC_BIAS, MACC_SCAL, MACC_ARBS, MACC_ARSL, MACC_PRBS, MACC_PRSL;
    
    double *MACC_EPBS, *MACC_EPSL, TT0, *DADBS, *DADSL;
    


double lagrange (double *y, int dim_y, int dim_x, double t, double *z);

double obs_vel (double jd, double utc, double *obs, int part, double *bmat);
void initsolvefor (double *xsm, double *x);
void updsolvefor (double *x);
void getsolvefor ();

double stidecs_k2(InfStruct *info, double k2, double *stcs, int *body, int nbody);

    void pt_orb (double ts_orb, double te_orb, double step_orb, int dim_eph);

    double obs_alt (double jd, double utc, double *obs, int part, double *bmat);
    double obs_dsn (double jd, double utc, double *obs, int part, double *bmat);

    int get_ephemeris (double tjd[2], int to, int from, double *x);

    double accel_pm_part (double *tjd, double *x, double *fnt, 
        int part, double *dadr);



    double accel_nb_part (double *tjd, double *xic, double *fnt, 
        int part, double *dadr);

    double f_bcrs (double *jd, double *xi, int exclude, 
                   double *acc, int part, double *dadr);

    double accel_sr_part (double *tjd, double *xic, double *acc,
        int part, double *dadr, double *dadpb, double *dadpt);

    double accel_gt_part (double *tjd, double *xic, double *ag, 
        int part, double *dadr, double *dadk2);

    double accel_ac_part (double *tjd, double *xic, double *acc);

    void readacc1 (char *infile);
    void readsca1 (char *infile);

    double cal_acc_1(int accid);

    double obs_mex (double utc, double *xsim, int part, double *xbmat);
    double openoteph (double tts, double tte);
    double openopt (char *f_pt, int nmax);

   
//    double stidecs_Anelastic(double *tjd, double *c_ie, int id_perm, double *stcs);
    double stidecs_Anelastic(InfStruct *info, int id_perm, double *stcs);


    double accel_gravtide (double *tjd, double *xic, 
            double *ag, double *at, double *ao);

    double stidecs(double *tjd, double *c_ie, int id_perm, double *stcs);

    double accel_pmiers (double *tjd, double *x, double *fnt, double *fgr);



    double earth_fullaccel (double jd, double tt, double *xic, double *fxic);

    double accel_point (double *tjd, double *x, double *fnt, double *fgr);

    double accel_slrad (double *tjd, double *xic, double *acc);

    double accel_nbody (double *tjd, double *xic, double *fnt, double *fgr);

    double force_bcrs (double *jd, double *xi, short int exclude, 
                   double *fnt, double *fgr);

    
    void fun_accel (int dim, double jd, double tt, double *xic, double *fxic);

    double opengrv (char file_grv[2][200], double *coef, int nmax, int mmax);

    double earth_pointmass (double jd, double tdbs, double *x, double *f);

    double accel_gravt (double *tjd, double *xic, double *a4);




    void geteop (double mjd, double *xp, double *yp, 
               double *ut1_utc, double *dx, double *dy);



double stidecs_old(double *tjd, double gma1, double k2, 
               double *c20, double *c21, double *s21, double *c22, double *s22);


double al_part (double *xc2, double *dodx, double *dodp);


double simula_altim (double tdb, double *calculable, short int part, double *bmat);

void azelev (double jd_ut1, double delta_t, short int accuracy,
             double x, double y, double *llh, double ra,
             double dec, double *zd, double *az);


double delta_iid (double *jd, double *xi, double *ii, double *id);

double delta_tdb (double *txice, double *txics, double *deltat);

double fun_fullaccel (double tdbs, double *xic, double *fxic);

double fun_fullstate (double tdbs, double *state, double *fstate);

double accel_bcrs (double *jd, double *xi, short int part, short int exclude, 
                   double *acc, double *dadr, double *dadp);

double accel_ntrel (double *tjd, double *xic, short int part, 
                    double *acc, double *dadr, double *dadp);

double accel_nonsp (double *tjd, double *xic, short int part, 
                    double *acc, double *dadr, double *dadp);

double accel_radpr (double *tjd, double *xic, short int part, 
                    double *acc, double *dadr, double *dadp);

void in2pa (double *jd, double *te);



double fun_pointmass (double tdbs, double *x, double *f);



double simula_phase (double utc3, double utc0, double *station3, 
                     short int uplink, double *station1, short int genrel,
                     double *calculable, double *azimuth, double *elevation, 
                     short int part, double *bmat);

double simula_dople (double utc3, double tc, double *station3, 
                     short int uplink, double *station1, short int genrel,
                     double *calculable, double *azimuth, double *elevation, 
                     short int part, double *bmat);


double simula_range (double utc3, double *station3, short int uplink, 
                     double *station1, short int genrel, 
                     double *calculable, double *azimuth, double *elevation, 
                     short int part, double *bmat);

double ltsolution (double utc_3, double *station3, short int uplink,  
                   double *station1, short int genrel, real128 *lt, 
                   double *azimuth, double *elevation, short int part,
                   double *bmat, double *txic);

real128 lt_form (double tdb3, double tdb2, double *re3, double *rp2, 
                   int genrel, real128 *rs3, real128 *rs2);

double lt_part (real128 *rs3, real128 *rs2, real128 *rs1, int uplink, 
                  double *dodx, double *dodp);

double opengravfile (double *cn0, double *cnm, double *snm, double *gma);





#endif


/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/

