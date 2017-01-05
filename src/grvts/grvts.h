/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/

#ifndef _GRVTS_H_
    #define _GRVTS_H_


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

    #ifndef _COORD_H_
        #include "coord.h"
    #endif

    #ifndef _NUMRS_H_
        #include "numrs.h"
    #endif


/*tides.c*/

        typedef struct
        {
            double ds;
            int argn[6];
            char name[5];
            int n;
            int m;
            double cp;
            double sp;
            double cm;
            double sm;
            double fnm;
        }OTSperturb;

        typedef struct 
        {
            double ds;
            int argn[6];
            double ang;
            double cosf;
            double sinf;
        }OTSconstit;


        extern int NMAJ_OTS, NALL_OTS;
        extern double *OTADM, *OPTM1, *OPTM2;
        extern OTSconstit *ALL_OTS;



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



/*sphcoe.c*/

        typedef struct CSStruct
        {
            int n;
            int m;
            int cs;
            double initv;
            double dadcsn[3];
            double dadcse[3];
            double dadcs[3];
        }CSStruct;


        extern CSStruct CSinfo[1000];
        extern int MGCS;



        double eig_open (char *grv_name, double dyear, 
                int nmax, int mmax, double *coef);

        double grv_open (char *grv_name, char *label, 
                int nmax, int mmax, double *coef);

        double cs2ada (double *llr, double *cs, double gm, double a, 
                int nmax, double *ae, 
                int part, double *dadre, int flagdadcs);

        double cs2ac (double *llr, double *cs, double gm, double a, 
                int nmax, double *acc);

        double cspt2gp (double *pt, double *cs, int nmax, 
                double *gpt);
    
        double cs2gp (double *llr, double *cs, double gm, double a, 
                int nmax, double *gpt);

        double cs2gp_pt (double *llr, double *cs, double gm, double a, 
                int nmax, double *gpt, double *pt);

        double lgdr(double t, int nmax, int m, double *pbar);

        double lgdr2i(double t, int nmax, int m, 
             double *pbar, double *pbar1, double *pbar2);

        double lgdr1i(double t, int nmax, int m, 
             double *pbar, double *pbar1);

        double lgdr2(double t, int nmax, int m, 
             double *pbar, double *pbar1, double *pbar2);


        double normfct (int n, int m);

        double addsubcs (double *cs, double *cs1, int nmax, int nmax1, 
            int flag);

        double zero0zero1 (double *cs, int nmax);

        double unit0zero1 (double *cs, int nmax);

        double csnm (double *cs, int nmax, int n, int m, 
            double *c, double *s);




#endif


/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/


