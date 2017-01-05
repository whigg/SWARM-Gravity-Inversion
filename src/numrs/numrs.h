/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/

#ifndef _NUMRS_H_
    #define _NUMRS_H_


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

//    #ifndef _GRVTS_H_
//        #include "grvts.h"
//    #endif



/*itgitp.c*/

        double lgr_order (double *y, int dim_y, int dim_x, 
            double t, double *z, int order);

        double intydt(int num, double dt, double *y, double *val);

        double rkf78_auto (double h, double t, double *x, int dim, double err, 
              double (*fun)(double,double *,double *), int autoadjust);

        double rkf78 (double jd, double t, double h, double *x, int dim, 
              void (*fun)(int, double, double,double *,double *));




/*matvec.c*/
        int brinv (double *a,int n);

        int brank(double *a, int m, int n);

        void choldc(double *a, int n, double p[]);

        void cholsl(double *a, int n, double p[], 
                double b[], double x[]);

        void solvels_chol(double *a, int n, 
            double *y, double *x, int nocov);

        void solvegaus(double *a, int n, double *y, double *x);

        double modvect (double *v);

        double dotvect (double *v1, double *v2);

        void crsvect (double *v1, double *v2, double *v);

        void mt (double *a, int m, int n, double *b);
    
        void brmul (double *a, double *b, int m,int n, int k,double *c);

        void rotmatx (double rad, double *matx, short int deri);

        void rotmaty (double rad, double *maty, short int deri);

        void rotmatz (double rad, double *matz, short int deri);

        int bssgj (double *a,int n);


/*stats*/

        void mgrns(double u, double g, double *r,int n,double a[]);

        double mgrn1(double u, double g, double *r);
    

        double mean (double * array, double N);

        double std_dev (double * array, double N);


/*dftflt*/


        int dft(long int length, double real_sample[], double imag_sample[]);

        int inverse_dft(long int length, double real_sample[], double imag_sample[]);
    
        double daydft (double *obs, double *fft, int num,  double dt, double tc);
    
        double getwinhp(int n, double *win, double kc, int m);

        double getwinbp(int n, double *win, double kc, double kh, int m);
 
        double getwinbp0(int n, double *win, double kc, int m);

        double dayfft (double *data, double *fft, int num,  double fs);

        void kkfft(double pr[], double pi[], int n, int k, double fr[],
           double fi[], int l, int il);

        double resfft(double *res, double *resflt, int num);

        void realft(double data[], int n, int isign);

        void four1(double data[], int nn, int isign);

        double k5filter (double x[], double y[], int size, double fln,
                 double fhn, double fs, int order);

        void firwin(int n, double fln, double fhn,double h[]);

        double hanningwin(int n,int i);

        double kaiserwin(int n,int i);


        double lsfl1c1d (double *xmax, double *x, double *y, int nmax, int n,
            double tpsec, double *ymaxflt, int nply, int ncpr, int nplycpr);

        double lsf_cpr_new (double *xmax, double *x, double *y, int nmax, int n,
            double tpsec, double *yflt, int order_poly, int order_cpr);


        int lsf_cpr_day ( double *x, double *y, int n, double tpsec, 
            double *yflt, int order_poly, int order_cpr);


        int lsf_cpr ( double *x, double *y, int n, double tpsec, 
            double *yflt, int order_poly, int order_cpr);

        int lsf_poly ( double *x, double *y, int n, double *a, int k);







#endif


/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/


