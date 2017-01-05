
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/

#ifndef _L1B2L1C_H_
    #define _L1B2L1C_H_

    #ifndef __OMP_H__
        #include <omp.h>
    #endif

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

    #ifndef _NOVAS_
       #include "novas.h"
    #endif

    #ifndef _EPHMAN_
       #include "eph_manager.h"
    #endif

#ifndef _GRVTS_H_

#include "grvts.h"
    #define _GRVTS_H_
#endif



#ifndef _COORD_H_

#include "coord.h"
    #define _COORD_H_
#endif


#ifndef _NUMRS_H_

#include "numrs.h"
    #define _NUMRS_H_
#endif

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/


    typedef struct GNV1C
    {
        int gps_time;
        char GRACE_id;
        char coord_ref;
        double xpos;
        double ypos;
        double zpos;
        double xvel;
        double yvel;
        double zvel;
        double pos[3];
        double vel[3];
        double lat;
        double lon;
        double r;
        double ele[6];
        double Tp;
    }GNV1C;



    typedef struct KBR1B
    {
        int gps_time;
        double biased_range;
        double range_rate;
        double range_accl;
        double iono_corr;
        double lighttime_corr;
        double lighttime_rate;
        double lighttime_accl;
        double ant_centr_corr;
        double ant_centr_rate;
        double ant_centr_accl;
        int K_A_SNR;
        int Ka_A_SNR;
        int K_B_SNR;
        int Ka_B_SNR;
        int qualflg;
        double range;
        double rate;
        double accl;
    }KBR1B;


    typedef struct ACC1B
    {
        int gps_time;
        char GRACE_id;
        double lin_accl_x;
        double lin_accl_y;
        double lin_accl_z;
        double ang_accl_x;
        double ang_accl_y;
        double ang_accl_z;
        double acl_x_res;
        double acl_y_res;
        double acl_z_res;
        int qualflg;
    }ACC1B;


    typedef struct SCA1B
    {
        int gps_time;
        char GRACE_id;
        int sca_id;
        double quatangle;
        double quaticoeff;
        double quatjcoeff;
        double quatkcoeff;
        double qual_rss;
        int qualflg;
    }SCA1B;




    



    typedef struct DATAL1C
    {
        int error;  //0: normal; 1: no orbit; 2: no kbrr; 3: no kbra; 9: outlier
        double t;
        
        double vl1cflt;
        double glosflt;

        double vl1c;
        double fitacc;
        double vpwfit;
        double vpwfitstd;
        double vpwfitres;
        double vmvfit;
        double vmvfitstd;
        double vmvfitres;
        double vc;      // L1C potential difference 

        double vl1cm;
        double vcm;
        double vpwfitm;

        double vcsr;      // Official L2 (n = 0 ~ 60)
        double vjpl;      // Official L2 (n = 0 ~ 60)
        double vgfz;      // Official L2 (n = 0 ~ 60)
        double vggm;      // Official L2 (n = 0 ~ 60)
   
        double acsr[3];
        double ajpl[3];
        double agfz[3];
        double glos;
        double gpwfit;
        double gpwfitstd;
        double gpwfitres;
        double gmvfit;
        double gmvfitstd;
        double gmvfitres;
        double gcsr;
        double gjpl;
        double ggfz;
        
        double vcorr;
        double dvcorr;
        double vref;
        double gref;
        double gc;
        double gcorr;
//        double losg;
//        double losb;
//        double dvtg;

//        double vg; 
//        double dvg;
//        double ag[3];  


        double vk;      // vk = vk1 + vk2;

        double vr;      // rotation term
        double vr3;
        double vr12;
        double vrm;


        double vrfl;
        double arfl[3];   
        double grfl; 
        double vrfh;  
        double arfh[3];
        double grfh;
        double vb;  
        double ab[3];
        double gb;
        double vbl;  
        double abl[3];
        double gbl;
        double dvab;
        double vab; 
        double dvrb;
        double vrb;
        double dvtb;
        double vtb;        
    

        double af[3];  
        double dvaf;
        double vaf; 
        double dvrf;     
        double vrf;    

        double vn;
        double an[3];
        double dvan;
        double van;
        double dvrn;     
        double vrn;
        double dvtn;
        double vtn;

        double vs;
        double as[3];
        double dvas;
        double vas;
        double dvrs;     
        double vrs;
        double dvts;
        double vts;

        double vo;
        double ao[3];
        double dvao;
        double vao;
        double dvro;     
        double vro;
        double dvto;
        double vto;

        double vh;
        double ah[3];
        double dvah;
        double vah;
        double dvrh;     
        double vrh;
        double dvth;
        double vth;


        double va;
        double aa[3];
        double dvaa;
        double vaa;
        double dvra;     
        double vra;
        double dvta;
        double vta;

        double vp;
        double ap[3];
        double dvap;
        double vap;
        double dvrp;
        double vrp;
        double dvtp;
        double vtp;

        double ag[3];
        double dvag;
        double vag;
        double dvrg;
        double vrg;


        double rp[3];
        double llr[3];
        double pos;
        double rv[3];
        double vel;
        double ra[3];

        double wr[3];

        double range;
        double rate;
        double accl;
        double Tp1; 
        double Tp2; 
   
    }DATAL1C;


    double *ACA_EPH, *ACB_EPH, *SCA_EPH, *SCB_EPH, *AOD_EPH, *OT_EPH, *AT_EPH, *OPTM1, *OPTM2;
    int DIM_ACA, DIM_ACB, DIM_SCA, DIM_SCB;

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/

    GNV1C *gnva, *gnvb;

    KBR1B *kbrx;

    ACC1B *acca, *accb;

    SCA1B *scaa, *scab;
    
    InfStruct *info;

    DATAL1C *sat1a, *sat2b, *sat12;

    int NDATA, DT, GPS_S, NFES, PERMT, STIDE, OTIDE, ATIDE, 
        MACC_BIAS, MACC_DTBS, MACC_SCAL, MACC_DTSL;
    double GMA[2], C20PERM;
    short int ACCURACY, INTERP;

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/

    double kbrvel2pos (int i);
    double kbrpos2vel (int i);

    double nbodyv (double *tjd, double *rp, double *vn, double *dvtn, double *an);

    double stidev (int num, double *llr1, double *vs, double *dvts, double *as);


    void reltivint ();
    double reltivdiff (double vb1, double vb2, double *vrel);
    void fitres (int order_poly, int order_cpr, int offsetcyc);

    void fitl1c_pw (int order_poly, int order_cpr);
    void fitlos_pw (int order_poly, int order_cpr);

    void fitl1c_mv (int nply, int ncpr, int nplycpr, int ndt);
    void fitlos_mv (int nply, int ncpr, int nplycpr, int ndt);
//    void fitres_new (int order_poly, int order_cpr);
    void fitres_new_overlap (int order_poly, int order_cpr, int overlap);


    double getinfo_1day(char *infile);


    double nbodypt (double *tjd, double *p, double *pt, double *a);

    double cspt2gp (double *pt, double *cs, int NMAX, double *gpt);
    double cs2pt (double *llr, double *cs, double gm, double a, int NMAX, double *gpt, double *pt);
    double cs2vdvdt (double *llr, double *cs, double gm, double a, int NMAX, double *v, double *dvdt, double *pt);
    double cs2acc (int num, double *llr, double *cs, double gm, double a, int nmax, 
                   double *v, double *dvdt, double *acc);


    double opengrav (char file_grv[2][200], double *coef, double *gma, 
        int nmax, int mmax, int head);


    double disse_a (double *v1, int i, double *ef1i, double *f1);
    double disse_b (double *v1, int i, double *ef1i, double *f1);


    int readkbr (char *infile, int *n);
    int readgnv (char *infile_a, char *infile_b, int *n_a, int *n_b);
    int readacc (char *infile_a, char *infile_b);
    int readsca (char *infile_a, char *infile_b);
    double cal_acc_01(void);
    double cal_acc_00(void);
    int calbias(char *infile_a, char *infile_b);
    int calbiaseph(char *infile_a, char *infile_b);
    double accia (double *tjd, int label, double *acc);
    double cal2_acc(void);
    double bias_acc(int order);
  

    double quat2mat_i2s (double *qvec, double *mat);

    double accgr (double *tjd, double *p, double *v, double *fgr);




#endif


/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/

