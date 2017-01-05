/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/*

  GRACE/GRAIL L2 solver (MPI, MKL, OMP modules)

  Version: June 2014

  Copyright (c) 2014 Kun Shang (shang.34@osu.edu) All Right Reserved 

*/
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
//#include <omp.h>
#include "mkl.h"
//#include "mkl_pblas.h"
#include "mpi.h"
//#include "mkl_scalapack.h"

#define LAPACK_ROW_MAJOR 101 
#define LAPACK_COL_MAJOR 102

#define max(A,B) ((A)>(B)?(A):(B))
#define min(A,B) ((A)>(B)?(B):(A))

   const double T0 = 2451545.00000000;
   const double TWOPI = 6.283185307179586476925287;
   const double ASEC2RAD = 4.848136811095359935899141e-6;
   const double DEG2RAD = 0.017453292519943296;
   const double RAD2DEG = 57.295779513082321;

        void cal_date (double tjd,
            short int *year, short int *month, short int *day,
            double *hour);

        double julian_date (short int year, short int month, short int day,
            double hour);

        double grv_open (char *grv_name, char *label, 
            int nmax, int mmax, double *coef);

        double cs2gp_pt (double *llr, double *cs, double gm, double a, int nmax,
            double *gpt, double *pt);

        double lgdr(double t, int nmax, int m, double *pbar);

        double zero0zero1 (double *cs, int nmax);

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




/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
int main (int argc, char *argv[])
{
    int ndata, ndata_all, ndata_avg, dt;
    double GMA[2];    
    FILE *fp_stdin, *fpin_l1c, *fpout_adj, *fpout_coef;  
    int i, n, nmax, mmax, slvls, l_old, ndat_max, factor_b1, factor_b2, 
        m, l, k, ind, num, *num_arc, *gps_arc, narc_max, *ind_m,
        LFIX = 0, ifix[4000], LSTK = 0, istk[4000], 
        nply, ncpr, nemp, ndt, narc, npd, *gpst, *gps_eph, 
        solvefor, solvegrv, solveemp;
    short int year, month, day;
    double *tp1, *csr, *gfz, *jpl, *yi,
        *obs, GMR2, r2R, hour,
        JD0, t, Tprd, 
        *pt1, *pt2, *ai, *coefaam, factor_n, 
        *lambda, *ksiupd, *k0, *aksi, *l1cf, 
        *factor_a, *factor_b, factor_p, factor_r,
        *ksi, *dksi, *dksiupd, r0[1], sigma0, pweight, std, *p0, *p0b0, 
        *knk0, *knk, *kksi, *id, *kfix, *kfixt,*kn, *nk, *nkknk,
        *atpa, *atpy, *atpai, *atpyi, *atpa_one, *cest, 
        omega, etpe, c, s, vfix[4000], vstk[4000],
        *llr1, *llr2, *l1c_eph, *l1ct, vb1, vb2;
    char line[500], card[20],
        f_aam[2][200]={"GIF48.GEO", "\0"},
         stdname[200], l1cfile[400], l2file[400];    
    time_t s0, s1, s2, s3, s4, s5, s6;
    long nl, solveforN, solvegrvN, solvefor2, ns, is, in;
    double tbeg, tend;
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/

    int ndata_s, ndata_e, npara_s, npara_e, npara_all, npara_avg, npara;

    int i_one = 1, i_negone = -1, i_zero = 0;
    double zero=0.0E+0, one=1.0E+0;

    double int_lim, mem_lim, int_lim1, mem_lim1, int_lim2, mem_lim2, alpha;

    int myrank, root_process, nprocs, ierr, num_intervals, N_s, N_e, iter;

    slvls = 4;

    ierr = MPI_Init(&argc, &argv);

    root_process = 0;

    myrank = 0;
    nprocs = 1;

    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    ierr = MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    printf ("(%4d) nprocs = %4d\n", myrank, nprocs);




    if (argc != 3)
    {
        printf ("usage: L2solver l1cfile nmax\n");
        exit(0);
    }
    else
    {
        sscanf (argv[1], "%s", l1cfile);
        sscanf (argv[2], "%d", &nmax);
    }

    if ( (fpin_l1c = fopen (l1cfile,"r")) == NULL)
    {
        printf ("Cannot open fpout_a file!\n");
        exit (0);
    }

    s0 = time(NULL);

    dt = 5;

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
// distribute ndata_all


    ndata_all = 0;
    while (1)
    {
        if (fgets(line,400, fpin_l1c) ==NULL) break;
        ndata_all ++;
    }
    rewind(fpin_l1c);

    ndata_avg = (int)(ndata_all / nprocs) + 1;

    ndata_s = myrank * ndata_avg;
    ndata_e = (myrank + 1) * ndata_avg - 1;
    if (ndata_e > ndata_all - 1)
        ndata_e = ndata_all - 1;
    ndata = ndata_e - ndata_s + 1;
    if (ndata < 0) ndata = 0;
    printf("(%4d) ndata_all = %d ndata_avg = %d ndata = %d ndata_s = %d ndata_e = %d\n",
        myrank, ndata_all, ndata_avg, ndata, ndata_s, ndata_e);


/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
// distribute solvefor

    mmax = nmax;

    npara_all = (nmax + 1) * (nmax + 1);
    printf("(%4d) nmax = %d npara_all = %d\n",
        myrank, nmax, npara_all);


/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
// index and memeory check

    int_lim1 = (double)npara_all * (double)ndata /1024.0/1024.0/1024.0;
    int_lim2 = (double)npara_all * (double)npara_all /1024.0/1024.0/1024.0;

    mem_lim1 = (double)npara_all * (double) ndata / 1024.0/1024.0/1024.0 * 8;
    mem_lim2 = (double)npara_all * (double) npara_all / 1024.0/1024.0/1024.0 * 8;

    mem_lim = mem_lim1 + mem_lim2;

    printf("(%4d) index  check: ", myrank);
    if (int_lim1 > 2)
        printf ("int_data ~= %9.3fG > INT_LIM (2G) : warning!!!\n", int_lim1);
    else
        printf ("int_data ~= %9.3fG < INT_LIM (2G) : safe~~~\n", int_lim1);
    printf("(%4d) index  check: ", myrank);
    if (int_lim2 > 2)
        printf ("int_para ~= %9.3fG > INT_LIM (2G) : warning!!!\n", int_lim2);
    else
        printf ("int_para ~= %9.3fG < INT_LIM (2G) : safe~~~\n", int_lim2);

    printf("(%4d) memory check: ", myrank);
    if (mem_lim1  > 4)
        printf ("mem_data ~= %9.3fG > MEM_LIM (4G) : warning!!!\n", mem_lim1);
    else
        printf ("mem_data ~= %9.3fG < MEM_LIM (4G) : safe~~~\n", mem_lim2);
    printf("(%4d) memory check: ", myrank);
    if (mem_lim2  > 4)
        printf ("mem_para ~= %9.3fG > MEM_LIM (4G) : warning!!!\n", mem_lim2);
    else
        printf ("mem_para ~= %9.3fG < MEM_LIM (4G) : safe~~~\n", mem_lim2);
    printf("(%4d) memory check: ", myrank);
    if (mem_lim  > 4)
        printf ("mem_all  ~= %9.3fG > MEM_LIM (4G) : warning!!!\n", mem_lim);
    else
        printf ("mem_all  ~= %9.3fG < MEM_LIM (4G) : safe~~~\n", mem_lim);


/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
// calloc local ndata






    l1c_eph = (double *) calloc (10 * ndata_all, sizeof(double));
    llr1 = (double *) calloc (3 * ndata, sizeof(double));
    llr2 = (double *) calloc (3 * ndata, sizeof(double));
    yi = (double *) calloc ( ndata, sizeof(double));



    i = 0;
    while (1)
    {
        if (fgets(line,400, fpin_l1c) ==NULL) break;
        sscanf (line, "%*d%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf",
            &l1c_eph[i * 10], &l1c_eph[i * 10 + 1], &l1c_eph[i * 10 + 2],
            &l1c_eph[i * 10 + 3], &l1c_eph[i * 10 + 4], &l1c_eph[i * 10 + 5],
            &l1c_eph[i * 10 + 6],
            &l1c_eph[i * 10 + 7], &l1c_eph[i * 10 + 8], &l1c_eph[i * 10 + 9]);
        i ++;
    }
    if (i != ndata_all)
    {
        printf ("(%4d) ndata_all (%d) != i (%d)\n", myrank, ndata_all, i);
        free (l1c_eph);
        fclose(fpin_l1c);
        exit(0);
    }

    fclose(fpin_l1c);




    n = 0;
    for (i = ndata_s; i <= ndata_e; i ++)
    {
        llr1[n * 3]     = l1c_eph[i * 10];
        llr1[n * 3 + 1] = l1c_eph[i * 10 + 1];
        llr1[n * 3 + 2] = l1c_eph[i * 10 + 2];
        llr2[n * 3]     = l1c_eph[i * 10 + 3];
        llr2[n * 3 + 1] = l1c_eph[i * 10 + 4];
        llr2[n * 3 + 2] = l1c_eph[i * 10 + 5];
        yi[n]         = l1c_eph[i * 10 + 6];
        n++;
    }
    if (n != ndata)
    {
        printf ("(%4d) ndata (%d) != n (%d)\n", myrank, ndata, i);
        free (l1c_eph);
        exit(0);
    }
    free (l1c_eph);



    ai   = (double *) calloc ( npara_all * ndata, sizeof(double));
    atpai = (double *) calloc ( npara_all * npara_all, sizeof(double));
    atpyi = (double *) calloc ( npara_all,         sizeof(double));
    ksi  = (double *) calloc ( npara_all,         sizeof(double));
    dksi = (double *) calloc ( npara_all,         sizeof(double));


    coefaam  = (double *) calloc (npara_all, sizeof(double));

    GMA[0] = 398600.44150E+09;
    GMA[1] = 6378136.3;

    grv_open (f_aam[0], f_aam[1], nmax, mmax, coefaam);
    zero0zero1 (coefaam, nmax);



    s1 = time(NULL);
    printf("(%4d) %5ld: seconds of initialization and scan data\n", myrank, s1-s0);
    fflush(stdout);


/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/

    pt1  = (double *) calloc ( npara_all, sizeof(double));
    pt2  = (double *) calloc ( npara_all, sizeof(double));

    for (i = 0; i < ndata; i ++)
    {
        cs2gp_pt (&llr1[i * 3], coefaam, GMA[0], GMA[1], nmax, &vb1, pt1);
        cs2gp_pt (&llr2[i * 3], coefaam, GMA[0], GMA[1], nmax, &vb2, pt2);
        for(n = 0; n < npara_all; n++)
        {
            ai[i * npara_all + n] = pt2[n] - pt1[n];
        }
    }

    free (llr1);
    free (llr2);
    free (pt1);
    free (pt2);


    s2 = time(NULL);
    printf("(%4d) %5ld: seconds of processing gravity and partial\n", myrank, s2-s1);
    fflush(stdout);

    std = 0.8e-3;
    pweight = 1.0/std/std;
    alpha = 1.0;
    alpha = pweight;

    cblas_dsyrk (CblasRowMajor, CblasUpper, CblasTrans, npara_all, ndata, alpha, ai, npara_all, 0.0, atpai, npara_all);
    cblas_dgemv (CblasRowMajor, CblasTrans, ndata, npara_all, alpha, ai, npara_all, yi, 1, 0.0, atpyi, 1);



    s3 = time(NULL);
    printf("(%4d) %5ld: seconds of finish ATPA&ATPY\n", myrank, s3-s2);
    fflush(stdout);

    free (ai);


    atpa = (double *) calloc ( npara_all * npara_all, sizeof(double));
    atpy = (double *) calloc ( npara_all,         sizeof(double));


    solvefor = npara_all;
    solvefor2 = npara_all * npara_all;
    MPI_Reduce(atpai, atpa, solvefor2, MPI_DOUBLE, MPI_SUM, root_process, MPI_COMM_WORLD);
    MPI_Reduce(atpyi, atpy, solvefor, MPI_DOUBLE, MPI_SUM, root_process, MPI_COMM_WORLD);

    printf("(%4d) Finish MPI_Reduce\n", myrank);
    fflush(stdout);


//    MPI_Barrier (MPI_COMM_WORLD);


    if (slvls == 4)
    {



    int nprow, npcol, context, ictxt, myrow, mycol, nb, lld, lld_distr, mp, nq, 
        info, desc_atpa[9], desc_atpy[9], desc_atpa_distr[9], desc_atpy_distr[9];
    double *atpa_distr, *atpy_distr;

//    Cblacs_pinfo(&myrank, &nprocs);

//    printf ("(%4d) in total processors of %4d\n", myrank, nprocs);


    nprow = (int)( sqrt((double)(nprocs) ) );
    npcol = nprocs / nprow;

//    nb = 3;
    nb = 64;

    n = solvefor;

    blacs_get_( &i_negone, &i_zero, &ictxt );
    blacs_gridinit_( &ictxt, "R", &nprow, &npcol );
    blacs_gridinfo_( &ictxt, &nprow, &npcol, &myrow, &mycol );

    printf ("(%4d) nprow = %4d npcol = %4d myrow = %4d mycol = %4d\n",
         myrank, nprow, npcol, myrow, mycol);


    mp = numroc_( &n, &nb, &myrow, &i_zero, &nprow );
    nq = numroc_( &n, &nb, &mycol, &i_zero, &npcol );
    atpa_distr = (double *)calloc( mp*nq, sizeof(double));
    atpy_distr = (double *)calloc( mp, sizeof(double));


    printf ("(%4d) mp = %4d nq = %4d \n",
         myrank, mp, nq);


    lld = max( numroc_( &n, &n, &myrow, &i_zero, &nprow ), 1 );
    lld_distr = max( mp, 1 );
    descinit_( desc_atpa, &n, &n, &n, &n, &i_zero, &i_zero, &ictxt, &lld, &info );
    descinit_( desc_atpa_distr, &n, &n, &nb, &nb, &i_zero, &i_zero, &ictxt, &lld_distr, &info );
    descinit_( desc_atpy, &n, &i_one, &n, &i_one, &i_zero, &i_zero, &ictxt, &lld, &info );
    descinit_( desc_atpy_distr, &n, &i_one, &nb, &nb, &i_zero, &i_zero, &ictxt, &lld_distr, &info );


    pdgeadd_( "N", &n, &n, &one, atpa, &i_one, &i_one, desc_atpa, &zero, atpa_distr, &i_one, &i_one, desc_atpa_distr );
    pdgeadd_( "N", &n, &i_one, &one, atpy, &i_one, &i_one, desc_atpy, &zero, atpy_distr, &i_one, &i_one, desc_atpy_distr );


/*
    for (i = 0; i < 9; i ++)
        printf ("(%4d) atpa[%d] = %4d atpa_distr[%d] = %4d atpy[%d] = %4d atpy_distr[%d] = %4d\n", 
            myrank, i, desc_atpa[i], i, desc_atpa_distr[i], i, desc_atpy[i], i, desc_atpy_distr[i] );

    for (i = 0; i < 9; i ++)
      printf ("(%4d) atpa[%d] = %f atpa_distr[%d] = %f atpy[%d] = %f atpy_distr[%d] = %f\n", 
            myrank, i, atpa[i], i, atpa_distr[i], i, atpy[i], i, atpy_distr[i] );

*/

//    pdgeqrf_(&m, &n, A_distr, &i_one, &i_one, descA_distr, tau, work, &lwork, &info);

/////////////////////////////////////////////////////////////////////////////////////////////////////////
    pdpotrf_("L", &n, atpa_distr, &i_one, &i_one, desc_atpa_distr, &info);

    pdpotrs_("L", &n , &i_one, atpa_distr, &i_one, &i_one, desc_atpa_distr, 
            atpy_distr, &i_one, &i_one, desc_atpy_distr, &info);

    pdpotri_("L", &n, atpa_distr, &i_one, &i_one, desc_atpa_distr, &info);



    pdgeadd_( "N", &n, &n, &one, atpa_distr, &i_one, &i_one, desc_atpa_distr, &zero, atpa, &i_one, &i_one, desc_atpa);
    pdgeadd_( "N", &n, &i_one, &one, atpy_distr, &i_one, &i_one, desc_atpy_distr, &zero, atpy, &i_one, &i_one, desc_atpy);


    blacs_barrier_( &ictxt, "A"); 


    }


if(myrank == root_process) 
{

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
    if (slvls == 3)
    {
        LAPACKE_dpotrf(LAPACK_ROW_MAJOR, 'U', solvefor, atpa, solvefor);
        LAPACKE_dpotrs(LAPACK_ROW_MAJOR, 'U', solvefor, 1, atpa, solvefor, atpy, 1);
        LAPACKE_dpotri(LAPACK_ROW_MAJOR, 'U', solvefor, atpa, solvefor);
        for (n = 0; n < solvefor; n++)
            ksi[n] = atpy[n];
    }

    if (slvls == 4)
    {
        for (n = 0; n < solvefor; n++)
            ksi[n] = atpy[n];
    }



    if (slvls == 1)
    {
        printf ("Cholesky Decomposition is used\n");
        solvels_chol(atpa, solvefor, atpy, ksi, 0);
    }
    if (slvls == 2)
    {
        printf ("Gauss Elimination is used\n");
        solvegaus(atpa, solvefor, atpy, ksi);
    }



    s4 = time(NULL);
    printf("(%4d) %5ld: seconds of inverting N\n", myrank, s4-s3);
    fflush(stdout);




    sigma0 = 1;

    
    sprintf (l2file, "%s.%d.NMA.L2", l1cfile, nmax);
  
//    strcat (l1cfile,".L2");

    printf ("OUTPUT: %s\n", l2file);

    if ( (fpout_coef = fopen (l2file,"w")) == NULL)
    {
        printf ("Cannot open fpout_coef file!\n");
        exit (0);
    }



    for (m = 0; m <= nmax; m ++)
    {
        l = nmax - m + 1;
        for (k = 0; k < l; k++)
        {
            if (m==0)
            {
                n = k + m;
                fprintf(fpout_coef, 
                    "%4d %4d %23.13e %23.13e %23.13e %23.13e\n", 
                    n, m, ksi[k], zero,
//                    sigma0 * sqrt(dksi[k * npara_all + k]), zero);
                    sigma0 * sqrt(atpa[k * npara_all + k]), zero);
//                    zero, zero);

            }
            else
            {
                ind = nmax + 1 + (2 * nmax - m + 2) * (m - 1);
                n = k + m;
                fprintf(fpout_coef, 
                    "%4d %4d %23.13e %23.13e %23.13e %23.13e\n", 
                    n, m, 
                    ksi[ind + n - m],
                    ksi[ind + n - m + l],
//                    sigma0 * sqrt(dksi[(ind + n - m) * solvefor + ind + n - m]),
//                    sigma0 * sqrt(dksi[(ind + n - m + l) * solvefor + ind + n - m + l]));
                    sigma0 * sqrt(atpa[(ind + n - m) * npara_all + ind + n - m]),
                    sigma0 * sqrt(atpa[(ind + n - m + l) * npara_all + ind + n - m + l]));
//                    zero, zero);

            }
        }
    }


    s5 = time(NULL);
    printf("(%4d) %5ld: seconds of output data\n", myrank,  s5-s4);
    fflush(stdout);



    free (ksi);
    free (dksi);

    free (atpa);
    free (atpy);
    fclose(fpout_coef);
   
    printf ("\nNormal end of ECHO!\n\npress any key to finish...\n");

}



//    blacs_gridexit_( &ictxt );
//    blacs_exit_( &i_zero );

      /* Close down this processes. */

    ierr = MPI_Finalize();


    exit(0);

}
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
 


double julian_date (short int year, short int month, short int day,
                    double hour)
{
   long int jd12h;

   double tjd;

   jd12h = (long) day - 32075L + 1461L * ((long) year + 4800L
      + ((long) month - 14L) / 12L) / 4L
      + 367L * ((long) month - 2L - ((long) month - 14L) / 12L * 12L)
      / 12L - 3L * (((long) year + 4900L + ((long) month - 14L) / 12L)
      / 100L) / 4L;
   tjd = (double) jd12h - 0.5 + hour / 24.0;

   return (tjd);
}


void cal_date (double tjd,
               short int *year, short int *month, short int *day,
               double *hour)
{
   long int jd, k, m, n;

   double djd;

   djd = tjd + 0.5;
   jd = (long int) djd;

   *hour = fmod (djd,1.0) * 24.0;

   k     = jd + 68569L;
   n     = 4L * k / 146097L;

   k     = k - (146097L * n + 3L) / 4L;
   m     = 4000L * (k + 1L) / 1461001L;
   k     = k - 1461L * m / 4L + 31L;

   *month = (short int) (80L * k / 2447L);
   *day   = (short int) (k - 2447L * (long int) *month / 80L);
   k      = (long int) *month / 11L;

   *month = (short int) ((long int) *month + 2L - 12L * k);
   *year  = (short int) (100L * (n - 49L) + m + k);

   return;
}



/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/*
*/
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
double grv_open (char *grv_name, char *label, 
        int nmax, int mmax, double *coef)
{
    FILE *fp_grv;
    double c,s, nd, md;
    int n=999,m=999, l, ind;
    char string[200], name[20];

    if ((fp_grv = fopen (grv_name,"r")) == NULL)
    {
        printf ("Cannot open gravity file?\n");
        exit (0);
    }

    memset (coef,0,(nmax+1)*(nmax+1));

//    coef[0] = 1;
    while (1)
    {
        if (fgets (string, 200, fp_grv) == NULL) break;

        if (strlen(label)==0)
        {
            sscanf (string, "%lf%lf%lf%lf", &nd, &md, &c, &s);	
            n = (int)nd;
            m = (int)md;
        }
        else 
        {
            sscanf (string, "%s", name);	
            if (strcmp (name, label) ==0)	
            {
                sscanf (string, "%*s%lf%lf%lf%lf", &nd, &md, &c, &s);	
                n = (int)nd;
                m = (int)md;
            }
        }

        if (n > nmax || m > mmax)
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
        }
    }

    fclose(fp_grv);
    return 0;
}











/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/*
*/
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
double cs2gp_pt (double *llr, double *cs, double gm, double a, int nmax, 
        double *gpt, double *pt)
{
    int n, m, k, l, ind;

    double sinf, cosf, *cosml, *sinml, *aprn, *pbar, lat, lon, r, 
        gpti, t, a2r;


/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/

    lat = llr[0];
    lon = llr[1];
    r = llr[2];

    sinf = sin(lat * DEG2RAD);
    cosf = cos(lat * DEG2RAD);

//    #pragma omp parallel private(cosml, sinml, aprn, pbar, pbar1, pbar2, n, m, k, l, ind, sinf, cosf)
    cosml = (double *) calloc ( nmax + 1, sizeof(double)); //cos(m*lamta)
    sinml = (double *) calloc ( nmax + 1, sizeof(double)); //sin(m*lamta)
    aprn = (double *) calloc ( nmax + 1, sizeof(double));  //sin(m*lamta)
    pbar = (double *) calloc ( nmax + 1, sizeof(double));

    cosml[0] = 1; sinml[0] = 0;
    cosml[1] = cos(lon * DEG2RAD); sinml[1] = sin(lon * DEG2RAD); 

    for (m = 2; m <= nmax; m++)
    {
//        cosml[m] = cos(m * lon * DEG2RAD);
//        sinml[m] = sin(m * lon * DEG2RAD);
        cosml[m] = 2.0 * cosml[1] * cosml[m-1] - cosml[m-2];
        sinml[m] = 2.0 * cosml[1] * sinml[m-1] - sinml[m-2];
    }

    aprn[0] = gm / r;
    a2r = a / r;

    for (n = 1; n <= nmax; n++)
    {
//        aprn[n] = pow (a / r, n) * gm / r;
        aprn[n] = aprn[n - 1] * a2r;
    }

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/

    t = sinf;
    for (m = 0; m <= nmax; m ++)
    {
        l = nmax - m + 1;
//        lgdr2(t, nmax, m, pbar, pbar1, pbar2);
        lgdr(t, nmax, m, pbar);
        for (k = 0; k < l; k++)
        {
            if (m==0)
            {
//                ind = 0;
                n = k + m;
                pt[k] = aprn[n] * pbar[k];
            }
            else
            {
                ind = nmax + 1 + (2 * nmax - m + 2) * (m - 1);
                n = k + m;
                pt[ind + n - m] = aprn[n] * pbar[k] * cosml[m];
                pt[ind + n - m + l] = aprn[n] * pbar[k] * sinml[m];
            }
        }
    }



//! ATPA
    gpti = 0;
    for(k = 0; k < (nmax + 1) * (nmax + 1); k++)
    {
        gpti = gpti + pt[k] * cs[k];
//        printf ("k = %d pt = %e\t cs = %e\t gp = %e \n", k, pt[k], cs[k], pt[k] * cs[k]);
//        pt[k] = pnmc[k];
    }


    *gpt = gpti; 


//    *gpt = gm/r + gm/r * a/r * a/r * (sqrt(5.0) * (1.5 * t * t - 0.5)) * cs[2];
    
//    printf ("cs[2] = %e\n", cs[2]);


    free (pbar);

    free (cosml);
    free (sinml);
    free (aprn);
    return 1;
}
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/




/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/

double lgdr(double t, int nmax, int m, double *pbar)

/*
! THIS CALCULATES THE FULLY NORMALIZED LEGENDRE FUNCTION WITH GIVEN ORDER(M),
! MAXIMUM DEGREE (NMAX), AND GIVEN EVALUATION POINT, T (COSINES OF COLATITUDE).
! THIS RETURNS ALL Pn,m, P'n,m, AND P''n,m (m=<n<=Nmax).
! THE RECURSION FORMULAR FOR THE FUNCTION ITSELF IS GIVEN IN JEKELI(1996).
! THE RECURSION FORMULAR FOR THE 1ST DERIVATIVE IS GIVEN IN TSCHERNING, ET AL(1983).
! THE FORMULAR FOR THE 2ND DERIVATIVE IS FROM THE ASSOCIATE LEGENDRE EQUATION.
! NOTE : EQUATIONS GIVEN IN TSCHERNING, ET AL(1983) HAVE ERRATA.
!
! S.C. Han, 1/24/01 (MODIFIED FOR CRAY T94 2/13/01)
!
*/
{
    int i;
//REAL*8 :: PBAR(NMAX-M+1),PBAR1(NMAX-M+1),PBAR2(NMAX-M+1),T,P00,P11,C,D
    double p00, p11, c, d;
//! THE FULLY NORMALIZED ASSOCIATED LEGENDRE FUNCTION
//! Pm,m : JEKEIL (A.3c) & (A.3d) , P'm,m : TSCHERNING (7)

    p00 = 1.0; 
    p11 = sqrt (3.0*(1.0-t*t));
    if (m>=1)
    {
        pbar[0] = p11; 

        for (i = 2; i <= m; i++)
        {
            pbar[0] = sqrt((2.0*i+1.0)/(2.0*i)*(1.0-t*t))*pbar[0];
        }
    }
    else 
    {
        pbar[0]=p00; 
    }

    if (nmax - m + 1 >= 2)
    {
        pbar[1] = sqrt(2.0*m +3.0) * t * pbar[0];
    }

    for(i = 3; i <= nmax-m+1; i++)
    {
        c=((2.0*m+2.0*i-3.0) * (2.0*m + 2.0*i-1.0)) / ((i-1.0)*(2.0*m+i-1.0));
        d=((2.0*m+2.0*i-1.0)*(2.0*m+i-2.0)*(i-2.0)) 
            / ((2.0*m+2.0*i-5.0)*(i-1.0)*(2.0*m+i-1.0));
        pbar[i-1] = sqrt(c)*t*pbar[i-2] - sqrt(d) * pbar[i-3];
    }

    return 0;
}





double zero0zero1 (double *cs, int nmax)
{
    int n, m, l, ind, ic, is;
    cs[0] = 0;
    cs[1] = 0;

    n = 1; m = 1;
    l = nmax - m + 1;
    ind = nmax + 1 + (2 * nmax - m + 2) * (m - 1);
    ic = ind + n - m;
    is = ind + n - m + l;
    cs[ic] = 0;
    cs[is] = 0;

    return 0;
}





/*
------------------------------------------------------------------------

    Purpose: some subroutines of matrix and vector operation
    Notes: 
    Programmer: Kun Shang @ 4.29.2014
    Functions:

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


------------------------------------------------------------------------
*/




/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/*
* brinv - 
* @param1: description of param1
* @param2: description of param2
*
* version: 20 Aug 2010
*/
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
int brinv (double *a,int n)
{ 
    int *is, *js, i, j, k, l, u, v;
    double d, p;

    is = (int *)malloc (n * sizeof(int));
    js = (int *)malloc (n * sizeof(int));
    for (k = 0; k <= n - 1; k++)
    { 
        d=0.0;

        for (i = k; i <= n - 1; i++)
        {
            for (j = k; j <= n - 1; j++)
            { 
                l = i * n + j; 
                p = fabs (a[l]);
                if (p > d) 
                { 
                    d = p; 
                    is[k] = i; 
                    js[k] = j;
                }
            }
        }

        if (d + 1.0 == 1.0)
        { 
//            rk = brank(a, n, n);
            printf ("warning: Matrix may be ill-conditioned\n");
//            printf ("rank = %d < dimension %d\n", rk, n);
//            free (is); 
//            free (js); 
//            exit(0);
        }

        if (is[k] != k)
        {
            for (j = 0; j <= n - 1; j++)
            { 
                u = k * n + j; 
                v = is[k] * n + j;
                p = a[u];
                a[u] = a[v];
                a[v] = p;
            }
        }

        if (js[k] != k)
        {
            for (i = 0; i <= n - 1; i++)
            { 
                u = i * n + k; 
                v = i * n + js[k];
                p = a[u];
                a[u] = a[v];
                a[v] = p;
            }
        }

        l = k * n + k;
        a[l] = 1.0 / a[l];
        
        for (j = 0; j <= n - 1; j++)
        {
            if (j != k)
            { 
                u = k * n + j; 
                a[u] = a[u] * a[l];
            }
        }
        
        for (i = 0; i <= n - 1; i++)
        {
            if (i != k)
                for (j = 0; j <= n - 1; j++)
                    if (j != k)
                    { 
                        u = i * n + j;
                        a[u] = a[u] - a[i * n + k] * a[k * n + j];
                    }
        }

        for (i = 0; i <= n - 1; i++)
        {
            if (i != k)
            { 
                u = i * n + k; 
                a[u] = - a[u] * a[l];
            }
        }
    }
    
    for (k = n - 1; k >= 0; k--)
    { 
        if (js[k] != k)
        {
            for (j = 0; j <= n - 1; j++)
            { 
                u = k * n + j; 
                v = js[k] * n + j;
                p = a[u]; 
                a[u] = a[v]; 
                a[v] = p;
            }
        }
        
        if (is[k] != k)
        {
            for (i = 0; i <= n - 1; i++)
            { 
                u = i * n + k; 
                v = i * n + is[k];
                p = a[u]; 
                a[u] = a[v]; 
                a[v] = p;
            }
        }
    }
    
    free(is); 
    free(js);
    return(1);
}







int brank(double *a, int m, int n)
  { int i,j,k,nn,is = 0,js = 0,l,ll,u,v;
    double q,d;
    nn=m;
    if (m>=n) nn=n;
    k=0;
    for (l=0; l<=nn-1; l++)
      { q=0.0;
        for (i=l; i<=m-1; i++)
        for (j=l; j<=n-1; j++)
          { ll=i*n+j; d=fabs(a[ll]);
        if (d>q) { q=d; is=i; js=j;}
          }
        if (q+1.0==1.0) return(k);
        k=k+1;
        if (is!=l)
          { for (j=l; j<=n-1; j++)
              { u=l*n+j; v=is*n+j;
                d=a[u]; a[u]=a[v]; a[v]=d;
              }
          }
        if (js!=l)
          { for (i=l; i<=m-1; i++)
              { u=i*n+js; v=i*n+l;
                d=a[u]; a[u]=a[v]; a[v]=d;
              }
          }
        ll=l*n+l;
        for (i=l+1; i<=n-1; i++)
          { d=a[i*n+l]/a[ll];
            for (j=l+1; j<=n-1; j++)
              { u=i*n+j;
                a[u]=a[u]-d*a[l*n+j];
              }
          }
      }
    return(k);
  }



void choldc(double *a, int n, double p[])
/*
 * Given a positive-definite symmetric matrix a[1..n][1..n], 
 * this routine constructs its Cholesky decomposition, A = L · LT . 
 * On input, only the upper triangle of a need be given; it is not modified. 
 * The Cholesky factor L is returned in the lower triangle of a, 
 * except for its diagonal elements which are returned in p[1..n].
 *
*/
{
    int i,j,k, rk;
    double sum;

    for (i=0;i<n;i++) 
    {
        for (j=i;j<n;j++) 
        {
            sum=a[i * n + j];
            for (k=i-1;k>=0;k--) 
                sum -= a[i * n + k]*a[j * n + k];
            if (i == j) 
            {
                if (sum <= 0.0)
                {
                    rk = brank(a, n, n);
                    printf("error: Matrix is not positive definite:\n");
                    printf("       i = %d\tj = %d\tsum = %e\n", i, j, sum);
                    printf("       rank = %d\n", rk);
                    exit(0);
                }
                p[i]=sqrt(sum);
            } 
            else a[j * n + i]=sum/p[i];
        }
    }


}


void cholsl(double *a, int n, double p[], double b[], double x[])
/*
 * Solves the set of n linear equations A · x = b, 
 * where a is a positive-definite symmetric matrix. 
 * a[1..n][1..n] and p[1..n] are input as the output of the routine choldc. 
 * Only the lower subdiagonal portion of a is accessed. 
 * b[1..n] is input as the right-hand side vector. 
 * The solution vector is returned in x[1..n]. 
 * a, n, and p are not modified and can be left in place for 
 * successive calls with different right-hand sides b. 
 * b is not modified unless you identify b and x in the calling sequence,
 * which is allowed.
 *
*/
{
    int i,k;
    double sum;

    for (i=0;i<n;i++) 
    {
        for (sum=b[i],k=i-1;k>=0;k--) 
            sum -= a[i * n + k]*x[k];
        x[i]=sum/p[i];
    }
    for (i=n-1;i>=0;i--) 
    {
        for (sum=x[i],k=i+1;k<n;k++) 
            sum -= a[k * n + i]*x[k];
        x[i]=sum/p[i];
    }


}






///////////********************////////////////
void solvels_chol(double *a, int n, double *y, double *x, int nocov)
{
    double *p, sum;
    int i, k, j;

    p = (double *)calloc (n, sizeof(double));

    choldc(a, n, p);
    cholsl(a, n, p, y, x);

    if (nocov == 1) 
    {
        free (p);
        return;
    }

    for (i=0;i<n;i++) 
    { 
        a[i * n + i]=1.0/p[i];
        for (j=i+1;j<n;j++) 
        {
            sum=0.0;
            for (k=i;k<j;k++) 
                sum -= a[j * n + k]*a[k * n + i]; 
            a[j * n + i]=sum/p[j];
        } 
    }
    
    for (i = 0; i <= n - 1; i++)
    {
        for (j = i; j <= n - 1; j++)
        {
            sum = 0.0;
            for (k = j; k <= n - 1; k++)
                sum = sum + a[k * n + i] * a[k * n + j];
            a[i * n + j] = sum;
        }
    }

    for (i = 0; i <= n - 1; i++)
    {
        for (j = 0; j <= i - 1; j++)
        {
            a[i * n + j] = a[j * n + i] ;
        }
    }

    free(p);
    return;
}    




void solvegaus(double *a, int n, double *y, double *x)
{
    brinv(a, n);
    brmul(a, y, n, n, 1, x);
}



double modvect (double *v)
{
    return  sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}


double dotvect (double *v1, double *v2)
{
    return  v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
}


void crsvect (double *v1, double *v2, double *v)
{
    v[0] = v1[1] * v2[2] - v1[2] * v2[1];
    v[1] = v1[2] * v2[0] - v1[0] * v2[2];
    v[2] = v1[0] * v2[1] - v1[1] * v2[0]; 
}






/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/*
* mt - 
* @param1: description of param1
* @param2: description of param2
*
* version: 20 Aug 2010
*/
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
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



