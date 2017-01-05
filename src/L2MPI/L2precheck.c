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


/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
int main (int argc, char *argv[])
{

    double int_lim, int_lim1, int_lim2, mem_lim, mem_lim1, mem_lim2;
    int nmax, nprocs, nodes, ndata_all, ndata, npara_all, npara,
        n, m, n_1m, m_120;
    double month, cpx_1m_120, cpx, cpx_d2n, cpx_inv, cpx_i, cput_i, cput, cput_1m_120, wallt;

    if (argc != 4)
    {
        printf ("usage: L2precheck nmax month_data nprocs\n");
        exit(0);
    }
    else
    {
        sscanf (argv[1], "%d", &nmax);
        sscanf (argv[2], "%lf", &month);
//        sscanf (argv[3], "%d", &nodes);
        sscanf (argv[3], "%d", &nprocs);
    }

//    printf ("nmax = %d month = %d nodes = %d\n", nmax, month, nodes);
    printf ("nmax = %d month = %.2f nodes = %d\n", nmax, month, nprocs);

    printf ("\n!!!!!!!!!!!!!! Design Matrix Accumulation !!!!!!!!!!!!!!\n");


//    nprocs = nodes * 12;
    
    ndata_all = (int)(17280 * 31 * month);
    ndata = (int)(ndata_all / nprocs) + 1;

    npara_all = (nmax + 1) * (nmax + 1);
    npara = (int)(npara_all / nprocs) + 1;


/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
// index and memeory check

    int_lim = (double)npara_all * (double)ndata /1024.0/1024.0/1024.0;
    mem_lim = ((double)npara_all * (double)(ndata + npara) 
                + (double) ndata_all * 10) / 1024.0/1024.0/1024.0 * 8;

    if (int_lim > 2) 
        printf ("int_lim ~= %9.3fG > INT_LIM (2G) : warning!!!\n", int_lim);
    else 
        printf ("int_lim ~= %9.3fG < INT_LIM (2G) : safe~~~\n", int_lim);

    if (mem_lim  > 4) 
        printf ("mem_lim ~= %9.3fG > MEM_LIM (4G) : warning!!!\n", mem_lim);
    else
        printf ("mem_lim ~= %9.3fG < MEM_LIM (4G) : safe~~~\n", mem_lim);


    n = ndata_all;

    m = npara_all;


    cpx_d2n = (2.0 * n - 1.0) * (m + 1.0) * m / 2.0;

    cpx_inv = (double)m * (double)m * (double)m / 3.0; 

    cpx = cpx_d2n + cpx_inv;

    n_1m = 17280 * 31;
    m_120 = 121 * 121;
    cpx_1m_120 = (2.0 * n_1m - 1.0) * (m_120 + 1.0) * m_120 / 2.0 + m_120 * m_120 * m_120 / 3.0;
    cput_1m_120 = 7.0 * 3600.0;
    
    cpx_i  = cpx_1m_120;
    cput_i = cput_1m_120; 
 
   
    cput = cpx / cpx_i * cput_i;

    wallt = cput / nprocs;

    
    printf ("wallt ~= %f min cput ~= %f hour cpx_d2n : cpx_inv ~= %f\n", 
            wallt / 60.0, cput / 3600.0, cpx_d2n / cpx_inv);
   

    printf ("\n!!!!!!!!!!!!!! Normal Matrix Accumulation !!!!!!!!!!!!!!\n");


//    nprocs = nodes * 12;
    
    ndata_all = (int)(17280 * 31 * month);
    ndata = (int)(ndata_all / nprocs) + 1;

    npara_all = (nmax + 1) * (nmax + 1);
    npara = (int)(npara_all / nprocs) + 1;


/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
// index and memeory check

    int_lim1 = (double)npara_all * (double)ndata /1024.0/1024.0/1024.0;
    int_lim2 = (double)npara_all * (double)npara_all /1024.0/1024.0/1024.0;
    
    mem_lim1 = (double)npara_all * (double) ndata / 1024.0/1024.0/1024.0 * 8;
    mem_lim2 = (double)npara_all * (double) npara_all / 1024.0/1024.0/1024.0 * 8;
    
    mem_lim = mem_lim1 + mem_lim2;


    if (int_lim1 > 2) 
        printf ("int_data ~= %9.3fG > INT_LIM (2G) : warning!!!\n", int_lim1);
    else 
        printf ("int_data ~= %9.3fG < INT_LIM (2G) : safe~~~\n", int_lim1);
    if (int_lim2 > 2) 
        printf ("int_para ~= %9.3fG > INT_LIM (2G) : warning!!!\n", int_lim2);
    else 
        printf ("int_para ~= %9.3fG < INT_LIM (2G) : safe~~~\n", int_lim2);


    if (mem_lim1  > 4) 
        printf ("mem_data ~= %9.3fG > MEM_LIM (4G) : warning!!!\n", mem_lim1);
    else
        printf ("mem_data ~= %9.3fG < MEM_LIM (4G) : safe~~~\n", mem_lim1);


    if (mem_lim2  > 4) 
        printf ("mem_para ~= %9.3fG > MEM_LIM (4G) : warning!!!\n", mem_lim2);
    else
        printf ("mem_para ~= %9.3fG < MEM_LIM (4G) : safe~~~\n", mem_lim2);


    if (mem_lim  > 4) 
        printf ("mem_all  ~= %9.3fG > MEM_LIM (4G) : warning!!!\n", mem_lim);
    else
        printf ("mem_all  ~= %9.3fG < MEM_LIM (4G) : safe~~~\n", mem_lim);


/*
    if (mem_lim  > 4) 
    {
        printf ("mem_lim ~= %fG > MEM_LIM (4G) : warning!!!\n", mem_lim);
        printf ("mem_data ~= %fG mem_para ~= %fG\n", mem_lim1, mem_lim2);
    }
    else
    {
        printf ("mem_lim ~= %fG < MEM_LIM (4G) : safe~~~\n", mem_lim);
        printf ("mem_data ~= %fG mem_para ~= %fG\n", mem_lim1, mem_lim2);
    }
*/
    n = ndata_all;

    m = npara_all;


    cpx_d2n = (2.0 * n - 1.0) * (m + 1.0) * m / 2.0;

    cpx_inv = (double)m * (double)m * (double)m / 3.0; 

    cpx = cpx_d2n + cpx_inv;

    n_1m = 17280 * 31;
    m_120 = 121 * 121;
    cpx_1m_120 = (2.0 * n_1m - 1.0) * (m_120 + 1.0) * m_120 / 2.0 + m_120 * m_120 * m_120 / 3.0;
    cput_1m_120 = 7.0 * 3600.0;
    
    cpx_i  = cpx_1m_120;
    cput_i = cput_1m_120; 
 
   
    cput = cpx / cpx_i * cput_i;

    wallt = cput / nprocs;

    
    printf ("wallt ~= %f min cput ~= %f hour cpx_d2n : cpx_inv ~= %f\n", wallt / 60.0, cput / 3600.0, cpx_d2n / cpx_inv);
   




}
