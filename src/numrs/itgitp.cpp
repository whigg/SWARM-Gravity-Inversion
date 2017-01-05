
/*
------------------------------------------------------------------------

    Purpose: 
    Notes: 
    Programmer: Kun Shang @ 4.29.2014
    Functions:


        double lgr_order (double *y, int dim_y, int dim_x, 
            double t, double *z, int order);

        double intydt(int num, double dt, double *y, double *val);

        double rkf78_auto (double h, double t, double *x, int dim, double err, 
              double (*fun)(double,double *,double *), int autoadjust)

        double rkf78 (double jd, double t, double h, double *x, int dim, 
              void (*fun)(int, double, double,double *,double *))

     Global variables:




------------------------------------------------------------------------
*/

#ifndef _NUMRS_H_

#include "numrs.h"
    #define _NUMRS_H_
#endif






double lgr_order (double *y, int dim_y, int dim_x, double t, double *z, int order)
{
    int i, j, k, m, dim;
    double s;

    if (order < 1)
    {
        printf ("error: order < 1 !\n");
        exit(0);
    }
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





double intydt(int num, double dt, double *y, double *val)
{
    int i;
    double cum;

    cum = 0;

    val[0] = 0;
    for (i = 1; i < num; i ++)
    {
        cum = cum + (y[i] + y[i-1]) * dt / 2.0;
        val[i] = cum;
    }

    return 0;
}












/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/****************************************************************************/
/*                                                                          */
/*		Functions for Runge-Kutta integrator                                */
/*                                                                          */
/*      Version:    2009-9-8                                                */
/*                                                                          */
/*      Copyright (c) 2009 shangkun@shao.ac.cn All Right Reserved           */
/*                                                                          */
/****************************************************************************/

/*
  Version: 2009-9-8 
  Version: 2009-9-13 integrate forwards & backwards
*/
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
double rkf78_auto (double h, double t, double *x, int dim, double err, 
              double (*fun)(double,double *,double *), int autoadjust)
/*
    purpose: auto-adjusted Runge-Kutta-Ful... integrator 
    input:  double h				integration step
            double t                integrate from t to t+h
            double *x               x(t)
            int dim	                dim(x)
            double err              tolerance of step control
            double (*fun)()			right(force) function
    output: double *x				x(t+h)
    return: h                       new step after adjustment
*/
{	
    int i, j, n, flag = 0;
    double *y, *k, *f, d = 0, tn;
    double a[13] = { 0, 2.0/27, 1.0/9, 1.0/6, 5.0/12, 1.0/2, 5.0/6, 1.0/6, 
        2.0/3, 1.0/3, 1.0, 0, 1.0 };
    double c[13] = { 0, 0, 0, 0, 0, 34.0/105, 9.0/35, 9.0/35, 9.0/280, 
        9.0/280, 0, 41.0/840, 41.0/840 };
    double b[13][12] = 
    {
        {0},
        {2.0/27},
        {1.0/36,1.0/12},
        {1.0/24,0,1.0/8},
        {5.0/12,0,-25.0/16,25.0/16},
        {1.0/20,0,0,1.0/4,1.0/5},
        {-25.0/108,0,0,125.0/108,-65.0/27,125.0/54},
        {31.0/300,0,0,0,61.0/225,-2.0/9,13.0/900},
        {2.0,0,0,-53.0/6,704.0/45,-107.0/9,67.0/90,3.0},
        {-91.0/108,0,0,23.0/108,-976.0/135,311.0/54,-19.0/60,17.0/6,-1.0/12},
        {2383.0/4100,0,0,-341.0/164,4496.0/1025,-301.0/82,2133.0/4100,
        45.0/82,45.0/164,18.0/41},
        {3.0/205,0,0,0,0,-6.0/41,-3.0/205,-3.0/41,3.0/41,6.0/41},
        {-1777.0/4100,0,0,-341.0/164,4496.0/1025,-289.0/82,2193.0/4100,
        51.0/82,33.0/164,12.0/41,0,1.0}
    };
    
    y = (double *) calloc (dim, sizeof(double));
    k = (double *) calloc (dim*13, sizeof(double));
    f = (double *) calloc (dim, sizeof(double));
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
    
    do
    {	
        for (i = 0; i <= 12; i++)
        {			
            tn = t + a[i] * h;
            for (n = 0; n <= dim - 1; n++)
            {	
                y[n] = x[n];
                for (j = 0; j <= i-1; j++)
                    y[n] = y[n] + h * b[i][j] * k[n*13+j];
            }
            fun (tn,y,f);
            for (n = 0; n <= dim - 1; n++)
            {
                k[n*13+i] = f[n];
            }
        }
        d = 0;
        for (n = 0; n <= dim - 1; n++)
        {
            d = d + fabs (41.0 / 840 * (k[n*13+0] + k[n*13+10] 
                - k[n*13+11] - k[n*13+12]) * h);
        }
		
        flag = 0;
        if (autoadjust == 1)
        {
            if (d > err) //adapting step h
            {
                h = h/2.0; 
                flag = 1; 	
            }
            if ( (d < err * 1e-4) && (h < 5e-3))
            {
                h = h*2.0; 
                flag = 2;
            }
        }
    }while (flag == 1);

    for (n = 0; n <= dim - 1; n++)
    {
        for (i = 0; i <= 12; i++)
            x[n] = x[n] + h * c[i] * k[n*13+i];	
    }
	
    free (y);
    free (f);
    free (k);
    return h;
}





/****************************************************************************/
/*                                                                          */
/*		Functions for Runge-Kutta integrator                                */
/*                                                                          */
/*      Version:    2009-9-8                                                */
/*                                                                          */
/*      Copyright (c) 2009 shangkun@shao.ac.cn All Right Reserved           */
/*                                                                          */
/****************************************************************************/

/*
  Version: 2009-9-8 
  Version: 2009-9-13 integrate forwards & backwards
*/
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
double rkf78 (double jd, double t, double h, double *x, int dim, 
              void (*fun)(int, double, double,double *,double *))
//              double (*fun)(double, double,double *,double *))
/*
    purpose: auto-adjusted Runge-Kutta-Ful... integrator 
    input:  double h				integration step
            double t                integrate from t to t+h
            double *x               x(t)
            int dim	                dim(x)
            double err              tolerance of step control
            double (*fun)()			right(force) function
    output: double *x				x(t+h)
    return: h                       new step after adjustment
*/
{	
    int i, j, n, flag = 0;
    double *y, *k, *f, d = 0, tn;
    static const double a[13] = { 0, 2.0/27, 1.0/9, 1.0/6, 5.0/12, 1.0/2, 5.0/6, 1.0/6, 
        2.0/3, 1.0/3, 1.0, 0, 1.0 };
    static const double c[13] = { 0, 0, 0, 0, 0, 34.0/105, 9.0/35, 9.0/35, 9.0/280, 
        9.0/280, 0, 41.0/840, 41.0/840 };
    static const double b[13][12] = 
    {
        {0},
        {2.0/27},
        {1.0/36,1.0/12},
        {1.0/24,0,1.0/8},
        {5.0/12,0,-25.0/16,25.0/16},
        {1.0/20,0,0,1.0/4,1.0/5},
        {-25.0/108,0,0,125.0/108,-65.0/27,125.0/54},
        {31.0/300,0,0,0,61.0/225,-2.0/9,13.0/900},
        {2.0,0,0,-53.0/6,704.0/45,-107.0/9,67.0/90,3.0},
        {-91.0/108,0,0,23.0/108,-976.0/135,311.0/54,-19.0/60,17.0/6,-1.0/12},
        {2383.0/4100,0,0,-341.0/164,4496.0/1025,-301.0/82,2133.0/4100,
        45.0/82,45.0/164,18.0/41},
        {3.0/205,0,0,0,0,-6.0/41,-3.0/205,-3.0/41,3.0/41,6.0/41},
        {-1777.0/4100,0,0,-341.0/164,4496.0/1025,-289.0/82,2193.0/4100,
        51.0/82,33.0/164,12.0/41,0,1.0}
    };
    
    y = (double *) calloc (dim, sizeof(double));
    k = (double *) calloc (dim*13, sizeof(double));
    f = (double *) calloc (dim, sizeof(double));
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
    
    do
    {	
        for (i = 0; i <= 12; i++)
        {			
            tn = t + a[i] * h;
            for (n = 0; n <= dim - 1; n++)
            {	
                y[n] = x[n];
                for (j = 0; j <= i-1; j++)
                    y[n] = y[n] + h * b[i][j] * k[n*13+j];
            }
            fun (dim, jd, tn, y, f);
//            fun (jd, tn, y, f);
            for (n = 0; n <= dim - 1; n++)
            {
                k[n*13+i] = f[n];
            }
        }
        d = 0;
        for (n = 0; n <= dim - 1; n++)
        {
            d = d + fabs (41.0 / 840 * (k[n*13+0] + k[n*13+10] 
                - k[n*13+11] - k[n*13+12]) * h);
        }
		
        flag = 0;
    }while (flag == 1);

    for (n = 0; n <= dim - 1; n++)
    {
        for (i = 0; i <= 12; i++)
            x[n] = x[n] + h * c[i] * k[n*13+i];	
    }
	
    free (y);
    free (f);
    free (k);
    return h;
}





