/*
------------------------------------------------------------------------

    Purpose: some subroutines of matrix and vector operation
    Notes: 
    Programmer: Kun Shang @ 4.29.2014
    Functions:
 
        void mgrns(double u, double g, double *r,int n,double a[]);

        double mgrn1(double u, double g, double *r);

        double mean (double * array, double N);

        double std_dev (double * array, double N);




*/



#ifndef _NUMRS_H_
    #include "numrs.h"
#endif




double mean (double * array, double N)
{
    int i;
    double sum = 0 ;
    for (i = 0; i < N; i++)
        sum = sum + array [i];
    return sum/N;
} // function calculating mean


double std_dev (double * array, double N)
{
    int i;
    double sum = 0;
    double STD_DEV = 0; // returning zero's

    for (i = 0; i < N; i++)
    {
        sum = sum + array [i];
        STD_DEV = STD_DEV + pow(array [i], 2);
    }
    return sqrt ((STD_DEV/N) - (pow(sum/N,2)));
} // function calculating standard deviation




  void mgrns(double u, double g, double *r,int n,double a[])
//  int n;
//  double u,g,*r,a[];
  { int i,k,m;
    double s,w,v,t;
    s=65536.0; w=2053.0; v=13849.0;
    for (k=0; k<=n-1; k++)
      { t=0.0;
        for (i=1; i<=12; i++)
          { *r=(*r)*w+v; m=(int)(*r/s);
            *r=*r-m*s; t=t+(*r)/s;
          }
        a[k]=u+g*(t-6.0);
      }
    return;
  }


  double mgrn1(double u, double g, double *r)
//  double u,g,*r;
  { int i,m;
    double s,w,v,t;
    s=65536.0; w=2053.0; v=13849.0;
    t=0.0;
    for (i=1; i<=12; i++)
      { *r=(*r)*w+v; m=(int)(*r/s);
        *r=*r-m*s; t=t+(*r)/s;
      }
    t=u+g*(t-6.0);
    return(t);
  }




