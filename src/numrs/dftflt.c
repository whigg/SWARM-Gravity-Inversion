/*
------------------------------------------------------------------------

    Purpose: some subroutines of matrix and vector operation
    Notes: 
    Programmer: Kun Shang @ 4.29.2014
    Functions:

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



*/


/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/

#ifndef _NUMRS_H_
    #include "numrs.h"
#endif









/*                               -*- Mode: C -*- 
 * Filename: dft.c
 * Copyright (C) Dan Noland 2003
 * Author: Dan Noland
 * Created: Thu Mar 11 03:19:19 2004
 *           By: Dan Noland
 * Last-Updated: Thu Mar 11 03:51:54 2004
 *     Update #: 12
 * Status: 
 */


#define TRUE 1
#define FALSE 0

// Take, eat, this is my code which is hacked for you...
// (In case you didn't grasp that one consider this public domain)
// nolandda Thu Mar 11 03:49:18 EST 2004

/*
  Discrete Fourier Transform
*/

int dft(long int length, double real_sample[], double imag_sample[])
{
  long int i, j;
  double arg;
  double cosarg,sinarg;
  double *temp_real=NULL,*temp_imag=NULL;

  temp_real = calloc(length, sizeof(double));
  temp_imag = calloc(length, sizeof(double));
  if (temp_real == NULL || temp_imag == NULL)
  {
    return(FALSE);
  }

  for(i=0; i<length; i+=1) 
  {
    temp_real[i] = 0;
    temp_imag[i] = 0;
    arg = -1.0 * TWOPI * (double)i / (double)length;
    for(j=0; j<length; j+=1) 
    {
      cosarg = cos(j * arg);
      sinarg = sin(j * arg);
      temp_real[i] += (real_sample[j] * cosarg - imag_sample[j] * sinarg);
      temp_imag[i] += (real_sample[j] * sinarg + imag_sample[j] * cosarg);
    }
  }

  /* Copy the data back */
  for (i=0; i<length; i+=1) 
  {
    real_sample[i] = temp_real[i];
    imag_sample[i] = temp_imag[i];
  }

  free(temp_real);
  free(temp_imag);
  return(TRUE);
}


/*
  Inverse Discrete Fourier Transform
*/

int inverse_dft(long int length, double real_sample[], double imag_sample[])
{
  long int i, j;
  double arg;
  double cosarg,sinarg;
  double *temp_real=NULL,*temp_imag=NULL;

  temp_real = calloc(length, sizeof(double));
  temp_imag = calloc(length, sizeof(double));
  if (temp_real == NULL || temp_imag == NULL)
  {
    return(FALSE);
  }

  for(i=0; i<length; i+=1) 
  {
    temp_real[i] = 0;
    temp_imag[i] = 0;
    arg = TWOPI * (double)i / (double)length;
    for(j=0; j<length; j+=1) 
    {
      cosarg = cos(j * arg);
      sinarg = sin(j * arg);
      temp_real[i] += (real_sample[j] * cosarg - imag_sample[j] * sinarg);
      temp_imag[i] += (real_sample[j] * sinarg + imag_sample[j] * cosarg);
    }
  }

  /* Copy the data back */
  for (i=0; i<length; i+=1) 
  {
    real_sample[i] = temp_real[i] / (double)length;
    imag_sample[i] = temp_imag[i] / (double)length;
  }

  free(temp_real);
  free(temp_imag);
  return(TRUE);
}



double daydft (double *obs, double *fft, int num,  double dt, double tc)
{
    double *xr, *xi, *yr, *yi, *win, fc, df, kc, th, fh, kh;
    int i, n, m, k;

    n = 1; k = 0;

    while(n<num)
    {
        n = n * 2;
        k = k + 1;
    }

//    n = num;

    xr = (double *) calloc (n, sizeof(double));
    xi = (double *) calloc (n, sizeof(double));
    win = (double *) calloc (n, sizeof(double));
    yr = (double *) calloc (n, sizeof(double));
    yi = (double *) calloc (n, sizeof(double));


    for (i = 0; i < n; i++)
    {
        if (i>=num)
        {
            xr[i] = 0;
            xi[i] = 0;
        }
        else
        {
            xr[i] = obs[i];
            xi[i] = 0;
        }
    }



//    dft(n, xr, xi);

    kkfft(xr,xi,n,k,yr,yi,0,1); 


/*       
    for (i = 0; i < num; i++)
    {
        if (i < 40 || (num - i) < 40)
        {
            xr[i] = 0;
            xi[i] = 0;
        }
    }
*/

    df = 1.0 / n / dt;

     
    fc = 1.0 / tc;
    kc = fc / df;

    th = 60; // seconds
    fh = 1.0 / th;
    kh = fh / df;

//    printf("fc = %f\t df = %f\t dt = %f\t kc = %f\t kh = %f\n", fc, df, dt, kc, kh);


//    kc = 120;

    m = 2;

    getwinhp(n, win, kc, m);

//    getwinbp(n, win, kc, kh, m);


    for (i = 0; i < n; i++)
    {
//        xr[i] = xr[i] * win[i];        xi[i] = xi[i] * win[i];
        yr[i] = yr[i] * win[i];        yi[i] = yi[i] * win[i];
    }

//    inverse_dft(n, xr, xi);

  

    kkfft(yr,yi,n,k,xr,xi,1,1); 


    for (i = 0; i < num; i++)
    {           
        fft[i] = xr[i];
//        fft[i] = sqrt(xr[i] * xr[i] + xi[i] * xi[i]);

    }

    free(xr);
    free(yr);
    free(xi);
    free(yi);
    free(win);


    return 0;
}



double getwinhp(int n, double *win, double kc, int m)
{
    int i;

    for (i = 0; i < n; i++)
    {
        if (i < n / 2 )
        {
            win[i] = 1 - 1/sqrt(1 + pow((i/kc), 2*m));
        }
        else if (i >= n / 2 )
        {
            win[i] = 1 - 1/sqrt(1 + pow(((n - i)/kc), 2*m));
        }

    }

    return 0;
}


double getwinbp(int n, double *win, double kc, double kh, int m)
{
    int i;

    for (i = 0; i < n; i++)
    {
        if (i < n / 2 )
        {
            win[i] = 1 - 1/sqrt(1 + pow((i/kc), 2*m));
        }
        else if (i >= n / 2 )
        {
            win[i] = 1 - 1/sqrt(1 + pow(((n - i)/kc), 2*m));
        }

    }

    for (i = 0; i < n; i++)
    {
        if (i < n / 2 )
        {
            win[i] = win[i] * 1/sqrt(1 + pow((i/kh), 2*m));
        }
        else if (i >= n / 2 )
        {
            win[i] = win[i] * 1/sqrt(1 + pow(((n - i)/kh), 2*m));
        }

    }


    return 0;
}

double getwinbp0(int n, double *win, double kc, int m)
{
    int i;

    for (i = 0; i < n; i++)
    {
        if (i < n / 4 )
        {
            win[i] = 1 - 1/sqrt(1 + pow((i/kc), 2*m));
        }
        else if (i >= n / 4 && i < n / 2)
        {
            win[i] = 1 - 1/sqrt(1 + pow(((n/2 - i)/kc), 2*m));
        }
        else if (i >= n / 2 && i < 3 * n / 4)
        {
            win[i] = 1 - 1/sqrt(1 + pow((( i - n/2 )/kc), 2*m));
        }
        else if (i >= 3 * n / 4)
        {
            win[i] = 1 - 1/sqrt(1 + pow((( n - i)/kc), 2*m));
        }

    }

    return 0;
}



double dayfft (double *data, double *fft, int num,  double fs)
/*
    ÓÃÍ¾£º  ¶Ôk5Êý¾Ý×öFFT£¬Êä³öFFT×î´óÊ±µÄÆµÂÊÖµ
    ÊäÈë£º  double data[]           Ô­Ê¼Êý¾Ý
            int res                 ·Ö±æÂÊ£¨Hz£©½øÐÐFFTµÄÊý¾Ý³¤¶È£¨¿ÉÒÔ´óÓÚÔ­Ê¼Êý¾Ý³¤¶È£¬²¹Áã£©
            int debug               µ÷ÊÔ¿ª¹Ø£¨=1Êä³öFFTµÄ½á¹ûÎÄ¼þ£¬=0²»Êä³ö£©£¬
            double fln              Êä³öÏÔÊ¾µÄµÍÆµ
            double fhn              Êä³öÏÔÊ¾µÄµÍÆµ
            double fs               ²ÉÑùÂÊ
            int loop                debugÄ£Ê½ÏÂÊä³öµÄÎÄ¼þÃû
    Êä³ö£º  double *f0              FFT×î´óÊ±µÄÆµÂÊÖµ
    ·µ»Ø£º  =0                      Õý³£
*/

{
    int i,  label = 0; 
    double *xr, *xi, *yr, *yi, fn, fln, fhn, df,dfn;
    int n = 1, k = 0;   //n ×ÊÁÏ³¤¶È
    double fl = 0;
    double fh = 0;
        
    while(n<num)
    {
        n = n * 2;
        k = k + 1;
    }
//  n = n / 2;
//  k = k - 1;

    xr = (double *) calloc (n, sizeof(double));
    xi = (double *) calloc (n, sizeof(double));

    

    for (i = 0; i < n; i++)
    {
        if (i>=num)
        {
            label = 1;
            xr[i] = 0;
            xi[i] = 0;
        }
        else
        {
            xr[i] = data[i];
            xi[i] = 0;
        }
    }

//  µÍÆµ£¬¸ßÆµ£¬
    fln = fl/fs;
    fhn = fh/fs;

    df = fs/n;
    dfn = df/fs;
    fn = (fhn - fln) / dfn;


    yr = (double *) calloc (n, sizeof(double));
    yi = (double *) calloc (n, sizeof(double));



    kkfft(xr,xi,n,k,yr,yi,0,1);     //xr,xiÊäÈëµÄÊµ²¿ºÍÐé²¿£¬·µ»ØÄ£ºÍ·ù½Ç£¨¶È£©, yr,yi·µ»Ø±ä»»ºóµÄÊµ²¿ºÍÐé²¿, 
                                    //n¸öµã, n=2^k, 0±íÊ¾¸µÁ¢Ò¶±ä»»£¨1ÎªÄæ±ä»»£©, 1±íÊ¾¼ÆËãÄ£ºÍ·ù½Ç£¨¶È£©
/*

    max = 0;
    for (f = fl; f < fh; f = f + df)
    {
        i = (int)(f/df);
        
        if (xr[i]>max)
        {
            max = xr[i];
            *f0 = f;
        }
    if (debug == 1)
        fprintf (fpout, "%e\t%e\n", f, xr[i]/sqrt(n));
    }
*/


    for (i = 0; i < n; i++)
    {
        if (i < 600 || (n - i) < 600)
        {
            yr[i] = 0;
            yi[i] = 0;
        }
    }

    kkfft(yr,yi,n,k,xr,xi,1,1);     //xr,xiÊäÈëµÄÊµ²¿ºÍÐé²¿£¬·µ»ØÄ£ºÍ·ù½Ç£¨¶È£©, yr,yi·µ»Ø±ä»»ºóµÄÊµ²¿ºÍÐé²¿, 
                                    //n¸öµã, n=2^k, 0±íÊ¾¸µÁ¢Ò¶±ä»»£¨1ÎªÄæ±ä»»£©, 1±íÊ¾¼ÆËãÄ£ºÍ·ù½Ç£¨¶È£©


    for (i = 0; i < num; i++)
    {           
        fft[i] = xr[i];
    }


    free(xr);
    free(yr);
    free(xi);
    free(yi);
    

    return 0;
}




void kkfft(double pr[], double pi[], int n, int k, double fr[], 
           double fi[], int l, int il)
  { int it,m,is,i,j,nv,l0;
    double p,q,s,vr,vi,poddr,poddi;
    for (it=0; it<=n-1; it++)
      { m=it; is=0;
        for (i=0; i<=k-1; i++)
          { j=m/2; is=2*is+(m-2*j); m=j;}
        fr[it]=pr[is]; fi[it]=pi[is];
      }
    pr[0]=1.0; pi[0]=0.0;
    p=6.283185306/(1.0*n);
    pr[1]=cos(p); pi[1]=-sin(p);
    if (l!=0) pi[1]=-pi[1];
    for (i=2; i<=n-1; i++)
      { p=pr[i-1]*pr[1]; q=pi[i-1]*pi[1];
        s=(pr[i-1]+pi[i-1])*(pr[1]+pi[1]);
        pr[i]=p-q; pi[i]=s-p-q;
      }
    for (it=0; it<=n-2; it=it+2)
      { vr=fr[it]; vi=fi[it];
        fr[it]=vr+fr[it+1]; fi[it]=vi+fi[it+1];
        fr[it+1]=vr-fr[it+1]; fi[it+1]=vi-fi[it+1];
      }
    m=n/2; nv=2;
    for (l0=k-2; l0>=0; l0--)
      { m=m/2; nv=2*nv;
        for (it=0; it<=(m-1)*nv; it=it+nv)
          for (j=0; j<=(nv/2)-1; j++)
            { p=pr[m*j]*fr[it+j+nv/2];
              q=pi[m*j]*fi[it+j+nv/2];
              s=pr[m*j]+pi[m*j];
              s=s*(fr[it+j+nv/2]+fi[it+j+nv/2]);
              poddr=p-q; poddi=s-p-q;
              fr[it+j+nv/2]=fr[it+j]-poddr;
              fi[it+j+nv/2]=fi[it+j]-poddi;
              fr[it+j]=fr[it+j]+poddr;
              fi[it+j]=fi[it+j]+poddi;
            }
      }
    if (l!=0)
      for (i=0; i<=n-1; i++)
        { fr[i]=fr[i]/(1.0*n);
          fi[i]=fi[i]/(1.0*n);
        }
    if (il!=0)
      for (i=0; i<=n-1; i++)
        { pr[i]=sqrt(fr[i]*fr[i]+fi[i]*fi[i]);
          if (fabs(fr[i])<0.000001*fabs(fi[i]))
            { if ((fi[i]*fr[i])>0) pi[i]=90.0;
              else pi[i]=-90.0;
            }
          else
            pi[i]=atan(fi[i]/fr[i])*360.0/6.283185306;
        }
    return;
  }






double resfft(double *res, double *resflt, int num)
{
    int n, i, isign;
    double *data;

    n = 1;
    while(n<num)
    {
        n = n * 2;
    }

    data = (double *) malloc (n * sizeof(double));

    for (i = 0; i < n; i++)
    {
        if (i>num)
        {
            data[i] = 0;
        }
        else
        {
            data[i] = res[i];
        }
    }
//  if (label == 1)
//      printf("²¹Áã!!!\n");

    isign = 1;

    realft(data,n,isign);

    for (i = 0; i < num; i++)
    {           
        resflt[i] = data[i];
    }


    return 0;

}

void realft(double data[], int n, int isign)
{
    int i,i1,i2,i3,i4,n2p3;
    double c1=0.5,c2,h1r,h1i,h2r,h2i;
    double wr,wi,wpr,wpi,wtemp,theta;

    theta=3.141592653589793/(double) n;
    if (isign == 1) {
        c2 = -0.5;
        four1(data,n,1);
    } else {
        c2=0.5;
        theta = -theta;
    }
    wtemp=sin(0.5*theta);
    wpr = -2.0*wtemp*wtemp;
    wpi=sin(theta);
    wr=1.0+wpr;
    wi=wpi;
    n2p3=2*n+3;
    for (i=2;i<=n/2;i++) {
        i4=1+(i3=n2p3-(i2=1+(i1=i+i-1)));
        h1r=c1*(data[i1]+data[i3]);
        h1i=c1*(data[i2]-data[i4]);
        h2r = -c2*(data[i2]+data[i4]);
        h2i=c2*(data[i1]-data[i3]);
        data[i1]=h1r+wr*h2r-wi*h2i;
        data[i2]=h1i+wr*h2i+wi*h2r;
        data[i3]=h1r-wr*h2r+wi*h2i;
        data[i4] = -h1i+wr*h2i+wi*h2r;
        wr=(wtemp=wr)*wpr-wi*wpi+wr;
        wi=wi*wpr+wtemp*wpi+wi;
    }
    if (isign == 1) {
        data[1] = (h1r=data[1])+data[2];
        data[2] = h1r-data[2];
    } else {
        data[1]=c1*((h1r=data[1])+data[2]);
        data[2]=c1*(h1r-data[2]);
        four1(data,n,-1);
    }
}


#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr

void four1(double data[], int nn, int isign)
{
    int n,mmax,m,j,istep,i;
    double wtemp,wr,wpr,wpi,wi,theta;
    double tempr,tempi;

    n=nn << 1;
    j=1;
    for (i=1;i<n;i+=2) {
        if (j > i) {
            SWAP(data[j],data[i]);
            SWAP(data[j+1],data[i+1]);
        }
        m=n >> 1;
        while (m >= 2 && j > m) {
            j -= m;
            m >>= 1;
        }
        j += m;
    }
    mmax=2;
    while (n > mmax) {
        istep=2*mmax;
        theta=6.28318530717959/(isign*mmax);
        wtemp=sin(0.5*theta);
        wpr = -2.0*wtemp*wtemp;
        wpi=sin(theta);
        wr=1.0;
        wi=0.0;
        for (m=1;m<mmax;m+=2) {
            for (i=m;i<=n;i+=istep) {
                j=i+mmax;
                tempr=wr*data[j]-wi*data[j+1];
                tempi=wr*data[j+1]+wi*data[j];
                data[j]=data[i]-tempr;
                data[j+1]=data[i+1]-tempi;
                data[i] += tempr;
                data[i+1] += tempi;
            }
            wr=(wtemp=wr)*wpr-wi*wpi+wr;
            wi=wi*wpr+wtemp*wpi+wi;
        }
        mmax=istep;
    }
}

#undef SWAP


/***************************************************************************************/
/*                                                                                     */
/*      Functions for k5vssp(32) filter                                                */
/*                                                                                     */
/*      Version:    2008-7-3                                                          */
/*                                                                                     */
/*      Copyright (c) 2008 shangkun@shao.ac.cn All Right Reserved                      */
/*                                                                                     */
/***************************************************************************************/

/*

  Version:    2008-6-17 
  Version:    2008-6-18     ÐÞ¸Äy[]ÀÛ¼ÓÎ´ÇåÁãµÄÖÂÃü´íÎó 
  Version:    2008-7-3      Ôö¼ÓÊäÈë²ÎÊýorderÂË²¨½×ÊýÐèÍâ²¿Ö¸¶¨ 

*/

double k5filter (double x[], double y[], int size, double fln, 
                 double fhn, double fs, int order)
/*
    ÓÃÍ¾£º  ¶Ôk5Ò»ÃëÊý¾ÝÂË²¨
    ÊäÈë£º  double x[]              Ô­Ê¼Ò»ÃëÊý¾Ý
            int size                Êý¾Ý³¤¶È
            double fln              µÍÆµ
            double fhn              ¸ßÆµ
            double fs                   ²ÉÑùÂÊ£¨Èç¹ûÊÇÒ»ÃëÊý¾Ý=size£©
    Êä³ö£º  double y[]              ÂË²¨ºóµÄÊý¾Ý£¬Ã»ÓÐ³Ë·Å´óÒò×Ó
    ·µ»Ø£º  =0                      Õý³£
*/
{
    double h[1000] = {0}, fl, fh;
    int n,k;

//  h = (double *) calloc ( order , sizeof(double));
    fl = fln/fs;
    fh = fhn/fs;

    firwin(order,fl,fh,h);
    for(n = 0 ; n < order/2 ; n++)
    {
        y[n] = 0;
        for(k = 0;k<order/2+n;k++)
        {
            y[n] = y[n] + x[n-k+order/2] * h[k];
        }
    }
    for(n = order/2 ; n < size - order/2 ; n++)
    {
        y[n] = 0;
        for(k=0;k<order;k++)
        {
            y[n] = y[n] + x[n-k+order/2] * h[k];
        }
    }
    for(n = size - order/2 ; n < size  ; n++)
    {
        y[n] = 0;
        for(k=order/2-size+n;k<order;k++)
        {
            y[n] = y[n] + x[n-k+order/2] * h[k];
        }
    }
    return 0;
}



void firwin(int n, double fln, double fhn,double h[])
{
    int i,n2,mid;
    double s,wc1,wc2,delay;
    double pis = 3.141592653589;

    if(0 == (n%2))
    {
        n2  = n/2-1;
        mid = 1;
    }
    else
    {
        n2=n/2;
        mid=0;
    }
    
    delay=n/2.0;
    
    wc1=2.0*pis*fln;
    wc2=2.0*pis*fhn;

    for(i=0;i<=n2;i++)
    {
        s=i-delay;
        h[i]=(sin(wc2*s)-sin(wc1*s))/(pis*s);
        h[i]=h[i]*hanningwin(n,i);
        h[n-i]=h[i];
    }
    
    if(mid==1)
        h[n/2]=(wc2-wc1)/pis;
}


double hanningwin(int n,int i)
{
    double w;
    double pis = 3.141592653589;
    w=0.5*(1.0-cos(2*i*pis/(n-1)));
    return(w);
}


double kaiserwin(int n,int i)
{
    double w;
    double pis = 3.141592653589;
    w=0.5*(1.0-cos(2*i*pis/(n-1)));
    return(w);
}








double lsfl1c1d (double *xmax, double *x, double *y, int nmax, int n, 
    double tpsec, double *ymaxflt, int nply, int ncpr, int nplycpr)
{

    int i, j, m, err, k, cpr;

    double T, *a, *f, *fmax, *ft, *fpf, *fpy, *yflt, *res, stdres;

    k = nply + ncpr * nplycpr;

    a = (double *) calloc ( (int)(k) , sizeof(double));
    f = (double *) calloc ( (int)(n*k) , sizeof(double));
    yflt = (double *) calloc ( (int)(n) , sizeof(double));
    res = (double *) calloc ( (int)(n) , sizeof(double));
    fmax = (double *) calloc ( (int)(nmax*k) , sizeof(double));
    ft = (double *) calloc ( (int)(n*k) , sizeof(double));
    fpf = (double *) calloc ( (int)(k*k) , sizeof(double));
    fpy = (double *) calloc ( (int)(k) , sizeof(double));




    for ( i = 0; i < n; i ++)
    {
        for (j = 0; j < nply; j ++)
        {
            f[i*k+j] = pow(x[i], j);
        }

        T = tpsec/86400.0;
        cpr = 1;
        for (j = 0; j < ncpr; j = j + 2)
        {
            for (m = 0; m < nplycpr; m ++)
            {
                f[i * k + nply + j * nplycpr + m ] = pow(x[i], m) * sin (TWOPI / T * x[i]);
                f[i * k + nply + j * nplycpr + nplycpr + m] = pow(x[i], m) *  cos (TWOPI / T * x[i]);
            }
            cpr ++;
//            T = T / cpr;
            T = tpsec/86400.0 / cpr;
        }
        
    }



    for ( i = 0; i < nmax; i ++)
    {
        for (j = 0; j < nply; j ++)
        {
            fmax[i*k+j] = pow(xmax[i], j);
        }

        T = tpsec/86400.0;
        cpr = 1;
        for (j = 0; j < ncpr; j = j + 2)
        {
            for (m = 0; m < nplycpr; m ++)
            {
                fmax[i * k + nply + j * nplycpr + m ] = pow(xmax[i], m) * sin (TWOPI / T * xmax[i]);
                fmax[i * k + nply + j * nplycpr + nplycpr + m] = pow(xmax[i], m) *  cos (TWOPI / T * xmax[i]);
            }
            cpr ++;
//            T = T / cpr;
            T = tpsec/86400.0 / cpr;
        }
        
    }



    mt (f, n, k, ft);
    brmul (ft, f, k, n, k, fpf);
    brmul (ft, y, k, n, 1, fpy);
    err = brinv(fpf,k);
    if (err == 0)
        return 1;
    brmul (fpf, fpy, k, k, 1, a);


    brmul (fmax, a, nmax, k, 1, ymaxflt);

    brmul (f, a, n, k, 1, yflt);
    for (j = 0; j < n; j ++)
        res[j] = yflt[j] - y[j];

    stdres = std_dev (res, n);

    free (f);
    free (ft);
    free (fpf);
    free (fpy);
    free (a);
    free (res);
    free (yflt);
    free (fmax);



    return stdres;
}

















double lsf_cpr_new (double *xmax, double *x, double *y, int nmax, int n, 
    double tpsec, double *yflt, int order_poly, int order_cpr)
{

    int i, j, err, k, cpr;

    double T, *a, *f, *fmax, *ft, *fpf, *fpy, *yf, *res, stdres;

    k = order_poly + 1 + order_cpr * 2;

    a = (double *) calloc ( (int)(k) , sizeof(double));
    f = (double *) calloc ( (int)(n*k) , sizeof(double));
    yf = (double *) calloc ( (int)(n) , sizeof(double));
    res = (double *) calloc ( (int)(n) , sizeof(double));
    fmax = (double *) calloc ( (int)(nmax*k) , sizeof(double));
    ft = (double *) calloc ( (int)(n*k) , sizeof(double));
    fpf = (double *) calloc ( (int)(k*k) , sizeof(double));
    fpy = (double *) calloc ( (int)(k) , sizeof(double));




    for ( i = 0; i < n; i ++)
    {
        for (j = 0; j < order_poly + 1; j ++)
        {
            f[i*k+j] = pow(x[i], j);
        }

        T = tpsec/86400.0;
        cpr = 1;
        for (j = order_poly + 1; j < k; j = j + 2)
        {
            f[i*k+j] = sin (TWOPI / T * x[i]);
            f[i*k+j + 1] = cos (TWOPI / T * x[i]);
            cpr ++;
//            T = T / cpr;
            T = tpsec/86400.0 / cpr;
        }
        
    }



    for ( i = 0; i < nmax; i ++)
    {
        for (j = 0; j < order_poly + 1; j ++)
        {
            fmax[i*k+j] = pow(xmax[i], j);
        }

        T = tpsec/86400.0;
        cpr = 1;
        for (j = order_poly + 1; j < k; j = j + 2)
        {
            fmax[i*k+j] = sin (TWOPI / T * xmax[i]);
            fmax[i*k+j + 1] = cos (TWOPI / T * xmax[i]);
            cpr ++;
//            T = T / cpr;
            T = tpsec/86400.0 / cpr;
        }
        
    }



    mt (f, n, k, ft);
    brmul (ft, f, k, n, k, fpf);
    brmul (ft, y, k, n, 1, fpy);
    err = brinv(fpf,k);
    if (err == 0)
        return 1;
    brmul (fpf, fpy, k, k, 1, a);


    brmul (fmax, a, nmax, k, 1, yflt);

    brmul (f, a, n, k, 1, yf);
    for (j = 0; j < n; j ++)
        res[j] = yf[j] - y[j];

    stdres = std_dev (res, n);

    free (f);
    free (ft);
    free (fpf);
    free (fpy);
    free (a);
    free (res);
    free (yf);
    free (fmax);



    return stdres;
}

















int lsf_cpr_day ( double *x, double *y, int n, double tpsec, double *yflt, int order_poly, int order_cpr)
{

    int i, j, err, k, cpr;

    double T, *a, *f, *ft, *fpf, *fpy;

    k = order_poly + 1 + order_cpr * 2 + 2;

    a = (double *) calloc ( (int)(k) , sizeof(double));
    f = (double *) calloc ( (int)(n*k) , sizeof(double));
    ft = (double *) calloc ( (int)(n*k) , sizeof(double));
    fpf = (double *) calloc ( (int)(k*k) , sizeof(double));
    fpy = (double *) calloc ( (int)(k) , sizeof(double));


    y[0] = y[4];
    y[1] = y[4];
    y[2] = y[4];


    for ( i = 0; i < n; i ++)
    {
        for (j = 0; j < order_poly + 1; j ++)
        {
            f[i*k+j] = pow(x[i], j);
        }

        T = tpsec/86400.0;
        cpr = 1;
        for (j = order_poly + 1; j < k-2; j = j + 2)
        {
            f[i*k+j] = sin (TWOPI / T * x[i]);
            f[i*k+j + 1] = cos (TWOPI / T * x[i]);
            cpr ++;
//            T = T / cpr;
            T = tpsec/86400.0 / cpr;
        }
        
        for (j = k-2; j < k; j = j + 2)
        {
            f[i*k+j] = sin (TWOPI  * x[i]);
            f[i*k+j + 1] = cos (TWOPI  * x[i]);
        }

    }

/*
    for ( i = 0; i < n; i ++)
    {
        for (j = 0; j < k; j ++)
        {
            printf("%f\t", f[i*k+j]);
        }

        printf("%f\t%f\n", y[i], x[i]);

    }
*/


    mt (f, n, k, ft);
    brmul (ft, f, k, n, k, fpf);
    brmul (ft, y, k, n, 1, fpy);
    err = brinv(fpf,k);
    if (err == 0)
        return 1;
    brmul (fpf, fpy, k, k, 1, a);

/*
    for ( i = 0; i < k; i ++)
    {
        for (j = 0; j < k; j ++)
        {
            printf("%f\t", fpf[i*k+j]);
        }

        printf("%f\n", fpy[i]);

    }

 
    for (j = 0; j < k; j ++)
    {
        printf("%f\t", a[j]);
    }
    printf("\n");

*/

    brmul (f, a, n, k, 1, yflt);

    free (f);
    free (ft);
    free (fpf);
    free (fpy);
    free (a);
    return 0;
}














int lsf_cpr ( double *x, double *y, int n, double tpsec, double *yflt, int order_poly, int order_cpr)
{

    int i, j, err, k, cpr;

    double T, *a, *f, *ft, *fpf, *fpy;

    k = order_poly + 1 + order_cpr * 2;

    a = (double *) calloc ( (int)(k) , sizeof(double));
    f = (double *) calloc ( (int)(n*k) , sizeof(double));
    ft = (double *) calloc ( (int)(n*k) , sizeof(double));
    fpf = (double *) calloc ( (int)(k*k) , sizeof(double));
    fpy = (double *) calloc ( (int)(k) , sizeof(double));


    y[0] = y[4];
    y[1] = y[4];
    y[2] = y[4];


    for ( i = 0; i < n; i ++)
    {
        for (j = 0; j < order_poly + 1; j ++)
        {
            f[i*k+j] = pow(x[i], j);
        }

        T = tpsec/86400.0;
        cpr = 1;
        for (j = order_poly + 1; j < k; j = j + 2)
        {
            f[i*k+j] = sin (TWOPI / T * x[i]);
            f[i*k+j + 1] = cos (TWOPI / T * x[i]);
            cpr ++;
//            T = T / cpr;
            T = tpsec/86400.0 / cpr;
        }
        
    }

/*
    for ( i = 0; i < n; i ++)
    {
        for (j = 0; j < k; j ++)
        {
            printf("%f\t", f[i*k+j]);
        }

        printf("%f\t%f\n", y[i], x[i]);

    }
*/


    mt (f, n, k, ft);
    brmul (ft, f, k, n, k, fpf);
    brmul (ft, y, k, n, 1, fpy);
    err = brinv(fpf,k);
    if (err == 0)
        return 1;
    brmul (fpf, fpy, k, k, 1, a);

/*
    for ( i = 0; i < k; i ++)
    {
        for (j = 0; j < k; j ++)
        {
            printf("%f\t", fpf[i*k+j]);
        }

        printf("%f\n", fpy[i]);

    }

 
    for (j = 0; j < k; j ++)
    {
        printf("%f\t", a[j]);
    }
    printf("\n");

*/

    brmul (f, a, n, k, 1, yflt);

    free (f);
    free (ft);
    free (fpf);
    free (fpy);
    free (a);
    return 0;
}


int lsf_poly ( double *x, double *y, int n, double *a, int k)
{

    int i, j, err;

    double *f, *ft, *fpf, *fpy;

    f = (double *) calloc ( (int)(n*k) , sizeof(double));
    ft = (double *) calloc ( (int)(n*k) , sizeof(double));
    fpf = (double *) calloc ( (int)(k*k) , sizeof(double));
    fpy = (double *) calloc ( (int)(k) , sizeof(double));

    
    for ( i = 0; i < n; i ++)
    {
        for (j = 0; j < k; j ++)
        {
            f[i*k+j] = pow(x[i], j); 
        }
    }

    mt (f, n, k, ft);           
    brmul (ft, f, k, n, k, fpf);    
    brmul (ft, y, k, n, 1, fpy);    
    err = brinv(fpf,k); 
    if (err == 0)
        return 1;   
    brmul (fpf, fpy, k, k, 1, a);   //¾ØÕóÏà³Ë, ¾ÍÊÇ¼ÆËã(FTPyF)-1FTPyFÕâ¸ö¾ØÕóÁË, ¾ÍÊÇËã³öÁË¹«Ê½1, ·µ»Ø¾ØÕóa¾ÍÊÇ¹«Ê½1ÀïµÄc, Ò²¾ÍËã³öÀ´ÁË

    free (f);
    free (ft);
    free (fpf);
    free (fpy);
    return 0;
}










