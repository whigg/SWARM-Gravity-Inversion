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

        void rotmatx (double rad, double *matx, short int deri);

        void rotmaty (double rad, double *maty, short int deri);

        void rotmatz (double rad, double *matz, short int deri);

        int bssgj (double *a,int n);

------------------------------------------------------------------------
*/




#ifndef _NUMRS_H_

#include "numrs.h"
    #define _NUMRS_H_
#endif


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








void rotmatx (double rad, double *matx, short int deri)
{
    double cosst, sinst;

    cosst = cos(rad);
    sinst = sin(rad);

    matx[0] = 1;
    matx[1] = 0;
    matx[2] = 0;
    matx[3] = 0;
    matx[4] = cosst;
    matx[5] = sinst;
    matx[6] = 0;
    matx[7] = -sinst;
    matx[8] = cosst;

    if (deri == 1)
    {
        matx[0] = 0;
        matx[1] = 0;
        matx[2] = 0;
        matx[3] = 0;
        matx[4] = -sinst;
        matx[5] = cosst;
        matx[6] = 0;
        matx[7] = -cosst;
        matx[8] = -sinst;
    }
}








/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/*
* rotmaty - 
* @param1: description of param1
* @param2: description of param2
*
* version: 20 Aug 2010
*/
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
void rotmaty (double rad, double *maty, short int deri)
{
    double cosst, sinst;

    cosst = cos(rad);
    sinst = sin(rad);

    maty[0] = cosst;
    maty[1] = 0;
    maty[2] = -sinst;
    maty[3] = 0;
    maty[4] = 1;
    maty[5] = 0;
    maty[6] = sinst;
    maty[7] = 0;
    maty[8] = cosst;

    if (deri == 1)
    {
        maty[0] = -sinst;
        maty[1] = 0;
        maty[2] = -cosst;
        maty[3] = 0;
        maty[4] = 0;
        maty[5] = 0;
        maty[6] = cosst;
        maty[7] = 0;
        maty[8] = -sinst;
    }

}



/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/*
* rotmatz - 
* @param1: description of param1
* @param2: description of param2
*
* version: 20 Aug 2010
*/
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
void rotmatz (double rad, double *matz, short int deri)
{
    double cosst, sinst;

    cosst = cos(rad);
    sinst = sin(rad);

    matz[0] = cosst;
    matz[1] = sinst;
    matz[2] = 0;
    matz[3] = -sinst;
    matz[4] = cosst;
    matz[5] = 0;
    matz[6] = 0;
    matz[7] = 0;
    matz[8] = 1;

    if (deri == 1)
    {
        matz[0] = -sinst;
        matz[1] = cosst;
        matz[2] = 0;
        matz[3] = -cosst;
        matz[4] = -sinst;
        matz[5] = 0;
        matz[6] = 0;
        matz[7] = 0;
        matz[8] = 0;
    }
}






/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/*
* bssgj - 
* @param1: description of param1
* @param2: description of param2
*
* version: 20 Aug 2010
*/
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
int bssgj (double *a,int n)
{
    int i, j, k, m;
    double w, g, *b;

    b = (double *)malloc (n * sizeof(double));
    for (k = 0; k <= n - 1; k++)
    {
        w = a[0];
        if (fabs (w) + 1.0 == 1.0)
        { 
            free (b); 
            printf ("fail\n"); 
            return (-2);
        }
        m = n - k - 1;
        for (i = 1; i <= n - 1; i++)
        { 
            g = a[i * n]; 
            b[i] = g / w;
            if (i <= m) 
                b[i] = - b[i];
            for (j = 1; j <= i; j++)
                a[(i - 1) * n + j - 1] = a[i * n + j] + g * b[j];
        }
        a[n * n - 1] = 1.0 / w;
        for (i = 1; i <= n - 1; i++)
            a[(n - 1) * n + i - 1] = b[i];
    }
	
    for (i=0; i<=n-2; i++)
    {
        for (j=i+1; j<=n-1; j++)
            a[i*n+j]=a[j*n+i];
    }
     
    free(b);
    return(2);
}



