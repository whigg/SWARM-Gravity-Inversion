
/*

------------------------------------------------------------------------

    Purpose: operation between spherical harmonics coefficients 
    Notes: 
    Programmer: Kun Shang @ 4.29.2014
    Functions:

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

        double lgdr2(double t, int nmax, int m, 
             double *pbar, double *pbar1, double *pbar2);

        double normfct (int n, int m);

        double addsubcs (double *cs, double *cs1, int nmax, int nmax1, 
            int flag);

        double zero0zero1 (double *cs, int nmax);
        double unit0zero1 (double *cs, int nmax);

        double csnm (double *cs, int nmax, int n, int m, 
            double *c, double *s);

     Global variables:

        extern int ;
        extern double ;



------------------------------------------------------------------------



*/




#ifndef _GRVTS_H_

#include "grvts.h"
    #define _GRVTS_H_
#endif



    CSStruct CSinfo[1000];
    
    int MGCS;



/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/*
*/
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
double eig_open (char *grv_name, double dyear, 
        int nmax, int mmax, double *coef)
{
    FILE *fp_grv;
    double c,s, nd,md, prd;
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
 
        sscanf (string, "%s", name);	
        if (strcmp (name,"gfc") ==0)	
        {
            sscanf (string, "%*s%lf%lf%lf%lf", &nd, &md, &c, &s);
            n = (int)nd;
            m = (int)md;
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


        if (strcmp (name,"gfct") ==0)
        {
            sscanf (string, "%*s%lf%lf%lf%lf", &nd, &md, &c, &s);
            n = (int)nd;
            m = (int)md;
            if (n > nmax || m > mmax)
                continue;
            else if (m == 0)
            {
                coef[n] += c;
            }
            else
            {
                l = nmax - m + 1;
                ind = nmax + 1 + (2 * nmax - m + 2) * (m - 1);
                coef[ind + n - m] += c;
                coef[ind + n - m + l] += s;
            }
        }


        if (strcmp (name,"trnd") ==0)
        {
            sscanf (string, "%*s%lf%lf%lf%lf", &nd, &md, &c, &s);
            n = (int)nd;
            m = (int)md;
            if (n > nmax || m > mmax)
                continue;
            else if (m == 0)
            {
                coef[n] += c * dyear;
            }
            else
            {
                l = nmax - m + 1;
                ind = nmax + 1 + (2 * nmax - m + 2) * (m - 1);
                coef[ind + n - m] += c * dyear;
                coef[ind + n - m + l] += s * dyear;
            }
        }


        if (strcmp (name,"asin") ==0)
        {
            sscanf (string, "%*s%lf%lf%lf%lf%*f%*f%lf", &nd, &md, &c, &s, &prd);
            n = (int)nd;
            m = (int)md;
            if (n > nmax || m > mmax)
                continue;
            else if (m == 0)
            {
                coef[n] += c * sin (TWOPI/prd*dyear);
            }
            else
            {
                l = nmax - m + 1;
                ind = nmax + 1 + (2 * nmax - m + 2) * (m - 1);
                coef[ind + n - m] += c * sin (TWOPI/prd*dyear);
                coef[ind + n - m + l] += s * sin (TWOPI/prd*dyear);
            }
        }


        if (strcmp (name,"acos") ==0)
        {
            sscanf (string, "%*s%lf%lf%lf%lf%*f%*f%lf", &nd, &md, &c, &s, &prd);
            n = (int)nd;
            m = (int)md;
            if (n > nmax || m > mmax)
                continue;
            else if (m == 0)
            {
                coef[n] += c * cos (TWOPI/prd*dyear);
            }
            else
            {
                l = nmax - m + 1;
                ind = nmax + 1 + (2 * nmax - m + 2) * (m - 1);
                coef[ind + n - m] += c * cos (TWOPI/prd*dyear);
                coef[ind + n - m + l] += s * cos (TWOPI/prd*dyear);
            }
        }

    }

    fclose(fp_grv);
    return 0;
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


























double cs2ada (double *llr, double *cs, double gm, double a, 
        int nmax, double *ae, 
        int part, double *dadre, int flagdadcs)
{
    int n, m, k, l, ind, ic, is, label;

    double a2r, slat, clat, slon, clon, sclt, cclt, *cosml, *sinml, 
        *aprn, *pbar, *pbar1, *pbar2, *pt, *pt1, *pt2, lat, lon, r, vi, t, 
        an[3], c_en[9], c_ne[9], dadrn[9], dadrne[9], 
        dcdrn[27], dcdre[27], dcdxe[9], dcdye[9], dcdze[9], 
        dadrecx[3], dadrecy[3], dadrecz[3], dadrec[9], dadrea[9];

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/

    lat = llr[0];
    lon = llr[1];
    r = llr[2];

    slat = sin(lat * DEG2RAD);
    clat = cos(lat * DEG2RAD);
    slon = sin(lon * DEG2RAD);
    clon = cos(lon * DEG2RAD);
    sclt = clat;
    cclt = slat;

    cosml = (double *) calloc ( nmax + 1, sizeof(double)); //cos(m*lamta)
    sinml = (double *) calloc ( nmax + 1, sizeof(double)); //sin(m*lamta)
    aprn = (double *) calloc ( nmax + 1, sizeof(double));  //sin(m*lamta)
    pbar = (double *) calloc ( nmax + 1, sizeof(double));
    pbar1 = (double *) calloc ( nmax + 1, sizeof(double));
    pbar2 = (double *) calloc ( nmax + 1, sizeof(double));
    pt  = (double *) calloc ( (nmax + 1) * (nmax + 1), sizeof(double));
    pt1  = (double *) calloc ( (nmax + 1) * (nmax + 1), sizeof(double));
    pt2  = (double *) calloc ( (nmax + 1) * (nmax + 1), sizeof(double));

/*
    cosml[0] = 1; sinml[0] = 0;
    for (m = 1; m <= nmax; m++)
    {
        cosml[m] = cos(m * lon * DEG2RAD);
        sinml[m] = sin(m * lon * DEG2RAD);
    }
    
    for (n = 0; n <= nmax; n++)
    {
        aprn[n] = pow (a / r, n) * gm / r;
    }
*/


//    printf ("!!\n");
    cosml[0] = 1; sinml[0] = 0;
    cosml[1] = cos(lon * DEG2RAD); sinml[1] = sin(lon * DEG2RAD);

    for (m = 2; m <= nmax; m++)
    {
        cosml[m] = 2.0 * cosml[1] * cosml[m-1] - cosml[m-2];
        sinml[m] = 2.0 * cosml[1] * sinml[m-1] - sinml[m-2];
    }

    aprn[0] = gm / r;
    a2r = a / r;

    for (n = 1; n <= nmax; n++)
    {
        aprn[n] = aprn[n - 1] * a2r;
    }

//    return 0;

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/

// from here


    an[0] = 0; an[1] = 0; an[2] = 0;
    t = cclt; vi = 0;

    for (m = 0; m <= nmax; m ++)
    {
        l = nmax - m + 1;

        if (part == 1)
            lgdr2i(t, nmax, m, pbar, pbar1, pbar2);
        else
            lgdr1i(t, nmax, m, pbar, pbar1);

        for (k = 0; k < l; k++)
        {
            if (m==0)
            {
                n = k;
                ic = n;
                is = 0;
            }

            else
            {
                ind = nmax + 1 + (2 * nmax - m + 2) * (m - 1);
                n = k + m;
                ic = ind + n - m;
                is = ind + n - m + l;
            }

            pt[ic] = aprn[n] * pbar[k] * cosml[m];
            pt[is] = aprn[n] * pbar[k] * sinml[m];
            pt1[ic] = aprn[n] * pbar1[k] * cosml[m];
            pt1[is] = aprn[n] * pbar1[k] * sinml[m];
            pt2[ic] = aprn[n] * pbar2[k] * cosml[m];
            pt2[is] = aprn[n] * pbar2[k] * sinml[m];

            vi = vi + pt[ic] * cs[ic] +  pt[is] * cs[is];
            an[0] = an[0] + pt1[ic] * cs[ic] + pt1[is] * cs[is];
            an[1] = an[1] - m *  pt[is] * cs[ic] + m *  pt[ic] * cs[is];          
            an[2] = an[2] + (n+1) * pt[ic] * cs[ic] + (n+1) * pt[is] * cs[is];
//            an[2] = an[2] + (n+1) * ( pt[ic] * cs[ic] +  pt[is] * cs[is]);

        }
    }


// to here

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
    free (pbar);
    free (pbar1);
    free (pbar2);

    free (cosml);
    free (sinml);
    free (aprn);


    an[0] = - an[0] / r;
    an[1] = + an[1] / r / sclt;
    an[2] = + an[2] / r;


    c_ne[0] = - slat * clon;      //nsys: North-East-Down
    c_ne[1] = - slon;
    c_ne[2] = - clat * clon;
    c_ne[3] = - slat * slon;
    c_ne[4] = clon;
    c_ne[5] = - clat * slon;
    c_ne[6] = clat;
    c_ne[7] = 0;
    c_ne[8] = - slat;



    brmul(c_ne, an, 3, 3, 1, ae);  //from n-sys to e-sys

    for (n = 0; n < 9; n++)
        dadre[n] = 0;

    if (part == 0)
    {
        free (pt);
        free (pt1);
        free (pt2);
        return vi;
    }



    for (n = 0; n < 9; n++)
        dadrn[n] = 0;
    for (m = 0; m <= nmax; m ++)
    {
        l = nmax - m + 1;
        for (k = 0; k < l; k++)
        {
            if (m==0)
            {
                n = k;
                ic = n;
                is = 0;
            }
            else
            {
                ind = nmax + 1 + (2 * nmax - m + 2) * (m - 1);
                n = k + m;
                ic = ind + n - m;
                is = ind + n - m + l;
            }
            
            dadrn[0] += pt2[ic] * cs[ic] + pt2[is] * cs[is];
            dadrn[3] += m * ( - pt[is] * cs[ic] + pt[ic] * cs[is]) * cclt / sclt / sclt
                      - m * ( - pt1[is] * cs[ic] + pt1[ic] * cs[is]) / sclt;
            dadrn[6] -= (n+1) * (pt1[ic] * cs[ic] + pt1[is] * cs[is]);

            dadrn[1] -= m * ( - pt1[is] * cs[ic] + pt1[ic] * cs[is]) / sclt;
            dadrn[4] -= m * m * (pt[ic] * cs[ic] + pt[is] * cs[is]) / sclt / sclt;
            dadrn[7] += m * (n+1) * ( - pt[is] * cs[ic] + pt[ic] * cs[is]) / sclt;

            dadrn[2] -= (n+2) * (pt1[ic] * cs[ic] + pt1[is] * cs[is]);
            dadrn[5] += m * (n+2) * ( - pt[is] * cs[ic] + pt[ic] * cs[is]) / sclt;
            dadrn[8] += (n+2) * (n+1)* (pt[ic] * cs[ic] + pt[is] * cs[is]);
        }
    }

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/

    for (n = 0; n < 9; n++)
        dadrn[n] = dadrn[n] / r / r;

    mt(c_ne, 3, 3, c_en);
    brmul(dadrn, c_en, 3, 3, 3, dadrne);			
    brmul(c_ne, dadrne, 3, 3, 3, dadrea);	




    for (n = 0; n < 27; n++)
        dcdrn[n] = 0;

    dcdrn[0]  = - sclt * clon;  dcdrn[1]  = cclt * slon / sclt;
    dcdrn[3]  = 0;              dcdrn[4]  = - clon / sclt;
    dcdrn[6]  = cclt * clon;/**/dcdrn[7]  = slon;

    dcdrn[9]  = - sclt * slon;  dcdrn[10] = - cclt * clon / sclt;
    dcdrn[12] = 0;              dcdrn[13] = - slon / sclt;
    dcdrn[15] = cclt * slon;    dcdrn[16] = - clon;

    dcdrn[18] = - cclt;         dcdrn[19] = 0;
    dcdrn[21] = 0;              dcdrn[22] = 0;
    dcdrn[24] = - sclt;         dcdrn[25] = 0;

    for (n = 0; n < 27; n++)
        dcdrn[n] = dcdrn[n] / r;

    brmul(dcdrn, c_en, 9, 3, 3, dcdre);		

    for (n = 0; n < 9; n++)	
    {
        dcdxe[n] = dcdre[n*3];
        dcdye[n] = dcdre[n*3+1];
        dcdze[n] = dcdre[n*3+2];
//        printf ("%e\n", dcdze[n]);
    }


    brmul(dcdxe, an, 3, 3, 1, dadrecx);	
    brmul(dcdye, an, 3, 3, 1, dadrecy);  
    brmul(dcdze, an, 3, 3, 1, dadrecz);  

    for (n = 0; n < 3; n++)	
    {
        dadrec[n*3] = dadrecx[n]; 
        dadrec[n*3+1] = dadrecy[n]; 
        dadrec[n*3+2] = dadrecz[n];
    }




    for (n = 0; n <= 8; n++)
    {
        dadre[n] = dadrec[n] + dadrea[n]; 
    }


    if (flagdadcs == 0)
    {
        free (pt);
        free (pt1);
        free (pt2);
        return vi;
    }




// for dadcs_nm
//



    for (k = 0; k < MGCS; k ++)
    {   
        n = CSinfo[k].n; m = CSinfo[k].m; label = CSinfo[k].cs;

        if (m == 0)
        {
            ic = n;
            CSinfo[k].dadcsn[0] = - pt1[ic] / r;
            CSinfo[k].dadcsn[1] = 0;          
            CSinfo[k].dadcsn[2] = (n+1) * pt[ic] / r;
        }
        else 
        {
            l = nmax - m + 1;
            ind = nmax + 1 + (2 * nmax - m + 2) * (m - 1);
            ic = ind + n - m;
            is = ind + n - m + l;
            if (label == 1)
            {
                CSinfo[k].dadcsn[0] = - pt1[ic] / r; 
                CSinfo[k].dadcsn[1] = - m *  pt[is] / r / sclt;          
                CSinfo[k].dadcsn[2] = (n+1) * pt[ic] / r;
            }

            if (label == -1)
            {
                CSinfo[k].dadcsn[0] = - pt1[is] / r;
                CSinfo[k].dadcsn[1] =   m *  pt[ic] / r / sclt;          
                CSinfo[k].dadcsn[2] = (n+1) * pt[is] / r;
            }
        }

        brmul(c_ne, CSinfo[k].dadcsn, 3, 3, 1, CSinfo[k].dadcse);	
    }





    free (pt);
    free (pt1);
    free (pt2);
    return vi;

    

}
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
































/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/*
*/
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
double cs2ac (double *llr, double *cs, double gm, double a, 
            int nmax, double *acc)
{
    int n, m, k, l, ind;

    double sinf, cosf, sinlon, coslon, sincolat, coscolat, *cosml, *sinml, 
        *aprn, *pbar, *pbar1, *pbar2, accn[3], c_en[9],
        *pt, *ptt, lat, lon, r, vi, dvdr, dvdcolat, dvdlon, t;


/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/

    lat = llr[0];
    lon = llr[1];
    r = llr[2];

    sinf = sin(lat * DEG2RAD);
    cosf = cos(lat * DEG2RAD);
    sinlon = sin(lon * DEG2RAD);
    coslon = cos(lon * DEG2RAD);
    sincolat = cosf;
    coscolat = sinf;

    cosml = (double *) calloc ( nmax + 1, sizeof(double)); //cos(m*lamta)
    sinml = (double *) calloc ( nmax + 1, sizeof(double)); //sin(m*lamta)
    aprn = (double *) calloc ( nmax + 1, sizeof(double));  //sin(m*lamta)
    pbar = (double *) calloc ( nmax + 1, sizeof(double));
    pbar1 = (double *) calloc ( nmax + 1, sizeof(double));
    pbar2 = (double *) calloc ( nmax + 1, sizeof(double));
    pt  = (double *) calloc ( (nmax + 1) * (nmax + 1), sizeof(double));
    ptt  = (double *) calloc ( (nmax + 1) * (nmax + 1), sizeof(double));

    for (m = 0; m <= nmax; m++)
    {
        cosml[m] = cos(m * lon * DEG2RAD);
        sinml[m] = sin(m * lon * DEG2RAD);
    }
    
    for (n = 0; n <= nmax; n++)
    {
        aprn[n] = pow (a / r, n) * gm / r;
    }

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/

    t = coscolat; vi = 0; dvdlon = 0; dvdcolat = 0; dvdr = 0;
    for (m = 0; m <= nmax; m ++)
    {
        l = nmax - m + 1;
        lgdr2(t, nmax, m, pbar, pbar1, pbar2);
//        lgdr(t, nmax, m, pbar);
        for (k = 0; k < l; k++)
        {
            if (m==0)
            {
//                ind = 0;
                n = k + m;
                pt[k] = aprn[n] * pbar[k];
                ptt[k] = aprn[n] * pbar1[k];
                vi = vi + pt[k] * cs[k];
            
//                if (n>=2)
                {
                    dvdr = dvdr + (n+1) * pt[k] * cs[k];
                    dvdcolat = dvdcolat + ptt[k] * cs[k];
                }
            }
            else
            {
                ind = nmax + 1 + (2 * nmax - m + 2) * (m - 1);
                n = k + m;
                pt[ind + n - m] = aprn[n] * pbar[k] * cosml[m];
                pt[ind + n - m + l] = aprn[n] * pbar[k] * sinml[m];
                ptt[ind + n - m] = aprn[n] * pbar1[k] * cosml[m];
                ptt[ind + n - m + l] = aprn[n] * pbar1[k] * sinml[m];
                vi = vi + pt[ind + n - m] * cs[ind + n - m];
                vi = vi + pt[ind + n - m + l] * cs[ind + n - m + l];

                dvdcolat = dvdcolat + ptt[ind + n - m] * cs[ind + n - m];
                dvdcolat = dvdcolat + ptt[ind + n - m + l] * cs[ind + n - m + l];

                dvdlon = dvdlon - m * pt[ind + n - m + l] * cs[ind + n - m];
                dvdlon = dvdlon + m * pt[ind + n - m] * cs[ind + n - m + l];          

                dvdr = dvdr + (n+1) * pt[ind + n - m] * cs[ind + n - m];
                dvdr = dvdr + (n+1) * pt[ind + n - m + l] * cs[ind + n - m + l];
            }
        }
    }

//    dvdcolat = - dvdcolat * sincolat; //tmd!!
    dvdcolat = dvdcolat;
    dvdlon = + dvdlon;
    dvdr = - dvdr / r;


    accn[0] = - dvdcolat / r;
    accn[1] = + dvdlon / r / sincolat;
    accn[2] = - dvdr;


    c_en[0] = - sinf * coslon;      //from fixed to up-east-north system: rmat
    c_en[1] = - sinlon;
    c_en[2] = - cosf * coslon;
    c_en[3] = - sinf * sinlon;
    c_en[4] = coslon;
    c_en[5] = - cosf * sinlon;
    c_en[6] = cosf;
    c_en[7] = 0;
    c_en[8] = - sinf;


    brmul(c_en, accn, 3, 3, 1, acc);  //from fixed acc to inertial acc

//    *v = vi;
//    *dvdt = - ANGVEL * dvdlon;

    free (pbar);
    free (pbar1);
    free (pbar2);
    free (pt);
    free (ptt);

    free (cosml);
    free (sinml);
    free (aprn);
    return 1;
}
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/













/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/*
*/
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
double cspt2gp (double *pt, double *cs, int nmax, 
        double *gpt)
{
    int k;
    double gpti;

    gpti = 0;
    for(k = 0; k < (nmax + 1) * (nmax + 1); k++)
    {
        gpti = gpti + pt[k] * cs[k];
    }

    *gpt = gpti; 

    return 1;
}
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/


double cs2gp (double *llr, double *cs, double gm, double a, 
        int nmax, double *gpt)
{
    int n, m, k, l, ind;

    double sinf, cosf, *cosml, *sinml, *aprn, *pbar, lat, lon, r, 
        gpti, t, *pt;


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
    pt  = (double *) calloc ( (nmax + 1) * (nmax + 1), sizeof(double));

    for (m = 0; m <= nmax; m++)
    {
        cosml[m] = cos(m * lon * DEG2RAD);
        sinml[m] = sin(m * lon * DEG2RAD);
    }
    
    for (n = 0; n <= nmax; n++)
    {
        aprn[n] = pow (a / r, n) * gm / r;
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
    free (pt);
    free (cosml);
    free (sinml);
    free (aprn);
    return 1;
}
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/














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






/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/

double lgdr2(double t, int nmax, int m, 
             double *pbar, double *pbar1, double *pbar2)

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
        pbar1[0] = sqrt(3.0) * t;

        for (i = 2; i <= m; i++)
        {
            pbar1[0] = sqrt((2.0*i+1.0)/(2.0*i))*(sqrt(1.0-t*t)*pbar1[0]+t*pbar[0]);
//            pbar1[0] = sqrt((2.0*i+1.0)/(2.0*i))*(sqrt(1.0-t*t)*pbar1[0]+t*pbar[0]/(-sqrt(1.0-t*t)));
            pbar[0] = sqrt((2.0*i+1.0)/(2.0*i)*(1.0-t*t))*pbar[0];
        }
    }
    else 
    {
        pbar[0]=p00; 
        pbar1[0]=0.0;
    }

//    ! Pm+1,m : JEKEIL (A.3b)

    if (nmax - m + 1 >= 2)
    {
        pbar[1] = sqrt(2.0*m +3.0) * t * pbar[0];
    }

//  ! Pn,m (n>=m+2) : JEKEIL (A.3a)

    for(i = 3; i <= nmax-m+1; i++)
    {
        c=((2.0*m+2.0*i-3.0) * (2.0*m + 2.0*i-1.0)) / ((i-1.0)*(2.0*m+i-1.0));
        d=((2.0*m+2.0*i-1.0)*(2.0*m+i-2.0)*(i-2.0))/((2.0*m+2.0*i-5.0)*(i-1.0)*(2.0*m+i-1.0));
        pbar[i-1] = sqrt(c)*t*pbar[i-2] - sqrt(d) * pbar[i-3];
    }

//  ! THE FULLY NORMALIZED ASSOCIATED LEGENDRE FUNCTION - 1ST DERIVATIVE
//  ! P'n,m (n>=m+1) : TSCHERNING (8)
    for (i=2; i<=nmax-m+1; i++)
    {
        c = 1.0/sqrt(1.0-t*t)*t*(m+i-1);
        d = 1.0/sqrt(1.0-t*t)*sqrt((((m+i-1)*(m+i-1)-m*m)*(2.0*(m+i-1)+1.0))/(2.0*(m+i-1)-1.0));
//!! found it different from TSCHERNING (8),dcl-2010-2-14
//!! Jianbin confirms code is correct, dcl-2010-2-15
//!!      D=1D0/SQRT(1D0-T**2)/SQRT((((M+I-1)**2-M**2)*(2D0*(M+I-1)+1D0))/(2D0*(M+I-1)-1D0))
        pbar1[i-1] = c * pbar[i-1] - d * pbar[i-2];
    }

//! THE FULLY NORMALIZED ASSOCIATED LEGENDRE FUNCTION - 2ND DERIVATIVE
//! P''n,m (n>=m) : ASSOCIATE LEGENDRE EQUATION (2ND ORDER DIFFERENTIAL EQN.)

    for (i=1;i<=nmax-m+1;i++)
    {
        pbar2[i-1] = (-t/sqrt(1.0-t*t)) * pbar1[i-1] 
            - ((m+i-1)*(m+i)-m*m/(1.0-t*t)) * pbar[i-1];
    }
    return 0;
}
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/





/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/

#define NMAXT 201

double lgdr2i (double t, int nmax, int m, 
             double *pbar, double *pbar1, double *pbar2)
{
    int i;
    double p00, p11, ts, sqrt3 = 1.73205080756888;
    
    static int flag = 0;
    static double ani[NMAXT], bni[NMAXT], ami[NMAXT][NMAXT], 
        bmi[NMAXT][NMAXT], cmi[NMAXT][NMAXT], dmi[NMAXT][NMAXT];
    
    int mi, ii, ni;

    if (flag != 9) 
    {
        flag = 9;

        for (ii = 0; ii < NMAXT; ii ++)
            ani[ii] = sqrt ((2.0 * ii + 1.0) / (2.0 * ii));

        for (ii = 0; ii < NMAXT; ii ++)
            bni[ii] = sqrt (2.0 * ii + 3.0);

        for (mi = 0; mi < NMAXT; mi ++)
        {
            for (ii = 3; ii < NMAXT - mi + 1; ii++)
            {
                ami[mi][ii] = sqrt(( (2.0 * mi + 2.0 * ii - 3.0) * (2.0 * mi + 2.0 * ii - 1.0) ) 
                    / ( (ii - 1.0) * (2.0 * mi + ii - 1.0)));

                bmi[mi][ii] = sqrt(( (2.0 * mi + 2.0 * ii - 1.0) * (2.0 * mi + ii - 2.0) * (ii - 2.0) ) 
                    / ( (2.0 * mi + 2.0 * ii - 5.0) * (ii - 1.0) * (2.0 * mi + ii - 1.0) ));
            }
        }

        for (mi = 0; mi < NMAXT; mi ++)
        {
            for (ii = 2; ii < NMAXT - m + 1; ii++)
            {
                ni = mi + ii - 1;
                cmi[mi][ii] = ni;
                dmi[mi][ii] = sqrt(((ni * ni - mi * mi) * (2.0 * ni + 1.0)) / (2.0 * ni - 1.0));
            }
        }
    }


    ts = sqrt (1.0 - t * t);

    p00 = 1.0; 
    p11 = sqrt3 * ts;

    if (m >= 1)
    {
        pbar[0]  = p11; 
        pbar1[0] = sqrt3 * t;

        for (i = 2; i <= m; i++)
        {
            pbar1[0] = ani[i] * ( ts * pbar1[0] + t * pbar[0]);
            pbar[0]  = ani[i] * ts * pbar[0];
        }
    }
    else 
    {
        pbar[0]  = p00; 
        pbar1[0] = 0.0;
    }


    if (nmax - m + 1 >= 2)
        pbar[1] = bni[m] * t * pbar[0];

    for (i = 3; i <= nmax - m + 1; i++)
        pbar[i-1] = ami[m][i] * t * pbar[i-2] - bmi[m][i] * pbar[i-3];

    for (i = 2; i <= nmax - m + 1; i++)
        pbar1[i-1] = cmi[m][i] * t / ts * pbar[i-1] - dmi[m][i] / ts * pbar[i-2];

    for (i = 1; i <= nmax - m + 1; i++)
        pbar2[i-1] = - t / ts * pbar1[i-1]  
            - ((m + i - 1) * (m + i) - m * m / ts / ts) * pbar[i-1];

    return 0;
}






double lgdr1i (double t, int nmax, int m, 
             double *pbar, double *pbar1)
{
    int i;
    double p00, p11, ts, sqrt3 = 1.73205080756888;
    
    static int flag = 0;
    static double ani[NMAXT], bni[NMAXT], ami[NMAXT][NMAXT], 
        bmi[NMAXT][NMAXT], cmi[NMAXT][NMAXT], dmi[NMAXT][NMAXT];
    
    int mi, ii, ni;

    if (flag != 9) 
    {
        flag = 9;

        for (ii = 0; ii < NMAXT; ii ++)
            ani[ii] = sqrt ((2.0 * ii + 1.0) / (2.0 * ii));

        for (ii = 0; ii < NMAXT; ii ++)
            bni[ii] = sqrt (2.0 * ii + 3.0);

        for (mi = 0; mi < NMAXT; mi ++)
        {
            for (ii = 3; ii < NMAXT - mi + 1; ii++)
            {
                ami[mi][ii] = sqrt(( (2.0 * mi + 2.0 * ii - 3.0) * (2.0 * mi + 2.0 * ii - 1.0) ) 
                    / ( (ii - 1.0) * (2.0 * mi + ii - 1.0)));

                bmi[mi][ii] = sqrt(( (2.0 * mi + 2.0 * ii - 1.0) * (2.0 * mi + ii - 2.0) * (ii - 2.0) ) 
                    / ( (2.0 * mi + 2.0 * ii - 5.0) * (ii - 1.0) * (2.0 * mi + ii - 1.0) ));
            }
        }

        for (mi = 0; mi < NMAXT; mi ++)
        {
            for (ii = 2; ii < NMAXT - m + 1; ii++)
            {
                ni = mi + ii - 1;
                cmi[mi][ii] = ni;
                dmi[mi][ii] = sqrt(((ni * ni - mi * mi) * (2.0 * ni + 1.0)) / (2.0 * ni - 1.0));
            }
        }
    }


    ts = sqrt (1.0 - t * t);

    p00 = 1.0; 
    p11 = sqrt3 * ts;

    if (m >= 1)
    {
        pbar[0]  = p11; 
        pbar1[0] = sqrt3 * t;

        for (i = 2; i <= m; i++)
        {
            pbar1[0] = ani[i] * ( ts * pbar1[0] + t * pbar[0]);
            pbar[0]  = ani[i] * ts * pbar[0];
        }
    }
    else 
    {
        pbar[0]  = p00; 
        pbar1[0] = 0.0;
    }


    if (nmax - m + 1 >= 2)
        pbar[1] = bni[m] * t * pbar[0];

    for (i = 3; i <= nmax - m + 1; i++)
        pbar[i-1] = ami[m][i] * t * pbar[i-2] - bmi[m][i] * pbar[i-3];

    for (i = 2; i <= nmax - m + 1; i++)
        pbar1[i-1] = cmi[m][i] * t / ts * pbar[i-1] - dmi[m][i] / ts * pbar[i-2];

//    for (i = 1; i <= nmax - m + 1; i++)
//        pbar2[i-1] = - t / ts * pbar1[i-1]  
//            - ((m + i - 1) * (m + i) - m * m / ts / ts) * pbar[i-1];

    return 0;
}



















double normfct (int n, int m)
{
    int i;
    double nmm, npm;

    if (m == 0)
        return sqrt(2.0 * n + 1.0);
    else
    {
        nmm = 1.0;
        npm = 1.0;
        for (i = 1; i <= n - m; i++)
            nmm = nmm * i;
        for (i = 1; i <= n + m; i++)
            npm = npm * i;

        return sqrt(2.0 * (2.0 * n + 1.0) * nmm / npm);
    }
}



/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/

double addsubcs (double *cs, double *cs1, int nmax, int nmax1, 
    int flag)
// cs = cs + cs1; nmax >= nmax1
{
    int n, m, k, l, ind, ic, is;
    double c, s;

    if (nmax < nmax1)
    {
        printf ("error in addsubcs: nmax(%d) < nmax1(%d)\n",
                nmax, nmax1);
        exit (0);
    }
    
    if (flag != 1 && flag != -1)
    {
        printf ("error in addsubcs: flag(%d) != +-1\n", flag);
        exit (0);
    }

    for (m = 0; m <= nmax; m ++)
    {
        l = nmax - m + 1;
        for (k = 0; k < l; k++)
        {
            if (m==0)
            {
                n = k;
                csnm (cs1, nmax1, n, m, &c, &s);
                if (flag ==  1) cs[n] += c;
                if (flag == -1) cs[n] -= c;
            }
            else
            {
                ind = nmax + 1 + (2 * nmax - m + 2) * (m - 1);
                n = k + m;
                csnm (cs1, nmax1, n, m, &c, &s);
                ic = ind + n - m;
                is = ind + n - m + l;
                if (flag ==  1) { cs[ic] += c; cs[is] += s; }
                if (flag == -1) { cs[ic] -= c; cs[is] -= s; }
            }
        }
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


double unit0zero1 (double *cs, int nmax)
{
    int n, m, l, ind, ic, is;
    cs[0] = 1;
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





/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/

double csnm (double *cs, int nmax, int n, int m, 
    double *c, double *s)
{
    int ind, l;

    if (n > nmax || m > n)
    {
        *c = 0;
        *s = 0;
    }
    else if (m == 0)
    {
        *c = cs[n];
        *s = 0;
    }
    else
    {
        l = nmax - m + 1;
        ind = nmax + 1 + (2 * nmax - m + 2) * (m - 1);
        *c = cs[ind + n - m];
        *s = cs[ind + n - m + l];
    }
     
    return 0;


}



