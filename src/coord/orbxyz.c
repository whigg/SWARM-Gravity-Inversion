/*
------------------------------------------------------------------------

    Purpose: transformation between xyz, aei, llr, rtn
    Notes: 
    Programmer: Kun Shang @ 4.29.2014
    Functions:

        int aei2xyz (double ele[6], double gm, 
                    double pos[3], double vel[3]);

        double kepler (double M,double e);
        
        double xyz2aei (double pos[3], double vel[3], double gm,
                        double ele[6]);

        void xyz2llr (double *vt, double *llr);
        
        void llr2xyz (double *llr, double *vt);

        double chosephase (double sinvalue, double cosvalue);

        void xyz2rtn (double *x, double *v, double *xyz, double *rtn);

------------------------------------------------------------------------
*/

#ifndef _COORD_H_

#include "coord.h"
    #define _COORD_H_
#endif





/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/

int aei2xyz (double ele[6], double gm, double pos[3], double vel[3])
{
    double a, e, i, omega, w, M, E,r, P[3], Q[3], n;
    int x;

    a = ele[0]; 
    e = ele[1];
    i = ele[2] * DEG2RAD;
    omega = ele[3] * DEG2RAD;
    w = ele[4] * DEG2RAD;
    M = ele[5] * DEG2RAD;
    
    if (a <= 0)
    {
        printf("error: a = %e <= 0 ! in aei2xyz \n", a);
        pos[0] = 0; pos[1] = 0; pos[2] = 0;
        vel[0] = 0; vel[1] = 0; vel[2] = 0;
        return 1;
    }

    n=sqrt(gm/(a*a*a));

    E=kepler(M,e);
    P[0]=cos(omega)*cos(w)-sin(omega)*sin(w)*cos(i);
    P[1]=sin(omega)*cos(w)+cos(omega)*sin(w)*cos(i);
    P[2]=sin(w)*sin(i);
    Q[0]=-cos(omega)*sin(w)-sin(omega)*cos(w)*cos(i);
    Q[1]=-sin(omega)*sin(w)+cos(omega)*cos(w)*cos(i);
    Q[2]=cos(w)*sin(i);
    for(x=0;x<3;x++)
    {
        pos[x]=a*(cos(E)-e)*P[x]+a*sqrt(1-e*e)*sin(E)*Q[x];
    }
    r = modvect (pos);

    if (r <= 0)
    {
        printf("error: r = %e <= 0 ! in aei2xyz \n", r);
        pos[0] = 0; pos[1] = 0; pos[2] = 0;
        vel[0] = 0; vel[1] = 0; vel[2] = 0;
        return 1;
    }

    for(x=0;x<3;x++)
    {
        vel[x]=-a*a*n/r*sin(E)*P[x]+a*a*n/r*sqrt(1-e*e)*cos(E)*Q[x];
    }
        
    return 0;
}





/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/

double kepler(double M,double e)
{
    double E0,E1=M;
    do
    {
        E0=E1;
        E1=M+e*sin(E0);
    }
    while(fabs(E0-E1)>=1e-10);
    return(E1);
}




/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/

double xyz2aei(double pos[3], double vel[3], double gm, double ele[6])
{
    double a, e, omega, i, w, E, M, r, v, h, HV[3], n, 
         Pz, Qz;
    
    r = modvect (pos);
    v = modvect (vel);
    
    if (r <= 0)
    {
        ele[0] = 0; ele[1] = 0; ele[2] = 0;
        ele[3] = 0; ele[4] = 0; ele[5] = 0;
        printf("error: r = %e <= 0 ! in xyz2aei \n", r);
        return 1;
    }

    a = 1.0 / (2.0 / r - v * v / gm);

    crsvect (pos, vel, HV);
    h = modvect(HV);
    e = sqrt (1.0 - h * h / gm / a);
    i = acos (HV[2] / h);     //unit: rad
    omega = chosephase (HV[0] / h / sin(i), 
            - HV[1] / h / sin(i));   //unit: rad

    if(a <= 0)
    {    
        ele[0] = 0; ele[1] = 0; ele[2] = 0;
        ele[3] = 0; ele[4] = 0; ele[5] = 0;
        printf("error: a = %e <= 0 ! in xyz2aei \n", a);
        return 1;
    }
    

    n = sqrt ( gm / (a*a*a) );

    E = chosephase ( dotvect(pos, vel) / (a * a * n * e), 
                    (1.0 - r / a) / e);  //unit: rad

    M = E - e * sin(E);        //unit: rad

    Pz = (cos(E) / r * pos[2] - sin(E) / n / a * vel[2]);
    Qz = (sin(E) / r / sqrt(1.0-e*e) * pos[2] 
         + (cos(E) - e) / n / a / sqrt(1.0-e*e) * vel[2]);

    w = chosephase ( Pz / sin(i), Qz /sin(i));  //unit: rad

    ele[0] = a;
    ele[1] = e;
    ele[2] = i * RAD2DEG;
    ele[3] = omega * RAD2DEG;
    ele[4] = w * RAD2DEG;
    ele[5] = M * RAD2DEG;
                          
    return TWOPI / n;                              
}
                


/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/


void xyz2llr (double *vt, double *llr)
{
    double r, cosphi, phi, costhe, sinthe;
    r = sqrt (vt[0] * vt[0] + vt[1] * vt[1] + vt[2] * vt[2]);

    cosphi = vt[2] / r;
    phi = acos(cosphi) ;
    costhe = vt[0] / r / sin(phi);
    sinthe = vt[1] / r / sin(phi);
    llr[2] = r;
    llr[1] = chosephase(sinthe, costhe) * RAD2DEG;
    llr[0] = 90.0 - phi * RAD2DEG;
}


void llr2xyz (double *llr, double *vt)
{
    double r, lon, lat;
    lat = llr[0] * DEG2RAD;
    lon = llr[1] * DEG2RAD;
    r = llr[2];
    vt[0] = r * cos(lat) * cos(lon);
    vt[1] = r * cos(lat) * sin(lon);
    vt[2] = r * sin(lat);

}







/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/

double chosephase (double sinvalue, double cosvalue)
{
    double sv = sinvalue, cv = cosvalue;
    if (sv >= 0 && cv >= 0) 
        return (asin (sv));
    if (sv > 0 && cv < 0) 
        return (acos (cv));
    if (sv < 0 && cv < 0) 
        return ( - asin (sv) + TWOPI / 2.0);
    else 
        return (asin (sv) + TWOPI);
}





/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/

void xyz2rtn(double *x, double *v, double *xyz, double *rtn)
{ 
    double scal_x, scal_v, vr[3], vn[3], vt[3];

    scal_x=sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
    scal_v=sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);

// c...unit vector in R direction
    vr[0]=x[0]/scal_x;
    vr[1]=x[1]/scal_x;
    vr[2]=x[2]/scal_x;
 
// c...unit direction in N direction
    vn[0]=(vr[1]*v[2]-vr[2]*v[1])/scal_v;
    vn[1]=(vr[2]*v[0]-vr[0]*v[2])/scal_v;
    vn[2]=(vr[0]*v[1]-vr[1]*v[0])/scal_v;
        
// c...unit direction in T direction
    vt[0]=(vn[1]*vr[2]-vn[2]*vr[1]);
    vt[1]=(vn[2]*vr[0]-vn[0]*vr[2]);
    vt[2]=(vn[0]*vr[1]-vn[1]*vr[0]);


    rtn[0]=xyz[0]*vr[0]+xyz[1]*vr[1]+xyz[2]*vr[2];
    rtn[1]=xyz[0]*vt[0]+xyz[1]*vt[1]+xyz[2]*vt[2];
    rtn[2]=xyz[0]*vn[0]+xyz[1]*vn[1]+xyz[2]*vn[2];

    return;
}


