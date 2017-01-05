/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/*

  GRACE L1B to L1C & L2

  Version: 20 Dec 2011 
           fixed and stochastic constraint solution

  Copyright (c) 2011 Kun Shang (shang.34@osu.edu) All Right Reserved 

*/
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/

#ifndef _L1B2L1C_H_
    #include "l1b2l1c.h"
#endif


#define MAXLINE 500








double accgr (double *tjd, double *p, double *v, double *fgr)
{
    double GME, GMS, J, Jv[3], beta, gamma, r, v2, pv, pJ, a, b,
        pxv[3], vxJ[3], ps[3], vs[3], rs, vsxps[3], vsxpsxv[3],
        term1[3], term2[3], term3[3];
    int n;

    GME = GMA[0]; //m^3/s^2
    r = modvect(p);

    J = 9.8e8; //m^2/s
    gamma = 1;
    beta = 1;
    GMS = 1.32712442076e20; //m^3/s^2

    Jv[0] = 0; Jv[1] = 0; Jv[2] = J;

    v2 = dotvect(v, v);
    pv = dotvect(p, v);
    pJ = dotvect(p, Jv);
    crsvect(p, v, pxv);
    crsvect(v, Jv, vxJ);

    planet_ephemeris (tjd, 2, 10, ps, vs);
    for (n = 0; n < 3; n++)
    {
        ps[n] = ps[n] * AU;
        vs[n] = vs[n] * AU / 86400.0;
    }

    rs = modvect(ps);
    crsvect(vs, ps, vsxps);
    crsvect(vsxps, v, vsxpsxv);

    a = 2 * (beta + gamma) * GME / r - gamma * v2;
    b = 2 * (1 + gamma) * pv;
    for (n = 0; n < 3; n++)
    {
        term1[n] = GME / C / C / r / r / r *
                ( a * p[n] + b * v[n] );
        term2[n] = GME / C / C / r / r / r * (1 + gamma) *
                ( 3/r/r * pxv[n] * pJ + vxJ[n] );
        term3[n] = - GMS / C / C / rs / rs / rs * (1 + 2 * gamma) *
                vsxpsxv[n];

        fgr[n] = term1[n]
                + term2[n] + term3[n];
    }


    return 0;


}








double kbrvel2pos (int i)
{
    int n, iter;
    double rou, rp[3], rv[3], rvp[3], rvpv[3], pos, vel, rpnew[3], ev[3],evt[3];


    for (n = 0; n < 3; n ++)
    {
        rp[n] = sat12[i].rp[n];
        rv[n] = sat12[i].rv[n];
    }

    iter = 0;
    while(1)
    {
        iter ++;            
        crsvect(rv, rp, rvp);        
        crsvect(rvp, rv, rvpv);
        pos = modvect(rp);
        vel = modvect(rv);

        if (kbrx[i].rate != 0)
            rou = sat12[i].rate * pos / vel;
        else
            rou = dotvect(rp, rv) / vel;

        for (n = 0; n < 3; n++)
        {
            evt[n] = rvpv[n] / modvect(rvpv);
            ev[n] = rv[n] / modvect(rv);
        }

        for (n = 0; n < 3; n ++)
        {
            rpnew[n] = rou * ev[n] + sqrt(pos * pos - rou * rou) * evt[n];
        }

//        printf ("iter = %d rx = %f ry = %f rz = %f dx = %f dy = %f dz = %f\n",
//            iter, rp[0], rp[1], rp[2], rp[0] - rpnew[0], rp[1] - rpnew[1], rp[2] - rpnew[2]);
            
        if (iter > 1)  break;

        for (n = 0; n < 3; n ++)
        {
            rp[n] = rpnew[n];
        }
    }
    
    
    for (n = 0; n < 3; n ++)
    {
        sat12[i].rp[n] = rpnew[n];
        sat2b[i].rp[n] = sat1a[i].rp[n] + sat12[i].rp[n];
    }
    return 0;

}




double kbrpos2vel (int i)
{
    int n, iter;
    double rou, vel, pos, rvnew[3], ep[3], ept[3], rp[3], rv[3], rpv[3], rpvp[3];

    rou = sat12[i].rate;
    for (n = 0; n < 3; n ++)
    {
        rp[n] = sat12[i].rp[n];
        rv[n] = sat12[i].rv[n];
    }

    iter = 0;
    while(1)
    {
        iter ++;            
        crsvect(rp, rv, rpv);        
        crsvect(rpv, rp, rpvp);
//        vel = modvect(rv);
        vel = sat12[i].vel;
        pos = modvect(rp);

        for (n = 0; n < 3; n++)
        {
            ept[n] = rpvp[n] / modvect(rpvp);
            ep[n] = rp[n] / modvect(rp);
        }

        for (n = 0; n < 3; n ++)
        {
            rvnew[n] = rou * ep[n] + sqrt(vel * vel - rou * rou) * ept[n];
        }                
//        printf ("iter = %d vx = %f vy = %f vz = %f dvx = %e dvy = %e dvz = %e\n",
//            iter, rv[0], rv[1], rv[2], rv[0] - rvnew[0], rv[1] - rvnew[1], rv[2] - rvnew[2]);
            
        if (iter > 1)  break;

        for (n = 0; n < 3; n ++)
        {
            rv[n] = rvnew[n];
        }
    }
    
    
    for (n = 0; n < 3; n ++)
    {
        sat12[i].rv[n] = rvnew[n];
        sat2b[i].rv[n] = sat1a[i].rv[n] + sat12[i].rv[n];
    }
    return 0;

}



















double stidev (int num, double *llr1, double *vs, double *dvts, double *as)
{
    double gp1, *stcs, dv1, acc1[3], p1e[3], p1[3], p1e_a[3], p1e_b[3], 
        llr1_b[3], llr1_a[3], *stcs_b, gp1_b, *stcs_a, gp1_a;
    int num_b, num_a, nmax;

    nmax = 4;
        
    stcs = (double *) calloc ( (nmax + 1) * (nmax + 1), sizeof(double));

//    id_perm = PERM;
//    stidecs(num, 1, stcs);
//    stidecs_Anelastic(num, 1, stcs);
            
    stidecs_earth(&info[num], stcs, C20PERM, 1,1,0);
//    stidecs_Anelastic(&info[num], 1, stcs);

    cs2acc (num, llr1, stcs, GMA[0], GMA[1], nmax, &gp1, &dv1, acc1);

    *vs = gp1;

    as[0] = acc1[0];
    as[1] = acc1[1];
    as[2] = acc1[2];
   

    if (num != 0)
    {
        num_b = num - 1;
    }
    else
    {
        num_b = 0;
    }

    if (num != NDATA - 1)
    {
        num_a = num + 1;
    }
    else 
    {
        num_a = NDATA - 1;
    }
    stcs_b = (double *) calloc ( (nmax + 1) * (nmax + 1), sizeof(double));
    stcs_a = (double *) calloc ( (nmax + 1) * (nmax + 1), sizeof(double));

//    stidecs(num_b, id_perm, stcs_b);
//    stidecs(num_a, id_perm, stcs_a);

//    stidecs_Anelastic(num_b, 1, stcs_b);
//    stidecs_Anelastic(num_a, 1, stcs_a);


    stidecs_earth(&info[num_b], stcs_b, C20PERM, 1,1,0);
    stidecs_earth(&info[num_a], stcs_a, C20PERM, 1,1,0);
//    stidecs_Anelastic(&info[num_b], 1, stcs_b);
//    stidecs_Anelastic(&info[num_a], 1, stcs_a);

    llr2xyz (llr1, p1e);
    brmul(info[num].c_ei, p1e, 3, 3, 1, p1);
    brmul(info[num_b].c_ie, p1, 3, 3, 1, p1e_b);
    brmul(info[num_a].c_ie, p1, 3, 3, 1, p1e_a);
    xyz2llr(p1e_b, llr1_b);
    xyz2llr(p1e_a, llr1_a);


    cs2acc (num_b, llr1_b, stcs_b, GMA[0], GMA[1], nmax, &gp1_b, &dv1, acc1);
    cs2acc (num_a, llr1_a, stcs_a, GMA[0], GMA[1], nmax, &gp1_a, &dv1, acc1);
    

    *dvts = (gp1_a - gp1_b)/DT/2.0;

    free(stcs);
    free(stcs_b);
    free(stcs_a);
    return 0;
}









double getinfo_1day(char *infile)
{
    double jd0, tjd[2];
    int mjd0, i;

    jd0 = GPS_S / 86400.0 + T0;
    mjd0 = (int)(jd0 - 2400000.5);

    eop_open(infile, mjd0, mjd0 + 1);

    for (i = 0; i < NDATA; i++)
    {
        tjd[0] = jd0;
        tjd[1] = gps_GRACE2tt(jd0,sat12[i].t) / 86400.0;
//        tjd[0] = T0;
//        tjd[1] = gps2tt(sat12[i].t) / 86400.0;

        getinfo(tjd, 2, &info[i]);

    }

    eop_close();


    return 0;


}









double nbodyv (double *tjd, double *rp, double *vn, double *dvtn, double *an)
{
    double tjd_s[2], tjd_e[2], pte, pts, pt,  acc[3], dt;
    dt = 5.0;

    nbodypt (tjd, rp, &pt, acc);

    *vn = pt;
    an[0] = acc[0];
    an[1] = acc[1];
    an[2] = acc[2];


    tjd_s[0] = tjd[0];
    tjd_s[1] = tjd[1] - dt / 86400.0;
    tjd_e[0] = tjd[0];
    tjd_e[1] = tjd[1] + dt / 86400.0;

    nbodypt (tjd_s, rp, &pts, acc);
    nbodypt (tjd_e, rp, &pte, acc);


    *dvtn = (pte - pts) / dt / 2.0;

    return 0;
}















double nbodypt (double *tjd, double *pi, double *ptt, double *acc)
{
    double pj[3],vj[3], pij[3],f[3], rij, ri, rj, gm[11], 
        ptt1, ptt2, coszt, zt;
    short int earth, j, n;
    
    gm[0] =   2.203208082807623e+13;
    gm[1] =      3.248586038641429e+14;
    gm[2] =     398600.44180E+09;
    gm[3] =     4.28283719012840e+13;
    gm[4] =      1.267127698227696e+17;
    gm[5] =     3.794062664949063e+16;
    gm[6] =      5.794549096929744e+15;
    gm[7] =     6.836534169987595e+15;
    gm[8] =    9.816009029289940e+11;
    gm[9] =      4.902801056E+12;
    gm[10] =      1.32712442076e20;


    earth = 2;
    ptt1 = 0; ptt2 = 0; f[0] = 0; f[1] = 0; f[2] = 0;
    for (j = 0; j <= 10; j++)
    {
        if (j == earth)
           continue;
        planet_ephemeris (tjd, j, earth, pj, vj);
        for (n = 0; n < 3; n++)
        {
            pj[n]  = pj[n] * AU;
            pij[n] = pj[n]  - pi[n];
        }
        rij= sqrt(pij[0] * pij[0] + pij[1] * pij[1] + pij[2] * pij[2]);
        ri = sqrt(pi[0] * pi[0] + pi[1] * pi[1] + pi[2] * pi[2]);
        rj = sqrt(pj[0] * pj[0] + pj[1] * pj[1] + pj[2] * pj[2]);

        coszt = (pi[0] * pj[0] + pi[1] * pj[1] + pi[2] * pj[2]) / ri/ rj; 
        
        zt = acos(coszt);

        ptt1 = ptt1 + gm[j] * (1 / rij - 1 / rj - 1 / rj / rj * ri * coszt);
//        ptt1 = ptt1 + gm[j] * (1 / rij);

        ptt2 = ptt2 + gm[j] * ( 
            + pow(ri,2)/pow(rj,3) * ( 0.75*cos(2.0*zt) + 0.25)
            + pow(ri,3)/pow(rj,4) * ( 5.0/8.0*cos(3.0*zt) + 3.0/8.0*cos(zt) ) ); 
        for (n = 0; n < 3; n++)
            f[n] = f[n]
            + gm[j] / (rij * rij * rij) * pij[n] 
            - gm[j] / (rj * rj * rj) * pj[n];
        
    }
    

    *ptt = ptt1;
    
    for (n = 0; n < 3; n++)    
        acc[n] = f[n];


    return 0;
}

















/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/

double disse_a (double *v1, int i, double *ef1i, double *f1)
{
    double mat1[9], mtt1[9], qvec1[4], acc1[3];

    if (scaa[i].gps_time == 0 || scab[i].gps_time == 0)
    {
        *ef1i = 0;
        return 1;
    }

    qvec1[0] = scaa[i].quatangle;
    qvec1[1] = scaa[i].quaticoeff;
    qvec1[2] = scaa[i].quatjcoeff;
    qvec1[3] = scaa[i].quatkcoeff;
    acc1[0] = acca[i].lin_accl_x;
    acc1[1] = acca[i].lin_accl_y;
    acc1[2] = acca[i].lin_accl_z;

    quat2mat_i2s (qvec1, mat1);
//    brmul(mat1, acc1, 3,3,1, f1);
    mt(mat1, 3, 3, mtt1);
    brmul(mtt1, acc1, 3,3,1, f1);

    *ef1i = (f1[0] * v1[0] + f1[1] * v1[1] + f1[2] * v1[2]);

//    *ef1i = (f1[0] * v1[0] + f1[1] * v1[1] + f1[2] * v1[2]) * DT;

//    brmul(mat1, v1, 3, 3, 1, v1s);

//    *ef1i = (acc1[0] * v1s[0] + acc1[1] * v1s[1] + acc1[2] * v1s[2]) * DT;

    return 0;
}
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/



double disse_b (double *v1, int i, double *ef1i, double *f1)
{
    double mat1[9], mtt1[9], qvec1[4], acc1[3];

    if (scaa[i].gps_time == 0 || scab[i].gps_time == 0)
    {
        *ef1i = 0;
        return 1;
    }

    qvec1[0] = scab[i].quatangle;
    qvec1[1] = scab[i].quaticoeff;
    qvec1[2] = scab[i].quatjcoeff;
    qvec1[3] = scab[i].quatkcoeff;
    acc1[0] = accb[i].lin_accl_x;
    acc1[1] = accb[i].lin_accl_y;
    acc1[2] = accb[i].lin_accl_z;

    quat2mat_i2s (qvec1, mat1);
//    brmul(mat1, acc1, 3,3,1, f1);
    mt(mat1, 3, 3, mtt1);
    brmul(mtt1, acc1, 3,3,1, f1);

    *ef1i = (f1[0] * v1[0] + f1[1] * v1[1] + f1[2] * v1[2]);

//    *ef1i = (f1[0] * v1[0] + f1[1] * v1[1] + f1[2] * v1[2]) * DT;
          
//    brmul(mat1, v1, 3, 3, 1, v1s);

//    *ef1i = (acc1[0] * v1s[0] + acc1[1] * v1s[1] + acc1[2] * v1s[2]) * DT;

    return 0;
}
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/


double accia (double *tjd, int label, double *acc)
{
    double tt, as[3], qvec[4], c_is[9], c_si[9];

    tt = tjd[1] * 86400.0;
    
    if (label == 1)
    {
//        lagrangelow (ACA_EPH, DIM_ACA, 4, tt, as);
//        lagrangelow (SCA_EPH, DIM_SCA, 5, tt, qvec);
        lgr_order (ACA_EPH, DIM_ACA, 4, tt, as, 4);
        lgr_order (SCA_EPH, DIM_SCA, 5, tt, qvec, 4);


    }

    if (label == 2)
    {
//        lagrangelow (ACB_EPH, DIM_ACB, 4, tt, as);
//        lagrangelow (SCB_EPH, DIM_SCB, 5, tt, qvec);
        lgr_order (ACB_EPH, DIM_ACB, 4, tt, as, 4);
        lgr_order (SCB_EPH, DIM_SCB, 5, tt, qvec, 4);
    }


    quat2mat_i2s (qvec, c_is);
    mt(c_is, 3, 3, c_si);


    brmul(c_si, as, 3,3,1, acc);

    return 0;

}



/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
double quat2mat_i2s (double *qvec, double *mat)
/*
 * SCA_1B data file containing edited quaternion 
 * rotating inertial frame to SRF (5-sec data interval)
 * Page 17 in
 * Algorithm Theoretical Basis Document for GRACE Level-1B Data Processing V1.2
 * Sien-Chong Wu Gerhard Kruizinga Willy Bertiger
 * May 9, 2006
*/
{
    double q[4];
    
    q[0]=qvec[0];
    q[1]=qvec[1]; 
    q[2]=qvec[2]; 
    q[3]=qvec[3]; 

    mat[0] = q[0] * q[0] + q[1] * q[1] - q[2] * q[2] - q[3] * q[3];
    mat[1] = 2.0 * (q[1]*q[2] + q[3]*q[0]);
    mat[2] = 2.0 * (q[1]*q[3] - q[2]*q[0]);
    mat[3] = 2.0 * (q[1]*q[2] - q[3]*q[0]);
    mat[4] = q[0] * q[0] + q[2] * q[2] - q[1] * q[1] - q[3] * q[3];
    mat[5] = 2.0 * (q[2]*q[3] + q[0]*q[1]);
    mat[6] = 2.0 * (q[1]*q[3] + q[2]*q[0]);
    mat[7] = 2.0 * (q[2]*q[3] - q[0]*q[1]);
    mat[8] = q[0] * q[0] + q[3] * q[3] - q[1] * q[1] - q[2] * q[2];
 
    return 0;
}
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/









void fitl1c_mv (int nply, int ncpr, int nplycpr, int ndt)
{
    int np, np2, is, ie, i, n, j, nlsf;
    double t0, x, *xlsf, *ylsf, *prd, tpmean;
    
    t0 = sat12[0].t;

    np = (int)(ndt / DT);
    np2 = (int) (np/2);

    xlsf = (double *) calloc ( np, sizeof(double));
    ylsf = (double *) calloc ( np, sizeof(double));
    prd = (double *) calloc ( np, sizeof(double));


    for (i = 0; i < NDATA; i ++)
    {
        t0 = sat12[i].t;

        x = (sat12[i].t - t0) /  86400.0;
        is = i - np2;
        ie = i + np2 - 1;
        if (is < 0)
        {
            is = 0;
            ie = np - 1;
        }
        if (ie > NDATA - 1)
        {
            ie = NDATA - 1;
            is = NDATA - np;
        }

        n = 0;
        for (j = is; j <= ie; j ++)
        {
            if (sat12[j].error == 9) continue;

            xlsf[n] = (sat12[j].t - t0) / 86400.0;
            ylsf[n] = sat12[j].vl1c;
            ylsf[n] = ylsf[n] - sat12[j].vref;  ///////////////!!!!!!!!!!!//////////
//            prd[n] = (gnva[i].Tp + gnvb[i].Tp) / 2.0;
            prd[n] = (gnva[j].Tp + gnvb[j].Tp) / 2.0;
            n ++;
        }

        if (n < np2)
        {
            sat12[i].vmvfit = 0;
            sat12[i].vmvfitres = 999;
            continue;
        }


        nlsf = n;
        tpmean = mean (prd, nlsf);
        sat12[i].Tp2 = tpmean;

//        printf ("n = %d\t i = %d\t tpmean = %f\n", n, i, tpmean);
        sat12[i].vmvfitres = lsfl1c1d (&x, xlsf, ylsf, 1, nlsf, tpmean, &sat12[i].vmvfit, nply, ncpr, nplycpr);
//        sat12[i].ystd = fabs(sat12[i].res - sat12[i].y - (sat12[i].vcsr + sat12[i].vgfz + sat12[i].vjpl - 3 * sat12[i].vbl) / 3.0);
    
    }
    
    return;

}






void fitlos_mv (int nply, int ncpr, int nplycpr, int ndt)
{
    int np, np2, is, ie, i, n, j, nlsf;
    double t0, x, *xlsf, *ylsf, *prd, tpmean;
    
    t0 = sat12[0].t;

    np = (int)(ndt / DT);
    np2 = (int) (np/2);

    xlsf = (double *) calloc ( np, sizeof(double));
    ylsf = (double *) calloc ( np, sizeof(double));
    prd = (double *) calloc ( np, sizeof(double));


    for (i = 0; i < NDATA; i ++)
    {
        x = (sat12[i].t - t0) /  86400.0;
        is = i - np2;
        ie = i + np2 - 1;
        if (is < 0)
        {
            is = 0;
            ie = np - 1;
        }
        if (ie > NDATA - 1)
        {
            ie = NDATA - 1;
            is = NDATA - np;
        }

        n = 0;
        for (j = is; j <= ie; j ++)
        {
            if (sat12[j].error == 9) continue;

            xlsf[n] = (sat12[j].t - t0) / 86400.0;
            ylsf[n] = sat12[j].glos;
            ylsf[n] = ylsf[n] - sat12[j].gref;  ///////////////!!!!!!!!!!!//////////
//            prd[n] = (gnva[i].Tp + gnvb[i].Tp) / 2.0;
            prd[n] = (gnva[j].Tp + gnvb[j].Tp) / 2.0;
            n ++;
        }

        if (n < np2)
        {
            sat12[i].gmvfit = 0;
            sat12[i].gmvfitres = 999;
            continue;
        }

        nlsf = n;
        tpmean = mean (prd, nlsf);

        sat12[i].gmvfitres = lsfl1c1d (&x, xlsf, ylsf, 1, nlsf, tpmean, &sat12[i].gmvfit, nply, ncpr, nplycpr);
//        sat12[i].ystd = fabs(sat12[i].res - sat12[i].y - (sat12[i].vcsr + sat12[i].vgfz + sat12[i].vjpl - 3 * sat12[i].vbl) / 3.0);
    
    }
    
    return;

}











void fitl1c_pw (int order_poly, int order_cpr)
{
    int i, n, in, nmax;
    double *xn, *xi, *res, *resm, *prd, *fit, *fitm, tpmean, t0, stdres, stdresm;

    nmax = 5400/DT;
//    nmax = 3600/DT;
            
    xn  = (double *) calloc ( nmax, sizeof(double));
    res = (double *) calloc ( nmax, sizeof(double));
    resm = (double *) calloc ( nmax, sizeof(double));
    prd = (double *) calloc ( nmax, sizeof(double));
    fit = (double *) calloc ( nmax, sizeof(double));
    fitm = (double *) calloc ( nmax, sizeof(double));
    xi  = (double *) calloc ( nmax, sizeof(double));
            
    n = 0;
    t0 = sat12[0].t;
    for (i = 0; i < NDATA; i ++)
    {
//        t0 = sat12[i].t;

        if (sat12[i].error == 0)
        {
            xn[n] = (sat12[i].t - t0) / 86400.0;
            prd[n] = (gnva[i].Tp + gnvb[i].Tp) / 2.0;
            res[n] = sat12[i].vl1c;
            res[n] = res[n] - sat12[i].vref;  ///////////////!!!!!!!!!!!//////////
            resm[n] = sat12[i].vl1cm;
            n++;
        }
        xi[(i)%nmax] = (sat12[i].t - t0) / 86400.0;

        if (i == 0) continue;
        if (((i+1)%(5400/DT)) == 0)   // possible useful for GPS orbit (orbit cut at 86399s)
//        if (((i)%(5400/DT)) == 0)
        {
            if (n < nmax / 2)
            {
                for (in = 0; in < nmax; in ++)
                {
                    sat12[i - (nmax-1) + in].vpwfit = 0;  
//                    sat12[i - (nmax-1) + in].error = 3; 
                }
                t0 = sat12[i+1].t;
                n = 0;
                continue;
            }
            tpmean = mean (prd, n);
//            tpmean = 5633.5;
//            printf ("n = %d\t i = %d\t tpmean = %f\n", n, i, tpmean);

            stdres = lsf_cpr_new (xi, xn, res, nmax, n, tpmean, fit, order_poly, order_cpr);
            stdresm =lsf_cpr_new (xi, xn, resm, nmax, n, tpmean, fitm, order_poly, order_cpr);
//            printf ("n = %d\t i = %d\t tpmean = %f\t, stdres = %f\n", n, i, tpmean, stdres);
    
            for (in = 0; in < nmax; in ++)
            {
                sat12[i - (nmax-1) + in].vpwfit = fit[in];   // possible useful for GPS orbit
                sat12[i - (nmax-1) + in].vpwfitm = fitm[in];   // possible useful for GPS orbit
                sat12[i - (nmax-1) + in].vpwfitres = stdres;
                sat12[i - (nmax-1) + in].Tp1 = tpmean;
//                if (stdres > 0.003)
//                    sat12[i - (nmax-1) + in].error = 2;
//                sat12[i - nmax + in].fit = fit[in];
            }

            t0 = sat12[i+1].t;
            n = 0;
        }

    }
        
    free(xi);
    free(xn);
    free(res);
    free(resm);
    free(prd);
    free(fit);
    free(fitm);

}







void fitlos_pw (int order_poly, int order_cpr)
{
    int i, n, in, nmax;
    double *xn, *xi, *res, *resm, *prd, *fit, *fitm, tpmean, t0, stdres;

    nmax = 5400/DT;
//    nmax = 3600/DT;
            
    xn  = (double *) calloc ( nmax, sizeof(double));
    res = (double *) calloc ( nmax, sizeof(double));
    resm = (double *) calloc ( nmax, sizeof(double));
    prd = (double *) calloc ( nmax, sizeof(double));
    fit = (double *) calloc ( nmax, sizeof(double));
    fitm = (double *) calloc ( nmax, sizeof(double));
    xi  = (double *) calloc ( nmax, sizeof(double));
            
    n = 0;
    t0 = sat12[0].t;
    for (i = 0; i < NDATA; i ++)
    {
        if (sat12[i].error == 0)
        {
            xn[n] = (sat12[i].t - t0) / 86400.0;
            prd[n] = (gnva[i].Tp + gnvb[i].Tp) / 2.0;
            res[n] = sat12[i].glos;
            res[n] = res[n] - sat12[i].gref;  ///////////////!!!!!!!!!!!//////////
//            resm[n] = sat12[i].resm;
            n++;
        }
        xi[(i)%nmax] = (sat12[i].t - t0) / 86400.0;

        if (i == 0) continue;
        if (((i+1)%(5400/DT)) == 0)   // possible useful for GPS orbit (orbit cut at 86399s)
//        if (((i)%(5400/DT)) == 0)
        {
            if (n < nmax / 2)
            {
                for (in = 0; in < nmax; in ++)
                {
                    sat12[i - (nmax-1) + in].gpwfit = 0;  
//                    sat12[i - (nmax-1) + in].error = 3; 
                }
                t0 = sat12[i+1].t;
                n = 0;
                continue;
            }
            tpmean = mean (prd, n);
//            printf ("n = %d\t i = %d\t tpmean = %f\n", n, i, tpmean);

            stdres = lsf_cpr_new (xi, xn, res, nmax, n, tpmean, fit, order_poly, order_cpr);
//            stdresm =lsf_cpr_new (xi, xn, resm, nmax, n, tpmean, fitm, order_poly, order_cpr);
//            printf ("n = %d\t i = %d\t tpmean = %f\t, stdres = %f\n", n, i, tpmean, stdres);
    
            for (in = 0; in < nmax; in ++)
            {
                sat12[i - (nmax-1) + in].gpwfit = fit[in];   // possible useful for GPS orbit
//                sat12[i - (nmax-1) + in].fitm = fitm[in];   // possible useful for GPS orbit
                sat12[i - (nmax-1) + in].gpwfitres = stdres;
//                if (stdres > 0.003)
//                    sat12[i - (nmax-1) + in].error = 2;
//                sat12[i - nmax + in].fit = fit[in];
            }

            t0 = sat12[i+1].t;
            n = 0;
        }

    }
        
    free(xi);
    free(xn);
    free(res);
    free(resm);
    free(prd);
    free(fit);
    free(fitm);

}










void fitres_new_overlap (int order_poly, int order_cpr, int overlap)
{
    int i, n, in, nmax;
    double *xn, *xi, *res, *resm, *prd, *fit, *fitm, tpmean, t0, stdres, stdresm;

    nmax = 5400/DT;
            
    xn  = (double *) calloc ( nmax, sizeof(double));
    res = (double *) calloc ( nmax, sizeof(double));
    resm = (double *) calloc ( nmax, sizeof(double));
    prd = (double *) calloc ( nmax, sizeof(double));
    fit = (double *) calloc ( nmax, sizeof(double));
    fitm = (double *) calloc ( nmax, sizeof(double));
    xi  = (double *) calloc ( nmax, sizeof(double));
            
    n = 0;
    t0 = sat12[0].t;
    for (i = 0; i < NDATA; i ++)
    {
        if (sat12[i].error == 0)
        {
            xn[n] = (sat12[i].t - t0) / 86400.0;
            prd[n] = (gnva[i].Tp + gnvb[i].Tp) / 2.0;
            res[n] = sat12[i].vl1c;
            resm[n] = sat12[i].vl1cm;
            n++;
        }
        xi[(i)%nmax] = (sat12[i].t - t0) / 86400.0;

        if (i == 0) continue;
        if (((i+1)%(5400/DT)) == 0)   // possible useful for GPS orbit (orbit cut at 86399s)
//        if (((i)%(5400/DT)) == 0)
        {
            if (n < nmax / 2)
            {
                for (in = 0; in < nmax; in ++)
                {
                    sat12[i - (nmax-1) + in].vpwfit = 0;  
                    sat12[i - (nmax-1) + in].error = 3; 
                }
                t0 = sat12[i+1].t;
                n = 0;
                continue;
            }
            tpmean = mean (prd, n);
//            printf ("n = %d\t i = %d\t tpmean = %f\n", n, i, tpmean);

            stdres = lsf_cpr_new (xi, xn, res, nmax, n, tpmean, fit, order_poly, order_cpr);
            stdresm =lsf_cpr_new (xi, xn, resm, nmax, n, tpmean, fitm, order_poly, order_cpr);
//            printf ("n = %d\t i = %d\t tpmean = %f\t, stdres = %f\n", n, i, tpmean, stdres);
    
            for (in = 0; in < nmax; in ++)
            {
                sat12[i - (nmax-1) + in].vpwfit = fit[in];   // possible useful for GPS orbit
                sat12[i - (nmax-1) + in].vpwfitm = fitm[in];   // possible useful for GPS orbit
                if (stdres > 0.003)
                    sat12[i - (nmax-1) + in].error = 2;
//                sat12[i - nmax + in].fit = fit[in];
            }

            t0 = sat12[i+1].t;
            n = 0;
        }

    }
        
    free(xi);
    free(xn);
    free(res);
    free(resm);
    free(prd);
    free(fit);
    free(fitm);

}




















void fitres (int order_poly, int order_cpr, int offsetcyc)
{
    int i, is, ie, n, cnum, ofs, ofe, offset;
    double *xt, *res, *prd, *fit, tpmean;


//    offset = 5400/DT;
    offset = offsetcyc;
    is = 0;
    for (i = 0; i < NDATA; i ++)
    {
        if (i == 0 || gnva[i].r == 0) continue;
//        if ( ((gnva[i+1].lat > gnva[i].lat) && (gnva[i-1].lat > gnva[i].lat)) 
        if ( ((i%(5400/DT)) == 0)
            || i == NDATA - 1)
        {
            ie = i + 1;
            cnum = ie - is;

            if (is - offset < 0)
            {
                ofs = 0;
                ofe = offset * 2;    
            }
            else if (ie + offset >=NDATA)
            {
                ofs = offset * 2;
                ofe = 0;
            }
            else 
            {
                ofs = offset;
                ofe = offset;
            }
            
            xt = (double *) calloc ( cnum + offset * 2, sizeof(double));
            res = (double *) calloc ( cnum + offset * 2, sizeof(double));
            prd = (double *) calloc ( cnum + offset * 2, sizeof(double));
            fit = (double *) calloc ( cnum + offset * 2, sizeof(double));
            
            for (n = 0; n < cnum + offset * 2; n ++)
            {
                xt[n] = (gnva[is - ofs + n].gps_time - gnva[is].gps_time) / 86400.0;
                prd[n] = gnva[is - ofs + n].Tp;
                res[n] = sat12[is - ofs + n].vl1c;
            }
            tpmean = mean (prd, cnum + offset * 2);
//            tpmean = 5633.13;
//            printf ("cnum = %d\t is = %d\t ie = %d\t tpmean = %f\n", cnum, is, ie, tpmean);

            lsf_cpr (xt, res, cnum + offset * 2, tpmean, fit, order_poly, order_cpr);
//            lsf_cpr_day (xt, res, cnum + offset * 2, tpmean, fit, order_poly, order_cpr);

            for (n = 0; n < cnum; n ++)
            {
                sat12[is + n].vpwfit = fit[n + ofs];
            }
            free(xt);
            free(res);
            free(prd);
            free(fit);
            is = ie;
        }
    }
        

}

















/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
int calbias(char *infile_a, char *infile_b)
{

    double x0a, y0a, z0a, x0b, y0b, z0b;
    FILE *fpACC_a, *fpACC_b;
    int i;
    char line[MAXLINE];

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
    if ( (fpACC_a = fopen (infile_a,"r")) == NULL)
    {
        printf ("Cannot open fpin_a file!\n");
        //getch();
        exit (0);
    }
    if ( (fpACC_b = fopen (infile_b,"r")) == NULL)
    {
        printf ("Cannot open fpin_b file!\n");
        //getch();
        exit (0);
    }


    fgets(line, 100, fpACC_a);
    sscanf (line, "%lf", &x0a);
    fgets(line, 100, fpACC_a);
    sscanf (line, "%lf", &y0a);
    fgets(line, 100, fpACC_a);
    sscanf (line, "%lf", &z0a);


    fgets(line, 100, fpACC_b);
    sscanf (line, "%lf", &x0b);
    fgets(line, 100, fpACC_b);
    sscanf (line, "%lf", &y0b);
    fgets(line, 100, fpACC_b);
    sscanf (line, "%lf", &z0b);


//    printf ("%e\n%e\n%e\n%e\n%e\n%e\n", x0a, y0a, z0a, x0b, y0b, z0b);


    for (i = 0; i < DIM_ACA; i++)
    {
        ACA_EPH[i * 4 + 1] = ACA_EPH[i * 4 + 1] + x0a;
        ACA_EPH[i * 4 + 2] = ACA_EPH[i * 4 + 2] + y0a;
        ACA_EPH[i * 4 + 3] = ACA_EPH[i * 4 + 3] + z0a;
    }
    

    for (i = 0; i < DIM_ACB; i++)
    {
        ACB_EPH[i * 4 + 1] = ACB_EPH[i * 4 + 1] + x0b;
        ACB_EPH[i * 4 + 2] = ACB_EPH[i * 4 + 2] + y0b;
        ACB_EPH[i * 4 + 3] = ACB_EPH[i * 4 + 3] + z0b;
    }
   
    return 0;

}
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/










int calbiaseph(char *infile_a, char *infile_b)
{

    double *MACC_EPBS = NULL, *MACC_EPSL = NULL, tt;
    FILE *fpACC_A, *fpACC_B;
    int i,k,m,n, lbs, lsl, MACC_ARBS = 0, MACC_ARSL = 0, 
        MACC_NOBS = 0, MACC_NOSL = 0, MACC_PRBS = 0, MACC_PRSL = 0;
//    char line[MAXLINE];

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
    if ( (fpACC_A = fopen (infile_a,"r")) == NULL)
    {
        printf ("Cannot open fpin_a file!\n");
        //getch();
        exit (0);
    }
    if ( (fpACC_B = fopen (infile_b,"r")) == NULL)
    {
        printf ("Cannot open fpin_b file!\n");
        //getch();
        exit (0);
    }

    if (MACC_BIAS != 0)
    {
        MACC_ARBS = (int)(86400 / MACC_DTBS);
        MACC_NOBS = 3 * MACC_BIAS;
        MACC_PRBS = MACC_NOBS * MACC_ARBS;
        MACC_EPBS = (double *) calloc (MACC_PRBS, sizeof(double));
    }

    if (MACC_SCAL != 0)
    {
        MACC_ARSL = (int)(86400 / MACC_DTSL);
        MACC_NOSL = 3 * MACC_SCAL;
        MACC_PRSL = MACC_NOSL * MACC_ARSL;
        MACC_EPSL = (double *) calloc (MACC_PRSL, sizeof(double));
    }

/*A*/   
 
    if (MACC_BIAS != 0)
    {
        for (i = 0; i < MACC_ARBS; i ++)
        {
            fscanf (fpACC_A, "%*s%*d");
            for (n = 0; n < MACC_NOBS; n++)
            {
                fscanf(fpACC_A, "%lf", &MACC_EPBS[i * MACC_NOBS + n]);
            }
        }
    }

    if (MACC_SCAL != 0)
    {
        for (i = 0; i < MACC_ARSL; i ++)
        {
            fscanf (fpACC_A, "%*s%*d");
            for (n = 0; n < MACC_NOSL; n++)
            {
                fscanf(fpACC_A, "%lf", &MACC_EPSL[i * MACC_NOSL + n]);
            }
        }
    }

    for (i = 0; i < DIM_ACA; i++)
    {
        tt = ACA_EPH[i * 4] / 86400.0;
//        tt = ( ACA_EPH[i * 4] - 51.184) / 86400.0;
    
        if (MACC_SCAL != 0)
        {
            lsl = (int)((ACA_EPH[i * 4] - 51.184)/MACC_DTSL);
            if (lsl < 0) lsl = 0;
            if (lsl > MACC_ARSL - 1) lsl = MACC_ARSL - 1;

            for (k = 0; k < 3; k ++)
            {
                ACA_EPH[i * 4 + k + 1] = ACA_EPH[i * 4 + k + 1] * MACC_EPSL[lsl * MACC_NOSL + k];
            }
        }   

        if (MACC_BIAS != 0)
        {
            lbs = (int)((ACA_EPH[i * 4] - 51.184)/MACC_DTBS);
            if (lbs < 0) lbs = 0;
            if (lbs > MACC_ARBS - 1) lbs = MACC_ARBS - 1;

            for (k = 0; k < 3; k ++)
            {
                for (m = 0; m < MACC_BIAS; m ++)
                {
                    ACA_EPH[i * 4 + k + 1] = ACA_EPH[i * 4 + k + 1] + MACC_EPBS[lbs * MACC_NOBS + k + 3 * m] * pow (tt, m);
                }
            }
        }
    }
    
/*B*/
    
    if (MACC_BIAS != 0)
    {
        for (i = 0; i < MACC_ARBS; i ++)
        {
            fscanf (fpACC_B, "%*s%*d");
            for (n = 0; n < MACC_NOBS; n++)
            {
                fscanf(fpACC_B, "%lf", &MACC_EPBS[i * MACC_NOBS + n]);
            }
        }
    }

    if (MACC_SCAL != 0)
    {
        for (i = 0; i < MACC_ARSL; i ++)
        {
            fscanf (fpACC_B, "%*s%*d");
            for (n = 0; n < MACC_NOSL; n++)
            {
                fscanf(fpACC_B, "%lf", &MACC_EPSL[i * MACC_NOSL + n]);
            }
        }
    }

    for (i = 0; i < DIM_ACB; i++)
    {
        tt = ACB_EPH[i * 4] / 86400.0;
    
        if (MACC_SCAL != 0)
        {
            lsl = (int)((ACB_EPH[i * 4] - 19 - 32.184)/MACC_DTSL);
            if (lsl < 0) lsl = 0;
            if (lsl > MACC_ARSL - 1) lsl = MACC_ARSL - 1;

            for (k = 0; k < 3; k ++)
            {
                ACB_EPH[i * 4 + k + 1] = ACB_EPH[i * 4 + k + 1] * MACC_EPSL[lsl * MACC_NOSL + k];
            }
        }   

        if (MACC_BIAS != 0)
        {
            lbs = (int)((ACB_EPH[i * 4] - 19 - 32.184)/MACC_DTBS);
            if (lbs < 0) lbs = 0;
            if (lbs > MACC_ARBS - 1) lbs = MACC_ARBS - 1;

            for (k = 0; k < 3; k ++)
            {
                for (m = 0; m < MACC_BIAS; m ++)
                {
                    ACB_EPH[i * 4 + k + 1] = ACB_EPH[i * 4 + k + 1] + MACC_EPBS[lbs * MACC_NOBS + k + 3 * m] * pow (tt, m);
                }
            }
        }
    }
    
   

    if (MACC_BIAS != 0)
        free (MACC_EPBS);
    if (MACC_SCAL != 0)
        free (MACC_EPSL);
    fclose (fpACC_A);
    fclose (fpACC_B);
    return 0;

}
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/






/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
double cal_acc_01(void)
{

    double scale_x_A, scale_y_A, scale_z_A, scale_x_B, scale_y_B, scale_z_B,
        bias_x_A, bias_y_A, bias_z_A, bias_x_B, bias_y_B, bias_z_B, 
        mjd, mjd0, mjdm, gps, gpss, mjds;
    int i;

    scale_x_A = 0.9595;
    scale_y_A = 0.9797;
    scale_z_A = 0.9485;
    scale_x_B = 0.9465;
    scale_y_B = 0.9842;
    scale_z_B = 0.9303;

    gpss = ACA_EPH[0] + GPS_S - 19 - 32.184;
    mjds = gpss / 86400.0 + T0 - 2400000.5;

    if (mjds > 55562) // Jan. 1, 2011
    {
        scale_x_A = scale_x_A * 0.98;
        scale_z_A = scale_z_A * 0.98;
        scale_x_B = scale_x_B * 0.98;
        scale_z_B = scale_z_B * 0.98;
    }


    mjdm = 52705; //March 7, 2003

    for (i = 0; i < DIM_ACA; i++)
    {
        gps = ACA_EPH[i * 4] + GPS_S - 19 - 32.184;
        mjd = gps / 86400.0 + T0 - 2400000.5;
        if ( mjd < mjdm)
        {
            mjd0 = 52532;
            bias_x_A = - 1.106 
                       + 2.233e-4 * ( mjd - mjd0 ) 
                       + 2.5e-7   * ( mjd - mjd0 ) * ( mjd - mjd0 );
            bias_y_A =   27.042 
                       + 4.46e-3  * ( mjd - mjd0 ) 
                       + 1.1e-6   * ( mjd - mjd0 ) * ( mjd - mjd0 );                       
            bias_z_A = - 0.5486 
                       - 1.139e-6 * ( mjd - mjd0 ) 
                       + 1.7e-7   * ( mjd - mjd0 ) * ( mjd - mjd0 );
        }
        else
        {
            mjd0 = 53736;
            bias_x_A = - 1.2095 
                       - 4.128e-5 * ( mjd - mjd0 ) 
                       + 9.7e-9   * ( mjd - mjd0 ) * ( mjd - mjd0 );
            bias_y_A =   29.3370 
                       + 6.515e-4 * ( mjd - mjd0 ) 
                       - 3.9e-7   * ( mjd - mjd0 ) * ( mjd - mjd0 );                       
            bias_z_A = - 0.5606 
                       - 2.352e-6 * ( mjd - mjd0 ) 
                       + 3.8e-9   * ( mjd - mjd0 ) * ( mjd - mjd0 );
        }

        ACA_EPH[i * 4 + 1] = bias_x_A * 1e-6 + scale_x_A * ACA_EPH[i * 4 + 1];
        ACA_EPH[i * 4 + 2] = bias_y_A * 1e-6 + scale_y_A * ACA_EPH[i * 4 + 2];
        ACA_EPH[i * 4 + 3] = bias_z_A * 1e-6 + scale_z_A * ACA_EPH[i * 4 + 3];

    }
    

    for (i = 0; i < DIM_ACB; i++)
    {
        gps = ACA_EPH[i * 4] + GPS_S - 19 - 32.184;
        mjd = gps / 86400.0 + T0 - 2400000.5;
        if ( mjd < mjdm)
        {
            mjd0 = 52532;
            bias_x_B = - 0.5647 
                       - 7.788e-5 * ( mjd - mjd0 ) 
                       + 2.4E-7   * ( mjd - mjd0 ) * ( mjd - mjd0 );
            bias_y_B =   7.5101 
                       + 7.495E-3 * ( mjd - mjd0 ) 
                       - 9.6E-6   * ( mjd - mjd0 ) * ( mjd - mjd0 );                       
            bias_z_B = - 0.8602 
                       + 1.399E-4 * ( mjd - mjd0 ) 
                       + 2.5E-7   * ( mjd - mjd0 ) * ( mjd - mjd0 );
        }
        else
        {
            mjd0 = 53736;
            bias_x_B = - 0.6049 
                       - 1.982E-5 * ( mjd - mjd0 ) 
                       + 3.5E-9   * ( mjd - mjd0 ) * ( mjd - mjd0 );
            bias_y_B =   10.6860 
                       + 1.159E-3 * ( mjd - mjd0 ) 
                       - 4.3E-7   * ( mjd - mjd0 ) * ( mjd - mjd0 );                       
            bias_z_B = - 0.7901 
                       + 4.783E-5 * ( mjd - mjd0 ) 
                       - 6.5E-9   * ( mjd - mjd0 ) * ( mjd - mjd0 );
        }

        ACB_EPH[i * 4 + 1] = bias_x_B * 1e-6 + scale_x_B * ACB_EPH[i * 4 + 1];
        ACB_EPH[i * 4 + 2] = bias_y_B * 1e-6 + scale_y_B * ACB_EPH[i * 4 + 2];
        ACB_EPH[i * 4 + 3] = bias_z_B * 1e-6 + scale_z_B * ACB_EPH[i * 4 + 3];

    }
   
    return 0;

}
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/








int readgnv (char *infile_a, char *infile_b, int *n_a, int *n_b)
{
    FILE *fpGNV_a, *fpGNV_b;
    int n_gnva, n_gnvb, i, gps_i;
    char line[MAXLINE];

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
    if ( (fpGNV_a = fopen (infile_a,"r")) == NULL)
    {
        printf ("Cannot open fpin_a file!\n");
        ////getch();
        exit (0);
    }
    if ( (fpGNV_b = fopen (infile_b,"r")) == NULL)
    {
        printf ("Cannot open fpin_b file!\n");
        ////getch();
        exit (0);
    }


    n_gnva = 0;
    while (1)
    {
        if (fgets(line,MAXLINE, fpGNV_a) ==NULL) break;
        n_gnva ++;
        sscanf (line, "%d", &gps_i);
//  printf ("%d \t %d \t %d \n",  GPS_S, gps_i, DT); 
        if ((gps_i - GPS_S) % DT == 0 && (gps_i - GPS_S) / DT < NDATA)
        {
            i = (gps_i - GPS_S) / DT;
            sscanf (line, "%d%s%s%lf%lf%lf%lf%lf%lf%lf", 
                &gnva[i].gps_time, &gnva[i].GRACE_id, &gnva[i].coord_ref, 
                &gnva[i].xpos, &gnva[i].ypos, &gnva[i].zpos,
                &gnva[i].xvel, &gnva[i].yvel, &gnva[i].zvel,
                &gnva[i].r);

            gnva[i].pos[0] = gnva[i].xpos;
            gnva[i].pos[1] = gnva[i].ypos;            
            gnva[i].pos[2] = gnva[i].zpos;
            gnva[i].vel[0] = gnva[i].xvel;
            gnva[i].vel[1] = gnva[i].yvel;            
            gnva[i].vel[2] = gnva[i].zvel;
            
            gnva[i].Tp = xyz2aei(gnva[i].pos, gnva[i].vel,GMA[0], gnva[i].ele);

        


//            printf ("%f\n",  gnva[i].Tp); 
        }

    }

    n_gnvb = 0;
    while (1)
    {
        if (fgets(line,MAXLINE, fpGNV_b) ==NULL) break;
        n_gnvb ++;
        sscanf (line, "%d", &gps_i);
        if ((gps_i - GPS_S) % DT == 0 && (gps_i - GPS_S) / DT < NDATA)
        {
            i = (gps_i - GPS_S) / DT;
            sscanf (line, "%d%s%s%lf%lf%lf%lf%lf%lf%lf", 
                &gnvb[i].gps_time, &gnvb[i].GRACE_id, &gnvb[i].coord_ref, 
                &gnvb[i].xpos, &gnvb[i].ypos, &gnvb[i].zpos,
                &gnvb[i].xvel, &gnvb[i].yvel, &gnvb[i].zvel,
                &gnvb[i].r);

            gnvb[i].pos[0] = gnvb[i].xpos;
            gnvb[i].pos[1] = gnvb[i].ypos;
            gnvb[i].pos[2] = gnvb[i].zpos;
            gnvb[i].vel[0] = gnvb[i].xvel;
            gnvb[i].vel[1] = gnvb[i].yvel;
            gnvb[i].vel[2] = gnvb[i].zvel;
            
            gnvb[i].Tp = xyz2aei(gnvb[i].pos, gnvb[i].vel,GMA[0], gnvb[i].ele);
//            printf ("%f\n",  gnvb[i].Tp); 
        }

    }

    fclose(fpGNV_a);
    fclose(fpGNV_b);
    *n_a = n_gnva;
    *n_b = n_gnvb;

    return 0;

}




int readkbr (char *infile, int *n)
{
    FILE *fpKBR;
    int i = 0, n_kbr, gps_i;
    char line[MAXLINE];    


/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
    if ( (fpKBR = fopen (infile,"r")) == NULL)
    {
        printf ("Cannot open fpKBR file!\n");
        ////getch();
        exit (0);
    }

    n_kbr = 0;
    while (1)
    {
        if (fgets(line,MAXLINE, fpKBR) ==NULL) break;
        n_kbr ++;
        sscanf (line, "%d", &gps_i);
        if ((gps_i - GPS_S) % DT == 0 && (gps_i - GPS_S) / DT < NDATA)
        {
            i = (gps_i - GPS_S) / DT;
            sscanf (line, "%d%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%d%d%d%d%d", 
                &kbrx[i].gps_time, &kbrx[i].biased_range, &kbrx[i].range_rate, 
                &kbrx[i].range_accl, &kbrx[i].iono_corr, 
                &kbrx[i].lighttime_corr, &kbrx[i].lighttime_rate, &kbrx[i].lighttime_accl,
                &kbrx[i].ant_centr_corr, &kbrx[i].ant_centr_rate, &kbrx[i].ant_centr_accl,
                &kbrx[i].K_A_SNR, &kbrx[i].Ka_A_SNR, &kbrx[i].K_B_SNR, &kbrx[i].Ka_B_SNR,
                &kbrx[i].qualflg);
        }
        kbrx[i].range = kbrx[i].biased_range + kbrx[i].lighttime_corr + kbrx[i].ant_centr_corr;
        kbrx[i].rate  = kbrx[i].range_rate   + kbrx[i].lighttime_rate + kbrx[i].ant_centr_rate;
        kbrx[i].accl  = kbrx[i].range_accl   + kbrx[i].lighttime_accl + kbrx[i].ant_centr_accl;

    }


    fclose(fpKBR);

    *n = n_kbr;
    return 0;

}



int readacc (char *infile_a, char *infile_b)
{
    FILE *fpACC_a, *fpACC_b;
    int i, gps;
    char line[MAXLINE];

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
    if ( (fpACC_a = fopen (infile_a,"r")) == NULL)
    {
        printf ("Cannot open fpin_a file!\n");
        //getch();
        exit (0);
    }
    if ( (fpACC_b = fopen (infile_b,"r")) == NULL)
    {
        printf ("Cannot open fpin_b file!\n");
        //getch();
        exit (0);
    }


    i = 0;
    while (1)
    {
        if (fgets(line,MAXLINE, fpACC_a) ==NULL) break;
//        sscanf (line, "%d%*s%lf%lf%lf",
        sscanf (line, "%d%lf%lf%lf",
            &gps, &ACA_EPH[i * 4 + 1], &ACA_EPH[i * 4 + 2], &ACA_EPH[i * 4 + 3]);
        ACA_EPH[i * 4] = gps - GPS_S + 19 + 32.184;
        i++;
    }

    DIM_ACA = i;

    i = 0;
    while (1)
    {
        if (fgets(line,MAXLINE, fpACC_b) ==NULL) break;
//        sscanf (line, "%d%*s%lf%lf%lf",
        sscanf (line, "%d%lf%lf%lf",
            &gps, &ACB_EPH[i * 4 + 1], &ACB_EPH[i * 4 + 2], &ACB_EPH[i * 4 + 3]);
        ACB_EPH[i * 4] = gps - GPS_S + 19 + 32.184;
        i++;
    }


    DIM_ACB = i;

    fclose(fpACC_a);
    fclose(fpACC_b);


    return 0;

}





/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
int readsca (char *infile_a, char *infile_b)
{
    FILE *fpSCA_a, *fpSCA_b;
    int i, gps;
    char line[MAXLINE];    


/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
    if ( (fpSCA_a = fopen (infile_a,"r")) == NULL)
    {
        printf ("Cannot open fpin_a file!\n");
        //getch();
        exit (0);
    }
    if ( (fpSCA_b = fopen (infile_b,"r")) == NULL)
    {
        printf ("Cannot open fpin_b file!\n");
        //getch();
        exit (0);
    }

    i = 0;
    while (1)
    {
        if (fgets(line,MAXLINE, fpSCA_a) ==NULL) break;
//        sscanf (line, "%d%*s%*d%lf%lf%lf%lf",
        sscanf (line, "%d%lf%lf%lf%lf",
            &gps, &SCA_EPH[i * 5 + 1], &SCA_EPH[i * 5 + 2],
            &SCA_EPH[i * 5 + 3], &SCA_EPH[i * 5 + 4]);
        SCA_EPH[i * 5] = gps - GPS_S + 19 + 32.184;
        i++;
    }

    DIM_SCA = i;


    i = 0;
    while (1)
    {
        if (fgets(line,MAXLINE, fpSCA_b) ==NULL) break;
//        sscanf (line, "%d%*s%*d%lf%lf%lf%lf",
        sscanf (line, "%d%lf%lf%lf%lf",
            &gps, &SCB_EPH[i * 5 + 1], &SCB_EPH[i * 5 + 2],
            &SCB_EPH[i * 5 + 3], &SCB_EPH[i * 5 + 4]);
        SCB_EPH[i * 5] = gps - GPS_S + 19 + 32.184;
        i++;
    }

    DIM_SCB = i;


    fclose(fpSCA_a);
    fclose(fpSCA_b);


    return 0;
}









/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/*
*/
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
double cs2acc (int num, double *llr, double *cs, double gm, double a, int nmax, 
               double *v, double *dvdt, double *acc)
{
    int n, m, k, l, ind;

    double a2r, sinf, cosf, sinlon, coslon, sincolat, coscolat, *cosml, *sinml, 
        *aprn, *pbar, *pbar1, *pbar2, accn[3], c_en[9], c_in[9],
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

//    #pragma omp parallel private(cosml, sinml, aprn, pbar, pbar1, pbar2, n, m, k, l, ind, sinf, cosf)
    cosml = (double *) calloc ( nmax + 1, sizeof(double)); //cos(m*lamta)
    sinml = (double *) calloc ( nmax + 1, sizeof(double)); //sin(m*lamta)
    aprn = (double *) calloc ( nmax + 1, sizeof(double));  //sin(m*lamta)
    pbar = (double *) calloc ( nmax + 1, sizeof(double));
    pbar1 = (double *) calloc ( nmax + 1, sizeof(double));
    pbar2 = (double *) calloc ( nmax + 1, sizeof(double));
    pt  = (double *) calloc ( (nmax + 1) * (nmax + 1), sizeof(double));
    ptt  = (double *) calloc ( (nmax + 1) * (nmax + 1), sizeof(double));

/*
    for (m = 0; m <= nmax; m++)
    {
        cosml[m] = cos(m * lon * DEG2RAD);
        sinml[m] = sin(m * lon * DEG2RAD);
    }
    
    for (n = 0; n <= nmax; n++)
    {
        aprn[n] = pow (a / r, n) * gm / r;
    }
*/

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



/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/

    t = coscolat; vi = 0; dvdlon = 0; dvdcolat = 0; dvdr = 0;
    for (m = 0; m <= nmax; m ++)
    {
        l = nmax - m + 1;
//        lgdr2(t, nmax, m, pbar, pbar1, pbar2);
        lgdr1i(t, nmax, m, pbar, pbar1);
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


//    mt(info[num].c_ie, 3, 3, info[num].c_ei);
    brmul (info[num].c_ei, c_en, 3, 3, 3, c_in);  //inertial to fixed matrix gmat = rmat*tbt

    brmul(c_in, accn, 3, 3, 1, acc);  //from fixed acc to inertial acc


    *v = vi;
//    *v = gm/r + gm/r * a/r * a/r * (sqrt(5.0) * (1.5 * t * t - 0.5)) * cs[2];
    *dvdt = - ANGVEL * dvdlon;
//    if (lon < 180) *dvdt = - ANGVEL * dvdlon - 3.08e-12 * dvdcolat;
//    if (lon >=180) *dvdt = - ANGVEL * dvdlon + 3.08e-12 * dvdcolat;

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
double cs2vdvdt (double *llr, double *cs, double gm, double a, int NMAX, double *v, double *dvdt, double *pt)
{
    int n, m, k, l, ind;

    double sinf, cosf, *cosml, *sinml, *aprn, *pbar, lat, lon, r, 
        vi, dvdti, t;


/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/

    lat = llr[0];
    lon = llr[1];
    r = llr[2];

    sinf = sin(lat * DEG2RAD);
    cosf = cos(lat * DEG2RAD);

//    #pragma omp parallel private(cosml, sinml, aprn, pbar, pbar1, pbar2, n, m, k, l, ind, sinf, cosf)
    cosml = (double *) calloc ( NMAX + 1, sizeof(double)); //cos(m*lamta)
    sinml = (double *) calloc ( NMAX + 1, sizeof(double)); //sin(m*lamta)
    aprn = (double *) calloc ( NMAX + 1, sizeof(double));  //sin(m*lamta)
    pbar = (double *) calloc ( NMAX + 1, sizeof(double));

    for (m = 0; m <= NMAX; m++)
    {
        cosml[m] = cos(m * lon * DEG2RAD);
        sinml[m] = sin(m * lon * DEG2RAD);
    }
    
    for (n = 0; n <= NMAX; n++)
    {
        aprn[n] = pow (a / r, n) * gm / r;
    }

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/

    t = sinf; vi = 0; dvdti = 0;
    for (m = 0; m <= NMAX; m ++)
    {
        l = NMAX - m + 1;
//        lgdr2(t, NMAX, m, pbar, pbar1, pbar2);
        lgdr(t, NMAX, m, pbar);
        for (k = 0; k < l; k++)
        {
            if (m==0)
            {
//                ind = 0;
                n = k + m;
                pt[k] = aprn[n] * pbar[k];
                vi = vi + pt[k] * cs[k];
            }
            else
            {
                ind = NMAX + 1 + (2 * NMAX - m + 2) * (m - 1);
                n = k + m;
                pt[ind + n - m] = aprn[n] * pbar[k] * cosml[m];
                pt[ind + n - m + l] = aprn[n] * pbar[k] * sinml[m];
                vi = vi + pt[ind + n - m] * cs[ind + n - m];
                vi = vi + pt[ind + n - m + l] * cs[ind + n - m + l];
                dvdti = dvdti + m * pt[ind + n - m + l] * cs[ind + n - m];
                dvdti = dvdti - m * pt[ind + n - m] * cs[ind + n - m + l];
            }
        }
    }



//! ATPA
//  gpti = 0;
//    for(k = 0; k < (NMAX + 1) * (NMAX + 1); k++)
//    {
//        gpti = gpti + pt[k] * cs[k];
//        pt[k] = pnmc[k];
//    }

    *v = vi;
    *dvdt = dvdti * ANGVEL;

    free (pbar);

    free (cosml);
    free (sinml);
    free (aprn);
    return 1;
}
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/





