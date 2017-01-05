/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/*

  GRACE L1B to L1C & L2

  Version: 20 Dec 2011 
           fixed and stochastic constraint solution

  Copyright (c) 2011 Kun Shang (shang.34@osu.edu) All Right Reserved 

*/
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/

#ifndef __OMP_H__
    #include <omp.h>
#endif

#include "l1b2l1c.h"
#include "grvts.h"
#include "coord.h"
#include "numrs.h"


/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
int main (int argc, char *argv[])
{
    FILE *fp_stdin, *fpout_sim, *fp_otbin;  
    int i, n, nmax, NOMAX, DIM_OT, ncut, mmax, mcut, 
        nply, ncpr, nplycpr, ndt, order_poly = 0, order_cpr = 0,
        n_gnva, n_gnvb, n_kbr, kbrr, kbra,
        solvefor, solvegrv, bias, 
        year, month, day, days;
    short int error, de_num;
    double *fkn, *fint, rpv[3], wi[3], we[3], tjd[2], 
        jd0, p1e[3], p2e[3], *coefcsr, 
        *coefjpl, *coefgfz, outlier, fltpar, endcut,
        *pt1, *pt2, *coefb, *coefbl, *coefrfl, *coefrfh, step_ot, *COEFO,
        jd_beg, jd_end;
    char line[500], card[20], f_gnv1c_a[200], f_gnv1c_b[200], 
        f_par_a[200], f_par_b[200], f_eph[200], 
        f_acc1b_a[200], f_acc1b_b[200], 
        f_sca1b_a[200], f_sca1b_b[200], f_kbr1b_x[200], f_aot[200],
        f_grv[2][200]={"\0", "\0"}, 
        f_ref[2][200]={"\0", "\0"}, 
        f_csr[2][200]={"\0", "\0"}, 
        f_jpl[2][200]={"\0", "\0"}, 
        f_gfz[2][200]={"\0", "\0"}, 
        f_eop[200], stdname[200];    
    time_t s0, s1, s2, s3, s4, s5;

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
    if (argc == 1)
    {
        printf ("input file name: ");
        scanf ("%s", stdname);
    }
    else if (argc == 2)
    {
        sscanf (argv[1], "%s", stdname);
    }
    else
    {
        printf ("input argument error!\n");
        exit (0);
    }
    if ( (fp_stdin = fopen (stdname,"r")) == NULL)
    {
        printf ("Cannot open stdin file!\n");
        exit (0);
    }

    while (feof(fp_stdin) == 0)
    {
        card[1] = '\0';
        if (fgets (line, 500, fp_stdin) == NULL) break;
        sscanf (line, "%s", card);  
        
        if (strcmp (card, "GNV1C") ==0)  
            sscanf (line, "%*s%s%s", f_gnv1c_a, f_gnv1c_b); 
        if (strcmp (card, "PAR1B") ==0)  
            sscanf (line, "%*s%s%s", f_par_a, f_par_b); 
        if (strcmp (card, "ACC1B") ==0)  
            sscanf (line, "%*s%s%s", f_acc1b_a, f_acc1b_b); 
        if (strcmp (card, "SCA1B") ==0)  
            sscanf (line, "%*s%s%s", f_sca1b_a, f_sca1b_b); 
        if (strcmp (card, "KBR1B") ==0)  
            sscanf (line, "%*s%s", f_kbr1b_x);  
        if (strcmp (card, "AOT") ==0)
            sscanf (line, "%*s%s%d%lf", f_aot, &NOMAX, &step_ot);
        if (strcmp (card, "STIDE") ==0)
            sscanf (line, "%*s%d%lf%d", &PERMT, &C20PERM, &STIDE);
        if (strcmp (card, "GRV") ==0)    
            sscanf (line, "%*s%s%s", f_grv[0], f_grv[1]);   
        if (strcmp (card, "REF") ==0)    
            sscanf (line, "%*s%s%s", f_ref[0], f_ref[1]);   
        if (strcmp (card, "L2MCSR") ==0)
            sscanf (line, "%*s%s%s", f_csr[0], f_csr[1]);
        if (strcmp (card, "L2MJPL") ==0)
            sscanf (line, "%*s%s%s", f_jpl[0], f_jpl[1]);
        if (strcmp (card, "L2MGFZ") ==0)
            sscanf (line, "%*s%s%s", f_gfz[0], f_gfz[1]);

        if (strcmp (card, "EPH") ==0)
            sscanf (line, "%*s%s", f_eph);
        if (strcmp (card, "NMAX") ==0)  
            sscanf (line, "%*s%d%d", &nmax, &mmax);
        if (strcmp (card, "NCUT") ==0)  
            sscanf (line, "%*s%d%d", &ncut, &mcut);
        if (strcmp (card, "DT") ==0)    
            sscanf (line, "%*s%d", &DT);
        if (strcmp (card, "YEAR") ==0)
            sscanf (line, "%*s%d", &year);
        if (strcmp (card, "MONTH") ==0)
            sscanf (line, "%*s%d", &month);
        if (strcmp (card, "DAY") ==0)
            sscanf (line, "%*s%d", &day);
        if (strcmp (card, "DAYS") ==0)
            sscanf (line, "%*s%d", &days);
        if (strcmp (card, "KBR") ==0)
            sscanf (line, "%*s%d%d", &kbrr, &kbra);
        if (strcmp (card, "FITRES") ==0)
        {
            sscanf (line, "%*s%d%d%d%d", &nply, &ncpr, &nplycpr, &ndt);
            order_poly = nply - 1;
            order_cpr = (int)(ncpr / 2);
        }

        if (strcmp (card, "FILTER") ==0)
            sscanf (line, "%*s%lf%lf", &fltpar, &endcut);
        if (strcmp (card, "OUTLIER") ==0)    
            sscanf (line, "%*s%lf", &outlier);
        if (strcmp (card, "EOP") ==0)    
            sscanf (line, "%*s%s", f_eop);  
        if (strcmp (card, "ACC_BIAS") ==0)
            sscanf (line, "%*s%d%d", &MACC_BIAS, &MACC_DTBS);
        if (strcmp (card, "ACC_SCAL") ==0)
            sscanf (line, "%*s%d%d", &MACC_SCAL, &MACC_DTSL);

    }


    jd0 = julian_date (year,month,day,0);
    NDATA = days * 86400 / DT;
    GPS_S = (int)((jd0 - T0) * 86400 + 0.5);
//    mjd0 = (int)(jd0 - 2400000.5);
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
// open output files 



    if ( (fpout_sim = fopen ("sim.asc","w")) == NULL)
    {
        printf ("Cannot open fpout_sim file!\n");
        exit (0);
    }


/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
// open ephemeris file

    if ((error = ephem_open (f_eph, &jd_beg,&jd_end,&de_num)) != 0)
    {
      if (error == 1)
         printf ("JPL ephemeris file not found.\n");
       else
         printf ("Error reading JPL ephemeris file header.\n");
      return (error);
    }

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
// allocate memory

    solvegrv = (nmax + 1) * (nmax + 1);
    solvefor = (nmax + 1) * (nmax + 1) + bias;

    kbrx = (KBR1B *) calloc ( NDATA, sizeof(KBR1B));
    gnva = (GNV1C *) calloc ( NDATA, sizeof(GNV1C));
    gnvb = (GNV1C *) calloc ( NDATA, sizeof(GNV1C));
    ACA_EPH  = (double *) calloc (86400 * 4 * days, sizeof(double));
    SCA_EPH  = (double *) calloc (17280 * 5 * days, sizeof(double));
    ACB_EPH  = (double *) calloc (86400 * 4 * days, sizeof(double));
    SCB_EPH  = (double *) calloc (17280 * 5 * days, sizeof(double));  
  
    info = (InfStruct *) calloc ( NDATA, sizeof(InfStruct));

    sat1a = (DATAL1C *) calloc ( NDATA, sizeof(DATAL1C));
    sat2b = (DATAL1C *) calloc ( NDATA, sizeof(DATAL1C));
    sat12 = (DATAL1C *) calloc ( NDATA, sizeof(DATAL1C));

    fkn = (double *) calloc ( NDATA, sizeof(double));
    fint = (double *) calloc ( NDATA, sizeof(double));


    coefbl  = (double *) calloc (solvegrv, sizeof(double));
    coefb  = (double *) calloc ((ncut + 1) * (ncut + 1), sizeof(double));
    coefrfl  = (double *) calloc (solvegrv, sizeof(double));
    coefrfh  = (double *) calloc ((ncut + 1) * (ncut + 1), sizeof(double));
    coefcsr  = (double *) calloc (solvegrv, sizeof(double));
    coefjpl  = (double *) calloc (solvegrv, sizeof(double));
    coefgfz  = (double *) calloc (solvegrv, sizeof(double));


/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
// read L1B data

    s0 = time(NULL);  

    for (i = 0; i < NDATA; i ++)
    {
        sat1a[i].t = GPS_S + i * DT;
        sat2b[i].t = GPS_S + i * DT;
        sat12[i].t = GPS_S + i * DT;

        sat1a[i].error = 0;
        sat2b[i].error = 0;
        sat12[i].error = 0;
    }

    GMA[0] = 398600.44150E+09;
    GMA[1] = 6378136.3;
    
    readkbr(f_kbr1b_x, &n_kbr);
    readgnv(f_gnv1c_a, f_gnv1c_b, &n_gnva, &n_gnvb);
    readsca(f_sca1b_a, f_sca1b_b);
    readacc(f_acc1b_a, f_acc1b_b);


/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/

    DIM_OT = (int)(86400/step_ot);
    OT_EPH = (double *) calloc ( ((NOMAX + 1) * (NOMAX + 1) + 1) * DIM_OT, 
            sizeof(double));
    if((fp_otbin=fopen(f_aot,"rb"))==NULL)
    {
        printf("Cannot write oteph.bin!\n");
        exit(0);
    }
    fread (OT_EPH, sizeof(double) * ((NOMAX + 1) * (NOMAX + 1) + 1) * DIM_OT,
            1, fp_otbin);

    getinfo_1day(f_eop);
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
// open gravity field file

    grv_open (f_grv[0], f_grv[1], nmax, mmax, coefbl); 
    zero0zero1 (coefbl, nmax);
    grv_open (f_grv[0], f_grv[1], ncut, mcut, coefb); 
    zero0zero1 (coefb, ncut);
    grv_open (f_ref[0], f_ref[1], nmax, mmax, coefrfl); 
    zero0zero1 (coefrfl, nmax);
    grv_open (f_ref[0], f_ref[1], ncut, mcut, coefrfh); 
    zero0zero1 (coefrfh, ncut);
    grv_open (f_csr[0], f_csr[1], nmax, mmax, coefcsr); 
    zero0zero1 (coefcsr, nmax);
    grv_open (f_jpl[0], f_jpl[1], nmax, mmax, coefjpl); 
    zero0zero1 (coefjpl, nmax);
    grv_open (f_gfz[0], f_gfz[1], nmax, mmax, coefgfz); 
    zero0zero1 (coefgfz, nmax);

//    opengrav (f_grv, coefbl, GMA, nmax, mmax, 0);
//    opengrav (f_grv, coefb, GMA, ncut, mcut, 0);
 //   opengrav (f_ref, coefrfl, GMA, nmax, mmax, 0);
//    opengrav (f_ref, coefrfh, GMA, ncut, mcut, 0);
//    opengrav (f_csr, coefcsr, GMA, nmax, mmax, 0);
//    opengrav (f_jpl, coefjpl, GMA, nmax, mmax, 0);
//    opengrav (f_gfz, coefgfz, GMA, nmax, mmax, 0);
    coefgfz[2] = coefgfz[2] + C20PERM;
//    printf ("PERM = %d\n", PERM);
//    if (PERMT == 0)
//    {
//        coefa[2] = coefa[2] + C20PERM;
//        coefe[2] = coefe[2] + C20PERM;
//    }


/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
// calibrate accelerometer data


    cal_acc_01();
//    cal_acc_00();
    calbiaseph(f_par_a, f_par_b);
//    bias_acc(4);

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/


    for (i = 0; i < NDATA; i ++)
    {    
        sat1a[i].rp[0] = gnva[i].xpos; 
        sat1a[i].rp[1] = gnva[i].ypos; 
        sat1a[i].rp[2] = gnva[i].zpos;
        sat1a[i].rv[0] = gnva[i].xvel; 
        sat1a[i].rv[1] = gnva[i].yvel; 
        sat1a[i].rv[2] = gnva[i].zvel;
        sat2b[i].rp[0] = gnvb[i].xpos; 
        sat2b[i].rp[1] = gnvb[i].ypos; 
        sat2b[i].rp[2] = gnvb[i].zpos;
        sat2b[i].rv[0] = gnvb[i].xvel; 
        sat2b[i].rv[1] = gnvb[i].yvel; 
        sat2b[i].rv[2] = gnvb[i].zvel;

        sat1a[i].pos = modvect(sat1a[i].rp);
        sat1a[i].vel = modvect(sat1a[i].rv);
        sat2b[i].pos = modvect(sat2b[i].rp);
        sat2b[i].vel = modvect(sat2b[i].rv);

        for (n = 0; n < 3; n ++)
        {
            sat12[i].rp[n] = sat2b[i].rp[n] - sat1a[i].rp[n];
            sat12[i].rv[n] = sat2b[i].rv[n] - sat1a[i].rv[n];
        }
        sat12[i].pos = modvect(sat12[i].rp);
        sat12[i].vel = modvect(sat12[i].rv);


        sat12[i].range = sat12[i].pos;
        sat12[i].rate = dotvect(sat12[i].rp, sat12[i].rv) / sat12[i].pos;

        if (kbrr != 0)
        {
            if (kbrx[i].gps_time != 0)
            {
                sat12[i].rate = kbrx[i].rate;
            }
            else
            {
                sat12[i].error = 9;
            }

            if (kbrr == 2)
            {
                kbrvel2pos (i);
            }
            else if (kbrr == 1)
            {
                kbrpos2vel (i);
            }
        }

        brmul(info[i].c_ie, sat1a[i].rp, 3, 3, 1, p1e);  
        brmul(info[i].c_ie, sat2b[i].rp, 3, 3, 1, p2e);

        xyz2llr(p1e, sat1a[i].llr);
        xyz2llr(p2e, sat2b[i].llr);

    }

    s1 = time(NULL);  
    printf("\n%5ld: seconds of reading L1B data\n", s1-s0);
    fflush(stdout);

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/

    
    we[0] = 0; we[1] = 0; we[2] = ANGVEL;

    for (i = 0; i < NDATA; i ++)
    {    

        tjd[0] = info[i].jd0;
        tjd[1] = info[i].tt/86400.0;

        wi[0] = info[i].wi[0];
        wi[1] = info[i].wi[1];
        wi[2] = info[i].wi[2];
    
        crsvect(wi, sat1a[i].rp, sat1a[i].wr);
        crsvect(wi, sat2b[i].rp, sat2b[i].wr);

/////////////////
//vk: 1/2*v^2 - GM/r
        sat1a[i].vk = dotvect (sat1a[i].rv, sat1a[i].rv) / 2.0 
                - GMA[0]/sat1a[i].pos;
        sat2b[i].vk = dotvect (sat2b[i].rv, sat2b[i].rv) / 2.0 
                - GMA[0]/sat2b[i].pos;
        sat12[i].vk = sat2b[i].vk - sat1a[i].vk;

////////////////
//vr:
        crsvect(sat1a[i].rp, sat1a[i].rv, rpv);
        sat1a[i].vr = dotvect(rpv, wi);
        sat1a[i].vr3 = dotvect(rpv, we);

        crsvect(sat2b[i].rp, sat2b[i].rv, rpv);
        sat2b[i].vr = dotvect(rpv, wi);
        sat2b[i].vr3 = dotvect(rpv, we);

        sat12[i].vr = sat2b[i].vr - sat1a[i].vr;
        sat12[i].vr3 = sat2b[i].vr3 - sat1a[i].vr3;

        sat1a[i].vr12 = sat1a[i].vr - sat1a[i].vr3;
        sat2b[i].vr12 = sat2b[i].vr - sat2b[i].vr3;
        sat12[i].vr12 = sat12[i].vr - sat12[i].vr3;

///////////////
//vf:
//        disse_a (sat1a[i].rv, i, &sat1a[i].dvaf, sat1a[i].af);
//        disse_b (sat2b[i].rv, i, &sat2b[i].dvaf, sat2b[i].af);

        accia (tjd, 1, sat1a[i].af);
        accia (tjd, 2, sat2b[i].af);
            
        sat1a[i].dvaf = dotvect(sat1a[i].af, sat1a[i].rv);
        sat2b[i].dvaf = dotvect(sat2b[i].af, sat2b[i].rv);
        sat12[i].dvaf = sat2b[i].dvaf - sat1a[i].dvaf;

        sat1a[i].dvrf =  dotvect(sat1a[i].wr, sat1a[i].af);
        sat2b[i].dvrf =  dotvect(sat2b[i].wr, sat2b[i].af);
        sat12[i].dvrf = sat2b[i].dvrf - sat1a[i].dvrf;

        for (n=0;n<3;n++)
            sat12[i].af[n] = sat2b[i].af[n] - sat1a[i].af[n];


////////////////
//vg:


        accgr (tjd, sat1a[i].rp, sat1a[i].rv, sat1a[i].ag);
        accgr (tjd, sat2b[i].rp, sat2b[i].rv, sat2b[i].ag);


        sat1a[i].dvag = dotvect(sat1a[i].ag, sat1a[i].rv);
        sat2b[i].dvag = dotvect(sat2b[i].ag, sat2b[i].rv);
        sat12[i].dvag = sat2b[i].dvag - sat1a[i].dvag;

        sat1a[i].dvrg =  dotvect(sat1a[i].wr, sat1a[i].ag);
        sat2b[i].dvrg =  dotvect(sat2b[i].wr, sat2b[i].ag);
        sat12[i].dvrg = sat2b[i].dvrg - sat1a[i].dvrg;

        for (n=0;n<3;n++)
            sat12[i].ag[n] = sat2b[i].ag[n] - sat1a[i].ag[n];



///////////////
//vn:
        nbodyv (tjd, sat1a[i].rp, &sat1a[i].vn, &sat1a[i].dvtn, sat1a[i].an);
        nbodyv (tjd, sat2b[i].rp, &sat2b[i].vn, &sat2b[i].dvtn, sat2b[i].an);

        sat12[i].vn = sat2b[i].vn - sat1a[i].vn;
        sat12[i].dvtn = sat2b[i].dvtn - sat1a[i].dvtn;

        sat1a[i].dvan = dotvect(sat1a[i].an, sat1a[i].rv);
        sat2b[i].dvan = dotvect(sat2b[i].an, sat2b[i].rv);
        sat12[i].dvan = sat2b[i].dvan - sat1a[i].dvan;

        sat1a[i].dvrn =  dotvect(sat1a[i].wr, sat1a[i].an);
        sat2b[i].dvrn =  dotvect(sat2b[i].wr, sat2b[i].an);
        sat12[i].dvrn = sat2b[i].dvrn - sat1a[i].dvrn;

        for (n=0;n<3;n++)
            sat12[i].an[n] = sat2b[i].an[n] - sat1a[i].an[n];



//////////////
//vs:
        stidev (i, sat1a[i].llr, &sat1a[i].vs, &sat1a[i].dvts, sat1a[i].as);
        stidev (i, sat2b[i].llr, &sat2b[i].vs, &sat2b[i].dvts, sat2b[i].as);

        sat12[i].vs = sat2b[i].vs - sat1a[i].vs;
        sat12[i].dvts = sat2b[i].dvts - sat1a[i].dvts;

        sat1a[i].dvas = dotvect(sat1a[i].as, sat1a[i].rv);
        sat2b[i].dvas = dotvect(sat2b[i].as, sat2b[i].rv);
        sat12[i].dvas = sat2b[i].dvas - sat1a[i].dvas;

        sat1a[i].dvrs =  dotvect(sat1a[i].wr, sat1a[i].as);
        sat2b[i].dvrs =  dotvect(sat2b[i].wr, sat2b[i].as);
        sat12[i].dvrs = sat2b[i].dvrs - sat1a[i].dvrs;

        for (n=0;n<3;n++)
            sat12[i].as[n] = sat2b[i].as[n] - sat1a[i].as[n];

    }

    s2 = time(NULL);  
    printf("\n%5ld: seconds of processing tides correction\n", s2-s1);
    fflush(stdout);


/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
//*simulation potential difference*/


#pragma omp parallel private(pt1, pt2, COEFO, i, n)
    {
        pt1  = (double *) calloc ( solvegrv, sizeof(double));
        pt2  = (double *) calloc ( solvegrv, sizeof(double));
        COEFO = (double *) calloc ( (NOMAX + 1) * (NOMAX + 1), sizeof(double));

        #pragma omp for
    for (i = 0; i < NDATA; i ++)
    {


//////////////
//vo:

        lgr_order (OT_EPH, DIM_OT, (NOMAX + 1) * (NOMAX + 1) + 1, 
                info[i].tt, COEFO, 4);

        cs2acc (i, sat1a[i].llr, COEFO, GMA[0], GMA[1], NOMAX, 
            &sat1a[i].vo, &sat1a[i].dvto, sat1a[i].ao);
        cs2acc (i, sat2b[i].llr, COEFO, GMA[0], GMA[1], NOMAX, 
            &sat2b[i].vo, &sat2b[i].dvto, sat2b[i].ao);
        sat12[i].vo = sat2b[i].vo - sat1a[i].vo;

        sat1a[i].dvro =  dotvect(sat1a[i].wr, sat1a[i].ao);
        sat2b[i].dvro =  dotvect(sat2b[i].wr, sat2b[i].ao);
        sat12[i].dvro = sat2b[i].dvro - sat1a[i].dvro;

        for (n=0;n<3;n++)
            sat12[i].ao[n] = sat2b[i].ao[n] - sat1a[i].ao[n]; 

        sat1a[i].dvao = dotvect(sat1a[i].ao, sat1a[i].rv);
        sat2b[i].dvao = dotvect(sat2b[i].ao, sat2b[i].rv);
        sat12[i].dvao = sat2b[i].dvao - sat1a[i].dvao;



/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
//reference l1c and los for RL05
//
        cs2acc (i, sat1a[i].llr, coefrfl, GMA[0], GMA[1], nmax, 
            &sat1a[i].vrfl, &sat1a[i].dvtb, sat1a[i].arfl);
        cs2acc (i, sat2b[i].llr, coefrfl, GMA[0], GMA[1], nmax, 
            &sat2b[i].vrfl, &sat2b[i].dvtb, sat2b[i].arfl);
        sat12[i].vrfl = sat2b[i].vrfl - sat1a[i].vrfl; 
        for (n=0;n<3;n++)
            sat12[i].arfl[n] = sat2b[i].arfl[n] - sat1a[i].arfl[n]; 
        sat12[i].grfl = dotvect (sat12[i].arfl, sat12[i].rp) / sat12[i].pos;


        cs2acc (i, sat1a[i].llr, coefrfh, GMA[0], GMA[1], ncut,
            &sat1a[i].vrfh, &sat1a[i].dvtb, sat1a[i].arfh);
        cs2acc (i, sat2b[i].llr, coefrfh, GMA[0], GMA[1], ncut,
            &sat2b[i].vrfh, &sat2b[i].dvtb, sat2b[i].arfh);
        sat12[i].vrfh = sat2b[i].vrfh - sat1a[i].vrfh;  
        for (n=0;n<3;n++)
            sat12[i].arfh[n] = sat2b[i].arfh[n] - sat1a[i].arfh[n];   
        sat12[i].grfh = dotvect (sat12[i].arfh, sat12[i].rp) / sat12[i].pos;



/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/


        cs2acc (i, sat1a[i].llr, coefcsr, GMA[0], GMA[1], nmax, 
            &sat1a[i].vcsr, &sat1a[i].dvtb, sat1a[i].acsr);
        cs2acc (i, sat2b[i].llr, coefcsr, GMA[0], GMA[1], nmax, 
            &sat2b[i].vcsr, &sat2b[i].dvtb, sat2b[i].acsr);
        sat12[i].vcsr = sat2b[i].vcsr - sat1a[i].vcsr;  
        for (n=0;n<3;n++)
            sat12[i].acsr[n] = sat2b[i].acsr[n] - sat1a[i].acsr[n];  
        sat12[i].gcsr = dotvect (sat12[i].acsr, sat12[i].rp) / sat12[i].pos;

        cs2acc (i, sat1a[i].llr, coefgfz, GMA[0], GMA[1], nmax, 
            &sat1a[i].vgfz, &sat1a[i].dvtb, sat1a[i].agfz);
        cs2acc (i, sat2b[i].llr, coefgfz, GMA[0], GMA[1], nmax, 
            &sat2b[i].vgfz, &sat2b[i].dvtb, sat2b[i].agfz);
        sat12[i].vgfz = sat2b[i].vgfz - sat1a[i].vgfz;  
        for (n=0;n<3;n++)
            sat12[i].agfz[n] = sat2b[i].agfz[n] - sat1a[i].agfz[n];  
        sat12[i].ggfz = dotvect (sat12[i].agfz, sat12[i].rp) / sat12[i].pos;

        cs2acc (i, sat1a[i].llr, coefjpl, GMA[0], GMA[1], nmax, 
            &sat1a[i].vjpl, &sat1a[i].dvtb, sat1a[i].ajpl);
        cs2acc (i, sat2b[i].llr, coefjpl, GMA[0], GMA[1], nmax, 
            &sat2b[i].vjpl, &sat2b[i].dvtb, sat2b[i].ajpl);
        sat12[i].vjpl = sat2b[i].vjpl - sat1a[i].vjpl;  
        for (n=0;n<3;n++)
            sat12[i].ajpl[n] = sat2b[i].ajpl[n] - sat1a[i].ajpl[n];  
        sat12[i].gjpl = dotvect (sat12[i].ajpl, sat12[i].rp) / sat12[i].pos;


/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
//background model of orbit 
//

        cs2acc (i, sat1a[i].llr, coefbl, GMA[0], GMA[1], nmax,
            &sat1a[i].vbl, &sat1a[i].dvtb, sat1a[i].abl);
        cs2acc (i, sat2b[i].llr, coefbl, GMA[0], GMA[1], nmax,
            &sat2b[i].vbl, &sat2b[i].dvtb, sat2b[i].abl);
        sat12[i].vbl = sat2b[i].vbl - sat1a[i].vbl;   
        for (n=0;n<3;n++)
            sat12[i].abl[n] = sat2b[i].abl[n] - sat1a[i].abl[n]; 
        sat12[i].gbl = dotvect (sat12[i].abl, sat12[i].rp) / sat12[i].pos;


        cs2acc (i, sat1a[i].llr, coefb, GMA[0], GMA[1], ncut, 
            &sat1a[i].vb, &sat1a[i].dvtb, sat1a[i].ab);
        cs2acc (i, sat2b[i].llr, coefb, GMA[0], GMA[1], ncut, 
            &sat2b[i].vb, &sat2b[i].dvtb, sat2b[i].ab);

        sat12[i].dvtb = sat2b[i].dvtb - sat1a[i].dvtb; 

        sat1a[i].dvrb =  dotvect(sat1a[i].wr, sat1a[i].ab);
        sat2b[i].dvrb =  dotvect(sat2b[i].wr, sat2b[i].ab);
        sat12[i].dvrb = sat2b[i].dvrb - sat1a[i].dvrb; 

        for (n=0;n<3;n++)
            sat12[i].ab[n] = sat2b[i].ab[n] - sat1a[i].ab[n];   

        sat12[i].vb = sat2b[i].vb - sat1a[i].vb; 

        sat1a[i].dvab = dotvect(sat1a[i].ab, sat1a[i].rv);
        sat2b[i].dvab = dotvect(sat2b[i].ab, sat2b[i].rv);
        sat12[i].dvab = sat2b[i].dvab - sat1a[i].dvab;

        
        for (n=0;n<3;n++)
        {
            sat12[i].ra[n] = - GMA[0]/pow(sat2b[i].pos,3) * sat2b[i].rp[n] 
                             + GMA[0]/pow(sat1a[i].pos,3) * sat1a[i].rp[n]
                           + sat12[i].ab[n] 
                           + sat12[i].af[n]
                           + sat12[i].an[n] 
                           + sat12[i].as[n]
                           + sat12[i].ao[n]
                           + sat12[i].ag[n];
        }
        sat12[i].dvcorr = 
                    - sat12[i].dvaf + sat12[i].dvrf
                    - sat12[i].dvan + sat12[i].dvrn 
                    - sat12[i].dvas + sat12[i].dvrs 
                    - sat12[i].dvao + sat12[i].dvro 
                    - sat12[i].dvag + sat12[i].dvrg;


        sat12[i].gb = dotvect (sat12[i].ab, sat12[i].rp) / sat12[i].pos;
        sat12[i].gc = dotvect (sat12[i].ra, sat12[i].rp) / sat12[i].pos;

        sat12[i].gcorr = sat12[i].gc - sat12[i].gb;


        sat12[i].accl = (sat12[i].vel * sat12[i].vel 
                        - sat12[i].rate * sat12[i].rate 
                        + dotvect(sat12[i].rp, sat12[i].ra)) / sat12[i].range;

        if (kbra != 0)
        {
            if (kbrx[i].gps_time != 0)
                sat12[i].accl = kbrx[i].accl;
//            if (kbrx[i].gps_time != 0 && i != 0 && i != NDATA - 1)
//                sat12[i].accl = (kbrx[i+1].rate - kbrx[i-1].rate)/2.0/DT;
            else
                sat12[i].error = 9;
        }


        sat12[i].glos = sat12[i].accl 
                        + (sat12[i].rate * sat12[i].rate 
                        - sat12[i].vel * sat12[i].vel) / sat12[i].pos
                        - sat12[i].gcorr -  sat12[i].gb;
        
    }
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/


        free(pt1);
        free(pt2);
        free(COEFO);

    }

    s3 = time(NULL);  
    printf("\n%5ld: seconds of processing gravity and partial\n", s3-s2);
    fflush(stdout);


///////////////////
//vint: all the potential needed to be intergrated
    for (i = 0; i < NDATA; i ++)
        fkn[i] = sat12[i].dvcorr;
    intydt(NDATA, DT, fkn, fint);
    for (i = 0; i < NDATA; i ++)
        sat12[i].vcorr = fint[i];



    s4 = time(NULL);  
    printf("\n%5ld: seconds of finshing\n", s4-s3);
    fflush(stdout);




    for (i = 0; i < NDATA; i ++)
    {
        sat12[i].vc = sat12[i].vk  - sat12[i].vr + sat12[i].vcorr;

        sat12[i].vl1c = sat12[i].vc - sat12[i].vb;

        sat12[i].vref = 0;
        sat12[i].gref = 0;

//        sat12[i].vref = (sat12[i].vcsr + sat12[i].vgfz + sat12[i].vjpl 
//                      - 3 * sat12[i].vbl) / 3.0;

//        sat12[i].gref = (sat12[i].gcsr + sat12[i].ggfz + sat12[i].gjpl 
//                      - 3 * sat12[i].gbl) / 3.0;

    }



    fflush(stdout);

    
    endcut = endcut / DT;

    for (i = 0; i < NDATA; i ++)
//        fkn[i] = sat12[i].vl1c - sat12[i].vref;
        fkn[i] = sat12[i].vl1c;
    daydft (fkn, fint, NDATA,  DT, 5400 / fltpar);
    for (i = endcut; i < NDATA - endcut; i ++)
//        sat12[i].vl1cflt = fint[i] + sat12[i].vref;
        sat12[i].vl1cflt = fint[i];

    for (i = 0; i < NDATA; i ++)
//        fkn[i] = sat12[i].glos - sat12[i].gref;
        fkn[i] = sat12[i].glos;
    daydft (fkn, fint, NDATA,  DT, 5400 / fltpar);
    for (i = endcut; i < NDATA - endcut; i ++)
//        sat12[i].glosflt = fint[i] + sat12[i].gref;
        sat12[i].glosflt = fint[i];



// HARD CODE for fitting!!!!
// HARD CODE for fitting!!!!
//    fitres (order_poly, order_cpr, overlapcyc);

//    if (order_poly != 0 || order_cpr !=0)
    fitl1c_pw (order_poly, order_cpr);
    fitlos_pw (order_poly, order_cpr);

    fitl1c_mv (nply, ncpr, nplycpr, ndt);
    fitlos_mv (nply, ncpr, nplycpr, ndt);

//    printf ("!!!!!!!!!\n");

    for (i = 0; i < NDATA; i ++)
    {

//        if (sat12[i].error != 0) continue;
// OUTPUT SIM
// OUTPUT SIM
        sat12[i].vmvfitstd = fabs(sat12[i].vl1c - sat12[i].vmvfit - sat12[i].vref);
        sat12[i].vpwfitstd = fabs(sat12[i].vl1c - sat12[i].vpwfit - sat12[i].vref);
        sat12[i].gmvfitstd = fabs(sat12[i].glos - sat12[i].gmvfit - sat12[i].gref);
        sat12[i].gpwfitstd = fabs(sat12[i].glos - sat12[i].gpwfit - sat12[i].gref);
    

        fprintf (fpout_sim, "%10d ", gnva[i].gps_time);                     //1: gps_time
        fprintf (fpout_sim, "%25.15f ", sat1a[i].llr[0]);                   //2
        fprintf (fpout_sim, "%25.15f ", sat1a[i].llr[1]);                   //3
        fprintf (fpout_sim, "%25.12f ", sat1a[i].llr[2]);                   //4
        fprintf (fpout_sim, "%25.15f ", sat2b[i].llr[0]);                   //5
        fprintf (fpout_sim, "%25.15f ", sat2b[i].llr[1]);                   //6
        fprintf (fpout_sim, "%25.12f ", sat2b[i].llr[2]);                   //7
        fprintf (fpout_sim, "%20.12e ", sat12[i].vl1c);                     //8
        fprintf (fpout_sim, "%20.12e ", sat12[i].glos);                     //9
        fprintf (fpout_sim, "%20.12e ", sat12[i].vl1c + sat12[i].vb - sat12[i].vrfh - sat12[i].vmvfit);   //10
        fprintf (fpout_sim, "%20.12e ", sat12[i].vl1c + sat12[i].vb - sat12[i].vrfh - sat12[i].vpwfit);   //11
        fprintf (fpout_sim, "%20.12e ", sat12[i].vl1c + sat12[i].vb - sat12[i].vrfh );   //11
//        fprintf (fpout_sim, "%20.12e ", sat12[i].vl1cflt + sat12[i].vb - sat12[i].vrfh);                  //12
        fprintf (fpout_sim, "%20.12e ", sat12[i].glos + sat12[i].gb - sat12[i].grfh - sat12[i].gmvfit);   //13
        fprintf (fpout_sim, "%20.12e ", sat12[i].glos + sat12[i].gb - sat12[i].grfh - sat12[i].gpwfit);   //14
        fprintf (fpout_sim, "%20.12e ", sat12[i].glos + sat12[i].gb - sat12[i].grfh );   //14
//        fprintf (fpout_sim, "%20.12e ", sat12[i].glosflt + sat12[i].gb - sat12[i].grfh);                  //15
        fprintf (fpout_sim, "%20.12e ", sat12[i].vcsr - sat12[i].vrfl);      //16
        fprintf (fpout_sim, "%20.12e ", sat12[i].vgfz - sat12[i].vrfl);      //17
        fprintf (fpout_sim, "%20.12e ", sat12[i].vjpl - sat12[i].vrfl);      //18
        fprintf (fpout_sim, "%20.12e ", sat12[i].gcsr - sat12[i].grfl);      //19
        fprintf (fpout_sim, "%20.12e ", sat12[i].ggfz - sat12[i].grfl);      //20
        fprintf (fpout_sim, "%20.12e ", sat12[i].gjpl - sat12[i].grfl);      //21

        fprintf (fpout_sim, "%d\t", sat12[i].error);                        //22
        fprintf (fpout_sim, "%20.12f ", sat12[i].vmvfitres);                //23
        fprintf (fpout_sim, "%20.12f ", sat12[i].vpwfitres);                //24
        fprintf (fpout_sim, "%20.12f ", sat12[i].gmvfitres);                //25
        fprintf (fpout_sim, "%20.12f ", sat12[i].gpwfitres);                //26
        fprintf (fpout_sim, "%20.12f ", sat12[i].vmvfitstd);                //27
        fprintf (fpout_sim, "%20.12f ", sat12[i].vpwfitstd);                //28
        fprintf (fpout_sim, "%20.12f ", sat12[i].gmvfitstd);                //29
        fprintf (fpout_sim, "%20.12f ", sat12[i].gpwfitstd);                //30

        fprintf (fpout_sim, "%20.12f ", sat12[i].Tp2);                      //31
        fprintf (fpout_sim, "%20.12f ", sat12[i].Tp1);                      //32
        fprintf (fpout_sim, "%20.12f ", gnva[i].Tp);                        //33
        fprintf (fpout_sim, "%20.12f ", gnvb[i].Tp);                        //34
        fprintf (fpout_sim, "%20.12f ", sat12[i].range);                      //35
        fprintf (fpout_sim, "%20.12f ", sat12[i].rate);                      //36
        fprintf (fpout_sim, "%20.12f ", sat12[i].accl);                      //37

        fprintf (fpout_sim, "\n");
    }



    s5 = time(NULL);  
    printf("\n%5ld: seconds of process L1C data\n", s5-s4);
    fflush(stdout);


    fflush(fpout_sim);
//        exit(1);




    free (coefb);
    free (coefbl);
    free (coefrfl);
    free (coefrfh);
    free (coefjpl);
    free (coefcsr);
    free (coefgfz);
    free (kbrx);
    free (gnva);
    free (gnvb);
    free (ACA_EPH);
    free (SCA_EPH);
    free (ACB_EPH);
    free (SCB_EPH);
    free (AOD_EPH);
    free (info);
    free (sat1a);
    free (sat2b);
    free (sat12);
    free (fkn);
    free (fint);
    fclose(fpout_sim);
   
    ephem_close();  /* remove this line for use with solsys version 2 */

    printf ("\nNormal end of ECHO!\n\npress any key to finish...\n");
    
    return 0;

}
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
 


