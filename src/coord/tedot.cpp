/*
------------------------------------------------------------------------

    Purpose: Compute the dot of transformation matrix from ITRF to ICRF
        a supplement code to ter2cel from novas lib
    Notes: method can only equal to 1 (equinox-based)
    Programmer: Kun Shang @ 3.21.2014    
    Functions:

        short int ter2cel_dot (double jd_ut_high, double jd_ut_low, 
                    double delta_t, short int method, 
                    short int accuracy, short int option,
                    double x, double y, double *vect,
                    double *vecc);

        short int sidereal_time_dot (double jd_high, double jd_low,
                    double delta_t,short int gst_type,
                    short int method, short int accuracy,
                    double *gst, double *gst_dot);

        void spin_dot (double angle, double angle_dot, double *pos1,
                    double *pos2);

------------------------------------------------------------------------
*/

#ifndef _COORD_H_

#include "coord.h"
    #define _COORD_H_
#endif


/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/

/********ter2cel_dot */

short int ter2cel_dot (double jd_ut_high, double jd_ut_low, 
                    double delta_t, short int method, 
                    short int accuracy, short int option,
                    double x, double y, double *vect,
                    double *vecc)
/*
------------------------------------------------------------------------
    Purpose: Compute the dot of ter2cel, 
        a supplement code to ter2cel from novas lib
    Notes: method can only equal to 1 (equinox-based)
    Programmer: Kun Shang @ 3.21.2014    
------------------------------------------------------------------------
*/

/*
------------------------------------------------------------------------


------------------------------------------------------------------------
*/
{
   short int error = 0;
   short int rs, j;

   double jd_ut1, jd_tt, dummy, secdiff, jd_tdb, gast, gast_dot, r_cio, theta,
   v1[3], v2[3], v3[3], v4[3], xx[3], yy[3], zz[3];

/*
   Invalid value of 'accuracy'.
*/

   if ((accuracy < 0) || (accuracy > 1))
      return (error = 1);

/*
   Compute the TT Julian date corresponding to the input UT1 Julian
   date.
*/

   jd_ut1 = jd_ut_high + jd_ut_low;
   jd_tt = jd_ut1 + (delta_t / 86400.0);

/*
   Compute the TDB Julian date corresponding to the input UT1 Julian
   date.
*/

   jd_tdb = jd_tt;
   tdb2tt (jd_tdb, &dummy,&secdiff);
   jd_tdb = jd_tt + secdiff / 86400.0;

   switch (method)
   {
      case (0):

        return (error +=10);

/*
   'CIO-TIO-THETA' method.

   See second reference, eq. (3) and (4).

   Apply polar motion, transforming the vector to the terrestrial
   intermediate system.
*/

         if ((x == 0.0) && (y == 0.0))
         {
            v1[0] = vect[0];
            v1[1] = vect[1];
            v1[2] = vect[2];
         }
          else
            wobble (jd_tdb,0, x,y,vect, v1);

/*
   Obtain the basis vectors, in the GCRS, of the celestial intermediate
   system.
*/

      if ((error = cio_location (jd_tdb,accuracy, &r_cio,&rs)) != 0)
         return (error += 10);

      if ((error = cio_basis (jd_tdb,r_cio,rs,accuracy, xx,yy,zz)) != 0)
         return (error += 20);

/*
   Compute and apply the Earth rotation angle, 'theta', transforming the
   vector to the celestial intermediate system.
*/

         theta = era (jd_ut_high,jd_ut_low);
         spin (-theta,v1, v2);

/*
   Transform the vector from the celestial intermediate system to the
   GCRS.
*/

         vecc[0] = xx[0] * v2[0] + yy[0] * v2[1] + zz[0] * v2[2];
         vecc[1] = xx[1] * v2[0] + yy[1] * v2[1] + zz[1] * v2[2];
         vecc[2] = xx[2] * v2[0] + yy[2] * v2[1] + zz[2] * v2[2];
         break;

      case (1):

/*
   Equinox mode.

   Apply polar motion.
*/

         if ((x == 0.0) && (y == 0.0))
         {
            for (j = 0; j < 3; j++)
            {
               v1[j] = vect[j];
            }
         }
          else
            wobble (jd_tdb,0,x,y,vect, v1);

/*
   Apply Earth rotation.
*/

         sidereal_time_dot (jd_ut_high,jd_ut_low,delta_t,1,1,
                            accuracy, &gast, &gast_dot);
         spin_dot (-gast * 15.0, -gast_dot, v1, v2);

/*
   'option' = 1 skips remaining transformations.
*/

         if (option == 1)
         {
            vecc[0] = v2[0];
            vecc[1] = v2[1];
            vecc[2] = v2[2];
         }
          else
         {

/*
   Apply precession, nutation, and frame tie.
*/

            nutation (jd_tdb,-1,accuracy,v2, v3);
            precession (jd_tdb,v3,T0, v4);
            frame_tie (v4,-1, vecc);
         }
         break;

/*
   Invalid value of 'method'.
*/

      default:
         error = 2;
         break;
   }

   return (error);
}







/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/

/********sidereal_time_dot */

short int sidereal_time_dot (double jd_high, double jd_low,
                         double delta_t,short int gst_type,
                         short int method, short int accuracy,
                         double *gst, double *gst_dot)
/*
------------------------------------------------------------------------
    Purpose: Compute the dot of Greenwich apparent sidereal time, 
        a supplement code to sidereal_time from novas lib
    Notes: method can only equal to 1 (equinox-based)
    Programmer: Kun Shang @ 3.21.2014    
------------------------------------------------------------------------
*/

/*
------------------------------------------------------------------------

   PURPOSE:
      Computes the Greenwich apparent sidereal time, at Julian date
      'jd_high' + 'jd_low'.

   REFERENCES:
      Kaplan, G. (2005), US Naval Observatory Circular 179.

   INPUT
   ARGUMENTS:
      jd_high (double)
         High-order part of UT1 Julian date.
      jd_low (double)
         Low-order part of UT1 Julian date.
      delta_t (double)
         Difference TT-UT1 at 'jd_high'+'jd_low', in seconds
         of time.
      gst_type (short int)
         = 0 ... compute Greenwich mean sidereal time
         = 1 ... compute Greenwich apparent sidereal time
      method (short int)
         Selection for method
            = 0 ... CIO-based method
            = 1 ... equinox-based method
      accuracy (short int)
         Selection for accuracy
            = 0 ... full accuracy
            = 1 ... reduced accuracy

   OUTPUT
   ARGUMENTS:
      *gst (double)
         Greenwich apparent sidereal time, in hours.

   RETURNED
   VALUE:
      (short int)
         = 0         ... everything OK
         = 1         ... invalid value of 'accuracy'
         = 2         ... invalid value of 'method'
         > 10, < 30  ... 10 + error from function 'cio_rai'

   GLOBALS
   USED:
      T0, RAD2DEG

   FUNCTIONS
   CALLED:
      era                novas.c
      tdb2tt             novas.c
      e_tilt             novas.c
      cio_location       novas.c
      cio_basis          novas.c
      nutation           novas.c
      precession         novas.c
      frame_tie          novas.c
      fabs               math.h
      fmod               math.h
      atan2              math.h

   VER./DATE/
   PROGRAMMER:
      V1.0/06-92/TKB (USNO/NRL Optical Interfer.) Translate Fortran.
      V1.1/08-93/WTH (USNO/AA) Update to C programing standards.
      V1.2/03-98/JAB (USNO/AA) Expand documentation.
      V1.3/08-98/JAB (USNO/AA) Match flow of the Fortran counterpart.
      V2.0/09-03/JAB (USNO/AA) Incorporate the 2nd-reference changes.
      V2.1/08-04/JAB (USNO/AA) Incorporate the 1st-reference changes.
      V2.2/12-05/WKP (USNO/AA) Updated error handling.
      V2.3/01-06/WKP (USNO/AA) Changed 'mode' to 'method' and 'accuracy'.
      V2.4/04-06/JAB (USNO/AA) Use precession-in-RA terms in mean
                               sidereal time from third reference.
      V2.5/07-06/JAB (USNO/AA) Implement 'cio_location' construct.
      V2.6/06-08/WKP (USNO/AA) Changed value of direction argument in
                               call to 'nutation' from 1 to -1 for
                               consistency.

   NOTES:
      1. The Julian date may be split at any point, but for highest
      precision, set 'jd_high' to be the integral part of the Julian
      date, and set 'jd_low' to be the fractional part.
      2. This function is the C version of NOVAS Fortran routine
      'sidtim'.

------------------------------------------------------------------------
*/
{
   short int error = 0;
   short int ref_sys;

   static double ee;
   static double jd_last = -99.0;
   double unitx[3] = {1.0, 0.0, 0.0};
   double jd_ut, jd_tt, jd_tdb, tt_temp, t, theta, a, b, c, d,
      ra_cio, x[3], y[3], z[3], w1[3], w2[3], eq[3], ha_eq, st, st_dot,
      secdiff, eqeq;

/*
   Invalid value of 'accuracy'.
*/

   if ((accuracy < 0) || (accuracy > 1))
      return (error = 1);

/*
   Time argument for precession and nutation components of sidereal
   time is TDB.  First approximation is TDB = TT, then refine.
*/

   jd_ut = jd_high + jd_low;
   jd_tt = jd_ut + (delta_t / 86400.0);
   jd_tdb = jd_tt;
   tdb2tt (jd_tdb, &tt_temp,&secdiff);
   jd_tdb = jd_tt + (secdiff / 86400.0);

   t = (jd_tdb - T0) / 36525.0;

/*
   Compute the Earth Rotation Angle.  Time argument is UT1.
*/

   theta = era (jd_high, jd_low);

/*
   Compute the equation of the equinoxes if needed, depending upon the
   input values of 'gst_type' and 'method'.  If not needed, set to zero.
*/

   if (((gst_type == 0) && (method == 0)) ||       /* GMST; CIO-TIO */
       ((gst_type == 1) && (method == 1)))         /* GAST; equinox */
   {
      if (fabs (jd_tdb - jd_last) > 1.0e-8)
      {
         e_tilt (jd_tdb,accuracy, &a,&b,&ee,&c,&d);
         jd_last = jd_tdb;
      }
      eqeq = ee * 15.0;
   }
    else
   {
      eqeq = 0.0;
   }

/*
   Compute Greenwich sidereal time depending upon input values of
   'method' and 'gst_type'.
*/

   switch (method)
   {
      case (0):

        return (error +=10);

/*
   Use 'CIO-TIO-theta' method.  See Circular 179, Section 6.5.4.
*/

/*
   Obtain the basis vectors, in the GCRS, of the celestial intermediate
   system.
*/

         if ((error = cio_location (jd_tdb,accuracy, &ra_cio,
            &ref_sys)) != 0)
         {
            *gst = 99.0;
            return (error += 10);
         }

         cio_basis (jd_tdb,ra_cio,ref_sys,accuracy, x,y,z);

/*
   Compute the direction of the true equinox in the GCRS.
*/

         nutation (jd_tdb,-1,accuracy,unitx, w1);
         precession (jd_tdb,w1,T0, w2);
         frame_tie (w2,-1, eq);

/*
   Compute the hour angle of the equinox wrt the TIO meridian
   (near Greenwich, but passes through the CIP and TIO).
*/

         ha_eq = theta - atan2 ((eq[0] * y[0] + eq[1] * y[1] +
            eq[2] * y[2]), (eq[0] * x[0] + eq[1] * x[1] +
            eq[2] * x[2])) * RAD2DEG;

/*
   For mean sidereal time, subtract the equation of the equinoxes.
*/

         ha_eq -= (eqeq / 240.0);

         ha_eq = fmod (ha_eq, 360.0) / 15.0;
         if (ha_eq < 0.0)
            ha_eq += 24.0;
         *gst = ha_eq;
         break;

      case (1):

/*
   Use equinox method.  See Circular 179, Section 2.6.2.
*/

/*
   Precession-in-RA terms in mean sidereal time taken from third
   reference, eq. (42), with coefficients in arcseconds.
*/

         st = eqeq + 0.014506 +
               (((( -    0.0000000368   * t
                    -    0.000029956  ) * t
                    -    0.00000044   ) * t
                    +    1.3915817    ) * t
                    + 4612.156534     ) * t;

/*
   Form the Greenwich sidereal time.
*/


         *gst = fmod ((st / 3600.0 + theta), 360.0) / 15.0;


//         st_dot = (((( - 5 * 0.0000000368 * t - 4 * 0.000029956 ) * t 
//         - 3 * 0.00000044 ) * t + 2 * 1.3915817 ) * t + 4612.156534 );

         st_dot = 
               (((( - 5 * 0.0000000368   * t
                    - 4 * 0.000029956  ) * t
                    - 4 * 0.00000044   ) * t
                    + 2 * 1.3915817    ) * t
                    +  4612.156534     );

         
         *gst_dot = 1.00273781191135448 * TWOPI / 86400.0 
            + st_dot / 206264.806 / 36525.0 / 86400.0;


         if (*gst < 0.0)
            *gst += 24.0;
         break;

/*
   Invalid value of 'method'.
*/

      default:
         *gst = 99.0;
         error = 2;
         break;
   }

   return (error);
}





/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/

/********spin_dot */

void spin_dot (double angle, double angle_dot, double *pos1,
           double *pos2)
/*
------------------------------------------------------------------------
    Purpose: Compute the dot of a transformation matrix along z-axis, 
        a supplement code to spin from novas lib
    Programmer: Kun Shang @ 3.21.2014    
------------------------------------------------------------------------
*/

/*
------------------------------------------------------------------------

   PURPOSE:
      This function transforms a vector from one coordinate system
      to another with same origin and axes rotated about the z-axis.


------------------------------------------------------------------------
*/
{
   static double ang_last = -999.0;
   static double xx, yx, zx, xy, yy, zy, xz, yz, zz;
   double angr, cosang, sinang;

   if (fabs (angle - ang_last) >= 1.0e-12)
   {
      angr = angle * DEG2RAD;
      cosang = cos (angr);
      sinang = sin (angr);

/*
   Rotation matrix follows.
*/

//      xx =  cosang;
//      yx =  sinang;
//      zx =  0.0;
//      xy =  -sinang;
//      yy =  cosang;
//      zy =  0.0;
//      xz =  0.0;
//      yz =  0.0;
//      zz =  1.0;

      xx =  -sinang;
      yx =  cosang;
      zx =  0.0;
      xy =  -cosang;
      yy =  -sinang;
      zy =  0.0;
      xz =  0.0;
      yz =  0.0;
      zz =  0.0;

      ang_last = angle;
   }

/*
   Perform rotation.
*/

   pos2[0] = (xx * pos1[0] + yx * pos1[1] + zx * pos1[2]) * angle_dot;
   pos2[1] = (xy * pos1[0] + yy * pos1[1] + zy * pos1[2]) * angle_dot;
   pos2[2] = (xz * pos1[0] + yz * pos1[1] + zz * pos1[2]) * angle_dot;

   return;
}


