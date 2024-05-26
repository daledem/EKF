//$Header$
//
// EKF_GEOS3
//
// Author: David Ledesma
// Created: 2024/05/23
//
//------------------------------------------------------------------------------
#include "../include/anglesdr.h"

extern Matrix eopdata;

//---------------------------------
// public methods
//---------------------------------

//------------------------------------------------------------------------------
// void anglesdr(Matrix& r2,Matrix& v2,double az1,double az2,double az3,
//                 double el1,double el2,double el3,double Mjd1,double Mjd2,
//                 double Mjd3,Matrix rsite1,Matrix rsite2,Matrix rsite3)
//------------------------------------------------------------------------------
/**
 *   This function solves the problem of orbit determination using three
 *       optical sightings.
 *
 * @param[out] <r2> ijk position vector at t2   [m]
 * @param[out] <v2> ijk velocity vector at t2   [m/s]
 * @param <az1> azimuth at t1   [rad]
 * @param <az2> azimuth at t2   [rad]
 * @param <az3> azimuth at t3   [rad]
 * @param <el1> elevation at t1   [rad]
 * @param <el2> elevation at t2   [rad]
 * @param <el3> elevation at t3   [rad]
 * @param <Mjd1> Modified julian date of t1
 * @param <Mjd2> Modified julian date of t2
 * @param <Mjd3> Modified julian date of t3
 * @param <rsite1> ijk site1 position vector   [m]
 * @param <rsite2> ijk site2 position vector   [m]
 * @param <rsite3> ijk site3 position vector   [m]
 *
 */
//------------------------------------------------------------------------------
void anglesdr(Matrix &r2, Matrix &v2, double az1, double az2, double az3, double el1, double el2, double el3, double Mjd1, double Mjd2, double Mjd3, Matrix rsite1, Matrix rsite2, Matrix rsite3) {
    double lon1, lat1, h1,lon2, lat2, h2,lon3, lat3, h3,Mjd_UTC,x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,
            TAI_UTC,UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC,Mjd_TT,Mjd_UT1,tol,pctchg,magr1,magr2, a,f1,magr1in,magr2in,
            t1,t3,magr1old,magr2old,magrsite1,magrsite2,magrsite3,cc1,cc2,ktr,f2,q1,deltae32,f,g,magr1o,deltar1,f1delr1,
            f2delr1,q2,pf1pr1,pf2pr1,magr2o,deltar2,f1delr2,f2delr2,q3,pf1pr2,pf2pr2,delta,delta1,delta2;
    std::string error;
    char direct;

    magr1in = 1.1*Const::R_Earth;
    magr2in = 1.11*Const::R_Earth;
    direct  = 'y';

    tol    = 1e-8*Const::R_Earth;
    pctchg = 0.005;

    t1 = (Mjd1 - Mjd2)*86400.0;
    t3 = (Mjd3 - Mjd2)*86400.0;

    double vLos1[] = {cos(el1)*sin(az1), cos(el1)*cos(az1), sin(el1)};
    Matrix los1(3,1,vLos1,3);

    double vLos2[] = {cos(el2)*sin(az2), cos(el2)*cos(az2), sin(el2)};
    Matrix los2(3,1,vLos2,3);

    double vLos3[] = {cos(el3)*sin(az3), cos(el3)*cos(az3), sin(el3)};
    Matrix los3(3,1,vLos3,3);


    Geodetic(lon1, lat1, h1,rsite1.trans());
    Geodetic(lon2, lat2, h2,rsite2.trans());
    Geodetic(lon3, lat3, h3,rsite3.trans());

    Matrix M1 = LTC(lon1, lat1);
    Matrix M2 = LTC(lon2, lat2);
    Matrix M3 = LTC(lon3, lat3);

    // body-fixed system
    los1 = M1.trans()*los1;
    los2 = M1.trans()*los2;
    los3 = M1.trans()*los3;

    // mean of date system (J2000)
    Mjd_UTC = Mjd1;
    IERS(x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC,eopdata,Mjd_UTC,'l');
    timediff(UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC,UT1_UTC,TAI_UTC);
    Mjd_TT = Mjd_UTC + TT_UTC/86400;
    Mjd_UT1 = Mjd_TT + (UT1_UTC-TT_UTC)/86400;

    Matrix P = PrecMatrix(Const::MJD_J2000,Mjd_TT);
    Matrix N = NutMatrix(Mjd_TT);
    Matrix T = N * P;
    Matrix E = PoleMatrix(x_pole,y_pole) * GHAMatrix(Mjd_UT1) * T;

    los1 = E.trans()*los1;
    rsite1 = E.trans()*rsite1;

    Mjd_UTC = Mjd2;
    IERS(x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC,eopdata,Mjd_UTC,'l');
    timediff(UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC,UT1_UTC,TAI_UTC);
    Mjd_TT = Mjd_UTC + TT_UTC/86400;
    Mjd_UT1 = Mjd_TT + (UT1_UTC-TT_UTC)/86400;

    P = PrecMatrix(Const::MJD_J2000,Mjd_TT);
    N = NutMatrix(Mjd_TT);
    T = N * P;
    E = PoleMatrix(x_pole,y_pole) * GHAMatrix(Mjd_UT1) * T;

    los2 = E.trans()*los2;
    rsite2 = E.trans()*rsite2;

    Mjd_UTC = Mjd3;
    IERS(x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC,eopdata,Mjd_UTC,'l');
    timediff(UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC,UT1_UTC,TAI_UTC);
    Mjd_TT = Mjd_UTC + TT_UTC/86400;
    Mjd_UT1 = Mjd_TT + (UT1_UTC-TT_UTC)/86400;

    P = PrecMatrix(Const::MJD_J2000,Mjd_TT);
    N = NutMatrix(Mjd_TT);
    T = N * P;
    E = PoleMatrix(x_pole,y_pole) * GHAMatrix(Mjd_UT1) * T;

    los3 = E.trans()*los3;
    rsite3 = E.trans()*rsite3;

    magr1old  = 99999999.9;
    magr2old  = 99999999.9;
    magrsite1 = Matrix::norm(rsite1);
    magrsite2 = Matrix::norm(rsite2);
    magrsite3 = Matrix::norm(rsite3);

    cc1 = 2.0*Matrix::dot(los1,rsite1);
    cc2 = 2.0*Matrix::dot(los2,rsite2);
    ktr = 0;

    Matrix r3;
    while (fabs(magr1in-magr1old) > tol || fabs(magr2in-magr2old) > tol) {
        ktr = ktr + 1;
        doubler(r2,r3,f1,f2,q1,magr1,magr2,a,deltae32, cc1,cc2,magrsite1,magrsite2,magr1in,magr2in,
                        los1,los2,los3,rsite1,rsite2,rsite3,t1,t3,direct);

        f  = 1.0 - a/magr2*(1.0-cos(deltae32));
        g  = t3 - sqrt(pow(a,3)/Const::GM_Earth)*(deltae32-sin(deltae32));
        v2 = (r3 - f*r2)/g;

        magr1o = magr1in;
        magr1in = (1.0+pctchg)*magr1in;
        deltar1 = pctchg*magr1in;
        doubler(r2,r3,f1delr1,f2delr1,q2,magr1,magr2,a,deltae32,cc1,cc2,magrsite1,magrsite2,magr1in,magr2in,
                               los1,los2,los3,rsite1,rsite2,rsite3,t1,t3,direct);
        pf1pr1 = (f1delr1-f1)/deltar1;
        pf2pr1 = (f2delr1-f2)/deltar1;

        magr1in = magr1o;
        deltar1 = pctchg*magr1in;
        magr2o = magr2in;
        magr2in = (1.0+pctchg)*magr2in;
        deltar2 = pctchg*magr2in;
        doubler(r2,r3,f1delr2,f2delr2,q3,magr1,magr2,a,deltae32,cc1,cc2,magrsite1,
                                magrsite2,magr1in,magr2in,los1,los2,los3,rsite1,rsite2,rsite3,t1,t3,direct);
        pf1pr2 = (f1delr2-f1)/deltar2;
        pf2pr2 = (f2delr2-f2)/deltar2;

        magr2in = magr2o;
        deltar2 = pctchg*magr2in;

        delta  = pf1pr1*pf2pr2 - pf2pr1*pf1pr2;
        delta1 = pf2pr2*f1 - pf1pr2*f2;
        delta2 = pf1pr1*f2 - pf2pr1*f1;

        deltar1 = -delta1/delta;
        deltar2 = -delta2/delta;

        magr1old = magr1in;
        magr2old = magr2in;

        magr1in = magr1in + deltar1;
        magr2in = magr2in + deltar2;
    }

    doubler(r2,r3,f1,f2,q1,magr1,magr2,a,deltae32, cc1,cc2,magrsite1,magrsite2,magr1in,magr2in,
												       los1,los2,los3,rsite1,rsite2,rsite3,t1,t3,direct);

    f  = 1.0 - a/magr2*(1.0-cos(deltae32));
    g  = t3 - sqrt(pow(a,3)/Const::GM_Earth)*(deltae32-sin(deltae32));
    v2 = (r3 - f*r2)/g;
}

