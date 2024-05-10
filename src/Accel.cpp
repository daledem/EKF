#include "../include/Accel.h"

extern struct auxParam AuxParam;
extern Matrix eopdata;

Matrix Accel(double x, const Matrix &Y) {
    double x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC,UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC,Mjd_UT1,
        Mjd_TT,MJD_TDB;

    IERS(x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC,eopdata,AuxParam.Mjd_UTC + x/86400,'l');
    timediff(UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC,UT1_UTC,TAI_UTC);
    Mjd_UT1 = AuxParam.Mjd_UTC + x/86400 + UT1_UTC/86400;
    Mjd_TT = AuxParam.Mjd_UTC + x/86400 + TT_UTC/86400;

    Matrix P = PrecMatrix(Const::MJD_J2000,Mjd_TT);
    Matrix N = NutMatrix(Mjd_TT);
    Matrix T = N * P;
    Matrix E = PoleMatrix(x_pole,y_pole) * GHAMatrix(Mjd_UT1) * T;

    MJD_TDB = Mjday_TDB(Mjd_TT);

    Matrix r_Mercury(3,1);
    Matrix r_Venus(3,1);
    Matrix r_Earth(3,1);
    Matrix r_Mars(3,1);
    Matrix r_Jupiter(3,1);
    Matrix r_Saturn(3,1);
    Matrix r_Uranus(3,1);
    Matrix r_Neptune(3,1);
    Matrix r_Pluto(3,1);
    Matrix r_Moon(3,1);
    Matrix r_Sun(3,1);

    JPL_Eph_DE430(r_Mercury,r_Venus,r_Earth,r_Mars,r_Jupiter,r_Saturn,r_Uranus,r_Neptune,r_Pluto,r_Moon,r_Sun,MJD_TDB);

    // Acceleration due to harmonic gravity field
    Matrix a = AccelHarmonic(Y.getColumnaByIndex(1,1,3), E, AuxParam.n, AuxParam.m);

    // Luni-solar perturbations
    if (AuxParam.sun)
        a = a + AccelPointMass(Y.getColumnaByIndex(1,1,3),r_Sun,Const::GM_Sun);

    if (AuxParam.moon)
        a = a + AccelPointMass(Y.getColumnaByIndex(1,1,3),r_Moon,Const::GM_Moon);

    // Planetary perturbations
    if (AuxParam.planets) {
        a = a + AccelPointMass(Y.getColumnaByIndex(1,1,3),r_Mercury,Const::GM_Mercury);
        a = a + AccelPointMass(Y.getColumnaByIndex(1,1,3),r_Venus,Const::GM_Venus);
        a = a + AccelPointMass(Y.getColumnaByIndex(1,1,3),r_Mars,Const::GM_Mars);
        a = a + AccelPointMass(Y.getColumnaByIndex(1,1,3),r_Jupiter,Const::GM_Jupiter);
        a = a + AccelPointMass(Y.getColumnaByIndex(1,1,3),r_Saturn,Const::GM_Saturn);
        a = a + AccelPointMass(Y.getColumnaByIndex(1,1,3),r_Uranus,Const::GM_Uranus);
        a = a + AccelPointMass(Y.getColumnaByIndex(1,1,3),r_Neptune,Const::GM_Neptune);
        a = a + AccelPointMass(Y.getColumnaByIndex(1,1,3),r_Pluto,Const::GM_Pluto);
    }

    return  Y.getColumnaByIndex(1,4,6).join(a);
}
