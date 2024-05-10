#include "../include/NutMatrix.h"

Matrix NutMatrix(double Mjd_TT) {
    double eps, dpsi, deps;
    // Mean obliquity of the ecliptic
    eps = MeanObliquity (Mjd_TT);

    // Nutation in longitude and obliquity
    NutAngles (dpsi, deps,Mjd_TT);

    // Transformation from mean to true equator and equinox
    return R_x(-eps-deps)*R_z(-dpsi)*R_x(+eps);

}

