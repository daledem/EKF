#include "../include/Position.h"

Matrix Position(double lon, double lat, double h){
    double R_equ = Const::R_Earth;
    double f     = Const::f_Earth;

    double e2     = f*(2.0-f);   // Square of eccentricity
    double CosLat = cos(lat);    // (Co)sine of geodetic latitude
    double SinLat = sin(lat);

    // Position vector
    double N = R_equ / sqrt(1.0-e2*SinLat*SinLat);

    Matrix r(1,3);

    r(1,1) =  (         N+h)*CosLat*cos(lon);
    r(1,2) =  (         N+h)*CosLat*sin(lon);
    r(1,3) =  ((1.0-e2)*N+h)*SinLat;

    return r;
}
