#include "../include/gmst.h"

double gmst(double Mjd_UT1) {
    double Secs,MJD_J2000,Mjd_0,UT1,T_0,T,gmst;

    Secs = 86400.0;                       // Seconds per day
    MJD_J2000 = 51544.5;

    Mjd_0 = floor(Mjd_UT1);
    UT1   = Secs*(Mjd_UT1-Mjd_0);         // [s]
    T_0   = (Mjd_0  -MJD_J2000)/36525.0;
    T     = (Mjd_UT1-MJD_J2000)/36525.0;

    gmst  = 24110.54841 + 8640184.812866*T_0 + 1.002737909350795*UT1 + (0.093104-6.2e-6*T)*T*T;    // [s]

    return 2*Const::pi*Frac(gmst/Secs);       // [rad], 0..2pi
}
