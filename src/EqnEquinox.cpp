#include "../include/EqnEquinox.h"

double EqnEquinox(double Mjd_TT) {
    // Nutation in longitude and obliquity
    double dpsi, deps;
    NutAngles (dpsi, deps,Mjd_TT);

    // Equation of the equinoxes
    return dpsi * cos ( MeanObliquity(Mjd_TT) );
}
