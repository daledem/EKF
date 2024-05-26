//$Header$
//
// EKF_GEOS3
//
// Author: David Ledesma
// Created: 2024/05/08
//
//------------------------------------------------------------------------------
#include "../include/EqnEquinox.h"

//---------------------------------
// public methods
//---------------------------------

//------------------------------------------------------------------------------
// double EqnEquinox(double Mjd_TT)
//------------------------------------------------------------------------------
/**
*   Computation of the equation of the equinoxes
*
* @param <Mjd_TT> Modified Julian Date (Terrestrial Time)
*
* @return Equation of the equinoxes
*
* @note The equation of the equinoxes dpsi*cos(eps) is the right ascension of
*   the mean equinox referred to the true equator and equinox and is equal
*   to the difference between apparent and mean sidereal time.
*/
//------------------------------------------------------------------------------
double EqnEquinox(double Mjd_TT) {
    // Nutation in longitude and obliquity
    double dpsi, deps;
    NutAngles (dpsi, deps,Mjd_TT);

    // Equation of the equinoxes
    return dpsi * cos ( MeanObliquity(Mjd_TT) );
}
