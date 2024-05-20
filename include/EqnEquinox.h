#ifndef EQNEQUINOX_H
#define EQNEQUINOX_H

#include <cmath>
#include "./NutAngles.h"
#include "./MeanObliquity.h"

//------------------------------------------------------------------------------
// EqnEquinox(double Mjd_TT)
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


double EqnEquinox(double Mjd_TT);

#endif //EQNEQUINOX_H
