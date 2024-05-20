#ifndef JPL_EPH_DE430_H
#define JPL_EPH_DE430_H

#include <cmath>
#include "./Matrix.h"
#include "./Cheb3D.h"

//------------------------------------------------------------------------------
// JPL_Eph_DE430(Matrix& r_Mercury,Matrix& r_Venus,Matrix& r_Earth, Matrix& r_Mars,Matrix& r_Jupiter,Matrix& r_Saturn,
//      Matrix &r_Uranus, Matrix&r_Neptune, Matrix&r_Pluto, Matrix&r_Moon, Matrix&r_Sun,double Mjd_TDB)
//------------------------------------------------------------------------------
/**
*   Computes the sun, moon, and nine major planets' equatorial
*                position using JPL Ephemerides
*
* @param[out] <r_Mercury> Geocentric equatorial position ([m]) referred to the International Celestial Reference Frame (ICRF)
* @param[out] <r_Venus> Geocentric equatorial position ([m]) referred to the International Celestial Reference Frame (ICRF)
* @param[out] <r_Earth> Geocentric equatorial position ([m]) referred to the International Celestial Reference Frame (ICRF)
*   (solar system barycenter (SSB))
* @param[out] <r_Mars> Geocentric equatorial position ([m]) referred to the International Celestial Reference Frame (ICRF)
* @param[out] <r_Jupiter> Geocentric equatorial position ([m]) referred to the International Celestial Reference Frame (ICRF)
* @param[out] <r_Saturn> Geocentric equatorial position ([m]) referred to the International Celestial Reference Frame (ICRF)
* @param[out] <r_Uranus> Geocentric equatorial position ([m]) referred to the International Celestial Reference Frame (ICRF)
* @param[out] <r_Neptune> Geocentric equatorial position ([m]) referred to the International Celestial Reference Frame (ICRF)
* @param[out] <r_Pluto> Geocentric equatorial position ([m]) referred to the International Celestial Reference Frame (ICRF)
* @param[out] <r_Moon> Geocentric equatorial position ([m]) referred to the International Celestial Reference Frame (ICRF)
* @param[out] <r_Sun> Geocentric equatorial position ([m]) referred to the International Celestial Reference Frame (ICRF)
* @param <Mjd_TDB> Modified julian date of TDB
*
* @note Light-time is already taken into account
*
*/
//------------------------------------------------------------------------------

void JPL_Eph_DE430(Matrix& r_Mercury,Matrix& r_Venus,Matrix& r_Earth, Matrix& r_Mars,Matrix& r_Jupiter,Matrix& r_Saturn,
    Matrix& r_Uranus,Matrix& r_Neptune,Matrix& r_Pluto,Matrix& r_Moon,Matrix& r_Sun, double Mjd_TDB);

#endif //JPL_EPH_DE430_H
