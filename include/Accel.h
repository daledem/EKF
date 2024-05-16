#ifndef ACCEL_H
#define ACCEL_H

#include "./Matrix.h"
#include "./SAT_Const.h"
#include "./Globals.h"
#include "./IERS.h"
#include "./timediff.h"
#include "./PrecMatrix.h"
#include "./NutMatrix.h"
#include "./PoleMatrix.h"
#include "./GHAMatrix.h"
#include "./Mjday_TDB.h"
#include "./JPL_Eph_DE430.h"
#include "./AccelHarmonic.h"
#include "./AccelPointMass.h"
//------------------------------------------------------------------------------
// Accel(double x,const Matrix& Y)
//------------------------------------------------------------------------------
/**
*   Computes the acceleration of an Earth orbiting satellite due to
*    - the Earth's harmonic gravity field,
*    - the gravitational perturbations of the Sun and Moon
*    - the solar radiation pressure and
*    - the atmospheric drag
*
* @param <Mjd_TT> Terrestrial Time (Modified Julian Date)
* @param <Y> Satellite state vector in the ICRF/EME2000 system
*
* @return Acceleration (a=d^2r/dt^2) in the ICRF/EME2000 system
*/
//------------------------------------------------------------------------------

Matrix Accel(double x,const Matrix& Y);

#endif //ACCEL_H
