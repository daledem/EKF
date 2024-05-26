//$Header$
//
// EKF_GEOS3
//
// Author: David Ledesma
// Created: 2024/04/28
//
//------------------------------------------------------------------------------
#include "../include/Frac.h"

//---------------------------------
// public methods
//---------------------------------

//------------------------------------------------------------------------------
// double Frac(double x)
//------------------------------------------------------------------------------
/**
 *   Fractional part of a number (y=x-[x])
 *
 * @param <x> Number from which the fractional part will be extracted
 *
 * @return Fractional part of x
 *
 */
//------------------------------------------------------------------------------
double Frac(double x) {
    return x-floor(x);
}
