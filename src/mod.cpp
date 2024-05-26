//$Header$
//
// EKF_GEOS3
//
// Author: David Ledesma
// Created: 2024/05/23
//
//------------------------------------------------------------------------------
#include "../include/mod.h"

//---------------------------------
// public methods
//---------------------------------

//------------------------------------------------------------------------------
// double Mjday_TDB(double Mjd_TT)
//------------------------------------------------------------------------------
/**
 *   Computes the absolute mod of two doubles
 *
 * @param <a> double
 * @param <b> double
 *
 * @return Mod of a / b
 *
 * @note This code was taken from https://stackoverflow.com/questions/28888619/modulo-function-in-c-that-behaves-like-mod-in-matlab
 *
 */
//------------------------------------------------------------------------------

double mod(double a, double b) {
    double result = fmod(a, b);
    return result >= 0 ? result : result + b;
}
