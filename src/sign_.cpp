//$Header$
//
// EKF_GEOS3
//
// Author: David Ledesma
// Created: 2024/04/19
//
//------------------------------------------------------------------------------
#include "../include/sign_.h"

//---------------------------------
// public methods
//---------------------------------

//------------------------------------------------------------------------------
// double sign_(double a,double b)
//------------------------------------------------------------------------------
/**
 *   returns absolute value of a with sign of b
 *
 * @param <a> Value whose absolute value will be taken
 * @param <b> Value whose sign will be applied to the absolute value of a
 *
 * @return Absolute value of a with sign of b
 *
 */
//------------------------------------------------------------------------------
double sign_(double a,double b){
    double result;
    if (b>=0.0)
        result = fabs(a);
    else
        result = - fabs(a);
    return result;
}
