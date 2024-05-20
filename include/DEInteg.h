#ifndef DEINTEG_H
#define DEINTEG_H

#include <cmath>
#include <algorithm>
#include "./Matrix.h"
#include "./sign_.h"
#include "./SAT_Const.h"

//------------------------------------------------------------------------------
// DEInteg(Matrix (*func)(double x,const Matrix& Y),double t,double tout,double relerr,double abserr,int n_eqn,Matrix y)
//------------------------------------------------------------------------------
/**
*   Numerical integration methods for ordinaray differential equations
*
*   This module provides implemenation of the variable order variable
*   stepsize multistep method of Shampine & Gordon.
*
* @param <func> function to compute derivative Y(x)
* @param <t> Input time
* @param <tout> Output time
* @param <relerr> Relative error
* @param <abserr> Absolute error
* @param <n_eqn>
* @param <y>
*
* @return y
*/
//------------------------------------------------------------------------------


Matrix DEInteg(Matrix (*func)(double x,const Matrix& Y),double t,double tout,double relerr,double abserr,int n_eqn,Matrix y);

#endif //DEINTEG_H
