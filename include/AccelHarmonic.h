#ifndef PROYECTO_ACCELHARMONIC_H
#define PROYECTO_ACCELHARMONIC_H

#include <cmath>
#include "./Matrix.h"
#include "./Legendre.h"
#include "./SAT_Const.h"

//------------------------------------------------------------------------------
// AccelHarmonic(const Matrix& r,const Matrix& E,int n_max,int m_max)
//------------------------------------------------------------------------------
/**
*   Computes the acceleration due to the harmonic gravity field of the
*   central body
*
* @param <r> Satellite position vector in the inertial system
* @param <E> Transformation matrix to body-fixed system
* @param <n_max> Maximum degree
* @param <m_max> Maximum order (m_max<=n_max; m_max=0 for zonals, only)
*
* @return Acceleration (a=d^2r/dt^2)
*/
//------------------------------------------------------------------------------

Matrix AccelHarmonic(const Matrix& r,const Matrix& E,int n_max,int m_max);

#endif //PROYECTO_ACCELHARMONIC_H
