#ifndef G_ACCELHARMONIC_H
#define G_ACCELHARMONIC_H

#include "./Matrix.h"
#include "./AccelHarmonic.h"

//------------------------------------------------------------------------------
// G_AccelHarmonic(const Matrix& r,const Matrix& U,int n_max,int m_max)
//------------------------------------------------------------------------------
/**
*   Computes the gradient of the Earth's harmonic gravity field
*
* @param <r> Satellite position vector in the true-of-date system
* @param <U> Transformation matrix to body-fixed system
* @param <n_max> Gravity model degree
* @param <m_max> Gravity model order
*
* @return Gradient (G=da/dr) in the true-of-date system
*
*/
//------------------------------------------------------------------------------

Matrix G_AccelHarmonic(const Matrix& r,const Matrix& U,int n_max,int m_max);

#endif //G_ACCELHARMONIC_H
