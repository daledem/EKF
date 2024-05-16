#ifndef CHEB3D_H
#define CHEB3D_H

#include "./Matrix.h"

//------------------------------------------------------------------------------
// Cheb3D(double t,int N,double Ta,double Tb,const Matrix& Cx,const Matrix& Cy,const Matrix& Cz)
//------------------------------------------------------------------------------
/**
*   Chebyshev approximation of 3-dimensional vectors
*
* @param <t> interval
* @param <N> Number of coefficients
* @param <Ta> Begin interval
* @param <Tb> End interval
* @param <Cx> Coefficients of Chebyshev polyomial (x-coordinate)
* @param <Cy> Coefficients of Chebyshev polyomial (y-coordinate)
* @param <Cz> Coefficients of Chebyshev polyomial (z-coordinate)
*
* @return ChebApp
*/
//------------------------------------------------------------------------------

Matrix Cheb3D(double t,int N,double Ta,double Tb,const Matrix& Cx,const Matrix& Cy,const Matrix& Cz);

#endif //CHEB3D_H
