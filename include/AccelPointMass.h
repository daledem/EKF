#ifndef ACCELPOINTMASS_H
#define ACCELPOINTMASS_H

#include <cmath>
#include "./Matrix.h"

//------------------------------------------------------------------------------
// AccelPointMass(const Matrix& r, Matrix& s, double GM)
//------------------------------------------------------------------------------
/**
*   Computes the perturbational acceleration due to a point
*	  mass
*
* @param <r> Satellite position vector
* @param <s> Point mass position vector
* @param <GM> Gravitational coefficient of point mass
*
* @return Acceleration (a=d^2r/dt^2)
*/
//------------------------------------------------------------------------------

Matrix AccelPointMass(const Matrix& r, Matrix& s, double GM);

#endif //ACCELPOINTMASS_H
