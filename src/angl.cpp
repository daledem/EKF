//$Header$
//
// EKF_GEOS3
//
// Author: David Ledesma
// Created: 2024/05/23
//
//------------------------------------------------------------------------------
#include "../include/angl.h"

//---------------------------------
// public methods
//---------------------------------

//------------------------------------------------------------------------------
// double angl(const Matrix &vec1, const Matrix &vec2)
//------------------------------------------------------------------------------
/**
 *   Computes the angle between the two vectors
 *
 * @param <vec1> vector 1
 * @param <vec2> vector 2
 *
 * @return Angle between the two vectors  -pi to pi
 */
//------------------------------------------------------------------------------
double angl(const Matrix &vec1, const Matrix &vec2) {
    double small, undefined, magv1,magv2,temp;

    small     = 0.00000001;
    undefined = 999999.1;

    magv1 = Matrix::norm(vec1);
    magv2 = Matrix::norm(vec2);

    if (magv1*magv2 > small^2) {
        temp= Matrix::dot(vec1,vec2) / (magv1*magv2);
        if (fabs( temp ) > 1.0) {
            temp= sign(temp) * 1.0;
        }
        return acos( temp );
    }

    return undefined;
}
