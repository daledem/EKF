//$Header$
//
// EKF_GEOS3
//
// Author: David Ledesma
// Created: 2024/04/24
//
//------------------------------------------------------------------------------
#include "../include/unit.h"

//---------------------------------
// public methods
//---------------------------------

//------------------------------------------------------------------------------
// Matrix unit(const Matrix& vec)
//------------------------------------------------------------------------------
/**
 *   this function calculates a unit vector given the original vector. if a
 *       zero vector is input, the vector is set to zero.
 *
 * @param <vec> vector
 *
 * @return unit vector
 *
 */
//------------------------------------------------------------------------------
Matrix unit(const Matrix& vec){
    double small = 0.000001;
    double magv = Matrix::norm(vec);
    Matrix outvec(1,3);

    if ( magv > small )
        for(int i = 1; i <= 3; i++)
            outvec(1,i) = vec(1,i)/magv;
    else
        for(int i = 1; i <= 3; i++)
            outvec(1,i) = 0.0;

    return outvec;
}