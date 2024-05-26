//$Header$
//
// EKF_GEOS3
//
// Author: David Ledesma
// Created: 2024/04/28
//
//------------------------------------------------------------------------------
#include "../include/AccelPointMass.h"

//---------------------------------
// public methods
//---------------------------------

//------------------------------------------------------------------------------
// Matrix AccelPointMass(const Matrix& r, Matrix& s, double GM)
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
Matrix AccelPointMass(const Matrix& r, Matrix& s, double GM) {

    Matrix d = r - s;

    // Acceleration 
    return -GM * ( d/(pow(Matrix::norm(d),3)) + s/(pow(Matrix::norm(s),3)) );
}


