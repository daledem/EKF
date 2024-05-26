//$Header$
//
// EKF_GEOS3
//
// Author: David Ledesma
// Created: 2024/05/09
//
//------------------------------------------------------------------------------
#include "../include/G_AccelHarmonic.h"

//---------------------------------
// public methods
//---------------------------------

//------------------------------------------------------------------------------
// Matrix G_AccelHarmonic(const Matrix& r,const Matrix& U,int n_max,int m_max)
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
Matrix G_AccelHarmonic(const Matrix &r, const Matrix &U, int n_max, int m_max) {
    double d = 1.0;   // Position increment [m]

    Matrix G(3,3);
    Matrix dr(3,1);

    // Gradient
    for(int i = 1; i <= 3; i++) {
        // Set offset in i-th component of the position vector
        dr = Matrix(3,1);
        dr(i,1) = d;
        // Acceleration difference
        Matrix da = AccelHarmonic ( r+dr/2,U, n_max, m_max ) -
             AccelHarmonic ( r-dr/2,U, n_max, m_max );
        // Derivative with respect to i-th axis
        for(int j = 1; j <= G.getFil(); j++)
            G(j,i) = (da/d)(j,1);
    }

    return G;
}

