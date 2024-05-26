//$Header$
//
// EKF_GEOS3
//
// Author: David Ledesma
// Created: 2024/04/28
//
//------------------------------------------------------------------------------
#include "../include/EccAnom.h"

//---------------------------------
// public methods
//---------------------------------

//------------------------------------------------------------------------------
// double EccAnom(double M, double e)
//------------------------------------------------------------------------------
/**
 *   Computes the eccentric anomaly for elliptic orbits
 *
 * @param <M> Mean anomaly in [rad]
 * @param <e> Eccentricity of the orbit [0,1]
 *
 * @return Eccentric anomaly in [rad]
 */
//------------------------------------------------------------------------------
double EccAnom(double M, double e) {
    int maxit = 15;
    double i = 1;

    // Starting value
    M = fmod(M,2.0*Const::pi);


    double E;
    if (e<0.8)
        E = M; 
    else
        E = Const::pi;

    double f = E - e*sin(E) - M;
    E = E - f / ( 1.0 - e*cos(E) );

    // Iteration
    while (fabs(f) > 1e2*Const::eps){
        f = E - e*sin(E) - M;
        E = E - f / ( 1.0 - e*cos(E) );
        i = i+1;
        if (i==maxit)
            throw "convergence problems in EccAnom";
    }
    return E;
}
