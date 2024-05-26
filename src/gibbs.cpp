//$Header$
//
// EKF_GEOS3
//
// Author: David Ledesma
// Created: 2024/05/23
//
//------------------------------------------------------------------------------
#include "../include/gibbs.h"

//---------------------------------
// public methods
//---------------------------------

//------------------------------------------------------------------------------
// void gibbs(Matrix &v2, double &theta, double &theta1, double &copa,
//              std::string &error, const Matrix &r1, const Matrix &r2,
//              const Matrix &r3)
//------------------------------------------------------------------------------
/**
 *  This function performs the gibbs method of orbit determination. this
 *  method determines the velocity at the middle point of the 3 given
 *  position vectors.
 *
 * @param[out] <v2> ijk velocity vector for r2 [m/s]
 * @param[out] <theta> angl between vectors [rad]
 * @param[out] <theta2> angl between vectors [rad]
 * @param[out] <copa> double ['ok',...]
 * @param[out] <error> flag indicating success
 * @param <r1> ijk position vector #1 [m]
 * @param <r2> ijk position vector #2 [m]
 * @param <r3> ijk position vector #3 [m]
 *
 */
//------------------------------------------------------------------------------
void gibbs(Matrix &v2, double &theta, double &theta1, double &copa, std::string &error, const Matrix &r1, const Matrix &r2, const Matrix &r3) {
    double small,magr1,magr2,magr3,magd,magn,r1mr2,r3mr1,r2mr3,l,tover2;

    small= 0.00000001;
    theta= 0.0;
    error = "          ok";
    theta1= 0.0;

    magr1 = Matrix::norm( r1 );
    magr2 = Matrix::norm( r2 );
    magr3 = Matrix::norm( r3 );
    for (int i = 1; i <= 3; i++)
        v2(i,1)= 0.0;

    Matrix p = Matrix::cross( r2,r3 );
    Matrix q = Matrix::cross( r3,r1 );
    Matrix w = Matrix::cross( r1,r2 );
    Matrix pn = unit( p.trans() );
    Matrix r1n = unit( r1.trans() );
    copa=  asin( Matrix::dot( pn,r1n ) );

    if ( fabs( Matrix::dot(r1n,pn) ) > 0.017452406 )
        error= "not coplanar";


    Matrix d = p + q + w;
    magd = Matrix::norm(d);
    Matrix n = magr1*p + magr2*q + magr3*w;
    magn = Matrix::norm(n);
    Matrix nn = unit( n.trans() );
    Matrix dn = unit( d.trans() );

    // -------------------------------------------------------------
    // determine if  the orbit is possible. both d and n must be in
    // the same direction, and non-zero.
    // -------------------------------------------------------------
    if ( ( fabs(magd)<small ) || ( fabs(magn)<small ) || ( Matrix::dot(nn,dn) < small ) ) {
        error= "  impossible";
    }else {
        theta  = angl( r1,r2 );
        theta1 = angl( r2,r3 );

        // ----------- perform gibbs method to find v2 -----------
        r1mr2= magr1-magr2;
        r3mr1= magr3-magr1;
        r2mr3= magr2-magr3;
        Matrix s  = r1mr2*r3 + r3mr1*r2 + r2mr3*r1;
        Matrix b  = Matrix::cross( d,r2 );
        l  = sqrt(Const::GM_Earth / (magd*magn) );
        tover2 = l / magr2;
        v2 = tover2 * b + l * s;
    }
}

