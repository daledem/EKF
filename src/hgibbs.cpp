#include "../include/hgibbs.h"

void hgibbs(Matrix &v2, double &theta, double &theta1, double &copa, std::string &error, const Matrix &r1, const Matrix &r2, const Matrix &r3, double Mjd1, double Mjd2, double Mjd3) {
    double magr1,magr2,magr3,tolangle,dt21,dt31,dt32,term1,term2,term3;

    error =  "          ok";
    theta = 0.0;
    theta1= 0.0;
    magr1 = Matrix::norm( r1 );
    magr2 = Matrix::norm( r2 );
    magr3 = Matrix::norm( r3 );

    for (int i = 1; i <= 3; i++)
        v2(i,1)= 0.0;

    tolangle= 0.01745329251994;
    dt21= (Mjd2-Mjd1)*86400.0;
    dt31= (Mjd3-Mjd1)*86400.0;
    dt32= (Mjd3-Mjd2)*86400.0;

    Matrix p = Matrix::cross( r2,r3 );
    Matrix pn = unit( p.trans() );
    Matrix r1n = unit( r1.trans() );
    copa=  asin( Matrix::dot( pn,r1n ) );

    if ( fabs( Matrix::dot(r1n,pn) ) > 0.017452406 )
        error= "not coplanar";

    theta  = angl( r1,r2 );
    theta1 = angl( r2,r3 );

    if ( (theta > tolangle) | (theta1 > tolangle) )
        error= "   angl > 1ø";

    term1= -dt32*( 1.0/(dt21*dt31) + Const::GM_Earth/(12.0*magr1*magr1*magr1) );
    term2= (dt32-dt21)*( 1.0/(dt21*dt32) + Const::GM_Earth/(12.0*magr2*magr2*magr2) );
    term3=  dt21*( 1.0/(dt32*dt31) + Const::GM_Earth/(12.0*magr3*magr3*magr3) );

    v2 =  term1*r1 + term2* r2 + term3* r3;
}

