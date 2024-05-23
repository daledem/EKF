#include "../include/anglesg.h"

extern Matrix eopdata;

void anglesg(Matrix &r2, Matrix &v2, double az1, double az2, double az3, double el1, double el2, double el3, double Mjd1, double Mjd2, double Mjd3,Matrix Rs1,Matrix Rs2,Matrix Rs3) {
    double lon1, lat1, h1,lon2, lat2, h2,lon3, lat3, h3,Mjd_UTC,x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,
            TAI_UTC,UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC,Mjd_TT,Mjd_UT1,tau1,tau3,a1,a3,b1,b3,d1s,d2s,Ccye,bigr2,u,
            C1,C2,C3,rho1,rho2,rho3,rhoold1,rhoold2,rhoold3,ll,magr1,magr2,magr3,theta,theta1,copa,p, a, e, i, Omega,
            omega, M,rdot,udot,tausqr,f1,g1,f3,g3,H1,H2,H3,G1,G2,G3,D1,D2,D3;
    std::string error;

    double vL1[] = {cos(el1)*sin(az1), cos(el1)*cos(az1), sin(el1)};
    Matrix L1(3,1,vL1,3);

    double vL2[] = {cos(el2)*sin(az2), cos(el2)*cos(az2), sin(el2)};
    Matrix L2(3,1,vL2,3);

    double vL3[] = {cos(el3)*sin(az3), cos(el3)*cos(az3), sin(el3)};
    Matrix L3(3,1,vL3,3);


    Geodetic(lon1, lat1, h1,Rs1.trans());
    Geodetic(lon2, lat2, h2,Rs2.trans());
    Geodetic(lon3, lat3, h3,Rs3.trans());

    Matrix M1 = LTC(lon1, lat1);
    Matrix M2 = LTC(lon2, lat2);
    Matrix M3 = LTC(lon3, lat3);

    // body-fixed system
    Matrix Lb1 = M1.trans()*L1;
    Matrix Lb2 = M1.trans()*L2;
    Matrix Lb3 = M1.trans()*L3;

    // mean of date system (J2000)
    Mjd_UTC = Mjd1;
    IERS(x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC,eopdata,Mjd_UTC,'l');
    timediff(UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC,UT1_UTC,TAI_UTC);
    Mjd_TT = Mjd_UTC + TT_UTC/86400;
    Mjd_UT1 = Mjd_TT + (UT1_UTC-TT_UTC)/86400;

    Matrix P = PrecMatrix(Const::MJD_J2000,Mjd_TT);
    Matrix N = NutMatrix(Mjd_TT);
    Matrix T = N * P;
    Matrix E = PoleMatrix(x_pole,y_pole) * GHAMatrix(Mjd_UT1) * T;

    Matrix Lm1 = E.trans()*Lb1;
    Rs1 = E.trans()*Rs1;

    Mjd_UTC = Mjd2;
    IERS(x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC,eopdata,Mjd_UTC,'l');
    timediff(UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC,UT1_UTC,TAI_UTC);
    Mjd_TT = Mjd_UTC + TT_UTC/86400;
    Mjd_UT1 = Mjd_TT + (UT1_UTC-TT_UTC)/86400;

    P = PrecMatrix(Const::MJD_J2000,Mjd_TT);
    N = NutMatrix(Mjd_TT);
    T = N * P;
    E = PoleMatrix(x_pole,y_pole) * GHAMatrix(Mjd_UT1) * T;

    Matrix Lm2 = E.trans()*Lb2;
    Rs2 = E.trans()*Rs2;

    Mjd_UTC = Mjd3;
    IERS(x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC,eopdata,Mjd_UTC,'l');
    timediff(UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC,UT1_UTC,TAI_UTC);
    Mjd_TT = Mjd_UTC + TT_UTC/86400;
    Mjd_UT1 = Mjd_TT + (UT1_UTC-TT_UTC)/86400.0;

    P = PrecMatrix(Const::MJD_J2000,Mjd_TT);
    N = NutMatrix(Mjd_TT);
    T = N * P;
    E = PoleMatrix(x_pole,y_pole) * GHAMatrix(Mjd_UT1) * T;

    Matrix Lm3 = E.trans()*Lb3;
    Rs3 = E.trans()*Rs3;

    // geocentric inertial position
    tau1 = (Mjd1-Mjd2)*86400;
    tau3 = (Mjd3-Mjd2)*86400;

    a1 = tau3/(tau3-tau1);
    a3 =-tau1/(tau3-tau1);

    b1 = tau3/(6*(tau3-tau1))*(pow((tau3-tau1),2)-pow(tau3,2));
    b3 =-tau1/(6*(tau3-tau1))*(pow((tau3-tau1),2)-pow(tau1,2));

    Matrix D = Lm1.append(Lm2).append(Lm3).inverse()*Rs1.append(Rs2).append(Rs2);

    d1s = D(2,1)*a1-D(2,2)+D(2,3)*a3;
    d2s = D(2,1)*b1+D(2,3)*b3;

    Ccye = 2*Matrix::dot(Lm2,Rs2);

    double poly[10], zeror[10], zeroi[10];
    int nr,degree =8;

    poly[0]=  1.0;  // R2^8... polynomial
    poly[1]=  0.0;
    poly[2]=  -(pow(d1s,2) + d1s*Ccye + pow((Matrix::norm(Rs2)),2));
    poly[3]=  0.0;
    poly[4]=  0.0;
    poly[5]=  -Const::GM_Earth*(d2s*Ccye + 2*d1s*d2s);
    poly[6]=  0.0;
    poly[7]=  0.0;
    poly[8]=  -pow(Const::GM_Earth,2)*pow(d2s,2);
    nr = real_poly_roots( poly,degree,zeror,zeroi );

    bigr2= -99999990.0;

    for (int j = 0; j < nr; j++) {
        if (( zeror[j] > bigr2 ) && ( fabs(zeroi[j]) < 10e-12)){
            bigr2 = zeror[j];
        }
    }

    u = Const::GM_Earth/(pow(bigr2,3));

    C1 = a1+b1*u;
    C2 = -1;
    C3 = a3+b3*u;

    double vaux[] = {C1,C2,C3};
    Matrix aux(1,3,vaux,3);
    Matrix temp = -D*aux.trans();
    rho1 = temp(1,1)/(a1+b1*u);
    rho2 = -temp(2,1);
    rho3 = temp(3,1)/(a3+b3*u);

    rhoold1 = rho1;
    rhoold2 = rho2;
    rhoold3 = rho3;

    rho2 = 99999999.9;
    ll   = 0;

    Matrix r1,r3;
    while ((fabs(rhoold2-rho2) > 1e-12) && (ll <= 0 )) {
        ll = ll + 1;
        rho2 = rhoold2;

        r1 = Rs1+rho1*Lm1;
        r2 = Rs2+rho2*Lm2;
        r3 = Rs3+rho3*Lm3;

        magr1 = Matrix::norm(r1);
        magr2 = Matrix::norm(r2);
        magr3 = Matrix::norm(r3);

        gibbs(v2, theta,theta1,copa,error,r1,r2,r3);

        if ( (error != "          ok") & (copa < Const::pi/180) )
            hgibbs(v2,theta,theta1,copa,error,r1,r2,r3,Mjd1,Mjd2,Mjd3);


        elements (p, a, e, i, Omega, omega, M,r2.append(v2));

        if ( ll <= 8 ) {
            u = Const::GM_Earth/pow(magr2,3);
            rdot= Matrix::dot(r2,v2)/magr2;
            udot= (-3*Const::GM_Earth*rdot)/(pow(magr2,4));

            tausqr= tau1*tau1;
            f1=  1 - 0.5*u*tausqr -(1./6)*udot*tausqr*tau1
                - (1./24) * u*u*tausqr*tausqr
                - (1./30)*u*udot*tausqr*tausqr*tau1;
            g1= tau1 - (1./6)*u*tau1*tausqr - (1./12) * udot*tausqr*tausqr
                - (1./120)*u*u*tausqr*tausqr*tau1
                - (1./120)*u*udot*tausqr*tausqr*tausqr;
            tausqr= tau3*tau3;
            f3=  1 - 0.5*u*tausqr -(1./6)*udot*tausqr*tau3
                - (1./24) * u*u*tausqr*tausqr
                - (1./30)*u*udot*tausqr*tausqr*tau3;
            g3= tau3 - (1./6)*u*tau3*tausqr - (1./12) * udot*tausqr*tausqr
                - (1./120)*u*u*tausqr*tausqr*tau3
                - (1./120)*u*udot*tausqr*tausqr*tausqr;
        }else {
            theta  = angl( r1,r2 );
            theta1 = angl( r2,r3 );

            f1= 1 - ( (magr1*(1 - cos(theta)) / p ) );
            g1= ( magr1*magr2*sin(-theta) ) / sqrt( p );
            f3= 1 - ( (magr3*(1 - cos(theta1)) / p ) );
            g3= ( magr3*magr2*sin(theta1) ) / sqrt( p );
        }

        C1 = g3/(f1*g3-f3*g1);
        C2 = -1;
        C3 =-g1/(f1*g3-f3*g1);

        H1 = Const::GM_Earth*tau3/12;
        H3 =-Const::GM_Earth*tau1/12;
        H2 = H1-H3;

        G1 = -tau3/(tau1*(tau3-tau1));
        G3 = -tau1/(tau3*(tau3-tau1));
        G2 = G1-G3;

        D1 = G1+H1/pow(magr1,3);
        D2 = G2+H2/pow(magr2,3);
        D3 = G3+H3/pow(magr3,3);

        /*
        temp = -[D1 D2 D3]*[C1 C2 C3]';
        rhoold1 = temp/(a1+b1*u);
        rhoold2 = -temp;
        rhoold3 = temp/(a3+b3*u);

        r1 = Rs1+rhoold1*Lm1;
        r2 = Rs2+rhoold2*Lm2;
        r3 = Rs3+rhoold3*Lm3;
        */
    }

    r1 = Rs1+rho1*Lm1;
    r2 = Rs2+rho2*Lm2;
    r3 = Rs3+rho3*Lm3;
}

