#include "../include/doubler.h"

void doubler(Matrix &r2, Matrix &r3, double &f1, double &f2, double &q1, double &magr1, double &magr2, double &a,
    double &deltae32, double cc1, double cc2, double magrsite1, double magrsite2, double magr1in, double magr2in,
    const Matrix &los1, const Matrix &los2, const Matrix &los3, const Matrix &rsite1, const Matrix &rsite2,
    const Matrix &rsite3, double t1, double t3, char direct) {
    double rho1,rho2,rho3,magr3,cosdv21,sindv21,dv21,cosdv31,sindv31,dv31,cosdv32,sindv32,dv32,c1,c3,p,ecosv1,ecosv2,
            ecosv3,esinv2,e,n,s,c,sinde32,cosde32,sinde21,cosde21,deltae21,deltam32,deltam12,sindh32,sindh21,deltah32,
            deltah21;

    rho1 = (-cc1 + sqrt(pow(cc1,2)-4*(pow(magrsite1,2)-pow(magr1in,2)))) / 2.0;
    rho2 = (-cc2 + sqrt(pow(cc2,2)-4*(pow(magrsite2,2)-pow(magr2in,2)))) / 2.0;

    Matrix r1 = rho1*los1 + rsite1;
    r2 = rho2*los2 + rsite2;

    magr1 = Matrix::norm(r1);
    magr2 = Matrix::norm(r2);

    Matrix w;
    if (direct == 'y')
        w = Matrix::cross(r1,r2)/(magr1*magr2);
    else
        w = -Matrix::cross(r1,r2)/(magr1*magr2);

    rho3 =  -Matrix::dot(rsite3,w)/Matrix::dot(los3,w);
    r3 = rho3*los3 + rsite3;
    magr3 = Matrix::norm(r3);

    cosdv21 = Matrix::dot(r2,r1)/(magr2*magr1);
    sindv21 = Matrix::norm(Matrix::cross(r2,r1))/(magr2*magr1);
    dv21 = atan2(sindv21,cosdv21);

    cosdv31 = Matrix::dot(r3,r1)/(magr3*magr1);
    sindv31 = sqrt(1.0 - pow(cosdv31,2));
    dv31 = atan2(sindv31,cosdv31);

    cosdv32 = Matrix::dot(r3,r2)/(magr3*magr2);
    sindv32 = Matrix::norm(Matrix::cross(r3,r2))/(magr3*magr2);
    dv32 = atan2(sindv32,cosdv32);

    if (dv31 > Const::pi) {
        c1 = (magr2*sindv32)/(magr1*sindv31);
        c3 = (magr2*sindv21)/(magr3*sindv31);
        p = (c1*magr1+c3*magr3-magr2)/(c1+c3-1);
    }else {
        c1 = (magr1*sindv31)/(magr2*sindv32);
        c3 = (magr1*sindv21)/(magr3*sindv32);
        p = (c3*magr3-c1*magr2+magr1)/(-c1+c3+1);
    }


    ecosv1 = p/magr1-1;
    ecosv2 = p/magr2-1;
    ecosv3 = p/magr3-1;

    if (dv21 != Const::pi)
        esinv2 = (-cosdv21*ecosv2+ecosv1)/sindv21;
    else
        esinv2 = (cosdv32*ecosv2-ecosv3)/sindv31;


    e = sqrt(pow(ecosv2,2)+pow(esinv2,2));
    a = p/(1-pow(e,2));

    if (e*e < 0.99) {
        n = sqrt((Const::GM_Earth)/pow(a,3));

        s = magr2/p*sqrt(1-pow(e,2))*esinv2;
        c = magr2/p*(pow(e,2)+ecosv2);

        sinde32 = magr3/sqrt(a*p)*sindv32-magr3/p*(1-cosdv32)*s;
        cosde32 = 1-magr2*magr3/(a*p)*(1-cosdv32);
        deltae32 = atan2(sinde32,cosde32);

        sinde21 = magr1/sqrt(a*p)*sindv21+magr1/p*(1-cosdv21)*s;
        cosde21 = 1-magr2*magr1/(a*p)*(1-cosdv21);
        deltae21 = atan2(sinde21,cosde21);

        deltam32 = deltae32+2*s*pow((sin(deltae32/2)),2)-c*sin(deltae32);
        deltam12 = -deltae21+2*s*pow((sin(deltae21/2)),2)+c*sin(deltae21);
    }else {
        n = sqrt((Const::GM_Earth)/pow(-a,3));

        s = magr2/p*sqrt(pow(e,2)-1)*esinv2;
        c = magr2/p*(pow(e,2)+ecosv2);

        sindh32 = magr3/sqrt(-a*p)*sindv32-magr3/p*(1-cosdv32)*s;
        sindh21 = magr1/sqrt(-a*p)*sindv21+magr1/p*(1-cosdv21)*s;

        deltah32 = log( sindh32 + sqrt(pow(sindh32,2) +1) );
        deltah21 = log( sindh21 + sqrt(pow(sindh21,2) +1) );

        deltam32 = -deltah32+2*s*pow((sinh(deltah32/2)),2)+c*sinh(deltah32);
        deltam12 = deltah21+2*s*pow((sinh(deltah21/2)),2)-c*sinh(deltah21);

        deltae32=deltah32;
    }

    f1 = t1-deltam12/n;
    f2 = t3-deltam32/n;

    q1 = sqrt(pow(f1,2)+pow(f2,2));
}

