#include "../include/JPL_Eph_DE430.h"

#include <cstdio>
#include <iostream>

extern Matrix PC;

void JPL_Eph_DE430(Matrix &r_Mercury, Matrix &r_Venus, Matrix &r_Earth, Matrix &r_Mars, Matrix &r_Jupiter, Matrix &r_Saturn, Matrix &r_Uranus, Matrix &r_Neptune, Matrix &r_Pluto, Matrix &r_Moon, Matrix &r_Sun, double Mjd_TDB) {
    double JD, t1, dt, Mjd0,EMRAT,EMRAT1;
    int i = 0,j;

    JD = Mjd_TDB + 2400000.5;

    // Apaño para sustituir el metodo find de matlab
    int k = 1;

    while(k < PC.getFil() && i == 0){
        if(PC(k,1) <= JD && JD <= PC(k,2)) {
            i = k;
        }
        k++;
    }

    Matrix PCtemp = PC.getFilaByIndex(i);

    t1 = PCtemp(1,1)-2400000.5; // MJD at start of interval

    dt = Mjd_TDB - t1;

    double vtemp1[] = {231,244,257,270};
    Matrix temp(1,4,vtemp1,4);
    Matrix Cx_Earth = PCtemp.getFilaByIndex(1,temp(1,1),temp(1,2)-1);
    Matrix Cy_Earth = PCtemp.getFilaByIndex(1,temp(1,2),temp(1,3)-1);
    Matrix Cz_Earth = PCtemp.getFilaByIndex(1,temp(1,3),temp(1,4)-1);
    temp = temp + 39.;
    Matrix Cx = PCtemp.getFilaByIndex(1,temp(1,1),temp(1,2)-1);
    Matrix Cy = PCtemp.getFilaByIndex(1,temp(1,2),temp(1,3)-1);
    Matrix Cz = PCtemp.getFilaByIndex(1,temp(1,3),temp(1,4)-1);
    Cx_Earth = Cx_Earth.append(Cx);
    Cy_Earth = Cy_Earth.append(Cy);
    Cz_Earth = Cz_Earth.append(Cz);
    if (0<=dt && dt<=16) {
        j=0;
        Mjd0 = t1;
    }else if(16<dt && dt<=32) {
        j=1;
        Mjd0 = t1+16*j;
    }
    r_Earth = 1e3*Cheb3D(Mjd_TDB, 13, Mjd0, Mjd0+16, Cx_Earth.getFilaByIndex(1,13*j+1,13*j+13),
                         Cy_Earth.getFilaByIndex(1,13*j+1,13*j+13), Cz_Earth.getFilaByIndex(1,13*j+1,13*j+13)).trans();

    double vtemp2[] = {441,454,467,480};
    temp = Matrix(1,4,vtemp2,4);
    Matrix Cx_Moon = PCtemp.getFilaByIndex(1,temp(1,1),temp(1,2)-1);
    Matrix Cy_Moon = PCtemp.getFilaByIndex(1,temp(1,2),temp(1,3)-1);
    Matrix Cz_Moon = PCtemp.getFilaByIndex(1,temp(1,3),temp(1,4)-1);
    for(i = 1; i <= 7; i++) {
        temp = temp+39;
        Cx = PCtemp.getFilaByIndex(1,temp(1,1),temp(1,2)-1);
        Cy = PCtemp.getFilaByIndex(1,temp(1,2),temp(1,3)-1);
        Cz = PCtemp.getFilaByIndex(1,temp(1,3),temp(1,4)-1);
        Cx_Moon = Cx_Moon.append(Cx);
        Cy_Moon = Cy_Moon.append(Cy);
        Cz_Moon = Cz_Moon.append(Cz);
    }
    if (0<=dt && dt<=4) {
        j=0;
        Mjd0 = t1;
    }else if(4<dt && dt<=8){
        j=1;
        Mjd0 = t1+4*j;
    }else if(8<dt && dt<=12){
        j=2;
        Mjd0 = t1+4*j;
    }else if(12<dt && dt<=16){
        j=3;
        Mjd0 = t1+4*j;
    }else if(16<dt && dt<=20){
        j=4;
        Mjd0 = t1+4*j;
    }else if(20<dt && dt<=24){
        j=5;
        Mjd0 = t1+4*j;
    }else if(24<dt && dt<=28){
        j=6;
        Mjd0 = t1+4*j;
    }else if(28<dt && dt<=32){
        j=7;
        Mjd0 = t1+4*j;
    }
    r_Moon = 1e3*Cheb3D(Mjd_TDB, 13, Mjd0, Mjd0+4, Cx_Moon.getFilaByIndex(1,13*j+1,13*j+13),
                        Cy_Moon.getFilaByIndex(1,13*j+1,13*j+13), Cz_Moon.getFilaByIndex(1,13*j+1,13*j+13)).trans();

    double vtemp3[] = {753,764,775,786};
    temp = Matrix(1,4,vtemp3,4);
    Matrix Cx_Sun = PCtemp.getFilaByIndex(1,temp(1,1),temp(1,2)-1);
    Matrix Cy_Sun = PCtemp.getFilaByIndex(1,temp(1,2),temp(1,3)-1);
    Matrix Cz_Sun = PCtemp.getFilaByIndex(1,temp(1,3),temp(1,4)-1);
    temp = temp+33;
    Cx = PCtemp.getFilaByIndex(1,temp(1,1),temp(1,2)-1);
    Cy = PCtemp.getFilaByIndex(1,temp(1,2),temp(1,3)-1);
    Cz = PCtemp.getFilaByIndex(1,temp(1,3),temp(1,4)-1);
    Cx_Sun = Cx_Sun.append(Cx);
    Cy_Sun = Cy_Sun.append(Cy);
    Cz_Sun = Cz_Sun.append(Cz);
    if (0<=dt && dt<=16){
        j=0;
        Mjd0 = t1;
    }else if(16<dt && dt<=32){
        j=1;
        Mjd0 = t1+16*j;
    }
    r_Sun = 1e3*Cheb3D(Mjd_TDB, 11, Mjd0, Mjd0+16, Cx_Sun.getFilaByIndex(1,11*j+1,11*j+11),
                       Cy_Sun.getFilaByIndex(1,11*j+1,11*j+11), Cz_Sun.getFilaByIndex(1,11*j+1,11*j+11)).trans();


    double vtemp4[] = {3,17,31,45};
    temp = Matrix(1,4,vtemp4,4);
    Matrix Cx_Mercury = PCtemp.getFilaByIndex(1,temp(1,1),temp(1,2)-1);
    Matrix Cy_Mercury = PCtemp.getFilaByIndex(1,temp(1,2),temp(1,3)-1);
    Matrix Cz_Mercury = PCtemp.getFilaByIndex(1,temp(1,3),temp(1,4)-1);
    for (i = 1; i <= 3; i++) {
        temp = temp+42;
        Cx = PCtemp.getFilaByIndex(1,temp(1,1),temp(1,2)-1);
        Cy = PCtemp.getFilaByIndex(1,temp(1,2),temp(1,3)-1);
        Cz = PCtemp.getFilaByIndex(1,temp(1,3),temp(1,4)-1);
        Cx_Mercury = Cx_Mercury.append(Cx);
        Cy_Mercury = Cy_Mercury.append(Cy);
        Cz_Mercury = Cz_Mercury.append(Cz);
    }
    if (0<=dt && dt<=8){
        j=0;
        Mjd0 = t1;
    }else if(8<dt && dt<=16){
        j=1;
        Mjd0 = t1+8*j;
    }else if (16<dt && dt<=24){
        j=2;
        Mjd0 = t1+8*j;
    }else if(24<dt && dt<=32){
        j=3;
        Mjd0 = t1+8*j;
    }
    r_Mercury = 1e3*Cheb3D(Mjd_TDB, 14, Mjd0, Mjd0+8, Cx_Mercury.getFilaByIndex(1,14*j+1,14*j+14),
                           Cy_Mercury.getFilaByIndex(1,14*j+1,14*j+14), Cz_Mercury.getFilaByIndex(1,14*j+1,14*j+14)).trans();

    double vtemp5[] = {171,181,191,201};
    temp = Matrix(1,4,vtemp5,4);
    Matrix Cx_Venus = PCtemp.getFilaByIndex(1,temp(1,1),temp(1,2)-1);
    Matrix Cy_Venus = PCtemp.getFilaByIndex(1,temp(1,2),temp(1,3)-1);
    Matrix Cz_Venus = PCtemp.getFilaByIndex(1,temp(1,3),temp(1,4)-1);
    temp = temp+30;
    Cx = PCtemp.getFilaByIndex(1,temp(1,1),temp(1,2)-1);
    Cy = PCtemp.getFilaByIndex(1,temp(1,2),temp(1,3)-1);
    Cz = PCtemp.getFilaByIndex(1,temp(1,3),temp(1,4)-1);
    Cx_Venus = Cx_Venus.append(Cx);
    Cy_Venus = Cy_Venus.append(Cy);
    Cz_Venus = Cz_Venus.append(Cz);
    if (0<=dt && dt<=16){
        j=0;
        Mjd0 = t1;
    }else if(16<dt && dt<=32){
        j=1;
        Mjd0 = t1+16*j;
    }
    r_Venus = 1e3*Cheb3D(Mjd_TDB, 10, Mjd0, Mjd0+16, Cx_Venus.getFilaByIndex(1,10*j+1,10*j+10),
                         Cy_Venus.getFilaByIndex(1,10*j+1,10*j+10), Cz_Venus.getFilaByIndex(1,10*j+1,10*j+10)).trans();

    double vtemp6[] = {309,320,331,342};
    temp = Matrix(1,4,vtemp6,4);
    Matrix Cx_Mars = PCtemp.getFilaByIndex(1,temp(1,1),temp(1,2)-1);
    Matrix Cy_Mars = PCtemp.getFilaByIndex(1,temp(1,2),temp(1,3)-1);
    Matrix Cz_Mars = PCtemp.getFilaByIndex(1,temp(1,3),temp(1,4)-1);
    j=0;
    Mjd0 = t1;
    r_Mars = 1e3*Cheb3D(Mjd_TDB, 11, Mjd0, Mjd0+32, Cx_Mars.getFilaByIndex(1,11*j+1,11*j+11),
                        Cy_Mars.getFilaByIndex(1,11*j+1,11*j+11), Cz_Mars.getFilaByIndex(1,11*j+1,11*j+11)).trans();

    double vtemp7[] = {342,350,358,366};
    temp = Matrix(1,4,vtemp7,4);
    Matrix Cx_Jupiter = PCtemp.getFilaByIndex(1,temp(1,1),temp(1,2)-1);
    Matrix Cy_Jupiter = PCtemp.getFilaByIndex(1,temp(1,2),temp(1,3)-1);
    Matrix Cz_Jupiter = PCtemp.getFilaByIndex(1,temp(1,3),temp(1,4)-1);
    j=0;
    Mjd0 = t1;
    r_Jupiter = 1e3*Cheb3D(Mjd_TDB, 8, Mjd0, Mjd0+32, Cx_Jupiter.getFilaByIndex(1,8*j+1,8*j+8),
                           Cy_Jupiter.getFilaByIndex(1,8*j+1,8*j+8), Cz_Jupiter.getFilaByIndex(1,8*j+1,8*j+8)).trans();

    double vtemp8[] = {366,373,380,387};
    temp = Matrix(1,4,vtemp8,4);
    Matrix Cx_Saturn = PCtemp.getFilaByIndex(1,temp(1,1),temp(1,2)-1);
    Matrix Cy_Saturn = PCtemp.getFilaByIndex(1,temp(1,2),temp(1,3)-1);
    Matrix Cz_Saturn = PCtemp.getFilaByIndex(1,temp(1,3),temp(1,4)-1);
    j=0;
    Mjd0 = t1;
    r_Saturn = 1e3*Cheb3D(Mjd_TDB, 7, Mjd0, Mjd0+32, Cx_Saturn.getFilaByIndex(1,7*j+1,7*j+7),
                          Cy_Saturn.getFilaByIndex(1,7*j+1,7*j+7), Cz_Saturn.getFilaByIndex(1,7*j+1,7*j+7)).trans();

    double vtemp9[] = {387,393,399,405};
    temp = Matrix(1,4,vtemp9,4);
    Matrix Cx_Uranus = PCtemp.getFilaByIndex(1,temp(1,1),temp(1,2)-1);
    Matrix Cy_Uranus = PCtemp.getFilaByIndex(1,temp(1,2),temp(1,3)-1);
    Matrix Cz_Uranus = PCtemp.getFilaByIndex(1,temp(1,3),temp(1,4)-1);
    j=0;
    Mjd0 = t1;
    r_Uranus = 1e3*Cheb3D(Mjd_TDB, 6, Mjd0, Mjd0+32, Cx_Uranus.getFilaByIndex(1,6*j+1,6*j+6),
                          Cy_Uranus.getFilaByIndex(1,6*j+1,6*j+6), Cz_Uranus.getFilaByIndex(1,6*j+1,6*j+6)).trans();

    double vtemp10[] = {405,411,417,423};
    temp = Matrix(1,4,vtemp10,4);
    Matrix Cx_Neptune = PCtemp.getFilaByIndex(1,temp(1,1),temp(1,2)-1);
    Matrix Cy_Neptune = PCtemp.getFilaByIndex(1,temp(1,2),temp(1,3)-1);
    Matrix Cz_Neptune = PCtemp.getFilaByIndex(1,temp(1,3),temp(1,4)-1);
    j=0;
    Mjd0 = t1;
    r_Neptune = 1e3*Cheb3D(Mjd_TDB, 6, Mjd0, Mjd0+32, Cx_Neptune.getFilaByIndex(1,6*j+1,6*j+6),
                           Cy_Neptune.getFilaByIndex(1,6*j+1,6*j+6), Cz_Neptune.getFilaByIndex(1,6*j+1,6*j+6)).trans();

    double vtemp11[] = {423,429,435,441};
    temp = Matrix(1,4,vtemp11,4);
    Matrix Cx_Pluto = PCtemp.getFilaByIndex(1,temp(1,1),temp(1,2)-1);
    Matrix Cy_Pluto = PCtemp.getFilaByIndex(1,temp(1,2),temp(1,3)-1);
    Matrix Cz_Pluto = PCtemp.getFilaByIndex(1,temp(1,3),temp(1,4)-1);
    j=0;
    Mjd0 = t1;
    r_Pluto = 1e3*Cheb3D(Mjd_TDB, 6, Mjd0, Mjd0+32, Cx_Pluto.getFilaByIndex(1,6*j+1,6*j+6),
                         Cy_Pluto.getFilaByIndex(1,6*j+1,6*j+6), Cz_Pluto.getFilaByIndex(1,6*j+1,6*j+6)).trans();

    double vtemp12[] = {819,829,839};
    temp = Matrix(1,3,vtemp12,3);
    Matrix Cx_Nutations = PCtemp.getFilaByIndex(1,temp(1,1),temp(1,2)-1);
    Matrix Cy_Nutations = PCtemp.getFilaByIndex(1,temp(1,2),temp(1,3)-1);
    for (i = 1; i <= 3; i++) {
        temp = temp+20;
        Cx = PCtemp.getFilaByIndex(1,temp(1,1),temp(1,2)-1);
        Cy = PCtemp.getFilaByIndex(1,temp(1,2),temp(1,3)-1);
        Cx_Nutations = Cx_Nutations.append(Cx);
        Cy_Nutations = Cy_Nutations.append(Cy);
    }
    if (0<=dt && dt<=8){
        j=0;
        Mjd0 = t1;
    }else if(8<dt && dt<=16){
        j=1;
        Mjd0 = t1+8*j;
    }else if (16<dt && dt<=24){
        j=2;
        Mjd0 = t1+8*j;
    }else if(24<dt && dt<=32){
        j=3;
        Mjd0 = t1+8*j;
    }
    Matrix Nutations = Cheb3D(Mjd_TDB, 10, Mjd0, Mjd0+8, Cx_Nutations.getFilaByIndex(1,10*j+1,10*j+10),
                       Cy_Nutations.getFilaByIndex(1,10*j+1,10*j+10),Matrix(10,1)).trans();

    double vtemp13[] = {899,909,919,929};
    temp = Matrix(1,4,vtemp13,4);
    Matrix Cx_Librations = PCtemp.getFilaByIndex(1,temp(1,1),temp(1,2)-1);
    Matrix Cy_Librations = PCtemp.getFilaByIndex(1,temp(1,2),temp(1,3)-1);
    Matrix Cz_Librations = PCtemp.getFilaByIndex(1,temp(1,3),temp(1,4)-1);
    for(i = 1; i <= 3; i++) {
        temp = temp+30;
        Cx = PCtemp.getFilaByIndex(1,temp(1,1),temp(1,2)-1);
        Cy = PCtemp.getFilaByIndex(1,temp(1,2),temp(1,3)-1);
        Cz = PCtemp.getFilaByIndex(1,temp(1,3),temp(1,4)-1);
        Cx_Librations = Cx_Librations.append(Cx);
        Cy_Librations = Cy_Librations.append(Cy);
        Cz_Librations = Cz_Librations.append(Cz);
    }
    if (0<=dt && dt<=8){
        j=0;
        Mjd0 = t1;
    }else if(8<dt && dt<=16){
        j=1;
        Mjd0 = t1+8*j;
    }else if (16<dt && dt<=24){
        j=2;
        Mjd0 = t1+8*j;
    }else if(24<dt && dt<=32){
        j=3;
        Mjd0 = t1+8*j;
    }
    Matrix Librations = Cheb3D(Mjd_TDB, 10, Mjd0, Mjd0+8, Cx_Librations.getFilaByIndex(1,10*j+1,10*j+10),
                        Cy_Librations.getFilaByIndex(1,10*j+1,10*j+10), Cz_Librations.getFilaByIndex(1,10*j+1,10*j+10)).trans();
    EMRAT = 81.30056907419062; // DE430
    EMRAT1 = 1/(1+EMRAT);
    r_Earth = r_Earth-EMRAT1*r_Moon;
    r_Mercury = (-r_Earth)+r_Mercury;
    r_Venus = (-r_Earth)+r_Venus;
    r_Mars = (-r_Earth)+r_Mars;
    r_Jupiter = (-r_Earth)+r_Jupiter;
    r_Saturn = (-r_Earth)+r_Saturn;
    r_Uranus = (-r_Earth)+r_Uranus;
    r_Neptune = (-r_Earth)+r_Neptune;
    r_Pluto = (-r_Earth)+r_Pluto;
    r_Sun = (-r_Earth)+r_Sun;
}

