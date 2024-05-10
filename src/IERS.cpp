#include "../include/IERS.h"

void IERS(double &x_pole, double &y_pole, double &UT1_UTC, double &LOD, double &dpsi, double &deps, double &dx_pole, double &dy_pole, double &TAI_UTC,const Matrix &eop, double Mjd_UTC, char interp) {
    double mjd,mfme,fixf;
    int i;

    if (interp =='l') {
        // linear interpolation
        mjd = (floor(Mjd_UTC));
        i = Matrix::find(mjd,eop,4);

        Matrix preeop = eop.getColumnaByIndex(i);
        Matrix nexteop = eop.getColumnaByIndex(i+1);
        mfme = 1440*(Mjd_UTC-floor(Mjd_UTC));
        fixf = mfme/1440;
        // Setting of IERS Earth rotation parameters
        // (UT1-UTC [s], TAI-UTC [s], x ["], y ["])
        x_pole  = preeop(5,1)+(nexteop(5,1)-preeop(5,1))*fixf;
        y_pole  = preeop(6,1)+(nexteop(6,1)-preeop(6,1))*fixf;
        UT1_UTC = preeop(7,1)+(nexteop(7,1)-preeop(7,1))*fixf;
        LOD     = preeop(8,1)+(nexteop(8,1)-preeop(8,1))*fixf;
        dpsi    = preeop(9,1)+(nexteop(9,1)-preeop(9,1))*fixf;
        deps    = preeop(10,1)+(nexteop(10,1)-preeop(10,1))*fixf;
        dx_pole = preeop(11,1)+(nexteop(11,1)-preeop(11,1))*fixf;
        dy_pole = preeop(12,1)+(nexteop(12,1)-preeop(12,1))*fixf;
        TAI_UTC = preeop(13,1);

        x_pole  = x_pole/Const::Arcs;  // Pole coordinate [rad]
        y_pole  = y_pole/Const::Arcs;  // Pole coordinate [rad]
        dpsi    = dpsi/Const::Arcs;
        deps    = deps/Const::Arcs;
        dx_pole = dx_pole/Const::Arcs; // Pole coordinate [rad]
        dy_pole = dy_pole/Const::Arcs; // Pole coordinate [rad]
    }else if (interp =='n') {
        mjd = (floor(Mjd_UTC));
        i = Matrix::find(mjd,eop,4);
        Matrix auxEop = eop.getColumnaByIndex(i);
        // Setting of IERS Earth rotation parameters
        // (UT1-UTC [s], TAI-UTC [s], x ["], y ["])
        x_pole  = auxEop(5,1)/Const::Arcs;  // Pole coordinate [rad]
        y_pole  = auxEop(6,1)/Const::Arcs;  // Pole coordinate [rad]
        UT1_UTC = auxEop(7,1);             // UT1-UTC time difference [s]
        LOD     = auxEop(8,1);             // Length of day [s]
        dpsi    = auxEop(9,1)/Const::Arcs;
        deps    = auxEop(10,1)/Const::Arcs;
        dx_pole = auxEop(11,1)/Const::Arcs; // Pole coordinate [rad]
        dy_pole = auxEop(12,1)/Const::Arcs; // Pole coordinate [rad]
        TAI_UTC = auxEop(13,1);            // TAI-UTC time difference [s]
    }

}
