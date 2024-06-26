#include <assert.h>
#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <math.h>

#include "./include/TimeUpdate.h"
#include "./include/VarEqn.h"
#include "./include/Globals.h"
#include "./include/Matrix.h"
#include "./include/R_x.h"
#include "./include/R_y.h"
#include "./include/R_z.h"
#include "./include/sign_.h"
#include "./include/unit.h"
#include "./include/timediff.h"
#include "./include/Position.h"
#include "./include/Mjday_TDB.h"
#include "./include/Mjday.h"
#include "./include/MeanObliquity.h"
#include "./include/Legendre.h"
#include "./include/Geodetic.h"
#include "./include/Frac.h"
#include "./include/EccAnom.h"
#include "./include/NutAngles.h"
#include "./include/AccelPointMass.h"
#include "./include/AzElPa.h"
#include "./include/Cheb3D.h"
#include "./include/IERS.h"
#include "./include/gmst.h"
#include "./include/EqnEquinox.h"
#include "./include/gast.h"
#include "./include/PrecMatrix.h"
#include "./include/PoleMatrix.h"
#include "./include/GHAMatrix.h"
#include "./include/AccelHarmonic.h"
#include "./include/JPL_Eph_DE430.h"
#include "./include/NutMatrix.h"
#include "./include/Accel.h"
#include "./include/LTC.h"
#include "./include/DEInteg.h"
#include "./include/G_AccelHarmonic.h"
#include "./include/MeasUpdate.h"
#include "./include/elements.h"
#include "./include/angl.h"
#include "./include/doubler.h"
#include "./include/hgibbs.h"
#include "./include/gibbs.h"
#include "./include/anglesg.h"
#include "./include/anglesdr.h"

#define TOL_ 10e-14

int tests_run = 0;

Matrix Snm(181,181);
Matrix Cnm(181,181);
Matrix PC(2285,1020);
Matrix eopdata(13,21413);
struct auxParam AuxParam;

#define FAIL() printf("\nfailure in %s() line %d\n", __func__, __LINE__)
#define _assert(test) do { if (!(test)) { FAIL(); return 1; } } while(0)
#define _verify(test) do { int r=test(); tests_run++; if(r) return r; } while(0)

using namespace std;

int proMat_01()
{
    double v1[] = {1.0, 2.0, 3.0, 4.0};
    double v2[] = {1.0, 0.0, 0.0, 1.0};
    Matrix m1(2, 2, v1, 4);
    Matrix m2(2, 2, v2, 4);
    Matrix sol(2, 2);

    sol = m1 * m2;

    m1.print();
    m2.print();
    sol.print();
    sol.trans().print();

    _assert(sol(1,1) == m1(1,1) && sol(1,2) == m1(1,2) && sol(2,1) == m1(2,1) && sol(2,2) == m1(2,2));

    cout << sol.determinant() << endl;
    sol.inverse().print();
    (sol * sol.inverse()).print();

    return 0;
}

int R_x_01(){
    double alpha = 1.0;
    Matrix sol(3,3);

    sol = R_x(alpha);

    _assert(fabs(sol(1,1) - 1) < TOL_ && fabs(sol(1,2)) < TOL_ && fabs(sol(1,3)) < TOL_);
    _assert(fabs(sol(2,1)) < TOL_ && fabs(sol(2,2) - 0.540302305868140 ) < TOL_ && fabs(sol(2,3) - 0.841470984807897 ) < TOL_);
    _assert(fabs(sol(3,1)) < TOL_ && fabs(sol(3,2) + 0.841470984807897 ) < TOL_ && fabs(sol(3,3) - 0.540302305868140 ) < TOL_);

    return 0;
}

int R_y_01(){
    double alpha = 1.0;
    Matrix sol(3,3);

    sol = R_y(alpha);

    _assert(fabs(sol(1,1) - 0.540302305868140 ) < TOL_ && fabs(sol(1,2)) < TOL_ && fabs(sol(1,3) + 0.841470984807897 ) < TOL_);
    _assert(fabs(sol(2,1)) < TOL_ && fabs(sol(2,2) - 1.0 ) < TOL_ && fabs(sol(2,3)) < TOL_);
    _assert(fabs(sol(3,1) - 0.841470984807897 ) < TOL_ && fabs(sol(3,2)) < TOL_ && fabs(sol(3,3) - 0.540302305868140 ) < TOL_);

    return 0;
}

int R_z_01(){
    double alpha = 1.0;
    Matrix sol(3,3);

    sol = R_z(alpha);

    _assert(fabs(sol(1,1) - 0.540302305868140 ) < TOL_ && fabs(sol(1,2) - 0.841470984807897 ) < TOL_ && fabs(sol(1,3)) < TOL_);
    _assert(fabs(sol(2,1) + 0.841470984807897 ) < TOL_ && fabs(sol(2,2) - 0.540302305868140 ) < TOL_ && fabs(sol(2,3)) < TOL_);
    _assert(fabs(sol(3,1)) < TOL_ && fabs(sol(3,2)) < TOL_ && fabs(sol(3,3) - 1 ) < TOL_);

    return 0;
}

int sign__01(){
    _assert(fabs(sign_(5,-1) + 5) < TOL_);
    _assert(fabs(sign_(5,1) - 5) < TOL_);
    _assert(fabs(sign_(123,-1) + 123) < TOL_);
    _assert(fabs(sign_(123,1) - 123) < TOL_);

    return 0;
}

int unit_01(){
    double vec [] = {1,2,3};
    Matrix entrada = Matrix(1,3,vec,3);

    Matrix sol(1,3);

    sol = unit(entrada);

    double vecRes [] = {0.267261241912424,0.534522483824849,0.801783725737273};
    Matrix res = Matrix(1,3,vecRes,3);

    _assert(sol.equals(res,TOL_));

    return 0;
}

int unit_02(){
    double vec[] = {1,2,3,5};
    Matrix entrada = Matrix(1,4,vec,4);

    Matrix sol(1,3);
    sol = unit(entrada);

    double vecRes [] = {0.160128153805087,0.320256307610174,0.480384461415261};
    Matrix res = Matrix(1,3,vecRes,3);

    _assert(sol.equals(res,TOL_));

    return 0;
}

int timediff_01(){
    double UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC;

    timediff(UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC,5,6);

    _assert(fabs(UT1_TAI + 1) < TOL_);
    _assert(fabs(UTC_GPS - 13) < TOL_);
    _assert(fabs(UT1_GPS - 18) < TOL_);
    _assert(fabs(TT_UTC - 38.183999999999997) < TOL_);
    _assert(fabs(GPS_UTC + 13) < TOL_);

    return 0;
}

int position_01(){
    Matrix sol(1,3);
    sol = Position(1,2,3);
    double matriz[] = {-1.438078785611559*1.0e06,-2.239675009373783*1.0e06,5.776810445003163*1.0e06};
    Matrix expected(1,3,matriz,3);

    _assert(sol.equals(expected,TOL_));


    return 0;
}

int Mjday_TDB_01(){
    double sol;
    sol = Mjday_TDB(123);

    _assert(fabs(sol - 1.230000000183996e+02) < TOL_);

    return 0;
}

int Mjday_01(){
    double sol;
    sol = Mjday(1,2,3);

    _assert(fabs(sol + 678557) < TOL_);

    return 0;
}

int Mjday_02(){
    double sol;
    sol = Mjday(1,2,3,4,5,6);

    _assert(fabs(sol + 6.785568297916667e+05) < TOL_);

    return 0;
}

int MeanObliquity_01(){
    double sol;
    sol = MeanObliquity(123);

    _assert(fabs(sol - 0.409412306065671) < TOL_);

    return 0;
}

int Legendre_01(){
    Matrix pnm(2,2);
    Matrix dpnm(2,2);
    Legendre(pnm,dpnm,1,1,5);

    double matriz0 [] = {1.0,0.0,-1.66090556432769,0.49131731740833};
    Matrix expected0(2,2,matriz0,4);

    double matriz1 [] = {0.0,0.0,0.49131731740833,1.66090556432769};
    Matrix expected1(2,2,matriz1,4);

    _assert(pnm.equals(expected0,TOL_));
    _assert(dpnm.equals(expected1,TOL_));

    return 0;
}


int Geodetic_01(){
    double lon,lat,h;
    double vr[]= {1,2,3};
    Matrix r(1,3,vr,3);
    Geodetic(lon,lat,h,r);

    _assert(fabs(lon - 1.107148717794090) < TOL_);
    _assert(fabs(lat - 1.570744136243924) < TOL_);
    _assert(fabs(h + 6.356748616533795e+06) < TOL_);

    return 0;
}

int Frac_01(){
    double sol;
    sol = Frac(1.26);

    _assert(fabs(sol - 0.26) < TOL_);

    return 0;
}

int EccAnom_01(){
    double sol;
    sol = EccAnom(12,0.5);

    _assert(fabs(sol - 5.300942369055246) < TOL_);

    return 0;
}


int NutAngles_01(){
    double dpsi,deps;
    NutAngles(dpsi,deps,12);

    _assert(fabs(dpsi - 3.145583548136575e-05) < TOL_);
    _assert(fabs(deps - 3.844741370508873e-05) < TOL_);

    return 0;
}


int AccelPointMass_01(){
    double v1[] = {1,2,3};
    Matrix r(1,3,v1,3);
    double v2[] = {4,5,6};
    Matrix s(1,3,v2,3);

    Matrix sol(1,3);
    sol = AccelPointMass(r,s,7);

    double matriz [] = {0.108243193501550,0.097883141096128,0.087523088690707};
    Matrix expected(1,3,matriz,3);

    _assert(sol.equals(expected,TOL_));


    return 0;
}


int AzElPa_01(){
    double Az,El;
    Matrix dAds(1,3);
    Matrix dEds(1,3);
    double v[] = {2159055.44810214,
        1212982.41102093,
        729669.13008018};
    Matrix s(3,1,v,3);

    AzElPa(Az,El,dAds,dEds,s);

    double v1[] = {1.97784562210404e-07,-3.52047839037874e-07,0};
    Matrix expected1(1,3,v1,3);
    double v2[] = {-9.54424040207097e-08,-5.36206503841478e-08,3.71546961476309e-07};
    Matrix expected2(1,3,v2,3);

    _assert(fabs(Az - 1.05892995381513) < TOL_);
    _assert(fabs(El - 0.286534142298276) < TOL_);
    _assert(dAds.equals(expected1,TOL_));
    _assert(dEds.equals(expected2,TOL_));

    return 0;
}


int Cheb3D_01() {
    double vx[] = {1,2,3};
    Matrix Cx(1,3,vx,3);

    double vy[] = {4,5,6};
    Matrix Cy(1,3,vy,3);

    double vz[] = {7,8,9};
    Matrix Cz(1,3,vz,3);

    Matrix sol(1,3);

    sol = Cheb3D(3,3,2,4,Cx,Cy,Cz);

    double v[] = {-2,-2,-2};
    Matrix expected(1,3,v,3);

    _assert(sol.equals(expected,TOL_));

    return 0;
}

int IERS_01() {
    double x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC;
    IERS(x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC,eopdata,37666.00000000000000,'l');

    _assert(fabs(x_pole + 7.708537529641623e-08) < TOL_);
    _assert(fabs(y_pole - 1.037986091255516e-06) < TOL_);
    _assert(fabs(UT1_UTC - 0.032054700000000) < TOL_);
    _assert(fabs(LOD - 0.001669000000000) < TOL_);
    _assert(fabs(dpsi - 3.101789450370700e-07) < TOL_);
    _assert(fabs(deps - 3.049478054178981e-08) < TOL_);
    _assert(fabs(dx_pole) < TOL_);
    _assert(fabs(dy_pole) < TOL_);
    _assert(fabs(TAI_UTC - 2) < TOL_);

    return 0;
}


int IERS_02() {
    double x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC;
    IERS(x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC,eopdata,37666);

    _assert(fabs(x_pole + 7.708537529641623e-08) < TOL_);
    _assert(fabs(y_pole - 1.037986091255516e-06) < TOL_);
    _assert(fabs(UT1_UTC - 0.032054700000000) < TOL_);
    _assert(fabs(LOD - 0.001669000000000) < TOL_);
    _assert(fabs(dpsi - 3.101789450370700e-07) < TOL_);
    _assert(fabs(deps - 3.049478054178981e-08) < TOL_);
    _assert(fabs(dx_pole) < TOL_);
    _assert(fabs(dy_pole) < TOL_);
    _assert(fabs(TAI_UTC - 2) < TOL_);

    return 0;
}


int gmst_01() {
    double sol;
    sol = gmst(37666);

    _assert(fabs(sol - 1.765465139678228) < TOL_);

    return 0;
}


int EqnEquinox_01() {
    double sol;
    sol = EqnEquinox(37668);

    _assert(fabs(sol + 4.84288292664905e-05) < TOL_);

    return 0;
}


int gast_01() {
    double sol;
    sol = gast(37669);

    _assert(fabs(sol - 1.81702562326181) < TOL_);

    return 0;
}


int PrecMatrix_01() {
    Matrix sol(3,3);

    sol = PrecMatrix(37666,87965);

    double vecRes [] = {0.999436212010292,-0.0307945515365567,-0.0133773584156897,
                          0.0307945508840588,0.999525714109532,-0.000206081727946424,
                          0.0133773599177334,-0.000205984202882939,0.999910497900757};
    Matrix res = Matrix(3,3,vecRes,9);

    _assert(sol.equals(res,TOL_));

    return 0;
}


int PoleMatrix_01() {
    Matrix sol(3,3);

    sol = PoleMatrix(37666,87965);

    double vecRes [] = {-0.124794009354021,-0.391576364226166,-0.91164379348996,
                          0,0.918826561815273,-0.394661562991286,
                          0.992182672308557,-0.0492513987836072,-0.114664050549898};
    Matrix res = Matrix(3,3,vecRes,9);

    _assert(sol.equals(res,TOL_));

    return 0;
}


int GHAMatrix_01() {
    Matrix sol(3,3);

    sol = GHAMatrix(37666);

    double vecRes [] = {-0.193393597637188,0.981121254684122,0,
                          -0.981121254684122,-0.193393597637188,0,
                          0,0,1};
    Matrix res = Matrix(3,3,vecRes,9);

    _assert(sol.equals(res,TOL_));

    return 0;
}


int AccelHarmoic_01() {
    double vecR [] = {5720694.2260585,2687728.41425142,3483000.08675422};
    Matrix r = Matrix(3,1,vecR,3);

    double vecE [] = {-0.976558757940107,0.215250556888025,-0.000435947096290288,
                              -0.2152505181354,-0.976558854525697,-0.000134498699131133,
                              -0.000454678916875739,-3.7508044211961e-05,0.999999895930109};
    Matrix E = Matrix(3,3,vecE,9);

    Matrix sol(3,1);

    sol = AccelHarmonic(r,E,20,20);

    double vecRes [] = {-6.06544113186609,-2.8497772091355,-3.70232504408661};
    Matrix res = Matrix(3,1,vecRes,3);

    _assert(sol.equals(res,TOL_));

    return 0;
}


int JPL_Eph_DE430_01() {
    Matrix r_Mercury(3,1);
    Matrix r_Venus(3,1);
    Matrix r_Earth(3,1);
    Matrix r_Mars(3,1);
    Matrix r_Jupiter(3,1);
    Matrix r_Saturn(3,1);
    Matrix r_Uranus(3,1);
    Matrix r_Neptune(3,1);
    Matrix r_Pluto(3,1);
    Matrix r_Moon(3,1);
    Matrix r_Sun(3,1);
    JPL_Eph_DE430(r_Mercury,r_Venus,r_Earth,r_Mars,r_Jupiter,r_Saturn,r_Uranus,r_Neptune,r_Pluto,r_Moon,r_Sun,37666);

    double vMerc [] = {0.735116603324514*1.0e+11,
          -1.691926538640355*1.0e+11,
          -0.816370690593120*1.0e+11};
    Matrix Merc(3,1,vMerc,3);
    double vVenus [] = {0.246563930605697*1.0e+11,
           -2.313576670076116*1.0e+11,
           -1.016066087671335*1.0e+11};
    Matrix Venus(3,1,vVenus,3);
    double vEarth [] = {-0.301765548935945*1.0e+11,
            1.330462426738211*1.0e+11,
            0.576994514956490*1.0e+11};
    Matrix Earth(3,1,vEarth,3);
    double vMars [] = { 0.428014158577189*1.0e+11,
           -3.285476270689033*1.0e+11,
           -1.477395954528360*1.0e+11};
    Matrix Mars(3,1,vMars,3);
    double vJupi [] = {5.806041993621722*1.0e+11,
           -6.039191097953667*1.0e+11,
           -2.729724763892779*1.0e+11};
    Matrix Jupi(3,1,vJupi,3);
    double vSatu [] = {0.824980894643224*1.0e+12,
           -1.289737921795669*1.0e+12,
           -0.569475935851322*1.0e+12};
    Matrix Satu(3,1,vSatu,3);
    double vUran [] = {-2.307757638486849*1.0e+12,
            1.171432960693235*1.0e+12,
            0.546761533865996*1.0e+12};
    Matrix Uran(3,1,vUran,3);
    double vNept [] = {-3.343537384563364*1.0e+12,
           -2.968732390487868*1.0e+12,
           -1.134405249457255*1.0e+12};
    Matrix Nept(3,1,vNept,3);
    double vPluto [] = {-4.510941275724197*1.0e+12,
             1.012707866105105*1.0e+12,
            1.667907922624643*1.0e+12};
    Matrix Pluto(3,1,vPluto,3);
    double vMoon [] = {-2.790184113338402*1.0e+08,
           -2.602453013405682*1.0e+08,
           -0.746811863302955*1.0e+08};
    Matrix Moon(3,1,vMoon,3);
    double vSun [] = {0.296997880114303*1.0e+11,
           -1.321769585597421*1.0e+11,
           -0.573185756659810*1.0e+11};
    Matrix Sun(3,1,vSun,3);

    _assert(r_Mercury.equals(Merc,10e-3));
    _assert(r_Venus.equals(Venus,10e-3));
    _assert(r_Earth.equals(Earth,10e-3));
    _assert(r_Mars.equals(Mars,10e-3));
    _assert(r_Jupiter.equals(Jupi,10e-3));
    _assert(r_Saturn.equals(Satu,10e-3));
    _assert(r_Uranus.equals(Uran,10e-3));
    _assert(r_Neptune.equals(Nept,10e-3));
    _assert(r_Pluto.equals(Pluto,10e-3));
    _assert(r_Moon.equals(Moon,10e-3));
    _assert(r_Sun.equals(Sun,10e-3));

    return 0;
}


int NutMatrix_01() {
    Matrix sol(3,3);

    sol = NutMatrix(37666);

    double vecRes [] = {0.999999998577198,0.000048940522818,0.000021223327599,
    -0.000048941303822,0.999999998125241,0.000036800373610,
    -0.000021221526530,-0.000036801412255,0.999999999097651};
    Matrix res = Matrix(3,3,vecRes,9);

    _assert(sol.equals(res,TOL_));

    return 0;
}


int Accel_01() {
    double vY[] = {5720694.22605799,
                    2687728.41425145,
                    3483000.08675447,
                    4371.83136151856,
                   -1905.47309295878,
                   -5698.58341611591};
    Matrix Y(6,1,vY,6);

    Matrix sol(6,1);

    sol = Accel(-543.476874884521,Y);

    double vecRes [] = {4371.83136151856,
                       -1905.47309295878,
                       -5698.58341611591,
                       -6.06544204261725,
                       -2.84977703178303,
                       -3.70232534578413};
    Matrix res = Matrix(6,1,vecRes,6);

    _assert(sol.equals(res,TOL_));

    return 0;
}


int LTC_01() {
    Matrix sol(3,3);

    sol = LTC(-2.76234307910694,0.376551295459273);

    double vecRes [] = {0.370223471399199,-0.928942722252092,0,
           0.341586711932422,0.136136938528208,0.929938305587722,
          -0.863859421119156,-0.344284987681776,0.367715580035218};
    Matrix res = Matrix(3,3,vecRes,9);

    _assert(sol.equals(res,TOL_));

    return 0;
}


int TimeUpdate_01() {
    double vP[] = {
          15877.8629547137,-5671.38396107528,8540.76030401962,48.5595994176006,-13.370054730257,22.4250145441674,
           -5671.38396107528,24449.542245505,-1552.86140802114,-4.6404293131461,59.5344880226286,-27.1265452888365,
            8540.76030401961,-1552.86140802112,6013.07747951391,26.7050171224438,-4.18241475345229,15.4937252972693,
            48.5595994176007,-4.64042931314612,26.7050171224438,0.178295892456546,-0.0189353069596782,0.0668704711111732,
            -13.370054730257,59.5344880226286,-4.18241475345233,-0.0189353069596781,0.152596870759033,-0.0751888003409192,
            22.4250145441673,-27.1265452888364,15.4937252972693,0.0668704711111729,-0.0751888003409188,0.0826395235292038};
    Matrix p(6,6,vP,36);

    double vPhi[] = {
          1.0000252553551,7.08259832338682e-06,1.91609345962692e-07,5.00001498082565,1.17818628680213e-05,2.68390167250762e-07,
         7.08259815193455e-06,0.999988040046621,3.53016112726607e-08,1.17818628283062e-05,4.99995293819717,4.93631355286403e-08,
         1.91608860896806e-07,3.53015287901517e-08,0.999986704774626,2.68389762643458e-07,4.93630677857456e-08,4.99995072081276,
         1.01043851885367e-05,2.82768356912433e-06,6.4413632355971e-08,1.00002526606745,7.05571115134612e-06,1.30455626093002e-07,
          2.8276833653238e-06,-4.78603729148054e-06,1.1852833714462e-08,7.05571100209626e-06,0.999988029832341,2.39619734956897e-08,
         6.44131450915386e-08,1.18527460879032e-08,-5.31820682447183e-06,1.3045513745088e-07,2.39618836715383e-08,0.999986704276552};
    Matrix phi(6,6,vPhi,36);

    Matrix sol(6,6);

    sol = TimeUpdate(p,phi);

    double vecRes [] = {
          16368.6544574661,-5761.69527096386,8788.17676914831,49.6007868070509,-13.3911936132518,22.7142859225539,
           -5761.69527096386,25048.0229177432,-1711.18164639471,-4.72250477776863,60.1612187344361,-27.4933417263728,
            8788.17676914831,-1711.18164639469,6169.9204142798,27.1227222044556,-4.52543359225436,15.8746551064867,
            49.6007868070509,-4.72250477776866,27.1227222044556,0.179264707099232,-0.0187400566610417,0.0668822749140359,
           -13.3911936132518,60.1612187344361,-4.5254335922544,-0.0187400566610416,0.151948217745369,-0.0749710842103104,
            22.7142859225538,-27.4933417263727,15.8746551064867,0.0668822749140356,-0.07497108421031,0.0824749546958939};
    Matrix res = Matrix(6,6,vecRes,36);

    _assert(sol.equals(res,10e-10));

    return 0;
}


int DEInteg_01() {
    double vY0_apr [] = {
          6221397.62857869,
            2867713.77965738,
            3006155.98509949,
            4645.04725161806,
           -2752.21591588204,
           -7507.99940987031
    };
    Matrix Y0_apr(6,1,vY0_apr,6);

    Matrix sol(6,1);

    sol = DEInteg(Accel,0,-134.999991953373,1e-13,1e-6,6,Y0_apr);

    double vecRes [] = {
          5542555.89427451,
            3213514.83814162,
            3990892.92789074,
            5394.06894044389,
           -2365.21290574021,
            -7061.8448137347,
    };
    Matrix res = Matrix(6,1,vecRes,6);

    _assert(sol.equals(res,10e-09));

    return 0;
}


int G_AccelHarmonic_01() {
    double vecR [] = {
          7101800.90695315,
            1293997.58115302,
             10114.014948955
    };
    Matrix r = Matrix(3,1,vecR,3);

    double vecU [] = {
          -0.984320311904791,0.17638970840918,-0.000440838949610109,
          -0.176389673507182,-0.984320409906027,-0.000117142904888635,
       -0.000454589578418276,-3.75467022865179e-05,0.999999895969275
    };
    Matrix U = Matrix(3,3,vecU,9);

    Matrix sol(3,1);

    sol = G_AccelHarmonic(r,U,20,20);

    double vecRes [] = {
          2.02233500257165e-06,5.61803303433805e-07,4.39856240319614e-09,
        5.61803301435404e-07,-9.58631634517815e-07,8.05634892131479e-10,
        4.39855909334375e-09,8.0563404905587e-10,-1.06370336962723e-06
    };
    Matrix res = Matrix(3,3,vecRes,9);

    _assert(sol.equals(res,TOL_));

    return 0;
}


int VarEqn_01() {
    double vyPhi[] = {7101576.98990384,
            1295199.87127754,
            12739.2823333893,
            576.004651192995,
           -3084.62203617269,
           -6736.02594582755,
             1.0000252553551,
        7.08259815193455e-06,
        1.91608860896806e-07,
        1.01043851885367e-05,
         2.8276833653238e-06,
        6.44131450915386e-08,
        7.08259832338682e-06,
           0.999988040046621,
        3.53015287901517e-08,
        2.82768356912433e-06,
       -4.78603729148054e-06,
        1.18527460879032e-08,
        1.91609345962692e-07,
        3.53016112726607e-08,
           0.999986704774626,
         6.4413632355971e-08,
         1.1852833714462e-08,
       -5.31820682447183e-06,
            5.00001498082565,
        1.17818628283062e-05,
        2.68389762643458e-07,
            1.00002526606745,
        7.05571100209626e-06,
         1.3045513745088e-07,
        1.17818628680213e-05,
            4.99995293819717,
        4.93630677857456e-08,
        7.05571115134612e-06,
           0.999988029832341,
        2.39618836715383e-08,
        2.68390167250762e-07,
        4.93631355286403e-08,
            4.99995072081276,
        1.30455626093002e-07,
        2.39619734956897e-08,
           0.999986704276552};
    Matrix yPhi(42,1,vyPhi,42);

    Matrix sol(42,1);

    sol = VarEqn(-543.476874884521,yPhi);

    double vecRes [] = {576.004651192995,
           -3084.62203617269,
           -6736.02594582755,
           -7.53466223591457,
           -1.37422019436626,
         -0.0135523187197559,
        1.01043851885367e-05,
         2.8276833653238e-06,
        6.44131450915386e-08,
        2.02219654070929e-06,
        5.62315203995695e-07,
        5.54306317728296e-09,
        2.82768356912433e-06,
       -4.78603729148054e-06,
        1.18527460879032e-08,
        5.62315388477347e-07,
       -9.58426101908334e-07,
        1.01508481987651e-09,
         6.4413632355971e-08,
         1.1852833714462e-08,
       -5.31820682447183e-06,
        5.54345424371089e-09,
        1.01515597613439e-09,
       -1.06368579634574e-06,
            1.00002526606745,
        7.05571100209626e-06,
         1.3045513745088e-07,
        1.01107443555988e-05,
        2.81153608575356e-06,
        2.77154085572794e-08,
        7.05571115134612e-06,
           0.999988029832341,
        2.39618836715383e-08,
        2.81153631885257e-06,
       -4.79215600623252e-06,
        5.07544128312951e-09,
        1.30455626093002e-07,
        2.39619734956897e-08,
           0.999986704276552,
        2.77159004655386e-08,
        5.07553139892566e-09,
       -5.31844727806426e-06};
    Matrix res = Matrix(42,1,vecRes,42);

    _assert(sol.equals(res,TOL_));

    return 0;
}


int MeasUpdate_01() {
    double vecY [] = {
              7101597.84250255,
              1295247.06100899,
              12762.8936334637,
              576.097627376785,
              -3084.51047032311,
              -6736.01185847831
    };
    Matrix Y(6,1,vecY,6);

    double vecG [] = {
              0.484876050827302,0.042004563746839,-0.873573598478432,0,0,0
    };
    Matrix G(1,6,vecG,6);

    double vecP [] = {
              15877.8629548516,-5671.383961132,8540.76030405896,48.5595994191812,-13.3700547302772,22.4250145443688,
              -5671.383961132,24449.5422452456,-1552.86140807498,-4.64042931270271,59.5344880218809,-27.1265452890941,
              8540.76030405897,-1552.86140807498,6013.07747951118,26.7050171232874,-4.18241475350498,15.4937252973352,
              48.5595994191812,-4.6404293127027,26.7050171232874,0.178295892465964,-0.0189353069583249,0.0668704711122519,
              -13.3700547302772,59.534488021881,-4.18241475350496,-0.018935306958325,0.152596870757018,-0.0751888003412653,
              22.4250145443688,-27.1265452890941,15.4937252973352,0.066870471112252,-0.0751888003412654,0.0826395235300218
    };
    Matrix P(6,6,vecP,36);

    Matrix K(6,1);

    MeasUpdate(K,Y,2653472,2653524.97225556,92.5,G,P,6);

    double vecResK [] = {
              -4.29770209254283e-05,
              -0.0382899545454505,
              -0.122992273850069,
              2.26462926762789e-06,
              -3.43254740338148e-05,
              -0.000397229919712573
    };
    Matrix resK = Matrix(6,1,vecResK,6);

    double vecResY [] = {
              7101597.84477914,
              1295249.08931425,
              12769.408811626,
              576.097507414265,
              -3084.50865202532,
              -6735.99081631349
    };
    Matrix resY = Matrix(6,1,vecResY,6);

    double vecResP [] = {
              15877.8629371778     ,    -5671.39970736124     ,     8540.70972513922     ,     48.5596003504796,         -13.3700688461671,          22.4248511889027,
              -5671.39970736124   ,       24435.5132970153    ,     -1597.92419968857    ,     -4.63959958158567,         59.5219116084608,         -27.2720852367213,
              8540.70972513922   ,      -1597.92419968857     ,      5868.3299792234     ,     26.7076823266803 ,        -4.22281181591431,          15.0262322092083,
              48.5596003504796   ,      -4.63959958158566     ,     26.7076823266803    ,      0.17829584339217 ,      -0.0189345631362376,        0.0668790789580673,
              -13.3700688461671   ,       59.5219116084608    ,      -4.2228118159143  ,     -0.0189345631362377 ,        0.152585596485282,       -0.0753192713160279,
              22.4248511889026    ,     -27.2720852367213     ,     15.0262322092083   ,     0.0668790789580674  ,      -0.075319271316028,        0.0811296543646693
    };
    Matrix resP = Matrix(6,6,vecResP,36);

    _assert(K.equals(resK,TOL_));
    _assert(Y.equals(resY,10e-09));
    _assert(P.equals(resP,10e-09));

    return 0;
}


int elements_01() {
    double p, a, e, i, Omega, omega, M;

    double vecy[] = {6221397.62857869,4645.04725161806,
      2867713.77965738,-2752.21591588204,
      3006155.98509949,-7507.99940987031};
    Matrix y(3,2,vecy,6);


    elements(p, a, e, i, Omega, omega, M,y);

    _assert(fabs(p - 12001693.597214) < 10e-7);
    _assert(fabs(a - 18943922.6607145) < 10e-7);
    _assert(fabs(e - 0.605361104987026) < TOL_);
    _assert(fabs(i - 2.02656295535017) < TOL_);
    _assert(fabs(Omega - 3.35671076650829) < TOL_);
    _assert(fabs(omega - 2.73757289772562) < TOL_);
    _assert(fabs(M - 6.27144693341967) < TOL_);

    return 0;
}


int angl_01() {
    double theta;

    double vR1[] = {5720303.71012986,
            3152426.6965331,
           3750056.80416402};
    Matrix r1(3,1,vR1,3);

    double vR2[] = {6221397.62857869,
           2867713.77965738,
           3006155.98509949};
    Matrix r2(3,1,vR2,3);

    theta = angl(r1,r2);

    _assert(fabs(theta - 0.125269502872995) < TOL_);

    return 0;
}


int doubler_01() {
    double f1,f2,q1,magr1,magr2,a,deltae32;

    double vLos1[] = {-0.0514407086877254,
          0.838593165124504,
           0.54232403311689};
    Matrix los1(3,1,vLos1,3);

    double vLos2[] = {0.185350449302478,
          0.924321654007447,
          0.333578612738146};
    Matrix los2(3,1,vLos2,3);

    double vLos3[] = {0.489992086298769,
           0.86577353456705,
          -0.10170517296508};
    Matrix los3(3,1,vLos3,3);

    double vRsite1[] = {5854667.14002558,
           962016.478096048,
           2333503.27946785};
    Matrix rsite1(3,1,vRsite1,3);

    double vRsite2[] = {5847642.25665158,
           1003838.15985758,
           2333501.56780477};
    Matrix rsite2(3,1,vRsite2,3);

    double vRsite3[] = {5839554.60139875,
           1049867.90406392,
           2333499.52243618};
    Matrix rsite3(3,1,vRsite3,3);

    Matrix r2(3,1);
    Matrix r3(3,1);

    doubler(r2,r3,f1,f2,q1,magr1,magr2,a,deltae32,3542174.25253467,5580277.36743298,6375565.90271882,6375565.90271882,
          7232409.19024033,7230729.66707348,los1,los2,los3,rsite1,rsite2,rsite3,-97.9999914765358,108.000017702579,'y');

    double vR2[] = {6.147303745912716*1.0e+06,
    2.498215996213244*1.0e+06,
    2.872807861665693*1.0e+06};
    Matrix r2Expc(3,1,vR2,3);

    double vR3[] = {6.515289823716909*1.0e+06,
    2.243833456195559*1.0e+06,
    2.193240588493806*1.0e+06};
    Matrix r3Expc(3,1,vR3,3);

    _assert(r2.equals(r2Expc,10e-9));
    _assert(r3.equals(r3Expc,10e-9));
    _assert(fabs(f1 + 2.063075044134166e-09) < 10e-9);
    _assert(fabs(f2 - 4.628583383237128e-08) < 10e-9);
    _assert(fabs(q1 - 4.633178921858744e-08) < 10e-9);
    _assert(fabs(magr1 - 7.232409190240329e+06) < 10e-9);
    _assert(fabs(magr2 - 7.230729667073472e+06) < 10e-9);
    _assert(fabs(a - 7.458444095976775e+06) < 10e-7);
    _assert(fabs(deltae32 - 0.109188763779170) < 10e-9);

    return 0;
}

int hgibbs_01() {
    double theta,theta1,copa;
    string error;

    double vR1[] = {5720303.71012986,
            3152426.6965331,
           3750056.80416402};
    Matrix r1(3,1,vR1,3);

    double vR2[] = {6221397.62857869,
           2867713.77965738,
           3006155.98509949};
    Matrix r2(3,1,vR2,3);

    double vR3[] = {6699811.80976796,
           2569867.80763881,
           2154940.29542389};
    Matrix r3(3,1,vR3,3);

    Matrix v2(3,1);

    hgibbs(v2,theta,theta1,copa,error,r1,r2,r3,49746.1101504629,49746.1112847221,49746.1125347223);

    double vV2[] = {4796.82507080883,
          -2839.41807616618,
          -7741.59421072284};
    Matrix v2Expc(3,1,vV2,3);

    _assert(v2.equals(v2Expc,10e-12));
    _assert(fabs(theta - 0.125269502872995) < TOL_);
    _assert(fabs(theta1 - 0.136454013492469) < TOL_);
    _assert(fabs(copa - 0.00509723347775616) < TOL_);
    _assert(error == "   angl > 1ø");

    return 0;
}


int gibbs_01() {
    double theta,theta1,copa;
    string error;

    double vR1[] = {5720303.71012986,
            3152426.6965331,
           3750056.80416402};
    Matrix r1(3,1,vR1,3);

    double vR2[] = {6221397.62857869,
           2867713.77965738,
           3006155.98509949};
    Matrix r2(3,1,vR2,3);

    double vR3[] = {6699811.80976796,
           2569867.80763881,
           2154940.29542389};
    Matrix r3(3,1,vR3,3);

    Matrix v2(3,1);

    gibbs(v2,theta,theta1,copa,error,r1,r2,r3);

    double vV2[] = {4645.04725161805,
          -2752.21591588211,
          -7507.99940987023};
    Matrix v2Expc(3,1,vV2,3);

    _assert(v2.equals(v2Expc,10e-10));
    _assert(fabs(theta - 0.125269502872995) < TOL_);
    _assert(fabs(theta1 - 0.136454013492469) < TOL_);
    _assert(fabs(copa - 0.00509723347775616) < TOL_);
    _assert(error == "          ok");

    return 0;
}


int anglesg_01() {

    double vRs[] = {-5512567.84003606811166,
    -2196994.44666933314875,
    2330804.96614688728005};
    Matrix Rs(3,1,vRs,3);

    Matrix r2(3,1);
    Matrix v2(3,1);

    anglesg(r2,v2,1.05590848949330,1.36310214580757,1.97615602688759,0.28262465643395,0.45343479433887,0.58642713801159,49746.110150462948,49746.111284722108,49746.112534722313,Rs,Rs,Rs);

    double vR2[] = {6221397.62857869,
           2867713.77965738,
           3006155.98509949};
    Matrix r2Expc(3,1,vR2,3);

    double vV2[] = {4645.04725161806,
          -2752.21591588204,
          -7507.99940987031};
    Matrix v2Expc(3,1,vV2,3);

    _assert(r2.equals(r2Expc,10e-7));
    _assert(v2.equals(v2Expc,10e-7));

    return 0;
}


int anglesdr_01() {

    double vRs[] = {-5512567.84003606811166,
    -2196994.44666933314875,
    2330804.96614688728005};
    Matrix Rs(3,1,vRs,3);

    Matrix r2(3,1);
    Matrix v2(3,1);

    anglesdr(r2,v2,1.05590848949330,1.36310214580757,1.97615602688759,0.28262465643395,0.45343479433887,0.58642713801159,49746.110150462948,49746.111284722108,49746.112534722313,Rs,Rs,Rs);

    double vR2[] = {6.147303745912718*1.0e+06,
     2.498215996213244*1.0e+06,
     2.872807861665694*1.0e+06};
    Matrix r2Expc(3,1,vR2,3);

    double vV2[] = {3.764629388579156*1.0e+03,
      -2.217845391101618*1.0e+03,
      -6.141471707171266*1.0e+03};
    Matrix v2Expc(3,1,vV2,3);

    _assert(r2.equals(r2Expc,10e-7));
    _assert(v2.equals(v2Expc,10e-7));

    return 0;
}


int all_tests()
{
    _verify(proMat_01);
    _verify(R_x_01);
    _verify(R_y_01);
    _verify(R_z_01);
    _verify(sign__01);
    _verify(unit_01);
    _verify(unit_02);
    _verify(timediff_01);
    _verify(position_01);
    _verify(Mjday_TDB_01);
    _verify(Mjday_01);
    _verify(Mjday_02);
    _verify(MeanObliquity_01);
    _verify(Legendre_01);
    _verify(Geodetic_01);
    _verify(Frac_01);
    _verify(EccAnom_01);
    _verify(NutAngles_01);
    _verify(AccelPointMass_01);
    _verify(AzElPa_01);
    _verify(Cheb3D_01);
    _verify(IERS_01);
    _verify(IERS_02);
    _verify(gmst_01);
    _verify(EqnEquinox_01);
    _verify(gast_01);
    _verify(PrecMatrix_01);
    _verify(PoleMatrix_01);
    _verify(GHAMatrix_01);
    _verify(AccelHarmoic_01);
    _verify(JPL_Eph_DE430_01);
    _verify(NutMatrix_01);
    _verify(Accel_01);
    _verify(LTC_01);
    _verify(TimeUpdate_01);
    _verify(DEInteg_01);
    _verify(G_AccelHarmonic_01);
    _verify(VarEqn_01);
    _verify(MeasUpdate_01);
    _verify(elements_01);
    _verify(angl_01);
    _verify(doubler_01);
    _verify(hgibbs_01);
    _verify(gibbs_01);
    _verify(anglesg_01);
    _verify(anglesdr_01);

    return 0;
}


int main() {

    // Cargamos PC con los datos de DE430Coeff
    FILE *fid = fopen("data/M_tab.txt","r");

    if(fid == NULL){
        printf("Fail to open DE430Coeff.txt");
        exit(EXIT_FAILURE);
    }

    for(int i = 1; i <= 2285; i++){
        for(int j = 1; j <= 1020; j++){
            fscanf(fid,"%lf",&(PC(i,j)));
        }
    }

    fclose(fid);


    // Cargamos Snm y Cnm con los datos de GGM03S
    fid = fopen("data/GGM03S.txt","r");
    int aux1,aux2;
    double aux3,aux4;
    if(fid == NULL){
        printf("Fail to open GGM03S.txt");
        exit(EXIT_FAILURE);
    }

    for(int i = 0; i <= 180; i++){
        for(int j = 0; j <= i; j++){
            fscanf(fid,"%d %d %lf %lf %lf %lf",&aux1,
                &aux2,&(Cnm(i+1,j+1)),&(Snm(i+1,j+1)),&aux3,
                &aux4);
        }
    }

    fclose(fid);

    // Model parameters
    AuxParam.Mjd_UTC= 49746.1163541665;
    AuxParam.n= 20;
    AuxParam.m= 20;
    AuxParam.sun= 1;
    AuxParam.moon= 1;
    AuxParam.planets= 1;
    AuxParam.Mjd_TT= 49746.1170623147;


    // Cargamos eopdata con los datos de eop19620101
    fid = fopen("data/eop19620101.txt","r");

    if(fid == NULL){
        printf("Fail to open eop19620101.txt");
        exit(EXIT_FAILURE);
    }

    for(int i = 1; i<=21413; i++){
        fscanf(fid,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",&(eopdata(1,i)),
            &(eopdata(2,i)),&(eopdata(3,i)),&(eopdata(4,i)),&(eopdata(5,i)),
            &(eopdata(6,i)),&(eopdata(7,i)),&(eopdata(8,i)),&(eopdata(9,i)),
            &(eopdata(10,i)),&(eopdata(11,i)),&(eopdata(12,i)),&(eopdata(13,i)));
    }

    fclose(fid);

    int result = all_tests();

    if (result == 0)
        printf("PASSED\n");

    printf("Tests run: %d\n", tests_run);

    return result != 0;
}


