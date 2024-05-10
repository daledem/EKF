#include <assert.h>
#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <math.h>

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

    sol.append(m1).append(m2).print();
    Matrix prueba = sol.append(m1).append(m2);
    prueba = (prueba.append(sol));
    prueba.print();
    sol.getFilaByIndex(1).print();

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
    double v[] = {1,2,3};
    Matrix s(1,3,v,3);

    AzElPa(Az,El,dAds,dEds,s);

    double v1[] = {0.400000000000000,-0.200000000000000,0};
    Matrix expected1(1,3,v1,3);
    double v2[] = {-0.095831484749991,-0.191662969499982,0.159719141249985};
    Matrix expected2(1,3,v2,3);

    _assert(fabs(Az - 0.463647609000806) < TOL_);
    _assert(fabs(El - 0.930274014115472) < TOL_);
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

    _assert(r_Mercury.equals(Merc,0.1));
    _assert(r_Venus.equals(Venus,0.1));
    _assert(r_Earth.equals(Earth,0.1));
    _assert(r_Mars.equals(Mars,0.1));
    _assert(r_Jupiter.equals(Jupi,0.1));
    _assert(r_Saturn.equals(Satu,0.1));
    _assert(r_Uranus.equals(Uran,0.1));
    _assert(r_Neptune.equals(Nept,0.1));
    _assert(r_Pluto.equals(Pluto,0.1));
    _assert(r_Moon.equals(Moon,0.1));
    _assert(r_Sun.equals(Sun,0.1));

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
 
    return 0;
}


int main() {

    // Cargamos PC con los datos de DE430Coeff
    FILE *fid = fopen("../data/M_tab.txt","r");

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
    fid = fopen("../data/GGM03S.txt","r");
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
    AuxParam.Mjd_UTC = 0;
    AuxParam.n      = 20;
    AuxParam.m      = 20;
    AuxParam.sun     = 1;
    AuxParam.moon    = 1;
    AuxParam.planets = 1;

    // Cargamos eopdata con los datos de eop19620101
    fid = fopen("../data/eop19620101.txt","r");

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


