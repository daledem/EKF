// #include "./include/Matrix.h"
// #include "./include/R_x.h"
// #include "./include/R_y.h"
// #include "./include/R_z.h"
// #include "./include/sign_.h"
// #include "./include/Globals.h"
// #include "./include/unit.h"
// #include "./include/timediff.h"
// #include "./include/Position.h"
// #include "./include/Mjday_TDB.h"
// #include "./include/Mjday.h"
// #include "./include/MeanObliquity.h"
// #include "./include/Legendre.h"
// #include "./include/Geodetic.h"
// #include "./include/Frac.h"
// #include "./include/EccAnom.h"
// #include "./include/NutAngles.h"
// #include "./include/AccelPointMass.h"
// #include "./include/AzElPa.h"
// #include "./include/Cheb3D.h"
// #include "./include/IERS.h"
// #include "./include/gmst.h"
// #include "./include/EqnEquinox.h"
// #include "./include/gast.h"
// #include "./include/PrecMatrix.h"
// #include "./include/PoleMatrix.h"
// #include "./include/GHAMatrix.h"
//
// #include "./include/SAT_Const.h"
//
// Matrix Snm(181,181);
// Matrix Cnm(181,181);
// Matrix PC(2285,1020);
// Matrix eopdata(13,21413);
// struct {
//     double Mjd_UTC;
//     int n;
//     int m;
//     int sun;
//     int moon;
//     int planets;
//     double Mjd_TT;
// }AuxParam;
//
// int main() {
//
//     // Cargamos PC con los datos de DE430Coeff
//     FILE *fid = fopen("../data/M_tab.txt","r");
//
//     if(fid == NULL){
//         printf("Fail to open DE430Coeff.txt");
//         exit(EXIT_FAILURE);
//     }
//
//     for(int i = 1; i <= 2285; i++){
//         for(int j = 1; j <= 1020; j++){
//             fscanf(fid,"%lf",&(PC(i,j)));
//
//         }
//     }
//
//     fclose(fid);
//
//
//     // Cargamos Snm y Cnm con los datos de GGM03S
//     fid = fopen("../data/GGM03S.txt","r");
//     int aux1,aux2;
//     double aux3,aux4;
//     if(fid == NULL){
//         printf("Fail to open GGM03S.txt");
//         exit(EXIT_FAILURE);
//     }
//
//     for(int i = 0; i <= 180; i++){
//         for(int j = 0; j <= i; j++){
//             fscanf(fid,"%d %d %lf %lf %lf %lf",&aux1,
//                    &aux2,&(Cnm(i+1,j+1)),&(Snm(i+1,j+1)),&aux3,
//                    &aux4);
//         }
//     }
//
//     fclose(fid);
//
//     // Model parameters
//     AuxParam.Mjd_UTC = 0;
//     AuxParam.n = 0;
//     AuxParam.m = 0;
//
//     // Cargamos eopdata con los datos de eop19620101
//     fid = fopen("../data/eop19620101.txt","r");
//
//     if(fid == NULL){
//         printf("Fail to open eop19620101.txt");
//         exit(EXIT_FAILURE);
//     }
//
//     for(int i = 1; i<=21413; i++){
//         fscanf(fid,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",&(eopdata(1,i)),
//                &(eopdata(2,i)),&(eopdata(3,i)),&(eopdata(4,i)),&(eopdata(5,i)),
//                &(eopdata(6,i)),&(eopdata(7,i)),&(eopdata(8,i)),&(eopdata(9,i)),
//                &(eopdata(10,i)),&(eopdata(11,i)),&(eopdata(12,i)),&(eopdata(13,i)));
//     }
//
//     fclose(fid);
//
//
//
//     return 0;
// }
