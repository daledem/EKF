#include "../include/Globals.h"

Matrix *Global::eopdata;

Matrix *Global::Cnm;
Matrix *Global::Snm;

Matrix *Global::PC;

void Global::GGM03S() {
    int aux1,aux2;
    double aux3,aux4;
    // read Earth orientation parameters
    Global::Cnm = new Matrix(181,181);
    Global::Snm = new Matrix(181,181);
    FILE *fid = fopen("../data/GGM03S.txt","r");

    if(fid == NULL){
        printf("Fail to open GGM03S.txt");
        exit(EXIT_FAILURE);
    }

    for(int i = 0; i <= 180; i++){
        for(int j = 0; j <= i; j++){
            fscanf(fid,"%d %d %lf %lf %lf %lf",&aux1,
                   &aux2,&((*Global::Cnm)(i+1,j+1)),&((*Global::Snm)(i+1,j+1)),&aux3,
                   &aux4);
        }
    }

    fclose(fid);
}

void Global::eop19620101(int c){
    // read Earth orientation parameters
    Global::eopdata = new Matrix(13,c);
    FILE *fid = fopen("../data/eop19620101.txt","r");

    if(fid == NULL){
        printf("Fail to open eop19620101.txt");
        exit(EXIT_FAILURE);
    }

    for(int i = 1; i<=c; i++){
        fscanf(fid,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",&((*Global::eopdata)(1,i)),
                &((*Global::eopdata)(2,i)),&((*Global::eopdata)(3,i)),&((*Global::eopdata)(4,i)),&((*Global::eopdata)(5,i)),
                &((*Global::eopdata)(6,i)),&((*Global::eopdata)(7,i)),&((*Global::eopdata)(8,i)),&((*Global::eopdata)(9,i)),
                &((*Global::eopdata)(10,i)),&((*Global::eopdata)(11,i)),&((*Global::eopdata)(12,i)),&((*Global::eopdata)(13,i)));
    }

    fclose(fid);
}

void Global::DE430Coeff() {
    Global::PC = new Matrix(2285,1020);

    FILE *fid = fopen("../data/M_tab.txt","r");

    if(fid == NULL){
        printf("Fail to open DE430Coeff.txt");
        exit(EXIT_FAILURE);
    }

    for(int i = 1; i <= 2285; i++){
        for(int j = 1; j <= 1020; j++){
            fscanf(fid,"%lf",&((*Global::PC)(i,j)));

        }
    }

    fclose(fid);
}