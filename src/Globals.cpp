#include "../include/Globals.h"

Matrix *Global::eopdata;

void Global::eop19620101(int c){
    // read Earth orientation parameters
    Global::eopdata = new Matrix(13,c);
    FILE *fid = fopen("../data/eop19620101.txt","r");

    if(fid == nullptr){
        printf("Error");
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