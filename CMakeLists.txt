cmake_minimum_required(VERSION 3.23)
project(Proyecto)

set(CMAKE_CXX_STANDARD 14)

include_directories(include)

add_executable(Proyecto
        include/Matrix.h
        src/Matrix.cpp
        EKF_Test.cpp include/R_x.h src/R_x.cpp include/R_y.h src/R_y.cpp include/R_z.h src/R_z.cpp include/SAT_Const.h include/sign_.h src/sign_.cpp include/Globals.h include/unit.h src/unit.cpp
        include/timediff.h
        src/timediff.cpp
        include/Position.h
        src/Position.cpp
        include/Mjday_TDB.h
        src/Mjday_TDB.cpp
        include/Mjday.h
        src/Mjday.cpp
        include/MeanObliquity.h
        src/MeanObliquity.cpp
        include/Legendre.h
        src/Legendre.cpp
        include/Geodetic.h
        src/Geodetic.cpp
        include/Frac.h
        src/Frac.cpp
        include/EccAnom.h
        src/EccAnom.cpp
        include/NutAngles.h
        src/NutAngles.cpp
        include/AccelPointMass.h
        src/AccelPointMass.cpp
        include/AzElPa.h
        src/AzElPa.cpp
        include/Cheb3D.h
        src/Cheb3D.cpp
        include/IERS.h
        src/IERS.cpp
        include/gmst.h
        src/gmst.cpp
        include/gast.h
        src/gast.cpp
        include/EqnEquinox.h
        src/EqnEquinox.cpp
        EKF_GEOS3.cpp
        include/PrecMatrix.h
        src/PrecMatrix.cpp
        include/PoleMatrix.h
        src/PoleMatrix.cpp
        include/GHAMatrix.h
        src/GHAMatrix.cpp include/AccelHarmonic.h
        src/AccelHarmonic.cpp
        include/JPL_Eph_DE430.h
        src/JPL_Eph_DE430.cpp
        include/Accel.h
        src/Accel.cpp
        include/NutMatrix.h
        src/NutMatrix.cpp
        include/LTC.h
        src/LTC.cpp
        include/TimeUpdate.h
        src/TimeUpdate.cpp
        include/DEInteg.h
        src/DEInteg.cpp)
