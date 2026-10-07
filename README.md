# EKF: Orbit Determination with Extended Kalman Filter

## Overview

This project is a C++ translation of a MATLAB algorithm for **Initial Orbit Determination (IOD) using Gauss’s method and the Extended Kalman Filter (EKF)**. It processes optical observations of a satellite (GEOS-3) to estimate its orbital state, combining a classical initial orbit determination step with a sequential EKF refinement.

The original algorithm is based on the MATLAB implementation by Meysam Mahooti, which uses three optical sightings and Gauss’s method to obtain a preliminary state vector, then optimizes it with an Extended Kalman Filter. This repository ports that workflow to C++ for improved performance and integration into C++-based astrodynamics pipelines.

## Purpose

The goal of this project is to provide a **fast, standalone C++ implementation** of a well-known orbit determination algorithm. It is intended for:

- Students and researchers learning astrodynamics and Kalman filtering.
- Engineers who need a lightweight EKF-based orbit estimator without MATLAB dependencies.
- Anyone interested in translating legacy MATLAB astrodynamics code into modern C++.

The project focuses on the **GEOS-3 satellite** as a test case, using real observation data and standard gravity/ Earth-orientation models.

## Features

- **Gauss’s method** for initial orbit determination from three optical sightings.
- **Extended Kalman Filter** for sequential state estimation.
- **Custom matrix algebra** implemented in C++ (no external linear algebra libraries).
- **JPL DE430 ephemerides** and **IERS Earth-orientation corrections** support.
- **Full astrodynamics modeling**: harmonic gravity field, point-mass perturbations, and more.
- **Test suite** (`EKF_Test.cpp`) to verify the implementation.

## Project Structure

```
EKF/
├── data/                 # Input data files
│   ├── GEOS3.txt         # Observations (date, time, azimuth, elevation, distance)
│   ├── GGM03S.txt        # Gravity field coefficients
│   ├── M_tab.txt         # JPL DE430 Chebyshev coefficients
│   ├── egm.txt           # Earth gravity model
│   └── eop19620101.txt   # Earth orientation parameters
├── include/              # Header files for all modules
├── src/                  # C++ source files
│   ├── Accel.cpp
│   ├── AccelHarmonic.cpp
│   ├── ...
│   └── anglesdr.cpp
├── EKF_GEOS3.cpp         # Main application
├── EKF_Test.cpp          # Unit / integration tests
├── compile               # Build script for main program
└── compileTests          # Build script for tests
```

## Prerequisites

- **C++ compiler** (g++ recommended, C++11 or later).
- No external libraries are required beyond the standard C++ library and `libm` (math library).

## Building and Running

### Compile the main program

Use the provided `compile` script:

```bash
./compile
```

This runs:

```bash
g++ EKF_GEOS3.cpp ./src/*.cpp -o yo
```

and then executes the resulting binary:

```bash
./yo
```

> **Note:** The main program expects the `data/` directory to be present in the working directory, as it loads the observation and model files from there.

### Compile and run the tests

Use the `compileTests` script:

```bash
./compileTests
```

Equivalent command:

```bash
g++ EKF_Test.cpp ./src/*.cpp -o yo
./yo
```

## Data Files

All input data is stored in the `data/` folder:

| File | Description |
|------|-------------|
| `GEOS3.txt` | Optical observations of GEOS-3 (year, month, day, hour, minute, second, azimuth, elevation, distance) |
| `GGM03S.txt` | Gravity field coefficients (Cnm, Snm) up to degree/order 180 |
| `M_tab.txt` | Chebyshev coefficients for JPL DE430 planetary ephemerides |
| `egm.txt` | Earth gravity model coefficients |
| `eop19620101.txt` | Earth Orientation Parameters (EOP) from 1962-01-01 onwards |

These files are read at runtime by the main program to initialize the force model and process observations.

## Algorithm Summary

1. **Initial Orbit Determination (Gauss’s Method)** – The first three optical sightings are used to compute a preliminary state vector (position and velocity) at the epoch of the first observation.
2. **Extended Kalman Filter** – The preliminary state is propagated through the full force model and updated sequentially with each subsequent observation, yielding an optimized orbit estimate.
3. **Force Model** – Includes Earth’s harmonic gravity field (up to degree 180), third-body perturbations (Sun, Moon, planets via DE430), and IERS-compliant Earth orientation.

## References

- Mahooti, M. (2024). *Initial Orbit Determination using Least Squares and Extended Kalman Filter*. MATLAB Central File Exchange.
- Vallado, D. A. *Fundamentals of Astrodynamics and Applications*, 4th Edition.
- JPL DE430 Planetary Ephemerides.
- IERS Conventions (2010).

## License

This project is a translation of the original MATLAB code. Please refer to the original MATLAB Central submission for licensing terms. The C++ translation is provided for educational and research purposes.
