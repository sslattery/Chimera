#!/bin/bash
##---------------------------------------------------------------------------##
## CONFIGURE Chimera ON BEAKER WITH CMAKE
##---------------------------------------------------------------------------##

EXTRA_ARGS=$@

rm -rf CMakeCache.txt

##---------------------------------------------------------------------------##

cmake \
    -D CMAKE_INSTALL_PREFIX:PATH=$PWD \
    -D CMAKE_BUILD_TYPE:STRING=DEBUG \
    -D CMAKE_VERBOSE_MAKEFILE:BOOL=ON \
    -D TPL_ENABLE_MPI:BOOL=ON \
    -D MPI_BASE_DIR:PATH=/local.hd/cnergg/stuart/builds/openmpi-1.4.6 \
    -D BLAS_LIBRARY_DIRS:PATH=/local.hd/cnergg/stuart/software/lapack-3.4.1 \
    -D BLAS_LIBRARY_NAMES:STRING="blas" \
    -D LAPACK_LIBRARY_DIRS:PATH=/local.hd/cnergg/stuart/software/lapack-3.4.1 \
    -D LAPACK_LIBRARY_NAMES:STRING="lapack" \
    -D TPL_ENABLE_Boost:BOOL=ON \
    -D Boost_INCLUDE_DIRS:PATH=/local.hd/cnergg/stuart/builds/boost_1_51_0/include \
    -D TPL_ENABLE_BoostLib:BOOL=ON \
    -D BoostLib_INCLUDE_DIRS:PATH=/local.hd/cnergg/stuart/builds/boost_1_51_0/include \
    -D BoostLib_LIBRARY_DIRS:PATH=/local.hd/cnergg/stuart/builds/boost_1_51_0/lib \
    -D Netcdf_LIBRARY_DIRS:PATH=/local.hd/cnergg/stuart/builds/netcdf-4.1.3/lib \
    -D Netcdf_INCLUDE_DIRS:PATH=/local.hd/cnergg/stuart/builds/netcdf-4.1.3/include \
    -D Trilinos_EXTRA_REPOSITORIES="Chimera" \
    -D Trilinos_ENABLE_EXPLICIT_INSTANTIATION:BOOL=ON \
    -D Trilinos_ENABLE_STK:BOOL=ON \
    -D STK_ENABLE_SEACASIoss:BOOL=ON \
    -D Panzer_ENABLE_Stratimikos:BOOL=ON \
    -D Panzer_ENABLE_Piro:BOOL=ON \
    -D Panzer_ENABLE_NOX:BOOL=ON \
    -D Panzer_ENABLE_Rythmos:BOOL=ON \
    -D Panzer_ENABLE_STK:BOOL=ON \
    -D Panzer_ENABLE_SEACASIoss:BOOL=ON \
    -D Panzer_ENABLE_SEACASExodus:BOOL=ON \
    -D Trilinos_ENABLE_Chimera:BOOL=ON \
    $EXTRA_ARGS \
    /local.hd/cnergg/stuart/software/Trilinos
