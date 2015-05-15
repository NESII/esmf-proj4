#!/usr/bin/env bash

export SRC=/home/ubuntu/src
export INSTALL_PREFIX=/home/ubuntu/sandbox
export CPU_COUNT=`nproc`
export CMAKE_INCLUDE_PATH=${INSTALL_PREFIX}/include
export CMAKE_LIBRARY_PATH=${INSTALL_PREFIX}/lib
export LD_LIBRARY_PATH=${INSTALL_PREFIX}/lib
#export DYLD_LIBRARY_PATH=${INSTALL_PREFIX}/lib

# install PROJ.4
cd ${SRC}
wget http://download.osgeo.org/proj/proj-4.9.1.tar.gz
tar -xzvf proj-4.9.1.tar.gz
cd proj-4.9.1/
wget http://download.osgeo.org/proj/proj-datumgrid-1.5.zip
unzip proj-datumgrid-1.5.zip -d ./nad/
./configure --prefix=${INSTALL_PREFIX}
make -j ${CPU_COUNT}
make check | tee check_proj.out
make install

# install CMake
cd ${SRC}
wget http://www.cmake.org/files/v3.2/cmake-3.2.2.tar.gz
tar -xzvf cmake-3.2.2.tar.gz
cd cmake-3.2.2
./configure --prefix=${INSTALL_PREFIX}
make -j ${CPU_COUNT}
make install

# install transform code
cd ${SRC}
git clone https://github.com/NESII/esmf-proj4.git
cd esmf-proj4
git checkout dev
mkdir build
cd build
${INSTALL_PREFIX}/bin/cmake .. -DCMAKE_INSTALL_PREFIX=${INSTALL_PREFIX}
make
make test | tee test_transform.out
make install
