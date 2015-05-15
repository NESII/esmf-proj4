#!/usr/bin/env bash

cd /tmp/build
sudo rm -r *
rm -r /tmp/install_transform/*
cmake ~/l/project/esmf-proj4/ -DCMAKE_INSTALL_PREFIX=/tmp/install_transform
make
make install
#export LD_LIBRARY_PATH=/home/benkoziol/anaconda/envs/esmf-proj4/lib
ldd /tmp/install_transform/lib/libtransform.so