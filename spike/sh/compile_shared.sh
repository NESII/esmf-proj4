#!/usr/bin/env bash

g++ -fPIC -shared transform.cpp -o /tmp/libtransform.so -I `pwd` -I /home/benkoziol/anaconda/envs/esmf-proj4/include/ \
    -L /home/benkoziol/anaconda/envs/esmf-proj4/lib/ -l proj
g++ test_transform.cpp -o /tmp/run_transform_tests -I `pwd` -I /home/benkoziol/anaconda/envs/esmf-proj4/include/ \
    -L /tmp -L /home/benkoziol/anaconda/envs/esmf-proj4/lib/ -l transform

export LD_LIBRARY_PATH=/tmp:/home/benkoziol/anaconda/envs/esmf-proj4/lib/
/tmp/run_transform_tests