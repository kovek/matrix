#!/bin/bash

# Install MPFR
curl -L http://www.mpfr.org/mpfr-current/mpfr-3.1.2.tar.gz > `pwd`/mpfr-3.1.2.tar.gz
tar -xvf mpfr-3.1.2.tar.gz
mv mpfr-3.1.2 mpfr

# Cleanup
rm ./mpfr-3.1.2.tar.gz

cd mpfr
./configure
make

# We now have src/.libs/libmpfr.a
cp src/.libs/libmpfr.a ../libs/libmpfr.a

# Get out of mpfr directory and into external
cd ..

