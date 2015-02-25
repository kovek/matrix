#!/bin/bash

# Install GMP
curl -L https://gmplib.org/download/gmp/gmp-6.0.0a.tar.lz > `pwd`/gmp-6.0.0a.tar.lz
lzip -d gmp-6.0.0a.tar.lz
tar -xvf gmp-6.0.0a.tar
mv gmp-6.0.0 gmp

# Cleanup
rm gmp-6.0.0a.tar

cd gmp
./configure
make
cp ./.libs/libgmp.a ../libs/libgmp.a
