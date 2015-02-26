#!/bin/bash
set -e

if [ "$(expr substr $(uname -s) 1 5)" == "Linux" ] || [ "$(uname)" == "Darwin" ]; then
	# Install glfw
	curl -L http://sourceforge.net/projects/plplot/files/plplot/5.10.0%20Source/plplot-5.10.0.tar.gz > ./plplot.tar.gz
	tar -xvf plplot.tar.gz

    # Cleanup
    rm plplot.tar.gz

    mv plplot-5.10.0 plplot
    cd plplot
    mkdir build
    cd build
    cmake -DBUILD_SHARED_LIBS=NO ..
    make

    # We now have glfw/build/src/libglfw3.a
    #cp src/libglfw3.a ../../libs/libglfw3.a
    make install

    # Get out of build and into the root of the glfw dir
    cd ..

    # ln -s include plplot
    # ln -s bindings/c++/plstream.h plplot/plstream.h

    # Get out of glfw dir and into external
    cd ..
fi
