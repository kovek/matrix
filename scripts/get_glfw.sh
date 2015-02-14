#!/bin/bash
set -e

if [ "$(expr substr $(uname -s) 1 5)" == "Linux" ] || [ "$(uname)" == "Darwin" ]; then
	# Install glfw
	curl -L http://downloads.sourceforge.net/project/glfw/glfw/3.1/glfw-3.1.tar.bz2 > ./glfw.tar.gz
	tar -xvf glfw.tar.gz

    # Cleanup
    rm glfw.tar.gz

	mv glfw-3.1 glfw
	cd glfw
	mkdir build
	cd build
	cmake -DBUILD_SHARED_LIBS=NO ..
	make
    # We now have glfw/build/src/libglfw3.a
    cp src/libglfw3.a ../../libs/libglfw3.a

	# Get out of build and into the root of the glfw dir
	cd ..
	# Get out of glfw dir and into external
    cd ..
fi
