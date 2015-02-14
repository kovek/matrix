#!/bin/bash

# Tools used during installation:
# 7z, git, cmake, make, tar, curl.

if [ "$(uname)" == "Darwin" ]; then
	# For mac osx
	ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
	brew install cmake
	brew install glew # Or should we link the archive instead?
	brew install p7zip

elif [ "$(expr substr $(uname -s) 1 5)" == "Linux" ]; then
	# Linux platform
	sudo apt-get install p7zip-full

elif [ "$(expr substr $(uname -s) 1 10)" == "MINGW32_NT" ]; then
	# For windows
	choco install cmake
	choco install git
	choco install 7zip
fi

# If we do not have the repository yet
git clone https://github.com/kovek/matrix
cd matrix

# Create the directory for all the external libraries
mkdir external
cd external

bash ../scripts/get_glfw.sh
bash ../scripts/get_mpfr.sh
bash ../scripts/get_glew.sh
bash ../scripts/get_glm.sh

# Get out of ./external/
cd .. # We go to root of project

# Compile the project
make matrix
