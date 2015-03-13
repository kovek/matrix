#!/bin/bash

curl -L http://sourceforge.net/projects/boost/files/boost/1.57.0/boost_1_57_0.zip/download > boost_1_57_0.zip

unzip boost_1_57_0.zip

mv boost_1_57_0 boost
cd boost
sh ./bootstraph.sh
mkdir result
./b2 install --prefix=`pwd`/result
