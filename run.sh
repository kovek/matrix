ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
brew install cmake
brew install glew
git clone https://github.com/kovek/matrix
cd matrix
curl -L http://downloads.sourceforge.net/project/ogl-math/glm-0.9.5.4/glm-0.9.5.4.zip\?r\=http%3A%2F%2Fglm.g-truc.net%2F0.9.5%2Findex.html\&ts\=1415828217\&use_mirror\=iweb > glm.tar.gz
tar -zxvf glm.tar.gz
curl -L http://downloads.sourceforge.net/project/glfw/glfw/3.0.4/glfw-3.0.4.zip\?r\=http%3A%2F%2Fwww.glfw.org%2Fdownload.html\&ts\=1415893559\&use_mirror\=softlayer-dal > `pwd`/glfw.tar.gz
tar -zxvf glfw.tar.gz
cd glfw-3.0.4
mkdir build
cd build
cmake -DBUILD_SHARED_LIBS=NO ..
make
cd ../..
make matrix
