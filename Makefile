all: build/main

run: build/main
	./build/main

build/main: src/main.o src/Scene.o src/initialize.o src/common/md5.o
	g++	-std=c++0x -o build/main -g \
	src/initialize.o src/main.o src/Scene.o src/common/md5.o \
	external/libs/libglfw3.a external/libs/libgmp.a external/libs/libmpfr.a \
   	-Iexternal/glm \
	-lmpfr -lglfw3 -lmpfr -lmpfr -lplplotcxxd -lsqlite3 -lboost_serialization \
	-framework Cocoa -framework CoreVideo -framework OpenGL -framework IOKit



src/main.o: src/main.cpp
	g++ -std=c++0x src/main.cpp -c -o src/main.o \
	-Isrc/mpreal \
	-Iexternal/plplot/include -Iexternal/plplot/build/include \
	-Isrc/common

src/Scene.o: src/Scene.cpp
	g++ -std=c++0x src/Scene.cpp -c -o src/Scene.o \
	-Isrc/mpreal \
	-Iexternal/glfw/include

src/initialize.o: src/initialize.cpp
	g++ -std=c++0x src/initialize.cpp -c -o src/initialize.o \
	-Isrc/mpreal

src/common/md5.o: src/common/md5.cpp
	g++ src/common/md5.cpp -c -o src/common/md5.o

.PHONY: build/main
