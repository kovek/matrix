all: build/main

run: build/main
	./build/main

build/main: src/main.cpp src/Scene.cpp
	g++ src/initialize.cpp src/main.cpp src/Scene.cpp src/common/md5.cpp external/libs/libglfw3.a external/libs/libgmp.a external/libs/libmpfr.a -Iexternal/glfw/include -Iexternal/glm -Iexternal/plplot/include -Iexternal/plplot/build/include -o build/main -Isrc/mpreal -Isrc/common -framework Cocoa -framework CoreVideo -framework OpenGL -framework IOKit -std=c++0x -lmpfr -lplplotcxxd -lglfw3 -lsqlite3 -lboost_serialization -g

.PHONY: build/main
