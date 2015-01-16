all: build/main
	
run: build/main
	./build/main

build/main: src/main.cpp src/Scene.cpp
	g++ src/main.cpp src/Scene.cpp -o build/main -Isrc/mpreal -framework Cocoa -framework CoreVideo -framework OpenGL -framework IOKit -std=c++0x -lmpfr -lplplotcxxd -lglfw3

.PHONY: build/main
