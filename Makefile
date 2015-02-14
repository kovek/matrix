all: build/main

run: build/main
	./build/main

build/main: src/main.cpp src/Scene.cpp
	g++ src/main.cpp src/Scene.cpp -o build/main -Isrc/mpreal -framework Cocoa -framework CoreVideo -framework OpenGL -framework IOKit -std=c++0x -lmpfr -lplplotcxxd -lglfw3 -lboost_serialization -g

build/linux_main: src/main.cpp src/Scene.cpp
	g++ src/main.cpp src/Scene.cpp -o build/linux_main -Isrc/mpreal -std=c++0x -lmpfr -lglfw3 -lboost_serialization -g

.PHONY: build/main
