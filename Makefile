all:
	cc main2.c -o main -Wall -std=c99 -DDEBUG -DW=512 -DH=512 -framework OpenGL -framework GLUT
run: all
	./main
step:
	g++ step.cpp -o step -I /usr/local/include/ -lglfw3 -framework Cocoa -framework CoreVideo -lGLEW -lGLUT -framework OpenGl  -framework IOKit; ./step
draw:
	g++ draw.cpp -o draw -I /usr/local/include/ -lglfw3 -framework Cocoa -framework CoreVideo -lGLEW -framework OpenGl  -lGLUT -framework IOKit; ./draw
test: shader.o
	g++ test.cpp shader.o -o test -I /usr/local/include/ -I./ -lglfw3 -framework Cocoa -framework CoreVideo -lGLEW -lGLUT -framework OpenGL -framework IOKit -std=c++0x ; ./test
shader.o:
	g++ -c common/shader.cpp -o shader.o
matrix:
	g++ matrix.cpp Scene.cpp -o main -I./ -lglfw3 -framework Cocoa -framework CoreVideo -framework OpenGL -framework IOKit -std=c++0x;
	./main
test_object: matrix.o Scene.o
	g++ matrix.o Scene.o -o main -lglfw3 -framework Cocoa -framework CoreVideo -lGLEW -lGLUT -framework OpenGL -framework IOKit;
matrix.o:
	g++ -c matrix.cpp Scene.cpp -I /usr/local/include/ -I./ -lglfw3 -framework Cocoa -framework CoreVideo -lGLEW -lGLUT -framework OpenGL -framework IOKit -stdlib=libc++ -std=c++0x;
Scene.o:
	g++ -c Scene.cpp -I /usr/local/include/ -I./ -lglfw3 -framework Cocoa -framework CoreVideo -lGLEW -lGLUT -framework OpenGL -framework IOKit -stdlib=libc++ -std=c++0x;


graph:
	g++ write.cpp -lplplotcxxd -o graph
	echo "done";
	./graph
.PHONY: step, draw, test, matrix.o, Scene.o


