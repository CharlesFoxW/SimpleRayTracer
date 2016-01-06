# Makefile - Yihao Wang - ECS 175 Project 4
all:
	g++ main.cpp elements.cpp elements.h AntTweakBar.h -std=c++11 -Wl,-rpath '-Wl,$$ORIGIN/lib' /lib/libAntTweakBar.so -o hw5.out  -lGL -lglut

clean:
	rm -f *.out *.o core *.core data_output.txt
