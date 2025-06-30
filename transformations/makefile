all:
	g++ -std=c++17 -g -I../util_headers -I./include -L../util_lib -L./bin -c bin/transformations.cpp -o bin/transformations.o -lutil

	ar rcs bin/libtransformations.a bin/transformations.o

	g++ -std=c++17 -g -I../util_headers -I./include -L../util_lib -L./bin main.cpp -o main -lutil -ltransformations