all:
	g++ -std=c++17 -g -I../util_headers -I./include -L../util_lib -L./bin -c bin/eigenvalues.cpp -o bin/eigenvalues.o -lutil

	ar rcs bin/libeigen.a bin/eigenvalues.o

	g++ -std=c++17 -g -I../util_headers -I./include -L../util_lib -L./bin main.cpp -o main -lutil -leigen