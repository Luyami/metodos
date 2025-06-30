#Util lib maker

all:
	g++ -std=c++17 -I./util_headers -c util_cpp/vec.cpp -o util_lib/vec.o
	g++ -std=c++17 -I./util_headers -c util_cpp/mat.cpp -o util_lib/mat.o

	ar rcs util_lib/libutil.a util_lib/vec.o util_lib/mat.o