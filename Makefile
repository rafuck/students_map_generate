all: map
CPP      = clang-omp++
OBJECTS  = main.o
PATH_LIB = -L/usr/local/lib
PATH_INC = -I/usr/local/include -I/usr/local/Cellar/libiomp/20150227/include/libiomp/
FLAGS    = -fopenmp
#libiomp5.dylib

map: $(OBJECTS)
	$(CPP) $(PATH_LIB) $^ $(FLAGS) -o $@
main.o: main.cpp
	$(CPP) $(PATH_INC) -O2 $(FLAGS) $< -c -o $@

clean:
	rm -rf *.o map
