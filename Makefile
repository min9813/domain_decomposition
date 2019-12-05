FNAME = none

# main: main.o keisan.o
	# g++ -o main main.o keisan.o

# main.o: main.cpp
	# g++ -c main.cpp -o main.o

keisan.o: keisan.cpp
	c++ -c -o keisan.o keisan.cpp

pde: src/pde.cpp
	c++ -o main src/pde.cpp src/lu.cpp -std=c++14

mp_test:mp_test.cpp
	c++ -o main mp_test.cpp -std=c++11 -lpthread

run: ${FNAME}
	./main

clean:
	rm -f *.o