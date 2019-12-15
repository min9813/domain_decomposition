FNAME = none
CONF_FILE = exp1
LPATH = /usr/local/lib

# main: main.o keisan.o
	# g++ -o main main.o keisan.o

# main.o: main.cpp
	# g++ -c main.cpp -o main.o

keisan.o: keisan.cpp
	c++ -c -o keisan.o keisan.cpp

pde: src/pde.cpp
	c++ -o main src/pde.cpp src/lu.cpp -std=c++17

pde_2dim: src/pde_2dim.cpp
	c++ -o main src/pde_2dim.cpp src/lu.cpp -std=c++17 -L${LPATH} -lyaml-cpp

mp_test:mp_test.cpp
	c++ -o main mp_test.cpp -std=c++17 -lpthread

run: ${FNAME}
	./main config/${CONF_FILE}.yaml result/${CONF_FILE}/${CONF_FILE}.txt

clean:
	rm -f *.o