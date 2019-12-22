FNAME = none
CONF_FILE = exp1
LPATH = /usr/local/lib



pde: src/pde.cpp
	c++ -o main src/pde.cpp src/lu.cpp -std=c++17

pde_2dim: src/pde_2dim.cpp
	c++ -o main src/pde_2dim.cpp src/lu.cpp src/exact_solution.cpp -std=c++17 -L${LPATH} -lyaml-cpp

mp_test:mp_test.cpp
	c++ -o main mp_test.cpp -std=c++17 -lpthread

main: src/main.cpp
	c++ -o main src/main.cpp src/pde.cpp src/pde_2dim.cpp src/lu.cpp src/exact_solution.cpp -std=c++17 \
	-L${LPATH} -lyaml-cpp

run: ${FNAME}
	./main config/${CONF_FILE}.yaml result/${CONF_FILE}/${CONF_FILE}

clean:
	rm -f *.o
	rm main