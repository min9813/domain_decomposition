CONF_FILE = exp1
LPATH = /usr/local/lib
LDFLAGS = -L/usr/local/lib
YAML_LIBRALY_FLAGS = -lyaml-cpp
CFLAGS = -std=c++17
LPTHREAD_FLAG = -lpthread
BIN_FOLDER = ./bin
SRC_FOLDER = ./src
TARGET = main



mp_test: ${SRC_FOLDER}/mp_test.cpp
	c++ -o main ${SRC_FOLDER}/mp_test.cpp -std=c++17 ${LPTHREAD_FLAG}

${BIN_FOLDER}/main.o: ${SRC_FOLDER}/main.cpp src/numerical.hpp
	c++ -c -o ${BIN_FOLDER}/main.o ${SRC_FOLDER}/main.cpp ${CFLAGS}

${BIN_FOLDER}/pde.o: ${SRC_FOLDER}/pde.cpp src/numerical.hpp
	c++ -c -o ${BIN_FOLDER}/pde.o ${SRC_FOLDER}/pde.cpp ${CFLAGS}

${BIN_FOLDER}/pde_2dim.o: ${SRC_FOLDER}/pde_2dim.cpp src/numerical.hpp
	c++ -c -o ${BIN_FOLDER}/pde_2dim.o ${SRC_FOLDER}/pde_2dim.cpp ${CFLAGS}

${BIN_FOLDER}/pde_2dim_mp.o: ${SRC_FOLDER}/pde_2dim_mp.cpp src/numerical.hpp
	c++ -c -o ${BIN_FOLDER}/pde_2dim_mp.o ${SRC_FOLDER}/pde_2dim_mp.cpp ${CFLAGS}

${BIN_FOLDER}/lu.o: ${SRC_FOLDER}/lu.cpp src/numerical.hpp
	c++ -c -o ${BIN_FOLDER}/lu.o ${SRC_FOLDER}/lu.cpp ${CFLAGS}

${BIN_FOLDER}/exact_solution.o: ${SRC_FOLDER}/exact_solution.cpp src/numerical.hpp
	c++ -c -o ${BIN_FOLDER}/exact_solution.o ${SRC_FOLDER}/exact_solution.cpp ${CFLAGS}

main: ${BIN_FOLDER}/main.o ${BIN_FOLDER}/pde.o ${BIN_FOLDER}/pde_2dim.o ${BIN_FOLDER}/lu.o ${BIN_FOLDER}/exact_solution.o ${BIN_FOLDER}/pde_2dim_mp.o
	c++ -o main ${BIN_FOLDER}/main.o ${BIN_FOLDER}/pde.o \
	${BIN_FOLDER}/pde_2dim.o ${BIN_FOLDER}/pde_2dim_mp.o ${BIN_FOLDER}/lu.o ${BIN_FOLDER}/exact_solution.o \
	${CFLAGS} ${LDFLAGS} ${YAML_LIBRALY_FLAGS} ${LPTHREAD_FLAG}

main_test: src/main.cpp
	c++ -o main_test src/main.cpp src/pde.cpp src/pde_2dim.cpp src/lu.cpp src/exact_solution.cpp -std=c++17 \
	-L${LPATH} -lyaml-cpp

run_test: main_test
	./main_test config/${CONF_FILE}.yaml result/${CONF_FILE}/${CONF_FILE}

run: ${TARGET}
	./main config/${CONF_FILE}.yaml result/${CONF_FILE}/${CONF_FILE}

clean:
	rm -f ${BIN_FOLDER}/*.o
	rm ./main