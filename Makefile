path_only = `dirname $(realpath $(lastword $(MAKEFILE_LIST)))`

all: heuristics

LIB_PATH = ${path_only}/lib
INC_PATH = ${path_only}/include
SRC_PATH = ${path_only}/src
OBJ_PATH = ${path_only}/

CPP = g++
PEDANTIC_PARANOID_FREAK =       -O0 -Wshadow -Wcast-align \
				-Waggregate-return -Wmissing-prototypes -Wmissing-declarations \
				-Wstrict-prototypes -Wmissing-prototypes -Wmissing-declarations \
				-Wmissing-noreturn -Wredundant-decls -Wnested-externs \
				-Wpointer-arith -Wwrite-strings -finline-functions
REASONABLY_CAREFUL_DUDE =	-Wall
NO_PRAYER_FOR_THE_WICKED =	-w -O2 
WARNINGS = $(REASONABLY_CAREFUL_DUDE)
CFLAGS = $(WARNINGS)  -O3 -g -DNOASSERT -std=c++14 -Xpreprocessor -fopenmp
INCLUDES = -I${INC_PATH}
DEFS = 
LDADD =
LIBS = ${LIB_PATH}/heuristics.a 



# all source files have associated object files

lib-io-tree-utils.o: src/lib-io-tree-utils.cpp include/lib-io-tree-utils.h
	$(CPP) $(INCLUDES) $(DEFS) $(CFLAGS) -o ${OBJ_PATH}/$@ -c $< -fPIC

cluster.o: src/cluster.cpp include/cluster.h lib-io-tree-utils.o
	$(CPP) $(INCLUDES) $(DEFS) $(CFLAGS) -o ${OBJ_PATH}/$@ -c $< -fPIC

lib-io-tree-liu-optimal.o: src/lib-io-tree-liu-optimal.cpp include/lib-io-tree-liu-optimal.h lib-io-tree-utils.o
	$(CPP) $(INCLUDES) $(DEFS) $(CFLAGS) -o ${OBJ_PATH}/$@ -c $< -fPIC

lib-io-tree-minmem.o: src/lib-io-tree-minmem.cpp include/lib-io-tree-minmem.h lib-io-tree-utils.o
	$(CPP) $(INCLUDES) $(DEFS) $(CFLAGS) -o ${OBJ_PATH}/$@ -c $< -fPIC

lib-io-tree.o: src/lib-io-tree.cpp include/lib-io-tree.h lib-io-tree-minmem.o lib-io-tree-utils.o lib-io-tree-liu-optimal.o
	$(CPP) $(INCLUDES) $(DEFS) $(CFLAGS) -o ${OBJ_PATH}/$@ -c $< -fPIC

heuristics.o: src/heuristics.cpp include/heuristics.h lib-io-tree-minmem.o lib-io-tree-liu-optimal.o lib-io-tree.o cluster.o
	$(CPP) $(INCLUDES) $(DEFS) $(CFLAGS) -o ${OBJ_PATH}/$@ -c $< -fPIC


heuristics:heuristics.o
	rm -f ${LIB_PATH}/$@.a ${LIB_PATH}/$@.so
	ar rcs ${LIB_PATH}/$@.a $(OBJ_PATH)/heuristics.o $(OBJ_PATH)/lib-io-tree.o $(OBJ_PATH)/lib-io-tree-utils.o $(OBJ_PATH)/lib-io-tree-minmem.o $(OBJ_PATH)/lib-io-tree-liu-optimal.o $(OBJ_PATH)/cluster.o

clean:
	rm -f ${OBJ_PATH}/*.o *~ 

.PHONY: clean all

.SUFFIXES:
.SECONDARY: