path_only = `dirname $(realpath $(lastword $(MAKEFILE_LIST)))`

all: heuristics
allmain: main

LIB_PATH = ${path_only}/lib
INC_PATH = ${path_only}/include
SRC_PATH = ${path_only}/src
OBJ_PATH = ${path_only}/
BIN_PATH = ${path_only}/
TEST_BIN_PATH = ${path_only}/test

CPP = g++
PEDANTIC_PARANOID_FREAK =       -O0 -Wshadow -Wcast-align \
				-Waggregate-return -Wmissing-prototypes -Wmissing-declarations \
				-Wstrict-prototypes -Wmissing-prototypes -Wmissing-declarations \
				-Wmissing-noreturn -Wredundant-decls -Wnested-externs \
				-Wpointer-arith -Wwrite-strings -finline-functions
REASONABLY_CAREFUL_DUDE =	-Wall
NO_PRAYER_FOR_THE_WICKED =	-w -O0
WARNINGS = $(REASONABLY_CAREFUL_DUDE)
CFLAGS = $(WARNINGS)  -O0 -g -DNOASSERT -std=c++14 -Xpreprocessor -fopenmp
INCLUDES = -I${INC_PATH}
DEFS = 
LDADD =
LIBS = ${LIB_PATH}/heuristics.a 



# all source files have associated object files

lib-io-tree-utils.o: src/tree.cpp include/tree.h
	$(CPP) $(INCLUDES) $(DEFS) $(CFLAGS) -o ${OBJ_PATH}/$@ -c $< -fPIC

cluster.o: src/cluster.cpp include/cluster.h lib-io-tree-utils.o
	$(CPP) $(INCLUDES) $(DEFS) $(CFLAGS) -o ${OBJ_PATH}/$@ -c $< -fPIC

lib-io-tree-minmem.o: src/lib-io-tree-minmem.cpp include/lib-io-tree-minmem.h lib-io-tree-utils.o
	$(CPP) $(INCLUDES) $(DEFS) $(CFLAGS) -o ${OBJ_PATH}/$@ -c $< -fPIC

lib-io-tree.o: src/lib-io-tree.cpp include/lib-io-tree.h lib-io-tree-minmem.o lib-io-tree-utils.o
	$(CPP) $(INCLUDES) $(DEFS) $(CFLAGS) -o ${OBJ_PATH}/$@ -c $< -fPIC

heuristics.o: src/heuristics.cpp include/heuristics.h lib-io-tree-minmem.o  lib-io-tree.o cluster.o
	$(CPP) $(INCLUDES) $(DEFS) $(CFLAGS) -o ${OBJ_PATH}/$@ -c $< -fPIC

heuristics:heuristics.o
	rm -f ${LIB_PATH}/$@.a ${LIB_PATH}/$@.so
	ar rcs ${LIB_PATH}/$@.a $(OBJ_PATH)/heuristics.o $(OBJ_PATH)/lib-io-tree.o $(OBJ_PATH)/lib-io-tree-utils.o $(OBJ_PATH)/lib-io-tree-minmem.o $(OBJ_PATH)/cluster.o

heuristicsnoclean:heuristics.o
	ar rcs ${LIB_PATH}/heuristics.a $(OBJ_PATH)/heuristics.o $(OBJ_PATH)/lib-io-tree.o $(OBJ_PATH)/lib-io-tree-utils.o $(OBJ_PATH)/lib-io-tree-minmem.o $(OBJ_PATH)/cluster.o

main: main.cpp heuristicsnoclean cluster.o
	$(CPP) $(INCLUDES) $(DEFS) $(CFLAGS) $< $(LIBS) $(LDADD) -o ${BIN_PATH}/$@

test: test/test.cpp heuristicsnoclean cluster.o
	$(CPP) $(INCLUDES) $(DEFS) $(CFLAGS) $< $(LIBS) $(LDADD) -o ${BIN_PATH}/$@


clean:
	rm -f ${OBJ_PATH}/*.o *~
	rm -f ${BIN_PATH}/main
	rm -f ${LIB_PATH}/heuristics.a


############################################################################################################################
#GTEST
# Points to the root of Google Test, relative to where this file is.
GTEST_DIR = ./googletest-main/googletest

# Where to find user code.
USER_DIR = ./test

# Flags passed to the preprocessor.
# Set Google Test's header directory as a system directory, such that
# the compiler doesn't generate warnings in Google Test headers.
CPPFLAGStest += -isystem $(GTEST_DIR)/include

# Flags passed to the C++ compiler.
CXXFLAGStest += -std=c++14 -g -Wall -Wextra -pthread

# All tests produced by t his Makefile.  Remember to add new tests you
# created to the list.
TESTS = cluster_unittest task_tree_unittest

# All Google Test headers.  Usually you shouldn't change this
# definition.
GTEST_HEADERS = $(GTEST_DIR)/include/gtest/*.h \
                $(GTEST_DIR)/include/gtest/internal/*.h

# House-keeping build targets.
alltests : $(TESTS)

cleantests :
	rm -f $(TESTS) gtest.a gtest_main.a *.o


GTEST_SRCS_ = $(GTEST_DIR)/src/*.cc $(GTEST_DIR)/src/*.h $(GTEST_HEADERS)

gtest-all.o : $(GTEST_SRCS_)
	$(CXX) $(CPPFLAGStest) -I$(GTEST_DIR) -I$(GTEST_DIR)/include $(CXXFLAGStest) -c \
            $(GTEST_DIR)/src/gtest-all.cc

gtest_main.o : $(GTEST_SRCS_)
	$(CXX) $(CPPFLAGStest) -I$(GTEST_DIR) -I$(GTEST_DIR)/include $(CXXFLAGStest) -c \
            $(GTEST_DIR)/src/gtest_main.cc

gtest.a : gtest-all.o
	$(AR) $(ARFLAGS) $@ $^

gtest_main.a : gtest-all.o gtest_main.o
	$(AR) $(ARFLAGS) $@ $^

# add new target here if new test files are created
sample1.o: $(USER_DIR)/test_tree.cpp $(GTEST_HEADERS)
	$(CXX) $(CPPFLAGStest) $(CXXFLAGStest) -c $(USER_DIR)/test_tree.cpp -o $@

sample2.o: $(USER_DIR)/test_cluster.cpp $(GTEST_HEADERS)
	$(CXX) $(CPPFLAGStest) $(CXXFLAGStest) -c $(USER_DIR)/test_cluster.cpp -o $@

task_tree_unittest : sample1.o gtest_main.a lib/heuristics.a
	$(CXX) $(CPPFLAGStest) $(CXXFLAGStest) $(LIBS) -lpthread $^ -o ${TEST_BIN_PATH}/$@

cluster_unittest : sample2.o gtest_main.a lib/heuristics.a
	$(CXX) $(CPPFLAGStest) $(CXXFLAGStest) $(LIBS) -lpthread $^ -o ${TEST_BIN_PATH}/$@

.PHONY: clean all

.SUFFIXES:
.SECONDARY: