CXX ?= g++
BUILD ?= debug

WARNINGS := -Wall -Wextra -Wpedantic
SANITIZERS := -fsanitize=address,undefined -fno-sanitize-recover=undefined

ifeq ($(BUILD),debug)
MODE_FLAGS := -O1 -g -fno-omit-frame-pointer $(SANITIZERS) -D_GLIBCXX_ASSERTIONS
MODE_LDFLAGS := $(SANITIZERS)
else ifeq ($(BUILD),release)
MODE_FLAGS := -O2 -DNDEBUG
MODE_LDFLAGS :=
else
$(error BUILD must be debug or release, got '$(BUILD)')
endif

CXXFLAGS ?= -std=c++20 $(WARNINGS) $(MODE_FLAGS)
LDFLAGS ?= $(MODE_LDFLAGS)
GTEST_CFLAGS := $(shell pkg-config --cflags gtest)
GTEST_LIBS := $(shell pkg-config --libs gtest)

BUILD_DIR := build/$(BUILD)
TEST_BINARY := $(BUILD_DIR)/grid_tests

.PHONY: all test clean

all: $(BUILD_DIR)/grid.o $(BUILD_DIR)/FdtdSolver.o $(BUILD_DIR)/output.o

$(BUILD_DIR)/grid.o: grid.cpp grid.h
	@mkdir -p $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c grid.cpp -o $@

$(BUILD_DIR)/FdtdSolver.o: FdtdSolver.cpp FdtdSolver.h grid.h
	@mkdir -p $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c FdtdSolver.cpp -o $@

$(BUILD_DIR)/output.o: output.cpp output.h FdtdSolver.h grid.h
	@mkdir -p $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c output.cpp -o $@

$(TEST_BINARY): grid.cpp grid.h tests/grid_tests.cpp
	@mkdir -p $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) $(GTEST_CFLAGS) -I. grid.cpp tests/grid_tests.cpp -o $@ $(LDFLAGS) $(GTEST_LIBS)

test: $(TEST_BINARY)
	./$(TEST_BINARY)

clean:
	rm -rf build
