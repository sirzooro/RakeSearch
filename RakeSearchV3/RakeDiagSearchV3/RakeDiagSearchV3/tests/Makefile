all: tests

#FLAGS = -DNO_SIMD=1
FLAGS = -msse2
#FLAGS = -mssse3
#FLAGS = -msse4.1
#FLAGS = -mavx
#FLAGS = -mavx2 -mbmi -mbmi2
#FLAGS = -march=skylake-avx512

CFLAGS = -O3 -ftree-vectorize -std=c++11 -g -MMD -MP -Wall -Wextra -Werror \
	-Iboinc -DUT_BUILD $(FLAGS)
CXX = g++

tests: TestSquare.o TestRakeSearch.o main.o
	$(CXX) -o $@ $^ $(CFLAGS)

%.o: %.cpp
	$(CXX) -c -o $@ $< $(CFLAGS)

-include TestSquare.d TestRakeSearch.d main.d
