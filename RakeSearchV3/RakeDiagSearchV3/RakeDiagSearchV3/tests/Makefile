all: tests

CFLAGS = -O3 -ftree-vectorize -std=c++11 -g -MMD -MP -Wall -Wextra -Werror \
	-Iboinc -DUT_BUILD
CXX = g++

tests: TestSquare.o TestRakeSearch.o main.o
	$(CXX) -o $@ $^ $(CFLAGS)

%.o: %.cpp
	$(CXX) -c -o $@ $< $(CFLAGS)

-include TestSquare.d TestRakeSearch.d main.d
