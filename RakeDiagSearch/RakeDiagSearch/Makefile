BOINC_DIR = ../../../../boinc-src
BOINC_API_DIR = $(BOINC_DIR)/api
BOINC_LIB_DIR = $(BOINC_DIR)/lib

CXX = g++

CXXFLAGS += -O3 -ftree-vectorize \
    -I$(BOINC_DIR) \
    -I$(BOINC_LIB_DIR) \
    -I$(BOINC_API_DIR) \
    -L$(BOINC_API_DIR) \
    -L$(BOINC_LIB_DIR) \

PROGRAM = rakesearch

ifeq ($(MinGW),1)
PROGRAM = rakesearch.exe
CXX = x86_64-w64-mingw32-g++
endif

ifeq ($(Native),1)
# Tune for machine where app is compiled
CXXFLAGS += -march=native -mtune=native
else
ifeq ($(AVX512),1)
# AVX512+BMI2, 64-bit
CXXFLAGS += -march=core2 -mtune=generic -msse4.2 -mpopcnt -mavx -mavx2 -mbmi -mbmi2 -mavx512f -mavx512cd -mavx512vl -mavx512bw -mavx512dq -m64
else
ifeq ($(AVX2),1)
# AVX2+BMI2, 64-bit
CXXFLAGS += -march=core2 -mtune=generic -msse4.2 -mpopcnt -mavx -mavx2 -mbmi -mbmi2 -m64
else
ifeq ($(AVX),1)
# AVX, 64-bit
CXXFLAGS += -march=core2 -mtune=generic -msse4.2 -mpopcnt -mavx -m64
else
ifeq ($(SSE2),1)
# SSE2, 64-bit
CXXFLAGS += -march=core2 -mtune=generic -m64
endif
endif
endif
endif
endif

all: $(PROGRAM)

clean:
	rm -f $(PROGRAM) *.o

$(PROGRAM): main.o Square.o Generator.o MovePairSearch.o $(BOINC_LIB_DIR)/libboinc.a $(BOINC_API_DIR)/libboinc_api.a
	$(CXX) $(CXXFLAGS) -o $(PROGRAM) main.o Square.o Generator.o MovePairSearch.o -pthread -lboinc_api -lboinc

Square.o: Square.cpp
	$(CXX) $(CXXFLAGS) -c Square.cpp

Generator.o: Generator.cpp
	$(CXX) $(CXXFLAGS) -c Generator.cpp
	
MovePairSearch.o: MovePairSearch.cpp
	$(CXX) $(CXXFLAGS) -c MovePairSearch.cpp 

main.o: main.cpp
	$(CXX) $(CXXFLAGS) -c main.cpp
