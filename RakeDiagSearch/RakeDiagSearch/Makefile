BOINC_DIR = ../../../../boinc-src
BOINC_API_DIR = $(BOINC_DIR)/api
BOINC_LIB_DIR = $(BOINC_DIR)/lib

CXX = g++

CXXFLAGS += -O3 -ftree-vectorize -std=c++11 -static-libgcc -static-libstdc++ \
    -I$(BOINC_DIR) \
    -I$(BOINC_LIB_DIR) \
    -I$(BOINC_API_DIR) \
    -L$(BOINC_API_DIR) \
    -L$(BOINC_LIB_DIR) \

PTHREAD = -pthread

PROGRAM = rakesearch

ifeq ($(MinGW),1)
$(info ===== Compiling MinGW app version =====)
PROGRAM = rakesearch.exe
CXX = x86_64-w64-mingw32-g++
PTHREAD = -static -pthread
endif

ifeq ($(Native),1)
# Tune for machine where app is compiled
$(info ===== Compiling native app version =====)
CXXFLAGS += -march=native -mtune=native
else ifeq ($(AVX512),1)
# AVX512+BMI2, 64-bit
$(info ===== Compiling AVX512+BMI2 app version =====)
CXXFLAGS += -march=skylake-avx512 -mprefer-vector-width=256 -m64
else ifeq ($(AVX2),1)
# AVX2+BMI2, 64-bit
$(info ===== Compiling AVX2+BMI2 app version =====)
CXXFLAGS += -march=core2 -mtune=haswell -mpopcnt -mavx2 -mbmi -mbmi2 -m64
else ifeq ($(AVX),1)
# AVX, 64-bit
$(info ===== Compiling AVX app version =====)
CXXFLAGS += -march=core2 -mtune=sandybridge -mpopcnt -mavx -m64
else ifeq ($(SSE2),1)
# SSE2, 64-bit
$(info ===== Compiling SSE2 app version =====)
CXXFLAGS += -march=core2 -mtune=generic -m64
else ifeq ($(NEON),1)
# NEON (ARM only; AARCH64 has NEON by default)
$(info ===== Compiling ARM NEON app version =====)
CXXFLAGS += -mfloat-abi=hard -mfpu=neon
else
$(info ===== Compiling default app version =====)
endif

all: $(PROGRAM)

clean:
	rm -f $(PROGRAM) *.o

$(PROGRAM): main.o Square.o Generator.o MovePairSearch.o $(BOINC_LIB_DIR)/libboinc.a $(BOINC_API_DIR)/libboinc_api.a
	$(CXX) $(CXXFLAGS) -o $(PROGRAM) main.o Square.o Generator.o MovePairSearch.o $(PTHREAD) -lboinc_api -lboinc

Square.o: Square.cpp
	$(CXX) $(CXXFLAGS) -c Square.cpp

Generator.o: Generator.cpp
	$(CXX) $(CXXFLAGS) -c Generator.cpp
	
MovePairSearch.o: MovePairSearch.cpp
	$(CXX) $(CXXFLAGS) -c MovePairSearch.cpp 

main.o: main.cpp
	$(CXX) $(CXXFLAGS) -c main.cpp
