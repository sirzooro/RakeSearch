BOINC_DIR = ../../../../boinc-src
BOINC_API_DIR = $(BOINC_DIR)/api
BOINC_LIB_DIR = $(BOINC_DIR)/lib

CXX = g++

CXXFLAGS += -O3 -ftree-vectorize -std=c++11 -Wall -static -static-libgcc -static-libstdc++ -pthread \
    -I$(BOINC_DIR) \
    -I$(BOINC_LIB_DIR) \
    -I$(BOINC_API_DIR) \
    -L$(BOINC_API_DIR) \
    -L$(BOINC_LIB_DIR) \

PROGRAM = rakesearch

ifeq ($(MinGW32),1)
$(info ===== Compiling MinGW 32-bit app version =====)
CXX = i686-w64-mingw32-g++
else ifeq ($(MinGW64),1)
$(info ===== Compiling MinGW 64-bit app version =====)
CXX = x86_64-w64-mingw32-g++
else ifeq ($(AARCH64),1)
$(info ===== Compiling AARCH64 app version =====)
CXX = aarch64-linux-gnu-g++
else ifeq ($(ARM),1)
ifeq ($(NEON),1)
$(info ===== Compiling ARM HF app version =====)
CXX = arm-linux-gnueabihf-g++
else
$(info ===== Compiling ARM SF app version =====)
CXX = arm-linux-gnueabi-g++
endif
else ifeq ($(M32),1)
$(info ===== Compiling 32-bit app version =====)
CXX = g++ -m32
else
$(info ===== Compiling default app version =====)
CXX = g++
endif

ifeq ($(Native),1)
# Tune for machine where app is compiled
$(info ===== Compiling native app version =====)
CXXFLAGS += -march=native -mtune=native
else ifeq ($(AVX512),1)
# AVX512+BMI2
$(info ===== Compiling AVX512+BMI2 app version =====)
CXXFLAGS += -march=skylake-avx512 -mprefer-vector-width=256
else ifeq ($(AVX2),1)
# AVX2+BMI2
$(info ===== Compiling AVX2+BMI2 app version =====)
CXXFLAGS += -march=core2 -mtune=haswell -mpopcnt -mavx2 -mbmi -mbmi2
else ifeq ($(AVX),1)
# AVX
$(info ===== Compiling AVX app version =====)
CXXFLAGS += -march=core2 -mtune=sandybridge -mpopcnt -mavx
else ifeq ($(SSSE3),1)
# SSSE3
$(info ===== Compiling SSSE3 app version =====)
CXXFLAGS += -mtune=core2 -mssse3
else ifeq ($(SSE2),1)
# SSE2
$(info ===== Compiling SSE2 app version =====)
CXXFLAGS += -mtune=core2 -msse2
else ifeq ($(NEON),1)
# NEON (ARM only; AARCH64 has NEON by default)
$(info ===== Compiling ARM NEON app version =====)
CXXFLAGS += -mfloat-abi=hard -mfpu=neon
else ifeq ($(NOSIMD),1)
$(info ===== Compiling No-SIMD app version =====)
CXXFLAGS += -DNO_SIMD=1
else
$(info ===== Compiling default app version =====)
endif

all: $(PROGRAM)

clean:
	rm -f $(PROGRAM) $(PROGRAM).exe *.o

$(PROGRAM): main.o Square.o Generator.o MovePairSearch.o $(BOINC_LIB_DIR)/libboinc.a $(BOINC_API_DIR)/libboinc_api.a
	$(CXX) $(CXXFLAGS) -o $(PROGRAM) main.o Square.o Generator.o MovePairSearch.o -lboinc_api -lboinc

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@
