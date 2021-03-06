LD_PTHREAD = -pthread

ifeq ($(MinGW32),1)
$(info ===== Compiling MinGW 32-bit app version =====)
CXX = i686-w64-mingw32-g++
BOINC_DIR = /boinc930/mingw32
else ifeq ($(MinGW64),1)
$(info ===== Compiling MinGW 64-bit app version =====)
CXX = x86_64-w64-mingw32-g++
BOINC_DIR = /boinc930/mingw64
else ifeq ($(AARCH64),1)
$(info ===== Compiling AARCH64 app version =====)
CXX = aarch64-unknown-linux-gnu-g++
BOINC_DIR = /boinc930/aarch64
else ifeq ($(ARM),1)
# Fix for "terminate called after throwing an instance of '__gnu_cxx::__concurrence_broadcast_error'" on ARM
LD_PTHREAD = -Wl,-whole-archive -lpthread -Wl,-no-whole-archive
ifeq ($(NEON),1)
$(info ===== Compiling ARMv7 HF app version =====)
CXX = arm-unknown-linux-gnueabihf-g++
BOINC_DIR = /boinc930/armhf
else
$(info ===== Compiling ARMv7 SF app version =====)
CXX = arm-unknown-linux-gnueabi-g++
BOINC_DIR = /boinc930/arm
# Fix for "terminate called after throwing an instance of '__gnu_cxx::__concurrence_broadcast_error'" on ARM
LD_PTHREAD = -Wl,-whole-archive -lpthread -Wl,-no-whole-archive
endif
else ifeq ($(ARMv6),1)
$(info ===== Compiling ARMv6 app version =====)
CXX = armv6-unknown-linux-gnueabi-g++
BOINC_DIR = /boinc930/armv6
else ifeq ($(M32),1)
$(info ===== Compiling 32-bit app version =====)
CXX = g++ -m32 -mtune=generic -march=pentiumpro
BOINC_DIR = /boinc930/linux32
else
$(info ===== Compiling default app version =====)
CXX = g++
BOINC_DIR = /boinc930/linux64
endif

ifeq ($(Native),1)
# Tune for machine where app is compiled
$(info ===== Compiling native app version =====)
TARGET_FLAGS = -march=native -mtune=native
else ifeq ($(AVX512),1)
# AVX512+BMI2
$(info ===== Compiling AVX512+BMI2 app version =====)
TARGET_FLAGS = -march=skylake-avx512 -mprefer-vector-width=256
else ifeq ($(AVX2),1)
# AVX2+BMI2
$(info ===== Compiling AVX2+BMI2 app version =====)
TARGET_FLAGS = -march=core2 -mtune=haswell -mpopcnt -mavx2 -mbmi -mbmi2
else ifeq ($(AVX),1)
# AVX
$(info ===== Compiling AVX app version =====)
TARGET_FLAGS = -march=core2 -mtune=sandybridge -mpopcnt -mavx -mprefer-vector-width=128
else ifeq ($(SSE41),1)
# SSE4.1
$(info ===== Compiling SSE4.1 app version =====)
TARGET_FLAGS = -mtune=core2 -msse4.1
else ifeq ($(SSSE3),1)
# SSSE3
$(info ===== Compiling SSSE3 app version =====)
TARGET_FLAGS = -mtune=core2 -mssse3
else ifeq ($(SSE2),1)
# SSE2
$(info ===== Compiling SSE2 app version =====)
TARGET_FLAGS = -mtune=core2 -msse2
else ifeq ($(NEON),1)
# NEON (ARM only; AARCH64 has NEON by default)
$(info ===== Compiling ARM NEON app version =====)
TARGET_FLAGS = -mfloat-abi=hard -mfpu=neon
else ifeq ($(NOSIMD),1)
$(info ===== Compiling No-SIMD app version =====)
TARGET_FLAGS = -DNO_SIMD=1
else
$(info ===== Compiling default app version =====)
endif

CXXFLAGS = $(TARGET_FLAGS) -O3 -ftree-vectorize -pthread -std=c++11 -Wall \
    -I$(BOINC_DIR)/include/boinc

LDFLAGS = $(TARGET_FLAGS) -O3 -ftree-vectorize -static -static-libgcc -static-libstdc++ $(LD_PTHREAD) -std=c++11 -Wall \
    -L$(BOINC_DIR)/lib

PROGRAM = rakesearch10

all: $(PROGRAM)

OBJ_FILES = main.o Square.o RakeSearch.o

clean:
	rm -f $(PROGRAM) $(PROGRAM).exe *.o

$(PROGRAM): $(OBJ_FILES)
	$(CXX) $(LDFLAGS) -o $(PROGRAM) $(OBJ_FILES) -lboinc_api -lboinc

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@
