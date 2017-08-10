
BOINC_DIR = /home/natalia/work/Conferences/2017_FAST/MSquares/boinc
BOINC_API_DIR = $(BOINC_DIR)/api
BOINC_LIB_DIR = $(BOINC_DIR)/lib

CXXFLAGS += -g \
    -I$(BOINC_DIR) \
    -I$(BOINC_LIB_DIR) \
    -I$(BOINC_API_DIR) \
    -L$(BOINC_API_DIR) \
    -L$(BOINC_LIB_DIR)

PROGS = main 

all: $(PROGS)

clean: distclean

distclean:
	rm -f $(PROGS) *.o

main: main.o Square.o Generator.o MovePairSearch.o $(BOINC_LIB_DIR)/libboinc.a $(BOINC_API_DIR)/libboinc_api.a
	g++ $(CXXFLAGS) $(CPPFLAGS) $(LDFLAGS) -o main main.o Square.o Generator.o MovePairSearch.o -pthread -lboinc_api -lboinc

Square.o: Square.cpp
	g++ -c Square.cpp

Generator.o: Generator.cpp
	g++ -c Generator.cpp
	
MovePairSearch.o: MovePairSearch.cpp
	g++ -c MovePairSearch.cpp

