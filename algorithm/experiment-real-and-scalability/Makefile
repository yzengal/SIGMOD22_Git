CC = gcc
CXX = g++
LIBS = 
# BOOSTLIB = -I /export/data/anonymous/boost -L /export/data/anonymous/boost/stage/lib -lboost_system
BOOSTLIB = -I /home/anonymous/boost -L /home/anonymous/boost/stage/lib -lboost_system
ANNLIB = -I /home/anonymous/ann/include -L /home/anonymous/ann/lib -lANN
CFLAGS = --std=c++11 -W -Wall -Wno-sign-compare -O3 -s -pipe -mmmx -msse -msse2 -msse3 # -mcmodel=medium
MEM = -D WATCH_MEM

all: DC DCnn DCsam FRT Bar96 Bar98

monitor.o: monitor.h monitor.cpp
	$(CXX) $(CFLAGS) -c monitor.cpp $(LIBS) $(MEM)
	
global.o: global.h global.cpp 
	$(CXX) $(CFLAGS) -c global.cpp $(LIBS) $(MEM)

HST_BASE.o: HST_BASE.h HST_BASE.cpp global.h global.cpp global.o monitor.h monitor.cpp monitor.o
	$(CXX) $(CFLAGS) -c global.o monitor.o HST_BASE.cpp $(LIBS) $(MEM)

FRT: FRT.cpp global.h global.cpp global.o HST_BASE.h HST_BASE.cpp HST_BASE.o monitor.h monitor.cpp monitor.o
	$(CXX) $(CFLAGS) -o FRT global.o monitor.o HST_BASE.o FRT.cpp $(LIBS) $(MEM)

HST_BARTAL.o: HST_BARTAL.h HST_BARTAL.cpp global.h global.cpp global.o monitor.h monitor.cpp monitor.o
	$(CXX) $(CFLAGS) -c global.o monitor.o HST_BARTAL.cpp $(LIBS) $(MEM)
	
Bar98: Bar98.cpp global.h global.cpp global.o HST_BARTAL.h HST_BARTAL.cpp HST_BARTAL.o monitor.h monitor.cpp monitor.o
	$(CXX) $(CFLAGS) -o Bar98 global.o monitor.o HST_BARTAL.o Bar98.cpp $(LIBS) $(MEM)
	
Bar96: Bar96.cpp global.h global.cpp global.o HST_BARTAL.h HST_BARTAL.cpp HST_BARTAL.o monitor.h monitor.cpp monitor.o
	$(CXX) $(CFLAGS) -o Bar96 global.o monitor.o HST_BARTAL.o Bar96.cpp $(LIBS) $(MEM)

HST_GDCO.o: HST_GDCO.h HST_GDCO.cpp global.h global.cpp global.o monitor.h monitor.cpp monitor.o
	$(CXX) $(CFLAGS) -c global.o monitor.o HST_GDCO.cpp $(LIBS) $(BOOSTLIB) $(MEM)

DC: DC.cpp global.h global.cpp global.o HST_GDCO.h HST_GDCO.cpp HST_GDCO.o monitor.h monitor.cpp monitor.o
	$(CXX) $(CFLAGS) -o DC global.o monitor.o HST_GDCO.o DC.cpp $(LIBS) $(BOOSTLIB) $(MEM)
			
DCnn: DCnn.cpp global.h global.cpp global.o HST_GDCO.h HST_GDCO.cpp HST_GDCO.o monitor.h monitor.cpp monitor.o
	$(CXX) $(CFLAGS) -o DCnn global.o monitor.o HST_GDCO.o DCnn.cpp $(LIBS) $(BOOSTLIB) $(MEM)
	
HST_PARS_Fast.o: HST_PARS_Fast.h HST_PARS_Fast.cpp global.h global.cpp global.o monitor.h monitor.cpp monitor.o
	$(CXX) $(CFLAGS) -c global.o monitor.o HST_PARS_Fast.cpp $(LIBS) $(BOOSTLIB) $(MEM)

DCsam: DCsam.cpp global.h global.cpp global.o HST_PARS_Fast.h HST_PARS_Fast.cpp HST_PARS_Fast.o monitor.h monitor.cpp monitor.o
	$(CXX) $(CFLAGS) -o DCsam global.o monitor.o HST_PARS_Fast.o DCsam.cpp $(LIBS) $(BOOSTLIB) $(MEM)
	
.PHONY: clean
clean:
		-@rm *.o *.gcno *~ 2> /dev/null || true
		-@rm DC DCnn DCsam FRT Bar96 Bar98 2> /dev/null || true
		
		