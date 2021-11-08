CC = gcc
CXX = g++
LIBS = 
BOOSTLIB = -I /home/anonymous/boost -L /home/anonymous/boost/stage/lib -lboost_system
ANNLIB = -I /home/data/1TB/yzengal/ann/include -L /home/data/1TB/yzengal/ann/lib -lANN
CFLAGS = --std=c++11 -W -Wall -Wno-sign-compare -O3 -s -pipe -mmmx -msse -msse2 -msse3 # -mcmodel=medium
MEM = -D WATCH_MEM

all: HSF HSFdc
monitor.o: monitor.h monitor.cpp
	$(CXX) $(CFLAGS) -c monitor.cpp $(LIBS) -D WATCH_MEM
	
global.o: global.h global.cpp
	$(CXX) $(CFLAGS) -c global.cpp $(LIBS) -D WATCH_MEM

HST_opt.o: HST_opt.h HST_opt.cpp global.h global.cpp global.o monitor.h monitor.cpp monitor.o
	$(CXX) $(CFLAGS) -c global.o monitor.o HST_opt.cpp $(LIBS) $(BOOSTLIB) $(MEM)
		
HSF: HSF.cpp global.h global.cpp global.o HST_opt.h HST_opt.cpp HST_opt.o monitor.h monitor.cpp monitor.o
	$(CXX) $(CFLAGS) -o HSF monitor.o global.o HST_opt.o HSF.cpp $(LIBS) $(BOOSTLIB) $(MEM)
	
global_DC.o: global_DC.h global_DC.cpp
	$(CXX) $(CFLAGS) -c global_DC.cpp $(LIBS) -D WATCH_MEM

HST_DC.o: HST_DC.h HST_DC.cpp global_DC.h global_DC.cpp global_DC.o monitor.h monitor.cpp monitor.o
	$(CXX) $(CFLAGS) -c global_DC.o monitor.o HST_DC.cpp $(LIBS) $(BOOSTLIB) $(MEM)
		
HSFdc: HSF_DC.cpp global_DC.h global_DC.cpp global_DC.o HST_DC.h HST_DC.cpp HST_DC.o monitor.h monitor.cpp monitor.o
	$(CXX) $(CFLAGS) -o HSFdc monitor.o global_DC.o HST_DC.o HSF_DC.cpp $(LIBS) $(BOOSTLIB) $(MEM)
		
.PHONY: clean
clean:
		-@rm *.o *.gcno *~ 2> /dev/null || true
		-@rm HSF HSFdc 2> /dev/null || true
		
