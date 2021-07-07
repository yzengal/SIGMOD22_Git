CC = gcc
CXX = g++
LIBS = 
BOOSTLIB = -I /home/anonymous/boost -L /home/anonymous/boost/stage/lib -lboost_system
ANNLIB = -I /home/anonymous/ann/include -L /home/anonymous/ann/lib -lANN
CFLAGS = --std=c++11 -W -Wall -Wno-sign-compare -O3 -s -pipe -mmmx -msse -msse2 -msse3 # -mcmodel=medium
MEM = -D WATCH_MEM

all: FRT-2 FRT-3 FRT-4 FRT-5 FRT-10 FRT-20 FRT-100 Bar98-2 Bar98-3 Bar98-4 Bar98-5 Bar98-10 Bar98-20 Bar98-100 Bar96-2 Bar96-3 Bar96-4 Bar96-5 Bar96-10 Bar96-20 Bar96-100 DCnn-2 DCnn-3 DCnn-4 DCnn-5 DCnn-10 DCnn-20 DCnn-100 DCsam-2 DCsam-3 DCsam-4 DCsam-5 DCsam-10 DCsam-20 DCsam-100
monitor.o: monitor.h monitor.cpp
	$(CXX) $(CFLAGS) -c monitor.cpp $(LIBS) $(MEM)
	
global-2.o: global.h global.cpp 
	$(CXX) $(CFLAGS) -o global-2.o -c global.cpp $(LIBS) $(MEM) -DDIM_V=2

global-3.o: global.h global.cpp 
	$(CXX) $(CFLAGS) -o global-3.o -c global.cpp $(LIBS) $(MEM) -DDIM_V=3

global-4.o: global.h global.cpp 
	$(CXX) $(CFLAGS) -o global-4.o -c global.cpp $(LIBS) $(MEM) -DDIM_V=4

global-5.o: global.h global.cpp 
	$(CXX) $(CFLAGS) -o global-5.o -c global.cpp $(LIBS) $(MEM) -DDIM_V=5

global-10.o: global.h global.cpp 
	$(CXX) $(CFLAGS) -o global-10.o -c global.cpp $(LIBS) $(MEM) -DDIM_V=10

global-20.o: global.h global.cpp 
	$(CXX) $(CFLAGS) -o global-20.o -c global.cpp $(LIBS) $(MEM) -DDIM_V=20

global-100.o: global.h global.cpp 
	$(CXX) $(CFLAGS) -o global-100.o -c global.cpp $(LIBS) $(MEM) -DDIM_V=100

HST_BASE-2.o: HST_BASE.h HST_BASE.cpp global.h global.cpp global-2.o monitor.h monitor.cpp monitor.o
	$(CXX) $(CFLAGS) -o HST_BASE-2.o -c global-2.o monitor.o HST_BASE.cpp $(LIBS) $(MEM) -DDIM_V=2

FRT-2: FRT.cpp global.h global.cpp global-2.o HST_BASE.h HST_BASE.cpp HST_BASE-2.o monitor.h monitor.cpp monitor.o
	$(CXX) $(CFLAGS) -o FRT-2 global-2.o monitor.o HST_BASE-2.o FRT.cpp $(LIBS) $(MEM) -DDIM_V=2

HST_BASE-3.o: HST_BASE.h HST_BASE.cpp global.h global.cpp global-3.o monitor.h monitor.cpp monitor.o
	$(CXX) $(CFLAGS) -o HST_BASE-3.o -c global-3.o monitor.o HST_BASE.cpp $(LIBS) $(MEM) -DDIM_V=3

FRT-3: FRT.cpp global.h global.cpp global-3.o HST_BASE.h HST_BASE.cpp HST_BASE-3.o monitor.h monitor.cpp monitor.o
	$(CXX) $(CFLAGS) -o FRT-3 global-3.o monitor.o HST_BASE-3.o FRT.cpp $(LIBS) $(MEM) -DDIM_V=3

HST_BASE-4.o: HST_BASE.h HST_BASE.cpp global.h global.cpp global-4.o monitor.h monitor.cpp monitor.o
	$(CXX) $(CFLAGS) -o HST_BASE-4.o -c global-4.o monitor.o HST_BASE.cpp $(LIBS) $(MEM) -DDIM_V=4

FRT-4: FRT.cpp global.h global.cpp global-4.o HST_BASE.h HST_BASE.cpp HST_BASE-4.o monitor.h monitor.cpp monitor.o
	$(CXX) $(CFLAGS) -o FRT-4 global-4.o monitor.o HST_BASE-4.o FRT.cpp $(LIBS) $(MEM) -DDIM_V=4

HST_BASE-5.o: HST_BASE.h HST_BASE.cpp global.h global.cpp global-5.o monitor.h monitor.cpp monitor.o
	$(CXX) $(CFLAGS) -o HST_BASE-5.o -c global-5.o monitor.o HST_BASE.cpp $(LIBS) $(MEM) -DDIM_V=5

FRT-5: FRT.cpp global.h global.cpp global-5.o HST_BASE.h HST_BASE.cpp HST_BASE-5.o monitor.h monitor.cpp monitor.o
	$(CXX) $(CFLAGS) -o FRT-5 global-5.o monitor.o HST_BASE-5.o FRT.cpp $(LIBS) $(MEM) -DDIM_V=5

HST_BASE-10.o: HST_BASE.h HST_BASE.cpp global.h global.cpp global-10.o monitor.h monitor.cpp monitor.o
	$(CXX) $(CFLAGS) -o HST_BASE-10.o -c global-10.o monitor.o HST_BASE.cpp $(LIBS) $(MEM) -DDIM_V=10

FRT-10: FRT.cpp global.h global.cpp global-10.o HST_BASE.h HST_BASE.cpp HST_BASE-10.o monitor.h monitor.cpp monitor.o
	$(CXX) $(CFLAGS) -o FRT-10 global-10.o monitor.o HST_BASE-10.o FRT.cpp $(LIBS) $(MEM) -DDIM_V=10

HST_BASE-20.o: HST_BASE.h HST_BASE.cpp global.h global.cpp global-20.o monitor.h monitor.cpp monitor.o
	$(CXX) $(CFLAGS) -o HST_BASE-20.o -c global-20.o monitor.o HST_BASE.cpp $(LIBS) $(MEM) -DDIM_V=20

FRT-20: FRT.cpp global.h global.cpp global-20.o HST_BASE.h HST_BASE.cpp HST_BASE-20.o monitor.h monitor.cpp monitor.o
	$(CXX) $(CFLAGS) -o FRT-20 global-20.o monitor.o HST_BASE-20.o FRT.cpp $(LIBS) $(MEM) -DDIM_V=20

HST_BASE-100.o: HST_BASE.h HST_BASE.cpp global.h global.cpp global-100.o monitor.h monitor.cpp monitor.o
	$(CXX) $(CFLAGS) -o HST_BASE-100.o -c global-100.o monitor.o HST_BASE.cpp $(LIBS) $(MEM) -DDIM_V=100

FRT-100: FRT.cpp global.h global.cpp global-100.o HST_BASE.h HST_BASE.cpp HST_BASE-100.o monitor.h monitor.cpp monitor.o
	$(CXX) $(CFLAGS) -o FRT-100 global-100.o monitor.o HST_BASE-100.o FRT.cpp $(LIBS) $(MEM) -DDIM_V=100

HST_BARTAL-2.o: HST_BARTAL.h HST_BARTAL.cpp global.h global.cpp global-2.o monitor.h monitor.cpp monitor.o
	$(CXX) $(CFLAGS) -o HST_BARTAL-2.o -c global-2.o monitor.o HST_BARTAL.cpp $(LIBS) $(MEM) -DDIM_V=2
	
Bar98-2: Bar98.cpp global.h global.cpp global-2.o HST_BARTAL.h HST_BARTAL.cpp HST_BARTAL-2.o monitor.h monitor.cpp monitor.o
	$(CXX) $(CFLAGS) -o Bar98-2 global-2.o monitor.o HST_BARTAL-2.o Bar98.cpp $(LIBS) $(MEM) -DDIM_V=2

HST_BARTAL-3.o: HST_BARTAL.h HST_BARTAL.cpp global.h global.cpp global-3.o monitor.h monitor.cpp monitor.o
	$(CXX) $(CFLAGS) -o HST_BARTAL-3.o -c global-3.o monitor.o HST_BARTAL.cpp $(LIBS) $(MEM) -DDIM_V=3
	
Bar98-3: Bar98.cpp global.h global.cpp global-3.o HST_BARTAL.h HST_BARTAL.cpp HST_BARTAL-3.o monitor.h monitor.cpp monitor.o
	$(CXX) $(CFLAGS) -o Bar98-3 global-3.o monitor.o HST_BARTAL-3.o Bar98.cpp $(LIBS) $(MEM) -DDIM_V=3

HST_BARTAL-4.o: HST_BARTAL.h HST_BARTAL.cpp global.h global.cpp global-4.o monitor.h monitor.cpp monitor.o
	$(CXX) $(CFLAGS) -o HST_BARTAL-4.o -c global-4.o monitor.o HST_BARTAL.cpp $(LIBS) $(MEM) -DDIM_V=4
	
Bar98-4: Bar98.cpp global.h global.cpp global-4.o HST_BARTAL.h HST_BARTAL.cpp HST_BARTAL-4.o monitor.h monitor.cpp monitor.o
	$(CXX) $(CFLAGS) -o Bar98-4 global-4.o monitor.o HST_BARTAL-4.o Bar98.cpp $(LIBS) $(MEM) -DDIM_V=4

HST_BARTAL-5.o: HST_BARTAL.h HST_BARTAL.cpp global.h global.cpp global-5.o monitor.h monitor.cpp monitor.o
	$(CXX) $(CFLAGS) -o HST_BARTAL-5.o -c global-5.o monitor.o HST_BARTAL.cpp $(LIBS) $(MEM) -DDIM_V=5
	
Bar98-5: Bar98.cpp global.h global.cpp global-5.o HST_BARTAL.h HST_BARTAL.cpp HST_BARTAL-5.o monitor.h monitor.cpp monitor.o
	$(CXX) $(CFLAGS) -o Bar98-5 global-5.o monitor.o HST_BARTAL-5.o Bar98.cpp $(LIBS) $(MEM) -DDIM_V=5

HST_BARTAL-10.o: HST_BARTAL.h HST_BARTAL.cpp global.h global.cpp global-10.o monitor.h monitor.cpp monitor.o
	$(CXX) $(CFLAGS) -o HST_BARTAL-10.o -c global-10.o monitor.o HST_BARTAL.cpp $(LIBS) $(MEM) -DDIM_V=10
	
Bar98-10: Bar98.cpp global.h global.cpp global-10.o HST_BARTAL.h HST_BARTAL.cpp HST_BARTAL-10.o monitor.h monitor.cpp monitor.o
	$(CXX) $(CFLAGS) -o Bar98-10 global-10.o monitor.o HST_BARTAL-10.o Bar98.cpp $(LIBS) $(MEM) -DDIM_V=10

HST_BARTAL-20.o: HST_BARTAL.h HST_BARTAL.cpp global.h global.cpp global-20.o monitor.h monitor.cpp monitor.o
	$(CXX) $(CFLAGS) -o HST_BARTAL-20.o -c global-20.o monitor.o HST_BARTAL.cpp $(LIBS) $(MEM) -DDIM_V=20
	
Bar98-20: Bar98.cpp global.h global.cpp global-20.o HST_BARTAL.h HST_BARTAL.cpp HST_BARTAL-20.o monitor.h monitor.cpp monitor.o
	$(CXX) $(CFLAGS) -o Bar98-20 global-20.o monitor.o HST_BARTAL-20.o Bar98.cpp $(LIBS) $(MEM) -DDIM_V=20

HST_BARTAL-100.o: HST_BARTAL.h HST_BARTAL.cpp global.h global.cpp global-100.o monitor.h monitor.cpp monitor.o
	$(CXX) $(CFLAGS) -o HST_BARTAL-100.o -c global-100.o monitor.o HST_BARTAL.cpp $(LIBS) $(MEM) -DDIM_V=100
	
Bar98-100: Bar98.cpp global.h global.cpp global-100.o HST_BARTAL.h HST_BARTAL.cpp HST_BARTAL-100.o monitor.h monitor.cpp monitor.o
	$(CXX) $(CFLAGS) -o Bar98-100 global-100.o monitor.o HST_BARTAL-100.o Bar98.cpp $(LIBS) $(MEM) -DDIM_V=100

Bar96-2: Bar96.cpp global.h global.cpp global-2.o HST_BARTAL.h HST_BARTAL.cpp HST_BARTAL-2.o monitor.h monitor.cpp monitor.o
	$(CXX) $(CFLAGS) -o Bar96-2 global-2.o monitor.o HST_BARTAL-2.o Bar96.cpp $(LIBS) $(MEM) -DDIM_V=2

Bar96-3: Bar96.cpp global.h global.cpp global-3.o HST_BARTAL.h HST_BARTAL.cpp HST_BARTAL-3.o monitor.h monitor.cpp monitor.o
	$(CXX) $(CFLAGS) -o Bar96-3 global-3.o monitor.o HST_BARTAL-3.o Bar96.cpp $(LIBS) $(MEM) -DDIM_V=3

Bar96-4: Bar96.cpp global.h global.cpp global-4.o HST_BARTAL.h HST_BARTAL.cpp HST_BARTAL-4.o monitor.h monitor.cpp monitor.o
	$(CXX) $(CFLAGS) -o Bar96-4 global-4.o monitor.o HST_BARTAL-4.o Bar96.cpp $(LIBS) $(MEM) -DDIM_V=4

Bar96-5: Bar96.cpp global.h global.cpp global-5.o HST_BARTAL.h HST_BARTAL.cpp HST_BARTAL-5.o monitor.h monitor.cpp monitor.o
	$(CXX) $(CFLAGS) -o Bar96-5 global-5.o monitor.o HST_BARTAL-5.o Bar96.cpp $(LIBS) $(MEM) -DDIM_V=5

Bar96-10: Bar96.cpp global.h global.cpp global-10.o HST_BARTAL.h HST_BARTAL.cpp HST_BARTAL-10.o monitor.h monitor.cpp monitor.o
	$(CXX) $(CFLAGS) -o Bar96-10 global-10.o monitor.o HST_BARTAL-10.o Bar96.cpp $(LIBS) $(MEM) -DDIM_V=10

Bar96-20: Bar96.cpp global.h global.cpp global-20.o HST_BARTAL.h HST_BARTAL.cpp HST_BARTAL-20.o monitor.h monitor.cpp monitor.o
	$(CXX) $(CFLAGS) -o Bar96-20 global-20.o monitor.o HST_BARTAL-20.o Bar96.cpp $(LIBS) $(MEM) -DDIM_V=20

Bar96-100: Bar96.cpp global.h global.cpp global-100.o HST_BARTAL.h HST_BARTAL.cpp HST_BARTAL-100.o monitor.h monitor.cpp monitor.o
	$(CXX) $(CFLAGS) -o Bar96-100 global-100.o monitor.o HST_BARTAL-100.o Bar96.cpp $(LIBS) $(MEM) -DDIM_V=100

HST_GDCO-2.o: HST_GDCO.h HST_GDCO.cpp global.h global.cpp global-2.o monitor.h monitor.cpp monitor.o 
	$(CXX) $(CFLAGS) -o HST_GDCO-2.o -c global-2.o monitor.o HST_GDCO.cpp $(LIBS) $(BOOSTLIB) $(MEM) -DDIM_V=2
			
DCnn-2: DCnn.cpp global.h global.cpp global-2.o HST_GDCO.h HST_GDCO.cpp HST_GDCO-2.o monitor.h monitor.cpp monitor.o 
	$(CXX) $(CFLAGS) -o DCnn-2 global-2.o monitor.o HST_GDCO-2.o DCnn.cpp $(LIBS) $(BOOSTLIB) $(MEM) -DDIM_V=2

HST_GDCO-3.o: HST_GDCO.h HST_GDCO.cpp global.h global.cpp global-3.o monitor.h monitor.cpp monitor.o 
	$(CXX) $(CFLAGS) -o HST_GDCO-3.o -c global-3.o monitor.o HST_GDCO.cpp $(LIBS) $(BOOSTLIB) $(MEM) -DDIM_V=3
			
DCnn-3: DCnn.cpp global.h global.cpp global-3.o HST_GDCO.h HST_GDCO.cpp HST_GDCO-3.o monitor.h monitor.cpp monitor.o 
	$(CXX) $(CFLAGS) -o DCnn-3 global-3.o monitor.o HST_GDCO-3.o DCnn.cpp $(LIBS) $(BOOSTLIB) $(MEM) -DDIM_V=3

HST_GDCO-4.o: HST_GDCO.h HST_GDCO.cpp global.h global.cpp global-4.o monitor.h monitor.cpp monitor.o 
	$(CXX) $(CFLAGS) -o HST_GDCO-4.o -c global-4.o monitor.o HST_GDCO.cpp $(LIBS) $(BOOSTLIB) $(MEM) -DDIM_V=4
			
DCnn-4: DCnn.cpp global.h global.cpp global-4.o HST_GDCO.h HST_GDCO.cpp HST_GDCO-4.o monitor.h monitor.cpp monitor.o 
	$(CXX) $(CFLAGS) -o DCnn-4 global-4.o monitor.o HST_GDCO-4.o DCnn.cpp $(LIBS) $(BOOSTLIB) $(MEM) -DDIM_V=4

HST_GDCO-5.o: HST_GDCO.h HST_GDCO.cpp global.h global.cpp global-5.o monitor.h monitor.cpp monitor.o 
	$(CXX) $(CFLAGS) -o HST_GDCO-5.o -c global-5.o monitor.o HST_GDCO.cpp $(LIBS) $(BOOSTLIB) $(MEM) -DDIM_V=5
			
DCnn-5: DCnn.cpp global.h global.cpp global-5.o HST_GDCO.h HST_GDCO.cpp HST_GDCO-5.o monitor.h monitor.cpp monitor.o 
	$(CXX) $(CFLAGS) -o DCnn-5 global-5.o monitor.o HST_GDCO-5.o DCnn.cpp $(LIBS) $(BOOSTLIB) $(MEM) -DDIM_V=5

HST_GDCO-10.o: HST_GDCO.h HST_GDCO.cpp global.h global.cpp global-10.o monitor.h monitor.cpp monitor.o 
	$(CXX) $(CFLAGS) -o HST_GDCO-10.o -c global-10.o monitor.o HST_GDCO.cpp $(LIBS) $(BOOSTLIB) $(MEM) -DDIM_V=10
			
DCnn-10: DCnn.cpp global.h global.cpp global-10.o HST_GDCO.h HST_GDCO.cpp HST_GDCO-10.o monitor.h monitor.cpp monitor.o 
	$(CXX) $(CFLAGS) -o DCnn-10 global-10.o monitor.o HST_GDCO-10.o DCnn.cpp $(LIBS) $(BOOSTLIB) $(MEM) -DDIM_V=10

HST_GDCO-20.o: HST_GDCO.h HST_GDCO.cpp global.h global.cpp global-20.o monitor.h monitor.cpp monitor.o 
	$(CXX) $(CFLAGS) -o HST_GDCO-20.o -c global-20.o monitor.o HST_GDCO.cpp $(LIBS) $(BOOSTLIB) $(MEM) -DDIM_V=20
			
DCnn-20: DCnn.cpp global.h global.cpp global-20.o HST_GDCO.h HST_GDCO.cpp HST_GDCO-20.o monitor.h monitor.cpp monitor.o 
	$(CXX) $(CFLAGS) -o DCnn-20 global-20.o monitor.o HST_GDCO-20.o DCnn.cpp $(LIBS) $(BOOSTLIB) $(MEM) -DDIM_V=20

HST_GDCO-100.o: HST_GDCO.h HST_GDCO.cpp global.h global.cpp global-100.o monitor.h monitor.cpp monitor.o 
	$(CXX) $(CFLAGS) -o HST_GDCO-100.o -c global-100.o monitor.o HST_GDCO.cpp $(LIBS) $(BOOSTLIB) $(MEM) -DDIM_V=100
			
DCnn-100: DCnn.cpp global.h global.cpp global-100.o HST_GDCO.h HST_GDCO.cpp HST_GDCO-100.o monitor.h monitor.cpp monitor.o 
	$(CXX) $(CFLAGS) -o DCnn-100 global-100.o monitor.o HST_GDCO-100.o DCnn.cpp $(LIBS) $(BOOSTLIB) $(MEM) -DDIM_V=100

HST_PARS_Fast-2.o: HST_PARS_Fast.h HST_PARS_Fast.cpp global.h global.cpp global-2.o monitor.h monitor.cpp monitor.o 
	$(CXX) $(CFLAGS) -o HST_PARS_Fast-2.o -c global-2.o monitor.o HST_PARS_Fast.cpp $(LIBS) $(BOOSTLIB) $(MEM) -DDIM_V=2

DCsam-2: DCsam_Rtree.cpp global.h global.cpp global-2.o HST_PARS_Fast.h HST_PARS_Fast.cpp HST_PARS_Fast-2.o monitor.h monitor.cpp monitor.o 
	$(CXX) $(CFLAGS) -o DCsam-2 global-2.o monitor.o HST_PARS_Fast-2.o DCsam_Rtree.cpp $(LIBS) $(BOOSTLIB) $(MEM) -DDIM_V=2

HST_PARSA_Fast-3.o: HST_PARSA_Fast.h HST_PARSA_Fast.cpp global.h global.cpp global-3.o monitor.h monitor.cpp monitor.o 
	$(CXX) $(CFLAGS) -o HST_PARSA_Fast-3.o -c global-3.o monitor.o HST_PARSA_Fast.cpp $(LIBS) $(ANNLIB) $(MEM) -DDIM_V=3

DCsam-3: DCsam_ANN.cpp global.h global.cpp global-3.o HST_PARSA_Fast.h HST_PARSA_Fast.cpp HST_PARSA_Fast-3.o monitor.h monitor.cpp monitor.o 
	$(CXX) $(CFLAGS) -o DCsam-3 global-3.o monitor.o HST_PARSA_Fast-3.o DCsam_ANN.cpp $(LIBS) $(ANNLIB) $(MEM) -DDIM_V=3

HST_PARSA_Fast-4.o: HST_PARSA_Fast.h HST_PARSA_Fast.cpp global.h global.cpp global-4.o monitor.h monitor.cpp monitor.o 
	$(CXX) $(CFLAGS) -o HST_PARSA_Fast-4.o -c global-4.o monitor.o HST_PARSA_Fast.cpp $(LIBS) $(ANNLIB) $(MEM) -DDIM_V=4

DCsam-4: DCsam_ANN.cpp global.h global.cpp global-4.o HST_PARSA_Fast.h HST_PARSA_Fast.cpp HST_PARSA_Fast-4.o monitor.h monitor.cpp monitor.o 
	$(CXX) $(CFLAGS) -o DCsam-4 global-4.o monitor.o HST_PARSA_Fast-4.o DCsam_ANN.cpp $(LIBS) $(ANNLIB) $(MEM) -DDIM_V=4

HST_PARSA_Fast-5.o: HST_PARSA_Fast.h HST_PARSA_Fast.cpp global.h global.cpp global-5.o monitor.h monitor.cpp monitor.o 
	$(CXX) $(CFLAGS) -o HST_PARSA_Fast-5.o -c global-5.o monitor.o HST_PARSA_Fast.cpp $(LIBS) $(ANNLIB) $(MEM) -DDIM_V=5

DCsam-5: DCsam_ANN.cpp global.h global.cpp global-5.o HST_PARSA_Fast.h HST_PARSA_Fast.cpp HST_PARSA_Fast-5.o monitor.h monitor.cpp monitor.o 
	$(CXX) $(CFLAGS) -o DCsam-5 global-5.o monitor.o HST_PARSA_Fast-5.o DCsam_ANN.cpp $(LIBS) $(ANNLIB) $(MEM) -DDIM_V=5

HST_PARSA_Fast-10.o: HST_PARSA_Fast.h HST_PARSA_Fast.cpp global.h global.cpp global-10.o monitor.h monitor.cpp monitor.o 
	$(CXX) $(CFLAGS) -o HST_PARSA_Fast-10.o -c global-10.o monitor.o HST_PARSA_Fast.cpp $(LIBS) $(ANNLIB) $(MEM) -DDIM_V=10

DCsam-10: DCsam_ANN.cpp global.h global.cpp global-10.o HST_PARSA_Fast.h HST_PARSA_Fast.cpp HST_PARSA_Fast-10.o monitor.h monitor.cpp monitor.o 
	$(CXX) $(CFLAGS) -o DCsam-10 global-10.o monitor.o HST_PARSA_Fast-10.o DCsam_ANN.cpp $(LIBS) $(ANNLIB) $(MEM) -DDIM_V=10

HST_PARSA_Fast-20.o: HST_PARSA_Fast.h HST_PARSA_Fast.cpp global.h global.cpp global-20.o monitor.h monitor.cpp monitor.o 
	$(CXX) $(CFLAGS) -o HST_PARSA_Fast-20.o -c global-20.o monitor.o HST_PARSA_Fast.cpp $(LIBS) $(ANNLIB) $(MEM) -DDIM_V=20

DCsam-20: DCsam_ANN.cpp global.h global.cpp global-20.o HST_PARSA_Fast.h HST_PARSA_Fast.cpp HST_PARSA_Fast-20.o monitor.h monitor.cpp monitor.o 
	$(CXX) $(CFLAGS) -o DCsam-20 global-20.o monitor.o HST_PARSA_Fast-20.o DCsam_ANN.cpp $(LIBS) $(ANNLIB) $(MEM) -DDIM_V=20

random.o: random.h random.cc def.h
	$(CXX) $(CFLAGS) -o random.o -c random.cc $(LIBS)
			
pri_queue.o: pri_queue.h pri_queue.cc def.h
	$(CXX) $(CFLAGS) -o pri_queue.o -c pri_queue.cc $(LIBS)
		
util.o: util.h util.cc def.h random.h random.cc random.o pri_queue.h pri_queue.cc pri_queue.o
	$(CXX) $(CFLAGS) -o util.o -c random.o pri_queue.o util.cc $(LIBS)
		
HST_PARSA_LSH-100.o: HST_PARSA_LSH.h HST_PARSA_LSH.cpp global.h global.cpp global-100.o monitor.h monitor.cpp monitor.o qalsh.h util.h util.cc util.o def.h random.h random.cc random.o pri_queue.h pri_queue.cc pri_queue.o
	$(CXX) $(CFLAGS) -o HST_PARSA_LSH-100.o -c global-100.o monitor.o util.o -c random.o pri_queue.o HST_PARSA_LSH.cpp $(LIBS) $(MEM) -DDIM_V=100

DCsam-100: DCsam_QALSH.cpp global.h global.cpp global-100.o HST_PARSA_LSH.h HST_PARSA_LSH.cpp HST_PARSA_LSH-100.o monitor.h monitor.cpp monitor.o qalsh.h util.h util.cc util.o def.h random.h random.cc random.o pri_queue.h pri_queue.cc pri_queue.o
	$(CXX) $(CFLAGS) -o DCsam-100 global-100.o monitor.o util.o random.o pri_queue.o HST_PARSA_LSH-100.o DCsam_QALSH.cpp $(LIBS) $(MEM) -DDIM_V=100

.PHONY: clean
clean:
		-@rm *.o *.gcno *~ 2> /dev/null || true
		-@rm FRT-2 FRT-3 FRT-4 FRT-5 FRT-10 FRT-20 FRT-100  2> /dev/null || true
		-@rm Bar98-2 Bar98-3 Bar98-4 Bar98-5 Bar98-10 Bar98-20 Bar98-100  2> /dev/null || true
		-@rm Bar96-2 Bar96-3 Bar96-4 Bar96-5 Bar96-10 Bar96-20 Bar96-100  2> /dev/null || true
		-@rm DCnn-2 DCnn-3 DCnn-4 DCnn-5 DCnn-10 DCnn-20 DCnn-100  2> /dev/null || true
		-@rm DCsam-2 DCsam-3 DCsam-4 DCsam-5 DCsam-10 DCsam-20 DCsam-100  2> /dev/null || true
