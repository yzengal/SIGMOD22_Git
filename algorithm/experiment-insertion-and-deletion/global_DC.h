#ifndef GLOBAL_DC_H
#define GLOBAL_DC_H

#include <bits/stdc++.h>
using namespace std;

// #define LOCAL_DEBUG

// typedef unorderd_map<int,int> hii;
typedef long long LL;
typedef pair<int,int> pii;
typedef pair<double,double> pdd;
typedef pair<int,double> pid;
typedef pair<double,int> pdi;
const int MAX_DIM = 10;

struct location_t {
	double x, y;

	location_t(double _x=0.0, double _y=0.0):
		x(_x), y(_y) {}
		
	location_t& operator=(const location_t &loc) { 
        if (this != &loc) {
			this->x = loc.x;
			this->y = loc.y;
        }
        return *this;	
	}
};

extern int nV;
extern int NUM_DIM;
extern const double speed;
extern const double EPS;
extern const double INF;
extern double usedTime;
extern long long usedMemory;
extern double timeLimit;

int dcmp(double x);  
double dist(const vector<location_t>& V, int x, int y);
double dist(location_t *V, int x, int y);                                                                                                                       double dist(location_t *V, int x, int y);
double dist(location_t& a, location_t& b);

#endif
