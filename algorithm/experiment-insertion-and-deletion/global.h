#ifndef GLOBAL_H
#define GLOBAL_H

#include <bits/stdc++.h>
using namespace std;

// #define LOCAL_DEBUG

// typedef unorderd_map<int,int> hii;
typedef long long LL;
typedef pair<int,int> pii;
typedef pair<double,double> pdd;
typedef pair<int,double> pid;
typedef pair<double,int> pdi;

const int MAX_DIM = 2;

struct location_t{
    vector<double> x;
	location_t() {x.resize(MAX_DIM, 0);}
	// location_t& operator=(const location_t &loc) { 
        // if (this != &loc) {
			// for (int i=0; i<MAX_DIM; ++i)
				// this->x[i] = loc.x[i];
        // }
        // return *this;	
	// }
};

extern int nV;
extern const double speed;
extern const double EPS;
extern const double INF;
extern double usedTime;
extern long long usedMemory;
extern double timeLimit;

int dcmp(double x);
double dist(const vector<location_t>& V, int x, int y);
double dist(location_t *V, int x, int y);
double dist(location_t& a, location_t& b);

#endif
