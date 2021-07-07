#ifndef GLOBAL_H
#define GLOBAL_H

#include <bits/stdc++.h>
using namespace std;

// #define DIM_V 4

typedef long long LL;
typedef pair<int,int> pii;
typedef pair<double,double> pdd;
typedef pair<int,double> pid;
typedef pair<double,int> pdi;
const int MAX_DIM = 10;

struct location_t {
    double x[DIM_V];
};

extern int nV;
extern location_t* V;
extern const double speed;
extern const double EPS;
extern const double INF;
extern double usedTime;
extern int usedMemory;
extern double timeLimit;

int dcmp(double x);  
double dist(location_t *V, int x, int y);                                                                                                                       double dist(location_t *V, int x, int y);
double dist(location_t& a, location_t& b);
void freeLocation();


#endif
