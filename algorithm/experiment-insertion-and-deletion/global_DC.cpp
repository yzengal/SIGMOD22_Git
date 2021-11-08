#include "global_DC.h"

int nV = 0;
const double speed = 1.0;
const double EPS = 1e-3;
const double INF = 1e20;
double usedTime = 0.0;
int NUM_DIM = 2;
long long usedMemory = 0;
double timeLimit = 12 * 60 * 60; // 12 hours

int dcmp(double x) {
	if (fabs(x) < EPS)
		return 0;
	return x>0 ? 1:-1;
}

double dist(const vector<location_t>& V, int x, int y) { 
	if (x == y) return 0.0;
	
	const location_t& a = V[x];
	const location_t& b = V[y];
	return sqrt(1.0*(a.x - b.x)*(a.x - b.x) + 1.0*(a.y - b.y)*(a.y - b.y));
} 

double dist(location_t *V, int x, int y) { 
	location_t& a = V[x];
	location_t& b = V[y];
	return sqrt(1.0*(a.x - b.x)*(a.x - b.x) + 1.0*(a.y - b.y)*(a.y - b.y));
} 

double dist(location_t &a, location_t &b) {
	return sqrt(1.0*(a.x - b.x)*(a.x - b.x) + 1.0*(a.y - b.y)*(a.y - b.y));
}
