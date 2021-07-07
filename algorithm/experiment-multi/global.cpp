#include "global.h"

int nV = 0;
location_t* V = NULL;
const double speed = 1.0;
const double EPS = 1e-3;
const double INF = 1e20;
double usedTime = 0.0;
int usedMemory = 0;
double timeLimit = 12 * 60 * 60;

int dcmp(double x) {
	if (fabs(x) < EPS)
		return 0;
	return x>0 ? 1:-1;
}

double dist(location_t *V, int x, int y) {
	if (x == y) return 0.0;

	location_t& a = V[x];
	location_t& b = V[y];
	double ret = 0.0;

	for (int i=0; i<DIM_V; ++i) {
		ret += (a.x[i]-b.x[i]) * (a.x[i]-b.x[i]);
	}

	return sqrt(ret);
}

double dist(location_t &a, location_t &b) {
	double ret = 0.0;

	for (int i=0; i<DIM_V; ++i) {
		ret += (a.x[i]-b.x[i])*(a.x[i]-b.x[i]);
	}

	return sqrt(ret);
}

void freeLocation() {
	delete[] V;
}
