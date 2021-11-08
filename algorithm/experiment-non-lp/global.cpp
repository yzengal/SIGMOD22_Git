#include "global.h"

int nV = 0;
location_t* V = NULL;
const double speed = 1.0;
const double EPS = 1e-3;
const double INF = 1e20;
double usedTime = 0.0;
int usedMemory = 0;
double timeLimit = 12 * 60 * 60;
double normalDistance = 1e7;

int dcmp(double x) {
	if (fabs(x) < EPS)
		return 0;
	return x>0 ? 1:-1;
}

double dist(location_t &a, location_t &b) {
	
#ifdef CHISQUARE
	return distChiSquare(a, b) * normalDistance;
#endif 

#ifdef HELLINGER
	return distHellinger(a, b) * normalDistance;
#else
	return distEuclid(a, b);
#endif 
}

double dist(location_t *V, int x, int y) {
	if (x == y) return 0.0;
	
#ifdef CHISQUARE
	return distChiSquare(V, x, y) * normalDistance;
#endif 

#ifdef HELLINGER
	return distHellinger(V, x, y) * normalDistance;
#else
	return distEuclid(V, x, y);
#endif 
}

void freeLocation() {
	delete[] V;
}

inline double distEuclid(location_t *V, int x, int y) {
	location_t& a = V[x];
	location_t& b = V[y];
	double ret = 0.0;

	for (int i=0; i<DIM_V; ++i) {
		ret += (a.x[i]-b.x[i]) * (a.x[i]-b.x[i]);
	}

	return sqrt(ret);
}

inline double distEuclid(location_t &a, location_t &b) {
	double ret = 0.0;

	for (int i=0; i<DIM_V; ++i) {
		ret += (a.x[i]-b.x[i])*(a.x[i]-b.x[i]);
	}

	return sqrt(ret);
}

inline double distHellinger(location_t *V, int x, int y) {
	location_t& a = V[x];
	location_t& b = V[y];
	static const int SQRT2 = sqrt(2.0);
	double ret = 0.0, tmp;

	for (int i=0; i<DIM_V; ++i) {
		tmp = sqrt(a.x[i]) - sqrt(b.x[i]);
		ret += tmp*tmp;
	}

	return sqrt(ret) / SQRT2;
}

inline double distHellinger(location_t &a, location_t &b) {
	static const int SQRT2 = sqrt(2.0);
	double ret = 0.0, tmp;

	for (int i=0; i<DIM_V; ++i) {
		tmp = sqrt(a.x[i]) - sqrt(b.x[i]);
		ret += tmp*tmp;
	}

	return sqrt(ret) / SQRT2;
}

inline double distChiSquare(location_t *V, int x, int y) {

	location_t& a = V[x];
	location_t& b = V[y];
	double ret = 0.0;

	for (int i=0; i<DIM_V; ++i) {
		if (a.x[i]+b.x[i] != 0)
			ret += (a.x[i]-b.x[i])*(a.x[i]-b.x[i])/(a.x[i]+b.x[i]);
	}

	return sqrt(ret * 0.5);
}

inline double distChiSquare(location_t &a, location_t &b) {
	double ret = 0.0;

	for (int i=0; i<DIM_V; ++i) {
		if (a.x[i]+b.x[i] != 0)
			ret += (a.x[i]-b.x[i])*(a.x[i]-b.x[i])/(a.x[i]+b.x[i]);
	}

	return sqrt(ret * 0.5);
}
