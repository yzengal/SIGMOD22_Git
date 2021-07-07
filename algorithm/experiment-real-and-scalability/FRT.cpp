#include <bits/stdc++.h>
using namespace std;
#include "global.h"
#include "HST_BASE.h"
#ifdef WATCH_MEM
#include "monitor.h"
#endif

string srcFileName("data_0.txt"), dataFileName, desFileName;

void FRT() {
	initVariable(dataFileName);
	initLocation(srcFileName);

	clock_t t0, t1;
	t0 = clock();
	
	_constructHST_ICDE_DPO(false, clock());
	
	t1 = clock();
	pair<double,double> tmp = getDistortion_fast();
	double timeCost = 1.0 * (t1-t0) / CLOCKS_PER_SEC;
	
	printf("FRT: maxDistor=%.4lf, avgDistor=%.4lf, time=%.4lf, memory=%.4lf\n", tmp.first, tmp.second, timeCost, usedMemory/1024.0);
	fflush(stdout);
}

int main(int argc, char **argv) {
	if (argc > 1)
		srcFileName = string(argv[1]);
	if (argc > 2)
		dataFileName = string(argv[2]);
	if (argc > 3) {
		desFileName = string(argv[3]);
		freopen(desFileName.c_str(), "w", stdout);
	}
	
	FRT();
	freeMemory_HST();
	
	return 0;
}
