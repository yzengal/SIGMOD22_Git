#ifndef HST_BASE_H
#define HST_BASE_H

#include <bits/stdc++.h>
using namespace std;
#include "global.h"
#ifdef WATCH_MEM
#include "monitor.h"
#endif

extern const int MAX_SAMPLE;
extern const int MAX_HEIGHT;
extern int H;
extern int alpha;
extern double* expks;
extern double* sumks;
extern int* pi;
extern int* reverse_pi;
extern double dmax;
extern double beta;


struct Node_t{
	int nid;		// the idx of the node
	int nomIdx; 	// the id of center point
	int lev;		// the level of the node before compression
	Node_t *far;	// the parent node
	int cid;		// the position of this node in its parent's child list.
	float wei;		// the weight of the edge from its parent 
	vector<Node_t*> child;
	
	Node_t(int nid_ = 0, int nomIdx_=0, int lev_=0, Node_t* far_=NULL, float wei_=0) {
		nid = nid_;
		nomIdx = nomIdx_;
		lev = lev_;
		far = far_;
		wei = wei_;
		cid = 0;
	}
};

extern Node_t* rt;
extern Node_t** leaves;

int countTree(Node_t* rt);
void dumpNode(Node_t* p);
void dumpTree(Node_t* rt);
void mergeChild(Node_t* far, Node_t* u);
void addChild(Node_t* far, Node_t* u);
void initLocation(string &fileName);
void initVariable(string &fileName);
void initMemory_HST(int n);
void freeMemory_HST();
void freeHST(Node_t* &rt);
pair<double,double> getDistortion();
pair<double,double> getDistortion_fast();
void randomization();
void calcDmax();
void calcDmaxPrune();
int whichLevel(double l);
double distAtLevel(int level);
double distOnHST(int u, int v);
double distOnHST(Node_t* u, Node_t* v);
pair<Node_t*, int> getLCA(int u, int v);
pair<Node_t*, int> getLCA(Node_t* u, Node_t* v);
int levelOfLCA(int u, int v);
int levelOfLCA(Node_t* u, Node_t* v);

void _constructHST_ICDE_DPO(bool load, clock_t startClock);

#endif
