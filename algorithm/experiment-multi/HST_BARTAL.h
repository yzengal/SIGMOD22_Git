#ifndef HST_BARTAL_H
#define HST_BARTAL_H

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
extern double beta;
extern double* expks;
extern double* sumks;
extern double dmax;


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
	
    Node_t &operator=(const Node_t &oth)
    {
        if (this != &oth)
        {
            nid = oth.nid;
            nomIdx = oth.nomIdx;
            lev = oth.lev;
            far = oth.far;
            cid = oth.cid;
            wei = oth.wei;
            child = oth.child;
        }

        return *this;
    }
};

struct cluster_t {
	int id, far;
	int centerId, lev;
	double diam;
	vector<int> VG;
	
	cluster_t() {
		id = far = centerId = lev = 0;
		diam = 0;
	}

    cluster_t &operator=(const cluster_t &oth)
    {
        if (this != &oth)
        {
            id = oth.id;
            far = oth.far;
			centerId = oth.centerId;
			lev = oth.lev;
			diam = oth.diam;
            VG = oth.VG;
        }

        return *this;
    }
	
	int calc_n_VG() {
		return VG.size();
	}
	
	double calc_diameter() {
		int sz = VG.size();
		if (sz <= 1) return 0.0;
		
		double ret = 0.0;
		for (int i=0; i<sz; ++i) {
			for (int j=i+1; j<sz; ++j) {
				ret = max(ret, dist(V, VG[i], VG[j]));
			}
		}
		ret += 1.0;
		diam = ret;
		
		return ret;
	}
};

extern Node_t* rt;
extern Node_t** leaves;

void dumpNode(Node_t* p);
void dumpTree(Node_t* rt);
int getTreeHeight(Node_t* rt);
void addChild(Node_t* far, Node_t* u);
void initVariable(string &fileName);
void initLocation(string &fileName);
void initMemory_HST(int n);
void freeMemory_HST();
void freeHST(Node_t* &rt);
pair<double,double> getDistortion();
void calcDmax();
double distOnHST(int u, int v);
double distOnHST(Node_t* u, Node_t* v);
int levelOfLCA(int u, int v);
int levelOfLCA(Node_t* u, Node_t* v);

double sampleRadius_focs96(int n, double r);
void _construct_probabilistic_partition(vector<cluster_t>& C, int clusterId, int& n_cluster);
void _constructHST_focs96(bool load, clock_t startClock);

double calc_g_stoc98(cluster_t& c);
double sampleRadius_stoc98(cluster_t& c);
void _constructHST_stoc98(bool load, clock_t startClock);
void _construct_HPM(cluster_t& _c, cluster_t& Hg1, cluster_t& Hg2, int& n_cluster);
void _construct_HPM(cluster_t& _c, double L, queue<pair<cluster_t,double> >& Q, int& n_cluster);


#endif
