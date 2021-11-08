#ifndef HST_DC_H
#define HST_DC_H


#include <bits/stdc++.h>
using namespace std;
#include "global_DC.h"
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/index/rtree.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/foreach.hpp>
#ifdef WATCH_MEM
#include "monitor.h"
#endif

//#define LOCAL_DEBUG

struct Node_t {
	int nid;		// the idx of the node
	int nomIdx; 	// the id of norminator point
	int lev;		// the level of the node before compression
	Node_t *far;	// the parent node
	float wei;		// the weight of the edge from its parent
	int cid;		// the position of this node in its parent's child list.
	vector<Node_t*> child;
	
	Node_t() {}
	Node_t(int nid_, int nomIdx_, int lev_, Node_t* far_, float wei_) {
		nid = nid_;
		nomIdx = nomIdx_;
		lev = lev_;
		far = far_;
		wei = wei_;
		cid = 0;
	}
};

struct partition_t {
	int lev;
	Node_t* ptr;
	vector<int> ids;
};

extern const int MAX_HEIGHT;
static const int RTREE_NODE_N = 20;
static const int PAGE_SIZE_N = 1;
namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;
typedef bg::model::point<double, 2, bg::cs::cartesian> rtree_point;
typedef std::pair<rtree_point, unsigned> value;
typedef bgi::rtree<value, bgi::quadratic<RTREE_NODE_N> > rtree_t;

struct Tree_t {
	int nV, dV;
	int H;
	int node_n, merge_n;
	double logAlpha;
	double dmax;
	double beta;
	double alpha;
	double distor, distor_, cdistor;
	double* expks;
	double* sumks;
	int num_partitions_LB;
	int num_points_LB;
	int num_points_UB;
	int num_points_SAMPLE;
	double res;
	double maxDistor;
	Node_t* rt;
	vector<Node_t*> leaves;
	vector<int> pi;
	vector<int> reverse_pi;
	vector<location_t> V;
	vector<int> ids;
	vector<int> mark;
	vector<bool> deleted;  // mark whether the points are deleted
	vector<vector<pair<int,double> > > nns; // store the nearest neighbours for each point
	vector<value> dataset;
	vector<rtree_point> RT_POINTS;
	vector<vector<double> > distorUB;
	
	int RMQ_top;
	vector<Node_t*> RMQ_V;
	vector<int> RMQ_D;
	vector<int> RMQ_beg;
	vector<vector<int> > RMQ_dp;

public:
	void initLocation(int nV, vector<location_t>& V, vector<int>& ids, double alpha);  // init the parameter of the tree 
	void constructHST(clock_t startClock); // startClock is used to control whether the algorithm is TLE.
	void freeMemory(); // free the memory usage of the structure
	pair<double,double> getDistortion(); // get the distortion of the tree						
	void initNeighbours(int treeId, Tree_t& oth); // init the neighbors of the points in the other tree.
	void updateNeighbours(int treeId, Tree_t& oth); // update the neighbors of the points in the other tree.
	void updateNeighbour(int vid, int treeId, Tree_t& oth); // update the neighbor of vid in the tree oth
	void loadLocation(int nV, vector<location_t>& V, vector<int>& ids, double alpha);
	double distOnHST(int u, int v); // query the distance on the tree
	pair<Node_t*, int> getLCA(int u, int v); // query the LCA of two boths
	void deleteOne(int idx);
	long long memoryCost();
	void calcDmax();
	pair<double,double> getDistortion_fast();
	
private:
	void initMemory(int n);
	void freeHST(Node_t*& rt);
	double distAtLevel(int level);
	double distOnHST(Node_t* u, Node_t* v);
	pair<Node_t*, int> getLCA(Node_t* u, Node_t* v);
	int levelOfLCA(int u, int v);
	int levelOfLCA(Node_t* u, Node_t* v);
	
	void initParameters(int flag, int m);
	double _calc_maxDistor_rtree0_UB(int centerId, partition_t& g, double UB);
	void _getParition(int& nid, int centerId, partition_t& g, vector<partition_t>& left);
	int _getCenter_maxDistor_complex(partition_t& g);
	int _getCenter_subspaceDistor_complex(partition_t& g);
	void _getSubspaces(int& nid, int centerId, partition_t& g, vector<partition_t>& left);
	double _calc_subspaceDistor_rtree0_UB2(int centerId, partition_t& g, double UB);
	
	double RMQ_distOnHST(Node_t* pu, Node_t* pv);
	Node_t* RMQ_query(int l, int r);
	void RMQ_init(int n);
	void RMQ_dfs(Node_t* rt, int d);
	void RMQ_resize();
	void RMQ_clear();
	
	inline double log2(double x) {
		return log(x) / logAlpha;
	}

	inline double pow2(int i) {
		return (i<0) ? (1.0/alpha) : expks[i];
	}
	
	inline int whichLevel(double l) {
		double radius = pow2(H)*beta;
		for (int lev=1; lev<=H+1; ++lev,radius/=alpha) {
			if (l > radius) return lev-1;
		}
		return 1;
	}
	
	inline int getLevel(double H, double dist) {
		if (dist < 1.0) return H+1;
		
		int k = ceil(log2(dist/beta));
		if (k>0 && expks[k-1]*beta>=dist) 
			--k;
		
		return H+1-k;
	}
	
	int countTree(Node_t* rt) {
		int ret = 0;
		if (rt != NULL) {
			ret = 1;
			for (int i=0; i<rt->child.size(); ++i)
				ret += countTree(rt->child[i]);
		}
		return ret;
	}
	
	inline void addChild(Node_t* far, Node_t* u) {
		if (far != NULL) {
			u->far = far;
			u->cid = far->child.size();
			far->child.push_back(u);
		}
	}

	inline void mergeChild(Node_t* far, Node_t* u) {
		if (far==NULL || far==rt || far->child.size()!=1) 
			return ;
		
		Node_t* gfar = far->far;
		int gcid = far->cid;
		
		u->far = far->far;
		u->cid = gcid;
		u->wei += far->wei;
		
		gfar->child[gcid] = u;
		delete far;
		merge_n++;
	}
};

#endif
