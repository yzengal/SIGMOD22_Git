#include "HST_PARSA_Fast.h"
#include "global.h"
#include <flann/flann.hpp>
#ifdef WATCH_MEM
#include "monitor.h"
#endif

const int MAX_HEIGHT = 55;
const int threshold_nV = 80;
int H = 0;
int alpha = 2;
double logAlpha = 1.0;
int* pi = NULL;
int* reverse_pi = NULL;
int* mark = NULL;
double dmax = -1.0;
double* expks = NULL;
double* sumks = NULL;
Node_t* rt = NULL;
Node_t** leaves = NULL;
double beta = 1.0;
static vector<vector<double> > distorUB;
int num_partitions_LB;
int num_points_LB;
int num_points_UB;
int num_points_SAMPLE;
static double res = 1.0;
double maxDistor = 1.0;
static double eps_ann = (DIM_V<=10) ? 4.0 : 4.5;
static int FLANN_CHECK = 128;
double* dataset = NULL;
static flann::Matrix<int> NNindices(new int[1], 1, 1);
static flann::Matrix<double> NNdists(new double[1], 1, 1);
static flann::Matrix<double> queryPoint(new double[DIM_V], 1, DIM_V);
	
inline double log2(double x) {
	return log(x) / logAlpha;
}

inline double pow2(int i) { 
	return (i<0) ? (1.0/alpha) : expks[i];
}

int whichLevel(double l) {
	double radius = pow2(H)*beta;
	for (int lev=1; lev<=H+1; ++lev,radius/=alpha) {
		if (l > radius) return lev-1;
	}
	return 1;
}

void dumpNode(Node_t* p) {
	for (int i=1; i<p->lev; ++i)
		putchar(' ');
	printf("nid=%d, wei=%.4lf, |child|=%d\n", p->nid, p->wei, (int)p->child.size());
}

void dumpTree(Node_t* rt) {
	if (rt != NULL) {
		dumpNode(rt);
		for (int i=0; i<rt->child.size(); ++i)
			dumpTree(rt->child[i]);
	}
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

inline int getLevel(double H, double dist) {
	if (dist < 1.0) return H+1;
	
	int k = ceil(log2(dist/beta));
	if (k>0 && expks[k-1]*beta>=dist) 
		--k;
	
	return H+1-k;
}

void initMemory_HST(int n) {
	nV = n;
	V = new location_t[nV];
	pi = new int[nV];
	reverse_pi = new int[nV];
	expks = new double[MAX_HEIGHT];
	expks[0] = 1.0;
	for (int i=1; i<MAX_HEIGHT; ++i)
		expks[i] = expks[i-1] * alpha;
	
	sumks = new double[MAX_HEIGHT+1];
	sumks[0] = 0.0;
	for (int i=1; i<=MAX_HEIGHT; ++i)
		sumks[i] = sumks[i-1] + expks[i];
	
	leaves = new Node_t*[nV];
	mark = new int[nV];
	
	logAlpha = log(alpha);
	distorUB.resize(nV);
	
	for (int i=0; i<nV; ++i) {
		distorUB[i].resize(MAX_HEIGHT);
		fill(distorUB[i].begin(), distorUB[i].end(), INF);
	}
	dataset = new double[nV*DIM_V];
}

void initVariable(string &fileName) {
	ifstream fin(fileName.c_str(), ios::in);   
    
	if (!fin.is_open()) {
		fprintf(stderr, "FILE %s IS INVALID.", fileName.c_str());
		exit(1);
	}
	
	fin >> beta >> alpha;
    fin.close();	
}

void freeMemory_HST() {
	freeLocation();
	delete[] pi;
	delete[] reverse_pi;
	delete[] expks;
	delete[] sumks;
    delete[] leaves;
	delete[] mark;
	freeHST(rt);
    delete[] NNindices.ptr();
    delete[] NNdists.ptr();
	delete[] queryPoint.ptr();
	delete[] dataset;
}

void freeHST(Node_t*& rt) {
	if (rt == NULL) return ;
	
	Node_t* chi;
	
	for (int i=0; i<rt->child.size(); ++i) {
		chi = rt->child[i];
		freeHST(chi);
	}
	
	delete rt;
	rt = NULL;
}

void initLocation(string &fileName) {
	ifstream fin(fileName.c_str(), ios::in);   
    
	if (!fin.is_open()) {
		fprintf(stderr, "FILE %s IS INVALID.", fileName.c_str());
		exit(1);
	}
	
	fin >> nV;
	initMemory_HST(nV);
    for (int i=0; i<nV; ++i) {
		for (int j=0; j<DIM_V; ++j)
			fin >> V[i].x[j];
    }
    fin.close();
}

void randomization() {
	for (int i=0; i<nV; ++i){
		pi[i] = i;
	}
	random_shuffle(pi, pi+nV);
	for (int i=0; i<nV; ++i){
		reverse_pi[pi[i]] = i;
	}
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
}

double distAtLevel(int lev) {
	if (lev >= H+1) return 0.0;
	if (lev <= 1) lev = 1;
	return sumks[H-lev+1] * beta * 2.0;

}

double distOnHST(int u, int v) {
	int level = levelOfLCA(u, v);
	return distAtLevel(level);
}

double distOnHST(Node_t* u, Node_t* v) {
	int level = levelOfLCA(u, v);
	return distAtLevel(level);
}

int levelOfLCA(int u, int v) {
	return levelOfLCA(leaves[u], leaves[v]);
}

int levelOfLCA(Node_t* u, Node_t* v) {
	if (u==NULL || v==NULL)
		return -1;
	
	while (u!=NULL && v!=NULL && u!=v) {
		if (u->lev > v->lev) {
			u = u->far;
		} else {
			v = v->far;
		}
	}
	
	if (u==NULL || v==NULL)
		return -1;
	
	return u->lev;
}

pair<Node_t*,int> getLCA(int u, int v) {
	return getLCA(leaves[u], leaves[v]);
}

pair<Node_t*,int> getLCA(Node_t* u, Node_t* v) {
	if (u==NULL || v==NULL)
            return make_pair(rt, -1);
	
	while (u!=NULL && v!=NULL && u!=v) {
		if (u->lev > v->lev) {
			u = u->far;
		} else {
			v = v->far;
		}
	}
	
	if (u==NULL || v==NULL)
		return make_pair(rt, -1);
	
	return make_pair(u, u->lev);
}

void calcDmax() {
	if (dmax >= 0.0) return ;

	dmax = 0.0;
	double dmax_ = 0.0;
	double tmp;
	for (int j=1; j<nV; ++j) {
		tmp = dist(V, pi[0], j);
		if (tmp >= dmax) {
			dmax_ = dmax;
			dmax = tmp;
		} else if (tmp > dmax_) {
			dmax_ = tmp;
		}
	}
	dmax += dmax_ + 1.0;
}

pair<double,double> getDistortion() {
	double cnt = 0;
	double maxDistor = 1.0;
	double avgDistor = 0.0;
	double d, dt, rat;
	
	for (int i=0; i<nV; ++i) {
		for (int j=i+1; j<nV; ++j) {
			d = dist(V, i, j);
			dt = distOnHST(leaves[i], leaves[j]);
			if (d > 0) {
				rat = max(dt/d, 1.0);
				maxDistor = max(maxDistor, rat);
				avgDistor += rat;
				cnt += 1.0;
			}
		}
	}
	if (cnt <= 0) {
		avgDistor = 1.0;
	} else {
		avgDistor /= cnt;
	}
	
	return make_pair(maxDistor, avgDistor);
}

double _calc_maxDistor_rtree0_UB(int centerId, partition_t& g, double UB) {
	double maxDistor = 0.0;
	vector<pair<double,int> > vp;
	int npos = g.ids.size();

	vector<int> levs;
	vector<int> cnt(H+2, 0);

	for (int i=0; i<npos; ++i) {
		double dis = dist(V, centerId, g.ids[i]);
		int lev = getLevel(H, dis);
		if (lev < g.lev) lev = g.lev;
		levs.push_back(lev);
		++cnt[lev];
	}
	for (int i=1; i<=H+1; ++i) cnt[i] += cnt[i-1];
	
	vector<int> U(npos, g.lev);
	for (int i=0; i<npos; ++i) {
		U[--cnt[levs[i]]] = i;
	}
	reverse(U.begin(), U.end());

#ifdef CHISQUARE
	flann::Index<flann::ChiSquareDistance<double> > index(flann::KDTreeIndexParams(1));
#elif defined(HELLINGER)
	flann::Index<flann::HellingerDistance<double> > index(flann::KDTreeIndexParams(1));
#else
	flann::Index<flann::L2_Simple<double> > index(flann::KDTreeIndexParams(1));
#endif 
	flann::SearchParams params(min(npos, FLANN_CHECK), eps_ann, false);	
	int nPts = 0, b = 0;
	
	for (int i=0,k=0,preLev=H+1,tid=0; i<npos&&maxDistor<UB; i=k,++tid) {
		while (k<npos && levs[U[k]]==levs[U[i]]) ++k;
		double dt = distAtLevel(levs[U[i]]);
		
		for (int j=i; j<k&&maxDistor<UB; ++j) {
			int pid = g.ids[U[j]], qid = 0;
			if (i==0 || nPts==0) continue;
			
			for (int idx=0; idx<DIM_V; ++idx)
				queryPoint[0][idx] = V[pid].x[idx];
			
			index.knnSearch(queryPoint, NNindices, NNdists, 1, params);
			qid = g.ids[U[b+NNindices[0][0]]];
			double d = dist(V, pid, qid);
			maxDistor = max(maxDistor, dt/d);
		}
		
		if (k > i) {
			nPts = k - i;
			b = i;
			int pos = 0;
			for (int j=i; j<k; ++j) {
				int pid = g.ids[U[j]];
				for (int idx=0; idx<DIM_V; ++idx)
					dataset[pos++] = V[pid].x[idx];
			}
			index.buildIndex(flann::Matrix<double>(dataset, nPts, DIM_V)); 
		}
		
		while (preLev >= levs[U[i]]) {
			distorUB[centerId][preLev] = min(distorUB[centerId][preLev], maxDistor);
			--preLev;
		}
	}
	
	return maxDistor;
}

void _getParition(int& nid, int centerId, partition_t& g, vector<partition_t>& left) {
	Node_t* far = g.ptr;
	partition_t gg;
	vector<pair<double,int> > vp;
	double radius = pow2(H-g.lev)*beta; 
	int npos = g.ids.size();
	
	vector<int> levs;
	vector<int> cnt(H+2, 0);
	for (int i=0; i<npos; ++i) {
		double dis = dist(V, centerId, g.ids[i]);
		int lev = getLevel(H, dis);
		if (lev < g.lev) lev = g.lev;
		levs.push_back(lev);
		++cnt[lev];
	}
	for (int i=1; i<=H+1; ++i) cnt[i] += cnt[i-1];
	
	vector<int> U(npos, g.lev);
	for (int i=0; i<npos; ++i) {
		U[--cnt[levs[i]]] = i;
	}
	
	for (int i=g.lev+1,idx=0; i<=H+1; ++i,radius/=alpha) {

		Node_t* child = new Node_t(nid++, centerId, i, far, radius);
		if (rt == NULL) rt = child; 
		addChild(far, child);

		if (i == H+1) {
			leaves[centerId] = child;
		}
		
		int _idx = idx;
		while (idx<npos && levs[U[idx]]<i) ++idx;
		if (idx > _idx) {
			gg.ids.clear();
			gg.ptr = far;
			gg.lev = i - 1;
			for (int j=_idx; j<idx; ++j) {
				gg.ids.push_back(g.ids[U[j]]);
			}
			left.push_back(gg);
		} else {
			mergeChild(far, child);
		}

		far = child;
	}
}

int _getCenter_maxDistor_complex(partition_t& g, FUNC_MAXDISTOR_UB _calc_maxDistor_UB) {
	int ret = -1;
	double mn = INF, tmp;
	int sz = g.ids.size();
	
	if (sz <= 1) {
		ret = g.ids[0];
	} else if (sz == nV) {
		int cnt_points_SAMPLE = max(1.0, log(sz)/log(2));
		for (int i=0; i<sz&&cnt_points_SAMPLE>0; ++i,--cnt_points_SAMPLE) {
			tmp = _calc_maxDistor_UB(g.ids[i], g, mn);
			if (tmp < mn) {
				mn = tmp;
				ret = g.ids[i];
			}
		}
		res = max(res, mn);
	} else if (sz < nV) {
		int cnt_points_SAMPLE = max(1.0, log(sz)/log(2));
		vector<pair<double,int> > UBs(sz);
		
		for (int i=0; i<sz; ++i)
			UBs[i] = make_pair(distorUB[g.ids[i]][g.lev+1], g.ids[i]);
		sort(UBs.begin(), UBs.end());
		
		for (int i=0; mn>res&&i<sz&&cnt_points_SAMPLE>0; ++i,--cnt_points_SAMPLE) {
			tmp = _calc_maxDistor_UB(UBs[i].second, g, mn);
			if (tmp < mn) {
				mn = tmp;
				ret = UBs[i].second;
			}
		}
		
		res = max(res, mn);
	}

	return ret;
}

void initParameters(int flag, int m) {	
	static double pc = max(1.0, log(m)/log(2)/6.0);
	static double pa = pow((double)nV, 0.5);
	pa = max(1.0, pa);
	pa = min((double)m, pa);
	if (flag == 0) {
		num_partitions_LB = pa;
		num_points_LB = pa;
		num_points_UB = pa;
		num_points_SAMPLE = min(m, (int)pc);
		if (num_points_SAMPLE == 0) num_points_SAMPLE = 1;
	} else if (flag == 1) {
		num_partitions_LB = pa;
		num_points_LB = pa / 2.0;
		num_points_UB = pa;
		num_points_SAMPLE = min(m, (int)pc);
		if (num_points_SAMPLE == 0) num_points_SAMPLE = 1;
	}
}

void _constructHST_partition_algo1(bool load, clock_t startClock) {
	_constructHST_partition_complex2(load, startClock, _calc_maxDistor_rtree0_UB, _calc_subspaceDistor_rtree0_UB2);
}

void _constructHST_partition_complex2(bool load, clock_t startClock, FUNC_MAXDISTOR_UB _calc_maxDistor_UB, FUNC_SUBDISTOR_UB _calc_subspaceDistor_UB) {
	res = maxDistor = 1.0;
	if (rt != NULL) {
		freeHST(rt);
		rt = NULL;
	}
	if (!load)
		randomization();
	calcDmax();

	H = ceil(log2(dmax+EPS)) + 1;

	int nid = 1;
	vector<partition_t> pars, _pars;
	{
		partition_t g;
		g.lev = 0;
		g.ptr = NULL;
		for (int i=0; i<nV; ++i)
			g.ids.push_back(i);
		_pars.push_back(g);
	}
	
	while (!_pars.empty()) {
		partition_t g = *_pars.rbegin();
		_pars.pop_back();
		int centerId = _getCenter_subspaceDistor_complex(g, _calc_subspaceDistor_UB);
		if (centerId == -1) {
			pars.push_back(g);
		} else {
			_getSubspaces(nid, centerId, g, _pars);
		}
	}
	
	while (!pars.empty()) {
		partition_t g = *pars.rbegin();
		pars.pop_back();
		int centerId = _getCenter_maxDistor_complex(g, _calc_maxDistor_UB);
		_getParition(nid, centerId, g, pars);
	}
	
#ifdef WATCH_MEM
	watchSolutionOnce(getpid(), usedMemory);
#endif
	int nid_ = countTree(rt);
	usedMemory = 1.0*nid_*sizeof(Node_t)/1024.0;
	printf("nV=%d, nid=%d, H=%d, mem=%.2lf\n", nV, nid_, H, 1.0*nid_*sizeof(Node_t)/1024.0);
	fflush(stdout);
}


int _getCenter_subspaceDistor_complex(partition_t& g, FUNC_SUBDISTOR_UB _calc_subspaceDistor_UB) {
	if (g.ids.size() <= num_points_UB) return -1;
	
	int ret = -1;
	double mn = INF, tmp;
	int sz = g.ids.size();
	
	if (sz <= 1) {
		ret = g.ids[0];
	} else if (sz == nV) {
		for (int i=0; i<sz; ++i) {
			mark[g.ids[i]] = -1;
		}
		
		vector<int> PI(sz, 0);
		for (int i=0; i<sz; ++i)
			PI[i] = g.ids[i];
		shuffle(PI.begin(), PI.end(), default_random_engine());
		int _centerId;
		
		initParameters(1, sz);
		int n_sample = num_points_SAMPLE;
		for (int i=0; i<sz&&n_sample>0; ++i) {
			int centerId = PI[i];
			if (mark[centerId] != -1) 
				continue;
			--n_sample;
			tmp = _calc_subspaceDistor_UB(centerId, g, mn);
			if (tmp < mn) {
				mn = tmp;
				ret = centerId;
			}
		}
		
		initParameters(0, sz);
		_centerId = ret;
		for (int i=0; i<sz; ++i) 
			mark[g.ids[i]] = -1;
		_calc_subspaceDistor_UB(_centerId, g, INF);
		PI.clear();
		for (int i=0; i<sz; ++i) {
			int centerId = g.ids[i];
			if (mark[centerId]==_centerId && centerId!=_centerId) {
				PI.push_back(centerId);
			}
		}
		
		shuffle(PI.begin(), PI.end(), default_random_engine());
		n_sample = num_points_SAMPLE;	
		for (int i=0; i<PI.size()&&n_sample>0; ++i) {
			int centerId = PI[i];
			--n_sample;
			tmp = _calc_subspaceDistor_UB(centerId, g, mn);
			if (tmp < mn) {
				mn = tmp;
				ret = centerId;
			}
		}
		
		res = max(res, mn);
	} else {
		
		for (int i=0; i<sz; ++i) {
			mark[g.ids[i]] = -1;
		}
		
		vector<pair<double,int> > UBs(sz);
		for (int i=0; i<sz; ++i)
			UBs[i] = make_pair(distorUB[g.ids[i]][g.lev+1], g.ids[i]);
		sort(UBs.begin(), UBs.end());
		
		vector<int> PI(sz, 0);
		for (int i=0; i<sz; ++i)
			PI[i] = UBs[i].second;
		int _centerId;
		
		initParameters(1, sz);
		int n_sample = num_points_SAMPLE;
		for (int i=0; mn>res&&i<sz&&n_sample>0; ++i) {
			int centerId = PI[i];
			if (mark[centerId] != -1) 
				continue;
			--n_sample;
			tmp = _calc_subspaceDistor_UB(centerId, g, mn);
			if (tmp < mn) {
				mn = tmp;
				ret = centerId;
			}
		}
		
		if (mn > res) {
			initParameters(0, sz);
			_centerId = ret;
			for (int i=0; i<sz; ++i) 
				mark[g.ids[i]] = -1;
			_calc_subspaceDistor_UB(_centerId, g, INF);
			PI.clear();
			for (int i=0; i<sz; ++i) {
				int centerId = UBs[i].second;
				if (mark[centerId]==_centerId && centerId!=_centerId) {
					PI.push_back(centerId);
				}
			}
			mark[_centerId] = _centerId;
			
			shuffle(PI.begin(), PI.end(), default_random_engine());
			n_sample = num_points_SAMPLE;
			for (int i=0; mn>res&&i<PI.size()&&n_sample>0; ++i) {
				int centerId = PI[i];
				--n_sample;
				tmp = _calc_subspaceDistor_UB(centerId, g, mn);
				if (tmp < mn) {
					mn = tmp;
					ret = centerId;
				}
			}
		}
		
		res = max(res, mn);
	}

	return ret;
}

void _getSubspaces(int& nid, int centerId, partition_t& g, vector<partition_t>& left) {
	if (g.ids.size() <= num_points_UB) return ;
	
	Node_t* far = g.ptr;
	partition_t gg;
	vector<pair<double,int> > vp;
	double radius = pow2(H-g.lev)*beta;
	int npos = g.ids.size();
	
	vector<int> levs;
	vector<int> cnt(H+2, 0);
	for (int i=0; i<npos; ++i) {
		double dis = dist(V, centerId, g.ids[i]);
		int lev = getLevel(H, dis);
		if (lev < g.lev) lev = g.lev;
		levs.push_back(lev);
		++cnt[lev];
	}
	for (int i=1; i<=H+1; ++i) cnt[i] += cnt[i-1];
	
	vector<int> U(npos, g.lev);
	for (int i=0; i<npos; ++i) {
		U[--cnt[levs[i]]] = i;
	}
	
	for (int i=g.lev+1,idx=0; i<=H+1; ++i,radius/=alpha) {

		Node_t* child = new Node_t(nid++, centerId, i, far, radius);
		if (rt == NULL) rt = child; 
		addChild(far, child);

		if (i == H+1) {
			leaves[centerId] = child;
		}

		int _idx = idx;
		while (idx<npos && levs[U[idx]]<i) ++idx;
		if (idx > _idx) {
			gg.ids.clear();
			gg.ptr = far;
			gg.lev = i - 1;
			for (int j=_idx; j<idx; ++j) {
				gg.ids.push_back(g.ids[U[j]]);
			}
			left.push_back(gg);
		} else {
			mergeChild(far, child);
		}

		if (npos-idx>0 && npos-idx<=num_points_UB) {
			gg.ids.clear();
			gg.ptr = far;
			gg.lev = i - 1;
			for (int j=idx; j<npos; ++j) {
				gg.ids.push_back(g.ids[U[j]]);
			}
			left.push_back(gg);
			break;
		}
		
		far = child;
	}
}

double _calc_subspaceDistor_rtree0_UB2(int centerId, partition_t& g, double UB) {
	double maxDistor = 0.0;
	int npos = g.ids.size();
	
	vector<int> levs;
	vector<int> cnt(H+2, 0);

	for (int i=0; i<npos; ++i) {
		double dis = dist(V, centerId, g.ids[i]);
		int lev = getLevel(H, dis);
		if (lev < g.lev) lev = g.lev;
		levs.push_back(lev);
		++cnt[lev];
	}
	for (int i=1; i<=H+1; ++i) cnt[i] += cnt[i-1];
	
	vector<int> U(npos, g.lev);
	for (int i=0; i<npos; ++i) {
		U[--cnt[levs[i]]] = i;
	}
	
#ifdef CHISQUARE
	flann::Index<flann::ChiSquareDistance<double> > index(flann::KDTreeIndexParams(1));
#elif defined(HELLINGER)
	flann::Index<flann::HellingerDistance<double> > index(flann::KDTreeIndexParams(1));
#else
	flann::Index<flann::L2_Simple<double> > index(flann::KDTreeIndexParams(1));
#endif 
	flann::SearchParams params(min(npos, FLANN_CHECK), eps_ann, false);	
	int nPts = 0, b = 0;
	
	for (int i=0,k=0; i<npos&&maxDistor<UB; i=k) {
		while (k<npos && levs[U[k]]==levs[U[i]]) ++k;
		double dt = distAtLevel(levs[U[i]]);
		
		if (k < npos) {
			nPts = 0;
			b = k;
			int pos = 0;
			for (int j=k; j<npos; ++j) {
				if (levs[U[j]] != levs[U[k]]) break;
				int pid = g.ids[U[j]];
				for (int dim_i=0; dim_i<DIM_V; ++dim_i)
					dataset[pos++] = V[pid].x[dim_i];
				++nPts;
			}
			index.buildIndex(flann::Matrix<double>(dataset, nPts, DIM_V));
			
			for (int j=i; j<k&&maxDistor<UB; ++j) {
				int pid = g.ids[U[j]], qid = 0;
				if (nPts == 0) continue;
				
				for (int dim_i=0; dim_i<DIM_V; ++dim_i)
					queryPoint[0][dim_i] = V[pid].x[dim_i];
				
				index.knnSearch(queryPoint, NNindices, NNdists, 1, params);	
				qid = g.ids[U[b+NNindices[0][0]]];
				double d = dist(V, pid, qid);
				maxDistor = max(maxDistor, dt/d);
			}
		}
		
		if (npos-k <= num_points_UB) break;
	}
	
	for (int i=1; i<=num_points_LB; ++i) {
		if (mark[U[npos-i]] == -1)
			mark[U[npos-i]] = centerId;
	}

	return maxDistor;
}

static int RMQ_top;
static vector<Node_t*> RMQ_V;
static vector<int> RMQ_D;
static vector<int> RMQ_beg;
static vector<vector<int> > RMQ_dp;

void RMQ_clear() {
	RMQ_V.clear();
	RMQ_D.clear();
	RMQ_beg.clear();
	RMQ_dp.clear();
}

void RMQ_resize() {
	RMQ_top = 0;
	RMQ_V.push_back((Node_t*)NULL);
	RMQ_D.push_back(0);
	
	int nid = 1;
	for (int i=0; i<nV; ++i)
		nid = max(nid, leaves[i]->nid);
	RMQ_beg.resize(nid+1, 0);
	
	int IN = 2*nid + 5;
	int JN = log(1.0*IN)/log(2.0)+4;
	RMQ_dp.resize(JN, vector<int>(IN, 0));
}

void RMQ_dfs(Node_t* rt, int d) {
	if (rt == NULL) return ;
	
	++RMQ_top;
	RMQ_V.push_back(rt);
	RMQ_D.push_back(d);
	RMQ_beg[rt->nid] = RMQ_top;
	
	for (int i=0; i<rt->child.size(); ++i) {
		RMQ_dfs(rt->child[i], d+1);
		++RMQ_top;
		RMQ_V.push_back(rt);
		RMQ_D.push_back(d);
	}
}	

void RMQ_init(int n) {
    int i, j;

    for (i=1; i<=n; ++i)
        RMQ_dp[0][i] = i;
    for (j=1; (1<<j)<=n; ++j)
        for (i=1; i+(1<<j)-1<=n; ++i)
            if (RMQ_D[RMQ_dp[j-1][i]] < RMQ_D[RMQ_dp[j-1][i+(1<<(j-1))]])
                RMQ_dp[j][i] = RMQ_dp[j-1][i];
            else
                RMQ_dp[j][i] = RMQ_dp[j-1][i+(1<<(j-1))];
}

Node_t* RMQ_query(int l, int r) {
    if (l > r)
        swap(l, r);

    int k = floor(log(r-l+1.0)/log(2.0) - 1.0);

    while (1<<(k+1) <= r-l+1)
        ++k;

    if (RMQ_D[RMQ_dp[k][l]] < RMQ_D[RMQ_dp[k][r-(1<<k)+1]])
        return RMQ_V[RMQ_dp[k][l]];
    else
        return RMQ_V[RMQ_dp[k][r-(1<<k)+1]];
}

double RMQ_distOnHST(Node_t* pu, Node_t* pv) {
	Node_t* lca = RMQ_query(RMQ_beg[pu->nid], RMQ_beg[pv->nid]);
	return distAtLevel(lca->lev);
}

pair<double,double> getDistortion_fast() {
	RMQ_resize();
	RMQ_dfs(rt, 1);
	RMQ_init(RMQ_top);
	
	double cnt = 0;
	double maxDistor = 1.0;
	double avgDistor = 0.0;
	double d, dt, rat;
	
	for (int i=0; i<nV; ++i) {
		for (int j=i+1; j<nV; ++j) {
			d = dist(V, i, j);
			dt = RMQ_distOnHST(leaves[i], leaves[j]);
			if (d > 0) {
				rat = max(dt/d, 1.0);
				maxDistor = max(maxDistor, rat);
				avgDistor += rat;
				cnt += 1.0;
			}
		}
	}
	if (cnt <= 0) {
		avgDistor = 1.0;
	} else {
		avgDistor /= cnt;
	}
	
	RMQ_clear();
	return make_pair(maxDistor, avgDistor);
}
