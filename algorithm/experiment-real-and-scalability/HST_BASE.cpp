#include "HST_BASE.h"
#include "global.h"
#ifdef WATCH_MEM
#include "monitor.h"
#endif

const int MAX_HEIGHT = 55;
int H = 0;
int alpha = 2;
double logAlpha = 1.0;
int* pi = NULL;
int* reverse_pi = NULL;
double dmax = -1.0;
double* expks = NULL;
double* sumks = NULL;
Node_t* rt = NULL;
Node_t** leaves = NULL;
double beta = 1.0;

inline double log2(double x) {
	return log(x) / logAlpha;
}

inline double pow2(int i) {
	return expks[i];
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
	
	if (expks[k]*beta == dist)
		++k;
	
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
	
	logAlpha = log(alpha);
}

void freeMemory_HST() {
	freeLocation();
	delete[] pi;
	delete[] reverse_pi;
	delete[] expks;
	delete[] sumks;
    delete[] leaves;
	freeHST(rt);
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
        fin >> V[i].x;
        fin >> V[i].y;
    }
    fin.close();
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

double sampleBeta_InJCSS04() {
	const double lnk = log10(alpha);
	const int base = 1000, beg = 1;
	double dx = (alpha - beg) * 1.0 / 1000;
	vector<double> P(base, 0.0);
	double sum = 0.0;
	
	for (int i=0; i<base; ++i) {
		P[i] = dx / (lnk * (beg + i*dx));
		sum += P[i];
	}
	
	for (int i=0; i<base; ++i) {
		P[i] /= sum;
		if (i>0) P[i] += P[i-1];
	}
	
	double randf = rand() % (base*10) / (base*10.0);
	
	double ret = alpha;
	for (int i=0; i<base; ++i) {
		if (randf <= P[i]) {
			ret = beg + i * dx;
			break;
		}
	}
	
	ret /= alpha;
	ret = min(ret, 1.0);
	ret = max(1.0/alpha, ret);	
	
	return ret;
}

double sampleBeta_InSTOC03() {
	int base = alpha*1000;
	beta = (rand()%base + 1) * 1.0	/ base;
	beta = min(beta, 1.0);
	beta = max(1.0/alpha, beta);
	return beta;
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
	if (far == rt || far->child.size()!=1) 
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

void calcDmaxPrune() {
	dmax = 0;
	for (int j=1; j<nV; ++j) {
		dmax = max(dmax, dist(V, pi[0], j));
	}
}

void calcDmax() {
	if (dmax >= 0.0) return ;

	dmax = 0.0;
	for (int i=0; i<nV; i++) {
		for (int j=i+1; j<nV; j++) {
			dmax = max(dmax, dist(V, i, j));
		}
	}
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

void _constructHST_ICDE_DPO(bool load, clock_t startClock) {
	int nid = 1;
	int *cnt = NULL, *cen = NULL;
	Node_t** nodes = NULL;
	Node_t** nodes_ = NULL; 
	double *cenl = NULL;
	int *p = NULL, *q = NULL, *far = NULL, *far_ = NULL;
	
	cnt = new int[nV];
	cen = new int[nV];
	cenl = new double[nV];
	p = new int[nV];
	q = new int[nV];
	far = new int[nV];
	far_ = new int[nV];
	nodes = new Node_t*[nV];
	nodes_ = new Node_t*[nV];
	
	if (rt != NULL)
		freeHST(rt);
	if (!load)
		randomization();

	dmax = 0;
	for (int i=0; i<nV; ++i) {
		cen[i] = 0;
		cenl[i] = dist(V, i, pi[cen[i]]);
		dmax = max(dmax, cenl[i]);
	}
	
	H = ceil(log2(dmax+EPS)) + 1;	
	double radius = pow2(H)*beta;
	
	rt = new Node_t(nid++, pi[0], 1, NULL, radius);
	for (int i=0; i<nV; ++i) {
        nodes_[i] = rt;
		far[i] = far_[i] = nid - 1;
		p[i] = i;
	}
	
	for (int k=2; k<=H+1; ++k) {
		radius /= alpha;
		for (int i=0; i<nV; ++i) {
			if (cenl[i] < radius)
				continue;
			
			int pid;
			while (cenl[i] >= radius) {
				 pid = pi[++cen[i]];
				 if (cen[pid] <= reverse_pi[i]) {
					 cenl[i] = dist(V, i, pi[cen[i]]);
				 }
			}
		}
		
		int bid = nid, i = 0, bi;
		memset(q, -1, sizeof(int)*nV);
		while (i < nV) {
			q[cen[p[i]]] = far[p[i]] = nid++;
			bi = i++;
			while (i<nV && far_[p[i]]==far_[p[i-1]]) {
				int pid = cen[p[i]];
				if (q[pid] == -1) {
					q[pid] = nid++;
				}
				far[p[i]] = q[pid];
				++i;
			}
			while (bi < i) {
				int pid = cen[p[bi]];
				q[pid] = -1;
				++bi;
			}
		}
		memset(cnt, 0, sizeof(int)*nV);
		for (int i=0; i<nV; ++i) {
			++cnt[far[i]-bid];
		}
		for (int j=1; j<nid-bid; ++j) {
			cnt[j] += cnt[j-1];
		}
		for (int i=nV-1; i>=0; --i) {
			int j = far[p[i]] - bid;
			q[--cnt[j]] = p[i];
		}
		
		for (int i=0,j=0; i<nV; i=j) {
			j = i;
			nodes[q[i]] = new Node_t(far[q[i]], pi[cen[q[i]]], k, nodes_[q[i]], radius);
			addChild(nodes_[q[i]], nodes[q[i]]);
			while (j<nV && far[q[j]]==far[q[i]]) {
				nodes[q[j]] = nodes[q[i]];
				++j;
			}
		}
		
		for (int i=0,j=0; i<nV; i=j) {
			j = i;
			mergeChild(nodes_[q[i]], nodes[q[i]]);
			while (j<nV && far[q[j]]==far[q[i]]) {
				++j;
			}
		}		
		
		if (k == H+1) {
			for (int i=0; i<nV; ++i) {
				leaves[q[i]] = nodes[q[i]];
			}
		}
		swap(far, far_);
		swap(p, q);
		swap(nodes, nodes_);
	}
	
#ifdef WATCH_MEM
	watchSolutionOnce(getpid(), usedMemory);
#endif

	int nid_ = countTree(rt);
	usedMemory = 1.0*nid_*sizeof(Node_t)/1024.0;
	printf("nV=%d, nid=%d, H=%d, mem=%.2lf\n", nV, nid_, H, 1.0*nid_*sizeof(Node_t)/1024.0);
	fflush(stdout);
	
	delete[] cnt;
	delete[] cen;
	delete[] cenl;
	delete[] p;
	delete[] q;
	delete[] far;
	delete[] far_;
	delete[] nodes;
	delete[] nodes_;
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
