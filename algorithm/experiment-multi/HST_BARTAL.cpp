#include "HST_BARTAL.h"
#include "global.h"
#ifdef WATCH_MEM
#include "monitor.h"
#endif

const int MAX_HEIGHT = 55;
int H = 0;
int alpha = 2;
double beta;
double logAlpha = 1.0;
double dmax = -1.0;
double* expks = NULL;
double* sumks = NULL;
Node_t* rt = NULL;
Node_t** leaves = NULL;

inline double log2(double x) {
	return log(x) / logAlpha;
}

inline double pow2(int i) {
	return expks[i];
}

void dumpNode(Node_t* p) {
	for (int i=1; i<p->lev; ++i)
		putchar(' ');
	printf("nid=%d, center=%d, wei=%.4lf, |child|=%d\n", p->nid, p->nomIdx, p->wei, (int)p->child.size());
}

void dumpTree(Node_t* rt) {
	if (rt != NULL) {
		dumpNode(rt);
		for (int i=0; i<rt->child.size(); ++i)
			dumpTree(rt->child[i]);
	}
}

int getTreeHeight(Node_t* rt) {
	int ret = 0, height;
	Node_t* p = rt;
	
	for (int i=0; i<nV; ++i) {
		p = leaves[i];
		height = 0;
		while (p!=NULL && p!=rt) {
			++height;
			p = p->far;
		}
		ret = max(ret, height);
	}
	
	return ret;
}

void initMemory_HST(int n) {
	srand(time(NULL));
	nV = n;
	V = new location_t[nV];
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

void initVariable(string &fileName) {
	ifstream fin(fileName.c_str(), ios::in);   
    
	if (!fin.is_open()) {
		fprintf(stderr, "FILE %s IS INVALID.", fileName.c_str());
		exit(1);
	}
	
	fin >> beta >> alpha;
    fin.close();	
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


double distOnHST(int u, int v) {
	return distOnHST(leaves[u], leaves[v]);
}

double distOnHST(Node_t* u, Node_t* v) {
	double ret = 0.0;
	
	while (u!=NULL && v!=NULL && u!=v) {
		if (u->lev > v->lev) {
			ret += u->wei;
			u = u->far;
		} else {
			ret += v->wei;
			v = v->far;
		}
	}
	
	if (u==NULL || v==NULL)
		ret = -1.0;
	
	return ret;
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

void calcDmax() {
	if (dmax >= 0.0) return ;

	dmax = 0.0;
	double tmp;
	for (int i=0; i<nV; i++) {
		for (int j=i+1; j<nV; j++) {
			tmp = dist(V, i, j);
			dmax = max(dmax, tmp);
		}
	}
	dmax += EPS;
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

double sampleRadius_focs96(int n, double r) {
	const int n_bucket = 20;
	double mx = r * log(n);
	double dx = mx / n_bucket;
	double CDF[n_bucket+1];
	
	CDF[0] = 0.0;
	for (int i=1; i<n_bucket; ++i) {
		double x = i * dx;
		CDF[i] = (1.0 - exp(0.0-x/r)) * n / (n-1);
	}
	CDF[n_bucket] = 1.0;
	
	double p = rand()%10000 / 10000.0;
	double ret = mx;
	for (int i=1; i<=n_bucket; ++i) {
		if (CDF[i] >= p) {
			ret = i * dx;
			ret = min(ret, mx);
			break;
		}
	}
	
	return ret;
}

void _construct_probabilistic_partition(vector<cluster_t>& C, int clusterId, int& n_cluster) {
	cluster_t c = C[clusterId], nc;
	if (c.calc_n_VG() <= 1) return ;
	
	while (true) {
		int n = c.calc_n_VG();
		if (n == 0) break;
		
		double rho = 2.0*log(n)+1, diam = c.calc_diameter();
		double r = floor(diam / (rho*alpha));
		double z = sampleRadius_focs96(n, r);
		
		if (diam < 2*r*log(n)) {
			nc = c;
			nc.id = n_cluster++;
			nc.far = c.id;
			nc.lev = c.lev + 1;
			C.push_back(nc);
			c.VG.clear();
		} else {
			int vt = c.VG[rand()%n];
			vector<bool> visit(n, false);
			
			nc.VG.clear();
			for (int i=0; i<n; ++i) {
				if (dist(V, vt, c.VG[i]) <= z) {
					visit[i] = true;
					nc.VG.push_back(c.VG[i]);
				}
			}
			
			double radius_LB = r / n;
			queue<int> Q;
			for (int i=0; i<nc.VG.size(); ++i) Q.push(nc.VG[i]);
			while (!Q.empty()) {
				int v = Q.front();
				Q.pop();
				for (int i=0; i<n; ++i) {
					if (!visit[i] && dist(V, v, c.VG[i])<=radius_LB) {
						visit[i] = true;
						nc.VG.push_back(c.VG[i]);
						Q.push(c.VG[i]);
					}
				}
			}

			nc.id = n_cluster++;
			nc.far = c.id;
			nc.centerId = vt;
			nc.lev = c.lev + 1;
			C.push_back(nc);
			
			for (int i=0,j=0; i<n; ++i) {
				if (!visit[i]) {
					c.VG[j++] = c.VG[i];
				}
				if (i == n-1) {
					c.VG.erase(c.VG.begin()+j, c.VG.end());
				}
			}
		}
	}
}

void _constructHST_focs96(bool load, clock_t startClock) {
	if (rt != NULL)
		freeHST(rt);
	
	vector<cluster_t> C;
	cluster_t c;
	c.id = 0, c.far = -1, c.centerId = rand()%nV, c.lev = 1;
	for (int i=0; i<nV; ++i) {
		c.VG.push_back(i);
	}
	C.push_back(c);
	
	int n_cluster = 1;
	int _C_sz = 0, C_sz;
	while (true) {
		C_sz = C.size();
		if (_C_sz == C_sz) break;
		while (_C_sz < C_sz) {
			_construct_probabilistic_partition(C, _C_sz, n_cluster);
			++_C_sz;
		}
	}
	
	int n_node = 1;
	vector<Node_t*> nodes;
	
	for (int i=0; i<C.size(); ++i) {
		Node_t* parent = (C[i].far==-1) ? NULL : nodes[C[i].far];
		double radius = (C[i].far==-1) ? 0.0 : C[C[i].far].calc_diameter()*(0.5+EPS);
		Node_t* child = new Node_t(n_node, C[i].centerId, C[i].lev, parent, radius);
		if (i == 0) {
			rt = child;
		} else {
			addChild(parent, child);
		}
		if (C[i].calc_n_VG() == 1) {
			int vid = C[i].VG[0];
			leaves[vid] = child;
		}
		nodes.push_back(child);
		n_node++;
	}
	
#ifdef WATCH_MEM
	watchSolutionOnce(getpid(), usedMemory);
#endif
	usedMemory = 1.0*n_node*sizeof(Node_t)/1024.0;
	printf("nV=%d, nid=%d, H=%d, sz=%.2lf\n", nV, n_node, getTreeHeight(rt), 1.0*n_node*sizeof(Node_t)/1024.0);
	fflush(stdout);
}

double calc_g_stoc98(cluster_t& c) {
	int n_sample = 20;
	double Delta = c.diam;
	int n = c.calc_n_VG();
	double T = 1.0 * n * n;
	double g;

	default_random_engine rng;
	uniform_real_distribution<double> uni(1.0, T*1.0-EPS);
	double mn = INF, gstar = 1.0;
	
	while (n_sample-- > 0) {
		g = uni(rng);
		if (g < 1.0) g = 1.0;
		if (g > T) g = T-EPS;
		double zg = floor(g) * (Delta*0.5) / T;
		int vt = c.VG[rand()%n];
		vector<int> I1, I2;
		for (int i=0; i<n; ++i) {
			if (dist(V, vt, c.VG[i]) <= zg) {
				I1.push_back(i);
			} else {
				I2.push_back(i);
			}
		}
		double tmp = 0;
		double le = Delta;
		for (int i=0; i<I1.size(); ++i) {
			for (int j=0; j<I2.size(); ++j) {
				double we = dist(V, c.VG[I1[i]], c.VG[I2[j]]);
				tmp += le / we;
			}
		}
		if (tmp < mn) {
			mn = tmp;
			gstar = g;
		}
	}
	
	return gstar;
}

double sampleRadius_stoc98(cluster_t& c) {
	double Delta = c.diam;
	long long n = c.calc_n_VG();
	long long T = n * n;
	double g = calc_g_stoc98(c);

	if (g > T) g = T;
	double zg = floor(g) * (Delta*0.5) / T;
	return zg;
}

void _construct_HPM(cluster_t& _c, cluster_t& Hg1, cluster_t& Hg2, int& n_cluster) {
	cluster_t c = _c;
	if (c.calc_n_VG() <= 1) return ;

	int n = c.calc_n_VG();
	double diam = c.calc_diameter();
	double T = 1.0 * n * n;
	double z = sampleRadius_stoc98(c);

	int vt = c.VG[rand()%n];
	vector<bool> visit(n, false);

	Hg1.VG.clear();
	for (int i=0; i<n; ++i) {
		if (dist(V, vt, c.VG[i]) <= z) {
			visit[i] = true;
			Hg1.VG.push_back(c.VG[i]);
		}
	}
	
	double radius_LB = diam / (2.0*T);
	queue<int> Q;
	for (int i=0; i<Hg1.VG.size(); ++i) Q.push(Hg1.VG[i]);
	while (!Q.empty()) {
		int v = Q.front();
		Q.pop();
		for (int i=0; i<n; ++i) {
			if (!visit[i] && dist(V, v, c.VG[i])<radius_LB) {
				visit[i] = true;
				Hg1.VG.push_back(c.VG[i]);
				Q.push(c.VG[i]);
			}
		}
	}
	Hg1.id = n_cluster++;
	Hg1.far = c.id;
	Hg1.centerId = vt;
	Hg1.lev = c.lev + 1;

	Hg2.VG.clear();
	for (int i=0; i<n; ++i) {
		if (!visit[i]) {
			Hg2.VG.push_back(c.VG[i]);
		}
	}
	if (Hg2.VG.size() > 0) {
		Hg2.id = n_cluster++;
		Hg2.far = c.id;
		Hg2.centerId = Hg2.VG[rand()%Hg2.calc_n_VG()];
		Hg2.lev = c.lev + 1;
	}
}

void _construct_HPM(cluster_t& _c, double L, queue<pair<cluster_t,double> >& Q, int& n_cluster) {
	if (_c.calc_n_VG() <= 1) return ;

	queue<cluster_t> C;
	C.push(_c);
	while (!C.empty()) {
		cluster_t c = C.front();
		C.pop();

		if (c.calc_diameter() < L/alpha) {
			cluster_t nc = c;
			nc.id = n_cluster++;
			nc.far = c.id;
			nc.centerId = c.VG[rand()%c.calc_n_VG()];
			nc.lev = c.lev + 1;
			Q.push(make_pair(nc, L/alpha));

		} else {

			int n = c.calc_n_VG();
			double diam = c.calc_diameter();
			double T = 1.0 * n * n;
			double z = sampleRadius_stoc98(c);

			int vt = c.VG[rand()%n];
			vector<bool> visit(n, false);
			cluster_t Hg1 = c, Hg2 = c;

			Hg1.VG.clear();
			for (int i=0; i<n; ++i) {
				if (dist(V, vt, c.VG[i]) <= z) {
					visit[i] = true;
					Hg1.VG.push_back(c.VG[i]);
				}
			}
			
			double radius_LB = diam / (2.0*T);
			queue<int> Q;
			for (int i=0; i<Hg1.VG.size(); ++i) Q.push(Hg1.VG[i]);
			while (!Q.empty()) {
				int v = Q.front();
				Q.pop();
				for (int i=0; i<n; ++i) {
					if (!visit[i] && dist(V, v, c.VG[i])<radius_LB) {
						visit[i] = true;
						Hg1.VG.push_back(c.VG[i]);
						Q.push(c.VG[i]);
					}
				}
			}
			Hg1.centerId = vt;
			C.push(Hg1);

			Hg2.VG.clear();
			for (int i=0; i<n; ++i) {
				if (!visit[i]) {
					Hg2.VG.push_back(c.VG[i]);
				}
			}
			if (Hg2.VG.size() > 0) {
				Hg2.centerId = Hg2.VG[rand()%Hg2.calc_n_VG()];
				C.push(Hg2);
			}
		}
	}
}
 
void _constructHST_stoc98(bool load, clock_t startClock) {
	if (rt != NULL)
		freeHST(rt);
	rt = NULL;
	calcDmax();
	
	double L = dmax + EPS;
	queue<pair<cluster_t,double> > C;
	cluster_t c;
	c.id = 0, c.far = -1, c.centerId = rand()%nV, c.lev = 1;
	for (int i=0; i<nV; ++i) {
		c.VG.push_back(i);
	}
	C.push(make_pair(c,L));

	int n_cluster = 1;
	int n_node = 1;
	vector<Node_t*> nodes;
	cluster_t nc, Hg1, Hg2;

	while (!C.empty()) {
		c = C.front().first;
		L = C.front().second;
		C.pop();

		_construct_HPM(c, L, C, n_cluster);
		Node_t* parent = (c.far==-1) ? NULL : nodes[c.far];
		double radius = (c.far==-1) ? 0.0 : L*alpha*(0.5+EPS);
		Node_t* child = new Node_t(n_node, c.centerId, c.lev, parent, radius);
		if (rt == NULL) {
			rt = child;
		} else {
			addChild(parent, child);
		}
		if (c.calc_n_VG() == 1) {
			int vid = c.VG[0];
			leaves[vid] = child;
		}
		nodes.push_back(child);
		n_node++;
	}
	
#ifdef WATCH_MEM
	watchSolutionOnce(getpid(), usedMemory);
#endif
	usedMemory = 1.0*n_node*sizeof(Node_t)/1024.0;
	printf("nV=%d, nid=%d, H=%d, sz=%.2lf\n", nV, n_node, getTreeHeight(rt), 1.0*n_node*sizeof(Node_t)/1024.0);
	fflush(stdout);
}
