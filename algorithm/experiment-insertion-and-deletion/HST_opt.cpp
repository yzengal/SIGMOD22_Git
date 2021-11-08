#include <bits/stdc++.h>
using namespace std;

#include "global.h"
#include "HST_opt.h"
#ifdef WATCH_MEM
#include "monitor.h"
#endif

#define MAXH 40


void Tree_t::initMemory(int n) {
	nV = n;
	expks = new double[MAXH];	
	expks[0] = 1.0;
	for (int i=1; i<MAXH; ++i)
		expks[i] = expks[i-1] * alpha;
	
	sumks = new double[MAXH+1];
	sumks[0] = expks[0];
	for (int i=1; i<=MAXH; ++i)
		sumks[i] = sumks[i-1] + expks[i];
	logAlpha = log(alpha);
	
	leaves.resize(nV);
	pi.resize(nV);
	reverse_pi.resize(nV);
	V.resize(nV);
	ids.resize(nV);
	deleted.resize(nV);
	fill(deleted.begin(), deleted.end(), false);
}

void Tree_t::freeMemory() {
	delete[] expks;
	delete[] sumks;
	freeHST(rt);
}

long long Tree_t::memoryCost() {
	long long ret = 0;
	
	ret += 1LL * node_n * sizeof(Node_t);
	ret += 1LL * deleted.size() * sizeof(bool);
	for (int i=0; i<nns.size(); ++i) {
		ret += 1LL * nns[i].size() * (sizeof(int)+sizeof(double));
	}
	
	return ret;	
}

void Tree_t::freeHST(Node_t*& rt) {
	if (rt == NULL) return ;
	
	Node_t* chi;
	
	for (int i=0; i<rt->child.size(); ++i) {
		chi = rt->child[i];
		freeHST(chi);
	}
	
	delete rt;
	rt = NULL;
}

void Tree_t::initLocation(int nV_, vector<location_t> &V_, vector<int>& ids_, double alpha_) {
	dV = 0, nV = nV_, alpha = alpha_;
	int base = alpha * 1000;
	beta = (rand()%base + 1000) * 1.0 / base;
	beta = min(beta, 1.0);
	beta = max(1.0/alpha, beta);
	
	//initMemory_HST(nV);
	initMemory(nV);
	for(int i=0; i<nV; ++i){
		V[i] = V_[i];
		ids[i] = ids_[i];
		pi[i] = i;
	}
	// generate the permutation pi
	random_shuffle(pi.begin(), pi.end());
	for (int i=0; i<nV; ++i){
		reverse_pi[pi[i]] = i;
	}
}

void Tree_t::constructHST(clock_t startClock) {
	int nid = 1;
	int *cnt = NULL, *cen = NULL;
	Node_t** nodes = NULL;
	Node_t** nodes_ = NULL; 
	double *cenl = NULL;
	int *p = NULL, *q = NULL, *far = NULL, *far_ = NULL;
	int *grp = NULL;
	distor = 1.0, distor_ = 0.0, cdistor = 0.0;
	node_n = merge_n = 0;
	
	cnt = new int[nV];
	cen = new int[nV];
	cenl = new double[nV];
	p = new int[nV];
	q = new int[nV];
	far = new int[nV];
	far_ = new int[nV];
	nodes = new Node_t*[nV];
	nodes_ = new Node_t*[nV];
	grp = new int[nV];
	
	if (rt != NULL)
		freeHST(rt);

	// Prunning 1: calcDmax()
	dmax = 0;
	for (int i=0; i<nV; ++i) {
		cen[i] = 0;
		cenl[i] = dist(V[i], V[pi[cen[i]]]);
		dmax = max(dmax, cenl[i]);
	}
	
	// initialization
	H = ceil(log2(dmax+EPS)) + 1;	
	double radius = pow2(H)*beta;
	
	// construct the root
	rt = new Node_t(0, pi[0], 1, NULL, 0);
	for (int i=0; i<nV; ++i) {
        nodes_[i] = rt;
		far[i] = far_[i] = 0;
		p[i] = i;
	}
	
	bool checkFlag = true;
	for (int k=2; k<=H+1; ++k) {
		radius /= alpha;
		for (int i=0; i<nV; ++i) {
			if (cenl[i] < radius)
				continue;
			
			int pid;
			while (cenl[i] >= radius) {
				// Prunning 2: skip j
				 pid = pi[++cen[i]];
				 if (cen[pid] <= reverse_pi[i]) {
					 cenl[i] = dist(V[i], V[pi[cen[i]]]);
				 }
			}
		}
		
		int bid = nid, i = 0, bi;
		memset(q, -1, sizeof(int)*nV);
		memset(cnt, 0, sizeof(int)*nV);
		while (i < nV) {
			int _nid = nid;
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
		// create the new node at the $k$-th level
		for (int i=0,j=0; i<nV; i=j) {
			j = i;
			nodes[q[i]] = new Node_t(far[q[i]], pi[cen[q[i]]], k, nodes_[q[i]], radius*alpha);
			addChild(nodes_[q[i]], nodes[q[i]]);
			while (j<nV && far[q[j]]==far[q[i]]) {
				nodes[q[j]] = nodes[q[i]];
				++j;
			}
		}
		// merge the new node with its parent 
		for (int i=0,j=0; i<nV; i=j) {
			j = i;
			mergeChild(nodes_[q[i]], nodes[q[i]]);
			while (j<nV && far[q[j]]==far[q[i]]) {
				++j;
			}
		}
		// fill in the leaves
		if (k == H+1) {
			for (int i=0; i<nV; ++i) {
				leaves[q[i]] = nodes[q[i]];
			}
		}
		swap(far, far_);
		swap(p, q);
		swap(nodes, nodes_);
	}
	// calculate distortion
	distor_ = (cdistor<=0.0) ? 1.0 : (distor_/cdistor);
	
	delete[] cnt;
	delete[] cen;
	delete[] cenl;
	delete[] p;
	delete[] q;
	delete[] far;
	delete[] far_;
	delete[] nodes;
	delete[] nodes_;	
	delete[] grp;
	
	node_n = nid - merge_n;
}

double Tree_t::distAtLevel(int lev) {
	if (lev >= H+1) return 0.0;
	return sumks[H-lev+1] * beta * 2.0;

}

double Tree_t::distOnHST(int u, int v) {
	int level = levelOfLCA(u, v);
	return distAtLevel(level);
}

double Tree_t::distOnHST(Node_t* u, Node_t* v) {
	int level = levelOfLCA(u, v);
	return distAtLevel(level);
}

int Tree_t::levelOfLCA(int u, int v) {
	return levelOfLCA(leaves[u], leaves[v]);
}

int Tree_t::levelOfLCA(Node_t* u, Node_t* v) {
	if (u==NULL || v==NULL)
		return -1;
	
	while (u!=NULL && v!=NULL && u!=v) {
		if (u->lev > v->lev) {
			u = u->far;
		} else if (u->lev < v->lev){
			v = v->far;
		} else {
			u = u->far;
			v = v->far;
		}
	}
	
	if (u==NULL || v==NULL)
		return -1;
	
	return u->lev;
}

pair<Node_t*,int> Tree_t::getLCA(int u, int v) {
	return getLCA(leaves[u], leaves[v]);
}

pair<Node_t*,int> Tree_t::getLCA(Node_t* u, Node_t* v) {
	if (u==NULL || v==NULL)
        return make_pair(rt, -1);
	
	while (u!=NULL && v!=NULL && u!=v) {
		if (u->lev > v->lev) {
			u = u->far;
		} else if (u->lev < v->lev){
			v = v->far;
		} else {
			u = u->far;
			v = v->far;
		}
	}
	
	if (u==NULL || v==NULL)
		return make_pair(rt, -1);
	
	return make_pair(u, u->lev);
}

pair<double,double> Tree_t::getDistortion() {
	double cnt = 0;
	double d, dt, rat;
	
	distor = 1.0;
	distor_ = 0.0;
	for (int i=0; i<nV; ++i) {
		if (deleted[i])
			continue;
		for (int j=i+1; j<nV; ++j) {
			if (deleted[j])
				continue;
			d = dist(V[i], V[j]);
			dt = distOnHST(leaves[i], leaves[j]);
			if (d > 0) {
				rat = max(dt/d, 1.0);
				distor = max(distor, rat);
				distor_ += rat;
				cnt += 1.0;
			}
		}
	}
	if (cnt <= 0) {
		distor_ = 1.0;
	} else {
		distor_ /= cnt;
	}
	
	return make_pair(distor, distor_);
}	
		
void Tree_t::initNeighbours(int idx, Tree_t& oth) {
	nns.push_back(vector<pair<int,double> >(nV, make_pair(-1,INF)));
	for (int i=0; i<nV; ++i) {
		if (deleted[i])
			continue;
		double dmin = INF, tmp;
		int id = -1;
		for (int j=0; j<oth.nV; ++j) {
			if (oth.deleted[j])
				continue;
			tmp = dist(V[i], oth.V[j]);
			if (tmp < dmin) {
				dmin = tmp;
				id = j;
			}
		}
		nns[idx][i] = make_pair(id, dmin);
	}
}

void Tree_t::updateNeighbours(int idx, Tree_t& oth) {
	for (int i=0; i<nV; ++i) {
		if (deleted[i]) continue;
		int nnId = nns[idx][i].first;
		if (idx==0 && nnId>=0)
			continue;
		double dmin = INF, tmp;
		int id = -1;
		for (int j=0; j<oth.nV; ++j) {
			if (oth.deleted[j])
				continue;
			tmp = dist(V[i], oth.V[j]);
			if (tmp < dmin) {
				dmin = tmp;
				id = j;
			}
		}
		nns[idx][i] = make_pair(id, dmin);
	}
}

void Tree_t::updateNeighbour(int i, int idx, Tree_t& oth) {
	int nnId = nns[idx][i].first;
	if (idx==0 && nnId>=0)
		return;

	double dmin = INF, tmp;
	int id = -1;
	for (int j=0; j<oth.nV; ++j) {
		if (oth.deleted[j])
			continue;
		tmp = dist(V[i], oth.V[j]);
		if (tmp < dmin) {
			dmin = tmp;
			id = j;
		}
	}
	nns[idx][i] = make_pair(id, dmin);
}

void Tree_t::deleteOne(int idx) {
	if (!deleted[idx]) {
		deleted[idx] = true;
		++dV;
	}
}

void Tree_t::loadLocation(int nV_, vector<location_t> &V_, vector<int>& ids_, double alpha_) {
	// 1. remove the deleted points
	int _nV = 0;
	for (int i=0; i<V.size(); ++i) {
		if (deleted[i])
			continue;
		V[_nV] = V[i], ids[_nV] = ids[i];
		for (int j=0; j<nns.size(); ++j) {
			nns[j][_nV] = nns[j][i];
		}
		++_nV;
	}
	dV=0, nV = _nV+nV_, alpha = alpha_;
	
	// 2. sample beta
	int base = alpha * 1000;
	beta = (rand()%base + 1000) * 1.0 / base;
	beta = min(beta, 1.0);
	beta = max(1.0/alpha, beta);
	logAlpha = log(alpha);
	
	// 3. clear the memory and resize the vector
	freeHST(rt);
	rt = NULL;
	leaves.resize(nV);
	pi.resize(nV);
	reverse_pi.resize(nV);
	V.resize(nV);
	ids.resize(nV);
	deleted.resize(nV);
	fill(leaves.begin(), leaves.end(), (Node_t*)NULL);
	fill(deleted.begin(), deleted.end(), false);
	for (int j=0; j<nns.size(); ++j) {
		nns[j].resize(nV);
		fill(nns[j].begin()+_nV, nns[j].end(), make_pair(-1,INF));
	}
	
	// 4. load the new locations
	for (int i=0; i<nV_; ++i) {
		V[_nV+i] = V_[i];
		ids[_nV+i] = ids_[i]; 
	}

	// 5. sample pi
	for(int i=0; i<nV; ++i){
		pi[i] = i;
	}
	random_shuffle(pi.begin(), pi.end());
	for (int i=0; i<nV; ++i){
		reverse_pi[pi[i]] = i;
	}
}
