#include <bits/stdc++.h>
using namespace std;

#include "global_DC.h"
#include "HST_DC.h"
#ifdef WATCH_MEM
#include "monitor.h"
#endif

// #define LOCAL_DEBUG

const int MAX_HEIGHT = 55;


void Tree_t::initMemory(int n) {
	nV = n;
	expks = new double[MAX_HEIGHT];
	expks[0] = 1.0;
	for (int i=1; i<MAX_HEIGHT; ++i)
		expks[i] = expks[i-1] * alpha;
	
	sumks = new double[MAX_HEIGHT+1];
	sumks[0] = 0.0;
	for (int i=1; i<=MAX_HEIGHT; ++i)
		sumks[i] = sumks[i-1] + expks[i];
	
	logAlpha = log(alpha);
	distorUB.resize(nV);
	for (int i=0; i<nV; ++i) {
		distorUB[i].resize(MAX_HEIGHT);
		fill(distorUB[i].begin(), distorUB[i].end(), INF);
	}
	dataset.resize(nV);
	
	mark.resize(nV);
	leaves.resize(nV);
	fill(leaves.begin(), leaves.end(), (Node_t*)NULL);
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
	
	initMemory(nV);
	for(int i=0; i<nV; ++i){
		V[i] = V_[i];
		ids[i] = ids_[i];
		pi[i] = i;
	}
	
	RT_POINTS.resize(nV);
	for (int i=0; i<nV; ++i) {
		RT_POINTS[i] = rtree_point(V[i].x, V[i].y);
	}
	
	// generate the permutation pi
	random_shuffle(pi.begin(), pi.end());
	for (int i=0; i<nV; ++i){
		reverse_pi[pi[i]] = i;
	}
}

void Tree_t::calcDmax() {
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

void Tree_t::constructHST(clock_t startClock) {
	distor = 1.0, distor_ = 0.0, cdistor = 0.0;
	node_n = merge_n = 0;

	clock_t t0 = clock();
	
	res = maxDistor = 1.0;
	if (rt != NULL) {
		freeHST(rt);
		rt = NULL;
	}
	calcDmax();
#ifdef LOCAL_DEBUG
	printf("Phase 0: time = %.6lf\n", 1.0 * (clock()-t0) / CLOCKS_PER_SEC);
	fflush(stdout);
#endif
	
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
		int centerId = _getCenter_subspaceDistor_complex(g);
#ifdef LOCAL_DEBUG
		printf("[BEGIN] centerId = %d\n", centerId);
		fflush(stdout);
#endif
		if (centerId == -1) {
			pars.push_back(g);
		} else {
			_getSubspaces(nid, centerId, g, _pars);
		}
#ifdef LOCAL_DEBUG
		printf("[END] centerId = %d\n", centerId);
		fflush(stdout);
#endif
	}
#ifdef LOCAL_DEBUG
	printf("Phase I: time = %.6lf\n", 1.0 * (clock()-t0) / CLOCKS_PER_SEC);
	fflush(stdout);
#endif
	t0 = clock();
	
	while (!pars.empty()) {
		partition_t g = *pars.rbegin();
		pars.pop_back();
		int centerId = _getCenter_maxDistor_complex(g);
		_getParition(nid, centerId, g, pars);
	}
#ifdef LOCAL_DEBUG
	printf("Phase II: time = %.6lf\n", 1.0 * (clock()-t0) / CLOCKS_PER_SEC);
	fflush(stdout);
#endif
	
	
	// calculate distortion
	distor_ = (cdistor<=0.0) ? 1.0 : (distor_/cdistor);
	node_n = nid - merge_n;
	usedMemory = 1.0*node_n*sizeof(Node_t)/1024.0;
#ifdef LOCAL_DEBUG
	printf("nV=%d, nid=%d, H=%d, mem=%.2lf\n", nV, node_n, H, 1.0*node_n*sizeof(Node_t)/1024.0);
	fflush(stdout);	
#endif
}

void Tree_t::initParameters(int flag, int m) {
	static double pc = max(1.0, log(m)/log(2.0));
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

int Tree_t::_getCenter_subspaceDistor_complex(partition_t& g) {
	if (g.ids.size() <= num_points_UB) return -1;
	if (g.ids.size() <= PAGE_SIZE_N)  return g.ids[rand()%g.ids.size()];
	
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
			tmp = _calc_subspaceDistor_rtree0_UB2(centerId, g, mn);
			if (tmp < mn) {
				mn = tmp;
				ret = centerId;
			}
		}
		
		initParameters(0, sz);
		_centerId = ret;
		for (int i=0; i<sz; ++i) 
			mark[g.ids[i]] = -1;
		_calc_subspaceDistor_rtree0_UB2(_centerId, g, INF);
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
			tmp = _calc_subspaceDistor_rtree0_UB2(centerId, g, mn);
			if (tmp < mn) {
				mn = tmp;
				ret = centerId;
			}
		}
		
		res = max(res, mn);
		maxDistor = max(maxDistor, mn);
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
			tmp = _calc_subspaceDistor_rtree0_UB2(centerId, g, mn);
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
			_calc_subspaceDistor_rtree0_UB2(_centerId, g, INF);
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
				tmp = _calc_subspaceDistor_rtree0_UB2(centerId, g, mn);
				if (tmp < mn) {
					mn = tmp;
					ret = centerId;
				}
			}
		}
		
		res = max(res, mn);
		maxDistor = max(maxDistor, mn);
	}

	return ret;
}

void Tree_t::_getSubspaces(int& nid, int centerId, partition_t& g, vector<partition_t>& left) {
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

double Tree_t::_calc_subspaceDistor_rtree0_UB2(int centerId, partition_t& g, double UB) {
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
	
	
	rtree_t rtree;
	int nPts = 0, b = 0;
	
	for (int i=0,k=0; i<npos&&maxDistor<UB; i=k) {
		while (k<npos && levs[U[k]]==levs[U[i]]) ++k;
		double dt = distAtLevel(levs[U[i]]);
		
		if (k < npos) {
			nPts = 0;
			b = k;
			for (int j=k,c=0; j<npos; ++j) {
				if (levs[U[j]] != levs[U[k]]) break;
				int pid = g.ids[U[j]];
				dataset[c++] = std::make_pair(RT_POINTS[pid], pid);
				++nPts;
			}		
			if (nPts > PAGE_SIZE_N) {
				rtree.clear();
				rtree.insert(dataset.begin(), dataset.begin()+nPts);
			}
			
			vector<value> nns;
			for (int j=i; j<k&&maxDistor<UB; ++j) {
				int pid = g.ids[U[j]], qid = 0;
				if (nPts == 0) continue;
				
				if (nPts > PAGE_SIZE_N) {
					rtree.query(bgi::nearest(rtree_point(V[pid].x, V[pid].y), 1), back_inserter(nns));
					qid = nns.rbegin()->second;
				} else {
					int kid = rand()%nPts;
					qid = g.ids[U[b+kid]];
				}
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

double Tree_t::_calc_maxDistor_rtree0_UB(int centerId, partition_t& g, double UB) {
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

	rtree_t rtree;
	int nPts = 0, b = 0, e = 0;
	
	for (int i=0,k=0,preLev=H+1,tid=0; i<npos&&maxDistor<UB; i=k,++tid) {
		while (k<npos && levs[U[k]]==levs[U[i]]) ++k;
		double dt = distAtLevel(levs[U[i]]);
	
		if (i>0 && nPts>0 && i<k) {
			if (nPts > PAGE_SIZE_N) {
				rtree.clear();
				for (int j=b,c=0; j<e; ++j) {
					int pid = g.ids[U[j]];
					dataset[c++] = std::make_pair(RT_POINTS[pid], pid);
				}
				rtree.insert(dataset.begin(), dataset.begin()+nPts);
			}
		}
	
		vector<value> nns;
		for (int j=i; j<k&&maxDistor<UB; ++j) {
			int pid = g.ids[U[j]], qid = 0;
			if (i==0 || nPts==0) continue;
			
			if (nPts > PAGE_SIZE_N) {
				rtree.query(bgi::nearest(rtree_point(V[pid].x, V[pid].y), 1), back_inserter(nns));
				qid = nns.rbegin()->second;
			} else {
				int kid = rand()%nPts;
				qid = g.ids[U[b+kid]];
			}
			double d = dist(V, pid, qid);
			maxDistor = max(maxDistor, dt/d);
		}
		
		if (k > i) {
			nPts = k - i;
			b = i, e = k;
		}
		
		while (preLev >= levs[U[i]]) {
			distorUB[centerId][preLev] = min(distorUB[centerId][preLev], maxDistor);
			--preLev;
		}
	}	
	
	return maxDistor;
}

void Tree_t::_getParition(int& nid, int centerId, partition_t& g, vector<partition_t>& left) {
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

int Tree_t::_getCenter_maxDistor_complex(partition_t& g) {
	if (g.ids.size() <= PAGE_SIZE_N)  return g.ids[rand()%g.ids.size()];
	
	int ret = -1;
	double mn = INF, tmp;
	int sz = g.ids.size();
	
	if (sz <= 1) {
		ret = g.ids[0];
	} else if (sz == nV) {
		int cnt_points_SAMPLE = max(1.0, log(sz)/log(2));
		for (int i=0; i<sz&&cnt_points_SAMPLE>0; ++i,--cnt_points_SAMPLE) {
			tmp = _calc_maxDistor_rtree0_UB(g.ids[i], g, mn);
			if (tmp < mn) {
				mn = tmp;
				ret = g.ids[i];
			}
		}
		res = max(res, mn);
		maxDistor = max(maxDistor, mn);
	} else if (sz < nV) {
		int cnt_points_SAMPLE = max(1.0, log(sz)/log(2));
		
		vector<pair<double,int> > UBs(sz);
		
		for (int i=0; i<sz; ++i)
			UBs[i] = make_pair(distorUB[g.ids[i]][g.lev+1], g.ids[i]);
		sort(UBs.begin(), UBs.end());
		
		for (int i=0; mn>res&&i<sz&&cnt_points_SAMPLE>0; ++i,--cnt_points_SAMPLE) {
			tmp = _calc_maxDistor_rtree0_UB(UBs[i].second, g, mn);
			if (tmp < mn) {
				mn = tmp;
				ret = UBs[i].second;
			}
		}
		
		res = max(res, mn);
		maxDistor = max(maxDistor, mn);
	}

	return ret;
}

double Tree_t::distAtLevel(int lev) {
	if (lev >= H+1) return 0.0;
	if (lev <= 1) lev = 1;
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
	mark.resize(nV);
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
	distorUB.resize(nV);
	for (int i=0; i<nV; ++i) {
		distorUB[i].resize(MAX_HEIGHT);
		fill(distorUB[i].begin(), distorUB[i].end(), INF);
	}
	dataset.resize(nV);
	RT_POINTS.resize(nV);
	
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

void Tree_t::RMQ_clear() {
	RMQ_V.clear();
	RMQ_D.clear();
	RMQ_beg.clear();
	RMQ_dp.clear();
}

void Tree_t::RMQ_resize() {
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

void Tree_t::RMQ_dfs(Node_t* rt, int d) {
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

void Tree_t::RMQ_init(int n) {
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

Node_t* Tree_t::RMQ_query(int l, int r) {
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

double Tree_t::RMQ_distOnHST(Node_t* pu, Node_t* pv) {
	Node_t* lca = RMQ_query(RMQ_beg[pu->nid], RMQ_beg[pv->nid]);
	return distAtLevel(lca->lev);
}

pair<double,double> Tree_t::getDistortion_fast() {
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