#include "HST_DC.h"

const int HSF_SIZE = 5;

Tree_t bs_trees[HSF_SIZE];
int* mp;
int* rmp;
bool* deleted;
int* whichTree;
int* posTree;
int loc_n;
int tot_n;
int insert_n = 0;
double alpha;
vector<location_t> locs;
vector<location_t> clocs;
vector<int> cids;
pair<double,double> obj(1.0, 0.0);
clock_t start_time, end_time, last_time;
double totUsedTime = 0.0;

pair<double,double> getPartDistortion() {
	double d, dt, rat;
	double ret = 1.0+EPS;
	double ret_ = 0.0;
	double cnt = 0.0;
	int u, v, uTreeId, vTreeId, uTreePos, vTreePos;
	pair<int,double> nnPair;
	int nnId;
	double nnDis;
	
	bool checkFlag = true;
	for (int i=0; i<cids.size(); ++i) {
		u = cids[i], uTreeId = whichTree[u], uTreePos = posTree[u];
		Tree_t& uTree = bs_trees[uTreeId];

		for (int j=0; j<locs.size()-clocs.size()+i; ++j) {
			v = rmp[j];
//			if (deleted[v])
//				continue;
			vTreeId = whichTree[v], vTreePos = posTree[v];
			Tree_t& vTree = bs_trees[vTreeId];

			d = dist(clocs[i], locs[j]);
			if (vTreeId < uTreeId) {
				nnPair = uTree.nns[vTreeId][uTreePos];
				nnId = nnPair.first, nnDis = nnPair.second;
				dt = vTree.distOnHST(vTreePos, nnId) + nnDis;
			} else if (vTreeId > uTreeId) {
				// the NN of vTreeId may be invalid since uTree has been updated.
				vTree.updateNeighbour(vTreePos, uTreeId, uTree);
				nnPair = vTree.nns[uTreeId][vTreePos];
				nnId = nnPair.first, nnDis = nnPair.second;
				dt = uTree.distOnHST(uTreePos, nnId) + nnDis;
			} else {// if (vTreeId == uTreeId)
				dt = uTree.distOnHST(uTreePos, vTreePos);
			}
			if (d > 0) {
				rat = max(dt/d, 1.0);
				ret = max(ret, rat);
				ret_ += rat;
				cnt += 1.0;
			}
			
			if (dt < d) {
				checkFlag = false;
			}
		}
	}
	if (cnt <= 0) {
		ret_ = 1.0;
	} else {
		ret_ /= cnt;
	}
	
	return make_pair(ret, ret_);	
}

int getNextTree() {
	int dmin = INT_MAX, ret = 1;
	
	for (int i=1; i<HSF_SIZE; ++i) {
		if (bs_trees[i].nV-bs_trees[i].dV < dmin) {
			dmin = bs_trees[i].nV-bs_trees[i].dV;
			ret = i;
		}
	}
	
	return ret;
}

void rebuild(int treeId) {
	last_time = clock();
	
	Tree_t& bs_tree = bs_trees[treeId];
	
	if (bs_tree.nV == 0) {// if the tree is null
		bs_tree.initLocation(clocs.size(), clocs, cids, alpha);
		bs_tree.constructHST(start_time);
		
		for (int i=cids.size()-1; i>=0; --i) {
			whichTree[cids[i]] = treeId;
			posTree[cids[i]] = i;
		}

		for (int i=0; i<treeId; ++i) {
			bs_tree.initNeighbours(i, bs_trees[i]);
		}
	} else {// we need to rebuild the tree
		int _nV = bs_tree.nV - bs_tree.dV;
		bs_tree.loadLocation(clocs.size(), clocs, cids, alpha);
		bs_tree.constructHST(start_time);

		for (int i=0; i<_nV; ++i) {
			// whichTree[bs_tree.ids[i]] = treeId;
			posTree[bs_tree.ids[i]] = i;			
		}
		for (int i=cids.size()-1; i>=0; --i) {
			whichTree[cids[i]] = treeId;
			posTree[cids[i]] = i + _nV;
		}
		for (int i=0; i<treeId; ++i) {
			bs_tree.updateNeighbours(i, bs_trees[i]);
		}
	}
	
	end_time = clock();
	pair<double,double> tmp = getPartDistortion();
	obj.first = max(obj.first, tmp.first);
	obj.second += tmp.first;
	
	double timeCost = (double)(end_time-last_time)/CLOCKS_PER_SEC;
	long long usedMemory_ = 0;
	for (int i=0; i<HSF_SIZE; ++i) {
		usedMemory_ += bs_trees[i].memoryCost();
	}
	printf("HSFdc %.4lf %.4lf %.4lf %lld\n", tmp.first, tmp.second, timeCost, usedMemory_);
	fflush(stdout);
	totUsedTime += timeCost;
	insert_n += 1;
}

inline void insert_loc(int id, location_t& loc) {
	mp[id] = locs.size();
	rmp[locs.size()] = id;
	cids.push_back( id );
	clocs.push_back( loc );
	locs.push_back( loc );
}

void delete_loc(int id) {
	int end = locs.size()-1, now = mp[id];
	int eid = rmp[end];
	
	// we can safely ignore ``now''
	locs[now] = locs[end];
	mp[eid] = now;
	rmp[now] = eid;
	locs.pop_back();
	deleted[id] = true;
	
	// we need to delete the point from HST
	int treeId = whichTree[id], treePos = posTree[id];
	Tree_t& bs_tree = bs_trees[treeId];
	bs_tree.deleteOne(treePos);
}

void static_input() {
	char op[10];
	int num;
	int id;
	location_t loc;
	
	scanf("%s %d", op, &num);
	for (int i = 0; i < num; ++i) {
		scanf("%d %lf %lf", &id, &loc.x, &loc.y);
		insert_loc(id, loc);
	}
	
	last_time = clock();
	Tree_t& bs_tree = bs_trees[0];
	// printf("[BEGIN] initLocation\n"); fflush(stdout);
	bs_tree.initLocation(clocs.size(), clocs, cids, alpha);
	// printf("[END] initLocation\n"); fflush(stdout);
	// update posTree and whichTree
	for (int i=cids.size()-1; i>=0; --i) {
		whichTree[cids[i]] = 0;
		posTree[cids[i]] = i;
	}
	// printf("[BEGIN] constructHST\n"); fflush(stdout);
	bs_tree.constructHST(start_time);
	// printf("[END] constructHST\n"); fflush(stdout);

	end_time = clock();
	pair<double,double> tmp = make_pair(bs_tree.distor, bs_tree.distor_);
	obj.first = max(obj.first, tmp.first);
	obj.second += tmp.first;
	
	double timeCost = (double)(end_time-last_time)/CLOCKS_PER_SEC;
	long long usedMemory_ = bs_tree.memoryCost();
	printf("HSFdc %.4lf %.4lf %.4lf %lld\n", tmp.first, tmp.second, timeCost, usedMemory_);
	fflush(stdout);
	totUsedTime += timeCost;
	insert_n += 1;
}

void dec_or_inc() {
	char op[10];
	int num;
	location_t loc;
	int id;
	
	scanf("%s %d",op,&num);
	if (op[0] == 'I') {
		clocs.clear();
		cids.clear();
		int treeId = getNextTree();
		for (int i = 0; i < num; ++i) {
			scanf("%d %lf %lf", &id, &loc.x, &loc.y);
			insert_loc(id, loc);
		}
		rebuild(treeId);
	}
	
	if (op[0] == 'D') {
		for (int i = 0; i < num; ++i) {
			scanf("%d", &id);
			delete_loc(id);
		}
	}
}

int main(int argc, char **argv) {
	string srcFileName;
	string desFileName;
	int nop;
	
	if (argc > 1)
		srcFileName = string(argv[1]);
	else
		srcFileName = "data_0.txt";
	if (argc > 2)
		desFileName = string(argv[2]);

	freopen(srcFileName.c_str(), "r", stdin);

	start_time = clock();
	
	scanf("%d %d %lf", &tot_n, &nop, &alpha);
	mp = new int[tot_n+5];
	rmp = new int[tot_n+5];
	deleted = new bool[tot_n+5];
	whichTree = new int[tot_n+5];
	posTree = new int[tot_n+5];
	memset(deleted, false, sizeof(bool)*(tot_n+5));
	
	static_input();
	
	/*
	* rebuild can only happen in the insert part
	* So I remove it to inside the dec_or_inc
	*/

	// rebuild(alpha,time_limit);
	for (int i = 2; i <= nop; ++i) { // the first op is static input...
		dec_or_inc();
		// rebuild(alpha,time_limit);
	}

#ifdef WATCH_MEM
    watchSolutionOnce(getpid(), usedMemory);
#endif
	printf("HSFdc %.4lf %.4lf %.4lf %lld\n", obj.first, obj.second/insert_n, totUsedTime, usedMemory);
	fclose(stdout);
	
	for (int i=0; i<HSF_SIZE; ++i) {
		bs_trees[i].freeMemory();
	}
	delete[] mp;
	delete[] rmp;
	delete[] deleted;
	delete[] whichTree;
	delete[] posTree;
	
	return 0;
}
