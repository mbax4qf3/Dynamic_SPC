#include "u_spc.h"

#include <algorithm>
#include <chrono>
#include <cinttypes>
#include <cstdint>
#include <functional>
#include <iostream>
#include <fstream>
#include <memory>
#include <numeric>
#include <queue>
#include <tuple>
#include <set>
#include <map>
#include <unordered_set>
#include <vector>
#include <string>
#include <thread>
#include <unistd.h>
#include <iomanip>
#include <omp.h>

#include "progressbar.h"
#include "macros.h"
#include "two_layer_queue.h"

namespace spc {

/*
*******************
****** Basic ******
*******************
*/

std::pair<uint32_t, uint64_t> USPCBasic::BFS_SPC(const Graph& const_graph, uint32_t s, uint32_t t) {
	std::vector<uint32_t> D(const_graph.size(), UINT32_MAX);
	std::vector<uint32_t> C(const_graph.size(), 0);

	std::queue<uint32_t> Q; Q.push(s);

	D[s] = 0; C[s] = 1;

	while (!Q.empty()) {
		uint32_t v = Q.front(); Q.pop();
		if (v == t) return std::make_pair(D[v], C[v]);

		for (const uint32_t w:const_graph[v]) {
			if (D[w] == UINT32_MAX) {
				D[w] = D[v] + 1;
				C[w] = C[v];
				Q.push(w);
			} else if (D[w] == D[v] + 1) {
				C[w] += C[v];
			}
		}
	}

	return std::make_pair(0,0);
}



/*
*******************
****** Index ******
*******************
*/

// construct an index for the graph "const_graph"
void USPCIndex::BuildIndex(const Graph& const_graph) {
	ASSERT(dL_.empty() && cL_.empty() && G_.empty());
	G_ = const_graph;
	n_ = G_.size();

	// initialization
	dL_.resize(n_);
	cL_.resize(n_);

	order_.resize(n_);
	rank_.resize(n_);

	// ordering
	(this->*of_[os_])(G_);
	OrderRank();
	
	// some auxiliary structures
	std::vector<uint32_t> dLu(n_, UINT32_MAX);
	std::vector<uint32_t> D(n_, UINT32_MAX);
	std::vector<uint32_t> C(n_, 0);

	progressbar bar(n_);

	for (size_t i = 0; i < n_; ++i) {
	bar.update();
	const uint32_t u = order_[i];

	// for fast distance computation
	for (const auto e : dL_[u]) dLu[LEExtractV(e)] = LEExtractD(e);

	// bfs
	std::vector<uint32_t> reset({u});
	std::queue<uint32_t> Q({u});
	D[u] = 0; C[u] = 1;

	while (!Q.empty()) {

		const uint32_t v = Q.front(); Q.pop();
		const uint32_t dSoFar = Distance(dLu, dL_[v]);
		if (D[v] > dSoFar) continue;

		// add a corresponding entry
		NormalD(D[v]); NormalC(C[v]);
		(D[v] < dSoFar? dL_[v] : cL_[v]).push_back(LEMerge(u, D[v], C[v]));

		// correct C[v]
		if (unlikely(kUBC <= C[v])) C[v] = kUBC;
		for (const uint32_t w : G_[v]) {
			if (rank_[w] <= rank_[u]) continue;
			if (UINT32_MAX == D[w]) {
				D[w] = D[v] + 1;
				C[w] = C[v];
				Q.push(w);
				reset.push_back(w);
			} else if (D[w] == D[v] + 1) {
				if (likely(kUBC - C[v] >= C[w])) C[w] += C[v];
				else C[w] = kUBC;
			}
		}
	}

	// clear
	for (const uint32_t v : reset) {
		D[v] = UINT32_MAX; C[v] = 0;
	}

	for (const auto e : dL_[u]) dLu[LEExtractV(e)] = UINT32_MAX;

	}
}


// merge dL_ and cL_ to cL_ at index stage
void USPCIndex::IndexMerge(){
	for (uint32_t i = 0; i < n_; ++i) {
	// merge
	std::vector<LabelEntry> mL;
	mL.reserve(dL_[i].size() + cL_[i].size());
	size_t di = 0, ci = 0;
	while (di < dL_[i].size() && ci < cL_[i].size()) {
		if (rank_[LEExtractV(dL_[i][di])] < rank_[LEExtractV(cL_[i][ci])]) {
			mL.push_back(dL_[i][di++]);
		} else {
			mL.push_back(cL_[i][ci++]);
		}
	}
	while (di < dL_[i].size()) mL.push_back(dL_[i][di++]);
	while (ci < cL_[i].size()) mL.push_back(cL_[i][ci++]);
		dL_[i].clear();
		cL_[i] = mL;
	}
	decltype(dL_)().swap(dL_);
}

 
// index write with graph, dL_, cL_, inverted label, and order 
// cL_ and dL_ are not merged
uint64_t USPCIndex::IndexWrite(const std::string& filename) {
	ASSERT(0 != n_);
	FILE* file = fopen(filename.c_str(), "wb");
	fwrite(&n_, sizeof(n_), 1, file);
	for (uint32_t u = 0; u < n_; ++u) {
		const uint32_t s = G_[u].size();
		fwrite(&s, sizeof(s), 1, file);
		fwrite(G_[u].data(), sizeof(G_[u].back()), s, file);
	}
	uint64_t num_labels = 0;
	for (uint32_t i = 0; i < n_; ++i) {
		// write canonical labels
		const uint32_t d_s = dL_[i].size();
		fwrite(&d_s, sizeof(d_s), 1, file);
		fwrite(dL_[i].data(), sizeof(dL_[i].back()), d_s, file);
		num_labels += d_s;

		// write non_canonical labels
		const uint32_t c_s = cL_[i].size();
		fwrite(&c_s, sizeof(c_s), 1, file);
		fwrite(cL_[i].data(), sizeof(cL_[i].back()), c_s, file);
		num_labels += c_s;

	}
	printf("total # of label entries: %" PRIu64 "\n", num_labels);

	// order information
	fwrite(order_.data(), sizeof(order_.back()), n_, file);

	fclose(file);
	return num_labels;
}

// fast computation
uint32_t USPCIndex::Distance(const std::vector<uint32_t>& dLu,
							 const std::vector<LabelEntry>& dLv) const {
	uint32_t d = UINT32_MAX;
	for (const auto e : dLv) {
		const uint32_t v = LEExtractV(e);
		if (UINT32_MAX == dLu[v]) continue;
		const uint32_t dd = dLu[v] + LEExtractD(e);
		if (dd < d) d = dd;
	}
	return d;
}

// degree order
void USPCIndex::DegreeOrder(const Graph& graph) {
	// std::cout << "degree-based ordering" << std::endl;
	std::vector<uint32_t> deg(n_);
	for (uint32_t i = 0; i < n_; ++i) {
		deg[i] = graph[i].size();
		order_[i] = i;
	}
	std::sort(order_.begin(), order_.end(),
			[&deg](const uint32_t v1, const uint32_t v2) {
				return deg[v1] > deg[v2];
			});
}


/*
*******************
****** Query ******
*******************
*/


// Query of Dis and Cnt
std::pair<uint32_t, uint64_t> USPCQuery::Count(uint32_t v1, uint32_t v2) const {
	ASSERT(v1 != v2);
	// count the # of shortest paths
	uint32_t sp_d = UINT32_MAX;
	uint64_t sp_c = 0;

	size_t p1 = 0, p2 = 0;
	while (p1 < cL_[v1].size() && p2 < cL_[v2].size()) {
		const uint32_t w1 = LEExtractV(cL_[v1][p1]);
		const uint32_t w2 = LEExtractV(cL_[v2][p2]);
		if (rank_[w1] < rank_[w2]) ++p1;
		else if (rank_[w1] > rank_[w2]) ++p2;
		else {
		const uint32_t d = LEExtractD(cL_[v1][p1]) +
			LEExtractD(cL_[v2][p2]);
		if (d < sp_d) {
			sp_d = d;
			sp_c = static_cast<uint64_t>(LEExtractC(cL_[v1][p1])) *
				LEExtractC(cL_[v2][p2]);
		} else if (d == sp_d) {
			uint64_t c = static_cast<uint64_t>(LEExtractC(cL_[v1][p1])) *
				LEExtractC(cL_[v2][p2]);
			sp_c += c;
		}
		++p1; ++p2;
		}
	}

	if (sp_d == UINT32_MAX || sp_c == 0) sp_d = 0;
	return std::make_pair(sp_d, sp_c);
}


// BiBFS Dis and Cnt
std::pair<uint32_t, uint64_t> USPCQuery::bi_BFS_Count(Graph& graph, uint32_t v1, uint32_t v2) {

	uint32_t sp_d = UINT32_MAX;
	uint64_t sp_c = 0;

    // bidirectional BFS
    std::vector<TwoLayerQueue> qque;
    std::vector<uint32_t> qdist[2];
    std::vector<uint32_t> qcnt[2];

    qdist[0].resize(n_, UINT32_MAX); qcnt[0].resize(n_, 0);
    qdist[1].resize(n_, UINT32_MAX); qcnt[1].resize(n_, 0);
    
    qque.push_back(TwoLayerQueue(n_)); qque.push_back(TwoLayerQueue(n_));
    
    int found_flag = 0;

    for (int dir = 0; dir < 2; dir++){
        int v = dir == 0 ? v1 : v2;
        qque[dir].clear();
        qque[dir].push(v);
        qque[dir].next();
        qdist[dir][v] = 0;
        qcnt[dir][v] = 1;
    }

    uint32_t res = UINT32_MAX, dis[2] = {0, 0};
    while (!qque[0].empty() && !qque[1].empty()) {
        int use = 0;
        use = (qque[0].size() <= qque[1].size()) ? 0 : 1;
        dis[use]++;

        while (!qque[use].empty()) {
            int v = qque[use].front();
            qque[use].pop();

			for (int w : graph[v]) {

				uint32_t &dst_d = qdist[1 - use][w];

				// first met, or meet again, accumulate counting
				if (qdist[use][w] == UINT32_MAX) {

					qdist[use][w] = qdist[use][v] + 1;
					qcnt[use][w] = qcnt[use][v];
					if (!found_flag)
						qque[use].push(w);

				} else if (qdist[use][w] == qdist[use][v] + 1) {

					if (!found_flag) 
						qcnt[use][w] += qcnt[use][v]; 
					else 
						qcnt[use][w] = qcnt[use][v];

				}

				if (dst_d != UINT32_MAX) {
					sp_d = qdist[use][w] + dst_d;
					sp_c += (qcnt[use][w] * qcnt[1 - use][w]);
					found_flag = 1;
				}
			}
        }

        if (found_flag) goto LOOP_END;
        qque[use].next();

    }
    LOOP_END:

	if (sp_d == UINT32_MAX || sp_c == 0) sp_d = 0;
	return std::make_pair(sp_d, sp_c);
	
}


// cL_ and dL_ are not merged when reading
void USPCQuery::IndexRead(const std::string& filename) {
	ASSERT(dL_.empty() && cL_.empty() && G_.empty()); // && inv_L_.empty()
	FILE* file = fopen(filename.c_str(), "rb");
	ASSERT(fread(&n_, sizeof(n_), 1, file) == 1);
	G_.resize(n_);
	// read graph
	for (uint32_t u = 0; u < n_; ++u) {
		uint32_t s = 0;
		ASSERT(fread(&s, sizeof(s), 1, file) == 1);
		G_[u].resize(s);
		ASSERT(fread(G_[u].data(), sizeof(G_[u].back()), s, file) == s);
	}
	// initialization
	dL_.resize(n_);
	cL_.resize(n_);

	// read labels
	uint64_t num_labels = 0;
	uint64_t num_dlabels = 0;
	uint64_t num_clabels = 0;
	uint32_t d_s, c_s, i_s;
	for (uint32_t i = 0; i < n_; ++i) {
		// read canonical labels
		ASSERT(fread(&d_s, sizeof(d_s), 1, file) == 1);
		dL_[i].resize(d_s);
		ASSERT(fread(dL_[i].data(), sizeof(dL_[i].back()), d_s, file) == d_s);
		num_labels += d_s;
		num_dlabels += d_s;

		// read non-canonical labels
		ASSERT(fread(&c_s, sizeof(c_s), 1, file) == 1);
		cL_[i].resize(c_s);
		ASSERT(fread(cL_[i].data(), sizeof(cL_[i].back()), c_s, file) == c_s);

		num_labels += c_s;
		num_clabels += c_s;
	}
	printf("total # of labels entries:\t%" PRIu64 "\n", num_labels);
	printf("total # of can-labels	:\t%" PRIu64 "\n", num_dlabels);
	printf("total # of non-can-labels:\t%" PRIu64 "\n", num_clabels);

	// order information
	order_.resize(n_);
	ASSERT(fread(order_.data(), sizeof(order_.back()), n_, file) == n_);
	rank_.resize(n_);
	OrderRank();
	// check
	for (uint32_t i = 0; i < n_; ++i) {
		for (size_t j = 1; j < dL_[i].size(); ++j) {
			ASSERT(rank_[LEExtractV(dL_[i][j])] >
					rank_[LEExtractV(dL_[i][j - 1])]);
		}
		for (size_t j = 1; j < cL_[i].size(); ++j) {
			ASSERT(rank_[LEExtractV(cL_[i][j])] >
					rank_[LEExtractV(cL_[i][j - 1])]);
		}
	}
	
	// reorganize the labels
	for (uint32_t i = 0; i < n_; ++i) {
		// merge
		std::vector<LabelEntry> mL;
		mL.reserve(dL_[i].size() + cL_[i].size());
		size_t di = 0, ci = 0;
		while (di < dL_[i].size() && ci < cL_[i].size()) {
			if (rank_[LEExtractV(dL_[i][di])] < rank_[LEExtractV(cL_[i][ci])]) {
				mL.push_back(dL_[i][di++]);
			} else {
				mL.push_back(cL_[i][ci++]);
			}
		}
		while (di < dL_[i].size()) mL.push_back(dL_[i][di++]);
		while (ci < cL_[i].size()) mL.push_back(cL_[i][ci++]);
		dL_[i].clear();
		cL_[i] = mL;
	}
	decltype(dL_)().swap(dL_);
	printf("labels merged; ");

	// check
	bool check = false;
	ASSERT(fread(&check, sizeof(check), 1, file) == 0);
	fclose(file);

	std::cout << "index read." << std::endl;
}

// cL_ and dL_ are merged
void USPCQuery::IndexRead_UPD(const std::string& filename) {
	ASSERT(cL_.empty() && G_.empty());
	FILE* file = fopen(filename.c_str(), "rb");
	ASSERT(fread(&n_, sizeof(n_), 1, file) == 1);
	G_.resize(n_);

	// read graph
	for (uint32_t u = 0; u < n_; ++u) {
		uint32_t s = 0;
		ASSERT(fread(&s, sizeof(s), 1, file) == 1);
		G_[u].resize(s);
		ASSERT(fread(G_[u].data(), sizeof(G_[u].back()), s, file) == s);
	}

	// initialization
	cL_.resize(n_);

	// read labels
	uint64_t num_labels = 0;
	uint32_t c_s, i_s;
	for (uint32_t i = 0; i < n_; ++i) {
		ASSERT(fread(&c_s, sizeof(c_s), 1, file) == 1);
		cL_[i].resize(c_s);
		ASSERT(fread(cL_[i].data(), sizeof(cL_[i].back()), c_s, file) == c_s);
		num_labels += c_s;

	}

	printf("total # of labels entries:\t%" PRIu64 "\n", num_labels);

	// order information
	order_.resize(n_);
	ASSERT(fread(order_.data(), sizeof(order_.back()), n_, file) == n_);
	rank_.resize(n_);
	OrderRank();

	fclose(file);
}

// update graph
void USPCQuery::UpdateGraph(Graph& graph, uint32_t v1, uint32_t v2, char upd_type) {
	ASSERT(v1 != v2);

	if (upd_type == 'i') {
		graph[v1].push_back(v2);
		graph[v2].push_back(v1);
	} else if (upd_type == 'd') {
		for (auto i = 0; i < graph[v1].size(); ++i) {
			if (graph[v1][i] == v2) {
				graph[v1].erase(graph[v1].begin() + i);
				break;
			}
		}

		for (auto i = 0; i < graph[v2].size(); ++i) {
			if (graph[v2][i] == v1) {
				graph[v2].erase(graph[v2].begin() + i);
				break;
			}
		}
	}
}

void USPCQuery::print_Label(uint32_t v) {
	std::ofstream lfile;
	std::string lfilename = "tmp_label/" + std::to_string(v) + "L.txt";
	lfile.open(lfilename.c_str());
	for (auto l:cL_[v]) {
		lfile << rank_[LEExtractV(l)] << ": " << LEExtractV(l) << "	" << LEExtractD(l) << "	" << LEExtractC(l) << "\n";
	}
	lfile.close();
}



/*
*******************
****** Update *****
*******************
*/

// cL_ and dL_ are not merged
void USPCUpdate::IndexRead(const std::string& filename) {
	ASSERT(dL_.empty() && cL_.empty() && G_.empty());
	FILE* file = fopen(filename.c_str(), "rb");
	ASSERT(fread(&n_, sizeof(n_), 1, file) == 1);
	G_.resize(n_);
	// read graph
	for (uint32_t u = 0; u < n_; ++u) {
		uint32_t s = 0;
		ASSERT(fread(&s, sizeof(s), 1, file) == 1);
		G_[u].resize(s);
		ASSERT(fread(G_[u].data(), sizeof(G_[u].back()), s, file) == s);
	}

	// initialization
	dL_.resize(n_);
	cL_.resize(n_);

	// read labels
	uint64_t num_labels = 0;
	uint64_t num_dlabels = 0;
	uint64_t num_clabels = 0;
	uint32_t d_s, c_s, i_s;
	for (uint32_t i = 0; i < n_; ++i) {
		// read canonical labels
		ASSERT(fread(&d_s, sizeof(d_s), 1, file) == 1);
		dL_[i].resize(d_s);
		ASSERT(fread(dL_[i].data(), sizeof(dL_[i].back()), d_s, file) == d_s);

		num_labels += d_s;
		num_dlabels += d_s;

		// read non-canonical labels
		ASSERT(fread(&c_s, sizeof(c_s), 1, file) == 1);
		cL_[i].resize(c_s);
		ASSERT(fread(cL_[i].data(), sizeof(cL_[i].back()), c_s, file) == c_s);
		
		num_labels += c_s;
		num_clabels += c_s;
	}

	printf("total # of labels entries:\t%" PRIu64 "\n", num_labels);
	printf("total # of can-labels	:\t%" PRIu64 "\n", num_dlabels);
	printf("total # of non-can-labels:\t%" PRIu64 "\n", num_clabels);

	// order information
	order_.resize(n_);
	ASSERT(fread(order_.data(), sizeof(order_.back()), n_, file) == n_);
	rank_.resize(n_);
	OrderRank();
	// check
	for (uint32_t i = 0; i < n_; ++i) {
		for (size_t j = 1; j < dL_[i].size(); ++j) {
			ASSERT(rank_[LEExtractV(dL_[i][j])] >
				rank_[LEExtractV(dL_[i][j - 1])]);
		}
		for (size_t j = 1; j < cL_[i].size(); ++j) {
			ASSERT(rank_[LEExtractV(cL_[i][j])] >
				rank_[LEExtractV(cL_[i][j - 1])]);
		}
	}

	// reorganize the labels
	for (uint32_t i = 0; i < n_; ++i) {
		// merge
		std::vector<LabelEntry> mL;
		mL.reserve(dL_[i].size() + cL_[i].size());
		size_t di = 0, ci = 0;
		while (di < dL_[i].size() && ci < cL_[i].size()) {
			if (rank_[LEExtractV(dL_[i][di])] < rank_[LEExtractV(cL_[i][ci])]) {
				mL.push_back(dL_[i][di++]);
			} else {
				mL.push_back(cL_[i][ci++]);
			}
		}
		while (di < dL_[i].size()) mL.push_back(dL_[i][di++]);
		while (ci < cL_[i].size()) mL.push_back(cL_[i][ci++]);
		dL_[i].clear();
		cL_[i] = mL;
	}
	decltype(dL_)().swap(dL_);
    
	printf("labels merged; ");

	// check
	bool check = false;
	ASSERT(fread(&check, sizeof(check), 1, file) == 0);
	fclose(file);

	std::cout << "index read." << std::endl;
}



// cL_ and dL_ are merged
uint64_t USPCUpdate::IndexWrite(const std::string& filename) {
	ASSERT(0 != n_);
	FILE* file = fopen(filename.c_str(), "wb");
	fwrite(&n_, sizeof(n_), 1, file);
	for (uint32_t u = 0; u < n_; ++u) {
		const uint32_t s = G_[u].size();
		fwrite(&s, sizeof(s), 1, file);
		fwrite(G_[u].data(), sizeof(G_[u].back()), s, file);
	}
	
	// write label
	uint64_t num_labels = 0;
	for (uint32_t i = 0; i < n_; ++i) {
		// write labels
		const uint32_t c_s = cL_[i].size();
		fwrite(&c_s, sizeof(c_s), 1, file);
		fwrite(cL_[i].data(), sizeof(cL_[i].back()), c_s, file);

		num_labels += c_s;
	}

	// order information
	fwrite(order_.data(), sizeof(order_.back()), n_, file);
	
	fclose(file);
	return num_labels;
}

// Query 
std::pair<uint32_t, uint64_t> USPCUpdate::Count(uint32_t v1, uint32_t v2) const {
	size_t p1 = 0, p2 = 0;
	uint32_t sp_d = UINT32_MAX;
	uint64_t sp_c = 0;
	while (p1 < cL_[v1].size() && p2 < cL_[v2].size()) {
		const uint32_t w1 = LEExtractV(cL_[v1][p1]);
		const uint32_t w2 = LEExtractV(cL_[v2][p2]);
		if (rank_[w1] < rank_[w2]) ++p1;
		else if (rank_[w1] > rank_[w2]) ++p2;
		else {
		const uint32_t d = LEExtractD(cL_[v1][p1]) +
			LEExtractD(cL_[v2][p2]);
		if (d < sp_d) {
			sp_d = d;
			sp_c = static_cast<uint64_t>(LEExtractC(cL_[v1][p1])) *
				LEExtractC(cL_[v2][p2]);
		} else if (d == sp_d) {
			uint64_t c = static_cast<uint64_t>(LEExtractC(cL_[v1][p1])) *
				LEExtractC(cL_[v2][p2]);
			sp_c += c;
		}
		++p1; ++p2;
		}
	}
	if (sp_d == UINT32_MAX) sp_d = 0;
	return std::make_pair(sp_d, sp_c);
}

/*
**************************
****Incremental update****
**************************
*/

// incremental update
std::tuple<uint32_t, size_t, uint32_t, size_t, uint32_t, uint32_t, uint32_t> USPCUpdate::Inc_SPC(uint32_t a, uint32_t b) {
	G_[a].push_back(b);
	G_[b].push_back(a);

	size_t aff_size = cL_[a].size() + cL_[b].size();
	std::vector<LabelEntry> aff_Labels;
	std::vector<int> aff_ab(aff_size, -1);

	size_t aL = 0, bL = 0;
	while (aL < cL_[a].size() && bL < cL_[b].size()) {
		auto aH = LEExtractV(cL_[a][aL]);
		auto bH = LEExtractV(cL_[b][bL]);
		if (rank_[aH] == rank_[bH]) {
			aff_Labels.push_back(cL_[a][aL]);
			aff_ab[aL+bL] = 0;
			++aL;
			aff_Labels.push_back(cL_[b][bL]);
			aff_ab[aL+bL] = 1;
			++bL;
		} else if (rank_[aH] < rank_[bH]) {
			aff_Labels.push_back(cL_[a][aL]);
			aff_ab[aL+bL] = 0;
			++aL;
		} else {
			aff_Labels.push_back(cL_[b][bL]);
			aff_ab[aL+bL] = 1;
			++bL;
		}
	}

	while (aL < cL_[a].size()) {
		aff_Labels.push_back(cL_[a][aL]);
		aff_ab[aL+bL] = 0;
		++aL;
	}
	while (bL < cL_[b].size()) {
		aff_Labels.push_back(cL_[b][bL]);
		aff_ab[aL+bL] = 1;
		++bL;
	}

	uint32_t renew_c = 0;
    uint32_t renew_d = 0;
	uint32_t new_ = 0;

	for (size_t i = 0; i < aff_size; ++i) {
		if (aff_ab[i] == 0 && rank_[LEExtractV(aff_Labels[i])] < rank_[b]) {
			auto ab_ = Inc_BFS(LEExtractV(aff_Labels[i]), b, LEExtractD(aff_Labels[i])+1, LEExtractC(aff_Labels[i]));
			renew_c += std::get<0>(ab_); renew_d += std::get<1>(ab_); new_ += std::get<2>(ab_);
		} else if (aff_ab[i] == 1 && rank_[LEExtractV(aff_Labels[i])] < rank_[a]) {
			auto ba_ = Inc_BFS(LEExtractV(aff_Labels[i]), a, LEExtractD(aff_Labels[i])+1, LEExtractC(aff_Labels[i]));
			renew_c += std::get<0>(ba_); renew_d += std::get<1>(ba_); new_ += std::get<2>(ba_);
		}
	}

    return std::make_tuple(a, cL_[a].size(), b, cL_[b].size(), renew_c, renew_d, new_);
}
 
// process of incremental update
std::tuple<uint32_t, uint32_t, uint32_t> USPCUpdate::Inc_BFS(uint32_t hub, uint32_t ab, uint32_t d, uint64_t c) {
	uint32_t num_add = 0;
	uint32_t num_renewc = 0;
    uint32_t num_renewd = 0;
	std::vector<uint32_t> D(n_, UINT32_MAX);
	std::vector<uint64_t> C(n_, 0);
	std::vector<uint32_t> hash_dist(n_, UINT32_MAX);
	std::queue<uint32_t> Q;
	Q.push(ab);
	D[ab] = d; C[ab] = c;

	for (const auto e : cL_[hub]) hash_dist[LEExtractV(e)] = LEExtractD(e);

	while (!Q.empty()) {
		auto v = Q.front(); Q.pop();
		auto previous = Distance(hash_dist, cL_[v], hub);

		uint64_t CC = C[v];
		if (D[v] > previous.first) continue;
		if (D[v] == previous.first && previous.first == LEExtractD(cL_[v][previous.second]) && LEExtractV(cL_[v][previous.second]) == hub) 
			CC += LEExtractC(cL_[v][previous.second]);

		if (LEExtractV(cL_[v][previous.second]) == hub) {

            if (LEExtractD(cL_[v][previous.second]) == D[v]) {
                num_renewc += 1;
            } else {
                num_renewd += 1;
            }

			cL_[v][previous.second] = LEMerge(hub,D[v],CC);

		} else {
			cL_[v].emplace(cL_[v].begin() + previous.second, LEMerge(hub,D[v],CC));
			num_add += 1;
		}

		for (auto nbr:G_[v]) {
			if (rank_[nbr] <= rank_[hub]) continue;
			if (D[nbr] == UINT32_MAX) {
				D[nbr] = D[v] + 1; C[nbr] = C[v];
				Q.push(nbr);
			} else if (D[nbr] == D[v] + 1) {
				C[nbr] += C[v];
			}
		}
	} // while

	return std::make_tuple(num_renewc,num_renewd,num_add);
}

// Calculate distance and also return the position of hub
std::pair<uint32_t, size_t> USPCUpdate::Distance(const std::vector<uint32_t>& dLu,
	const std::vector<LabelEntry>& dLv, uint32_t hub) const {
	size_t pos = 0;
	uint64_t d = UINT64_MAX;
	for (pos = 0; pos < dLv.size(); ++pos) {
		auto v = LEExtractV(dLv[pos]);
		if (rank_[v] > rank_[hub]) break;

		if (dLu[v] == UINT32_MAX) continue;

		auto dd = dLu[v] + LEExtractD(dLv[pos]);
		if (dd < d) d = dd;

		if (rank_[v] == rank_[hub]) break;
	}
	return std::make_pair(d, pos);
}
 
// Fast calculate distance
uint32_t USPCUpdate::FastDistance(const std::vector<uint32_t>& dLu,
	const std::vector<LabelEntry>& dLv) const {
		uint32_t d = UINT32_MAX;
		for (const auto e : dLv) {
			const uint32_t v = LEExtractV(e);
			// if (1023 < dLu[v]) continue;
			if (dLu[v] == UINT32_MAX) continue;
			const uint32_t dd = dLu[v] + LEExtractD(e);
			if (dd < d) d = dd;
		}
		return d;
}

// Fast calculate distance and counting
std::pair<uint32_t, uint32_t> USPCUpdate::FastDistanceCount(const std::vector<std::pair<uint32_t, uint32_t>>& dLu,
	const std::vector<LabelEntry>& dLv) const {
		uint32_t d = UINT32_MAX;
		uint32_t c = 0;
		for (const auto e : dLv) {
			const uint32_t v = LEExtractV(e);
			if (dLu[v].first == UINT32_MAX) continue;
			const uint32_t dd = dLu[v].first + LEExtractD(e);
			if (dd < d) { 
				d = dd; 
				c = 0;
			} else if (dd == d) 
				c += (dLu[v].second * LEExtractC(e));
		}
		return std::make_pair(d, c);
}

/*
**************************
****Decremental update****
**************************
*/

std::tuple<uint32_t,uint32_t,size_t,size_t,size_t,size_t,uint32_t,uint32_t,uint32_t,uint32_t> 
USPCUpdate::Dec_SPC(uint32_t a, uint32_t b) {

	std::cout << "Degree: " << G_[a].size() << " - " << G_[b].size() << "\n";

    const auto start = std::chrono::steady_clock::now();

	std::vector<int> hubList_a(n_, 0);
	std::vector<int> hubList_b(n_, 0);

	for (auto cL_a:cL_[a]) {
		hubList_a[LEExtractV(cL_a)] = 1;
	}

	for (auto cL_b:cL_[b]) {
		hubList_b[LEExtractV(cL_b)] = 1;
	}

	int recA = 0, recB = 0;
	std::vector<int> Aff_a_flag(n_, 0); // -1 for receriver, 0 for unkown, 1 for affected
	std::vector<int> Aff_b_flag(n_, 0);
	std::vector<uint32_t> Aff_a; // store affected A
	std::vector<uint32_t> Aff_b;
	std::vector<uint32_t> Rec_a; // store receiver A
	std::vector<uint32_t> Rec_b;

	// Find Affected A
	std::vector<uint32_t> Da(n_, UINT32_MAX);
	std::vector<uint32_t> Ca(n_, 0);
	Da[a] = 0; Ca[a] = 1;

	std::queue<uint32_t> Qa({a});

	while (!Qa.empty()) {
        
        auto current = std::chrono::steady_clock::now();

		auto u = Qa.front(); Qa.pop();
		auto dc_u_b = Count(u, b);

		if (Da[u] + 1 != dc_u_b.first) continue;

		if (Ca[u] < dc_u_b.second && (hubList_a[u] == 0 || hubList_b[u] == 0)) {

			Aff_a_flag[rank_[u]] = -1; 
			Rec_a.push_back(u);

		} else {

			// else, an affected, both sender and receriver
			Aff_a_flag[rank_[u]] = 1;
			Aff_a.push_back(rank_[u]);

		}

		for (auto nbr:G_[u]) {
			if (Da[nbr] == UINT32_MAX) {
				Da[nbr] = Da[u] + 1;
				Ca[nbr] = Ca[u];
				Qa.push(nbr);
			} else if (Da[nbr] == Da[u] + 1) {
				Ca[nbr] += Ca[u];
			}
		}
	}

	// Find Affected B
	std::vector<uint32_t> Db(n_, UINT32_MAX);
	std::vector<uint32_t> Cb(n_, 0);
	Db[b] = 0; Cb[b] = 1;

	std::queue<uint32_t> Qb({b});

	while (!Qb.empty()) {

        auto current = std::chrono::steady_clock::now();

		auto u = Qb.front(); Qb.pop();
		auto dc_u_a = Count(u, a);

		if (Db[u] + 1 != dc_u_a.first) continue;

		if (Cb[u] < dc_u_a.second && (hubList_a[u] == 0 || hubList_b[u] == 0)) {

			Aff_b_flag[rank_[u]] = -1; 
			Rec_b.push_back(u);

		} else {

			Aff_b_flag[rank_[u]] = 1;
			Aff_b.push_back(rank_[u]);

		}

		for (auto nbr:G_[u]) {
			if (Db[nbr] == UINT32_MAX) {
				Db[nbr] = Db[u] + 1;
				Cb[nbr] = Cb[u];
				Qb.push(nbr);
			} else if (Db[nbr] == Db[u] + 1) {
				Cb[nbr] += Cb[u];
			}
		}
	}
    
    uint32_t renew_C = 0, renew_D = 0, insert = 0, remove = 0;

	// fast update for some special cases
	auto fast_res = Fast_update(a, b, Aff_a_flag, Aff_a, Aff_b, Rec_a, Rec_b);

	if (std::get<0>(fast_res) == 1) {

		return std::make_tuple(rank_[a], rank_[b], Aff_a.size(), Aff_b.size(), Rec_a.size(), Rec_b.size(), 
        std::get<1>(fast_res), std::get<2>(fast_res), std::get<3>(fast_res), std::get<4>(fast_res));

	}
	
	// remove b from G_[a] and a from G_[b]
	G_[a].erase(std::remove(G_[a].begin(), G_[a].end(), b), G_[a].end());
	G_[b].erase(std::remove(G_[b].begin(), G_[b].end(), a), G_[b].end());

	// Update labels here

	std::sort(Aff_a.begin(), Aff_a.end());
	std::sort(Aff_b.begin(), Aff_b.end());

	for (uint32_t ai = 0, bi = 0; ai < Aff_a.size() || bi < Aff_b.size(); ) {

        auto current = std::chrono::steady_clock::now();

		if (bi == Aff_b.size() || ((ai < Aff_a.size() && bi < Aff_b.size()) && (Aff_a[ai] < Aff_b[bi]))) {

			// update hub order_[ai]
			int is_hub = (hubList_a[order_[Aff_a[ai]]] == 1 && hubList_b[order_[Aff_a[ai]]] == 1) ? 1 : 0;

			auto res_ = Update_hub(order_[Aff_a[ai]], Aff_b_flag, Aff_b, Rec_b, is_hub);

            renew_C += std::get<0>(res_); renew_D += std::get<1>(res_); 
            insert += std::get<2>(res_); remove += std::get<3>(res_); 

			++ai;

		} else if (ai == Aff_a.size() || ((ai < Aff_a.size() && bi < Aff_b.size()) && (Aff_a[ai] > Aff_b[bi]))) {

			//update hub order_[bi]
			int is_hub = (hubList_a[order_[Aff_b[bi]]] == 1 && hubList_b[order_[Aff_b[bi]]] == 1) ? 1 : 0;

			auto res_ = Update_hub(order_[Aff_b[bi]], Aff_a_flag, Aff_a, Rec_a, is_hub);

            renew_C += std::get<0>(res_); renew_D += std::get<1>(res_); 
            insert += std::get<2>(res_); remove += std::get<3>(res_); 

			++bi;

		}
		
	}

    return std::make_tuple(rank_[a], rank_[b], Aff_a.size(), Aff_b.size(), Rec_a.size(), Rec_b.size(), 
    renew_C, renew_D, insert, remove);
}

// Dec_Update: renew_only_C, renew_D, insert, remove
std::tuple<uint32_t, uint32_t, uint32_t, uint32_t> USPCUpdate::Update_hub(uint32_t hub, 
const std::vector<int>& Aff_list, const std::vector<uint32_t>& Affs,
const std::vector<uint32_t>& Recs, int is_hub) {

    uint32_t renew_C = 0, renew_D = 0, insert = 0, remove = 0;

	std::vector<int> updated_list(n_, 0);
	std::vector<uint32_t> D(n_, UINT32_MAX);
	std::vector<uint32_t> C(n_, 0);
	D[hub] = 0; C[hub] = 1;

	std::queue<uint32_t> Q({hub});

	while (!Q.empty()) {
		auto v = Q.front(); Q.pop();
		
        if (v != hub) {

			if (Aff_list[rank_[v]] == 0) {

				auto dis_h_v = Query_Distance(hub, v);
				if (dis_h_v < D[v]) continue;

			} else {
				auto dcp_soFar = Query_Search(hub, v);
				auto d_over = std::get<0>(dcp_soFar);
				auto c_over = std::get<1>(dcp_soFar);
				auto d_h = std::get<2>(dcp_soFar);
				auto c_h = std::get<3>(dcp_soFar);
				auto pos = std::get<4>(dcp_soFar);

				if (D[v] > d_over) {
                    continue; // DIS PRUNER
                }

				if (d_h == UINT32_MAX) {

						cL_[v].insert(cL_[v].begin() + pos, LEMerge(hub, D[v], C[v]));
                        ++insert;
						updated_list[v] = 1;

				} else {

						if (d_h != D[v] || c_h != C[v]) {

							cL_[v][pos] = LEMerge(hub, D[v], C[v]);
							updated_list[v] = 1;

                            if (d_h == D[v]) ++renew_C;
                            else ++renew_D;

						} else {
							updated_list[v] = 1;
						}
				}
			}
		}

		for (auto nbr:G_[v]) {
			if (rank_[nbr] <= rank_[hub]) continue;
			if (D[nbr] == UINT32_MAX) {
				D[nbr] = D[v] + 1;
				C[nbr] = C[v];
				Q.push(nbr);
			} else if (D[nbr] == D[v] + 1) {
				C[nbr] += C[v];
			}
		}
	}

	// if hub is the common hub of a, b, then there are potential erased labels
	if (is_hub == 1) {
		for (auto affV:Affs) {

			if (affV <= rank_[hub]) continue;

			auto cur_v = order_[affV];
			if (updated_list[cur_v] == 0) {

				for (size_t di = 0; di < cL_[cur_v].size(); ++di) {
					if (LEExtractV(cL_[cur_v][di]) == hub) {

						cL_[cur_v].erase(cL_[cur_v].begin() + di);
						++remove;
						updated_list[cur_v] = 1;

						break;
					}
				}
			}
		}

		for (auto recv:Recs) {

			if (rank_[recv] <= rank_[hub]) continue;

			if (updated_list[recv] == 0) {

				for (size_t di = 0; di < cL_[recv].size(); ++di) {
					if (LEExtractV(cL_[recv][di]) == hub) {

						cL_[recv].erase(cL_[recv].begin() + di);
						++remove;
						updated_list[recv] = 1;

						break;
					}
				}
			}
		}
	}

    return std::make_tuple(renew_C, renew_D, insert, remove);
}

// Isolated vertex optimization
std::tuple<int, uint32_t, uint32_t, uint32_t, uint32_t> USPCUpdate::Fast_update(uint32_t a, uint32_t b, 
const std::vector<int>& Aff_list, const std::vector<uint32_t>& AffA, const std::vector<uint32_t>& AffB,
const std::vector<uint32_t>& RecA, std::vector<uint32_t> RecB) {
	
	// disconnect a node from the rest graph
	if (RecA.size() == 0 && RecB.size() == 0) {

		if (AffA.size() == 1 && G_[a].size() == 1) {
			if (rank_[a] > rank_[b]) {
				uint32_t erase_cnt = cL_[a].size() - 1;
				std::vector <spc::LabelEntry>().swap(cL_[a]); 
				cL_[a].push_back(LEMerge(a, 0, 1));
				return std::make_tuple(1,0,0,0,erase_cnt);
			}

		} else if (AffB.size() == 1 && G_[b].size() == 1) {

			if (rank_[b] > rank_[a]) {
				uint32_t erase_cnt = cL_[b].size() - 1;
				std::vector <spc::LabelEntry>().swap(cL_[b]); 
				cL_[b].push_back(LEMerge(b, 0, 1));
				return std::make_tuple(1,0,0,0,erase_cnt);
			}

		}
	}

	return std::make_tuple(0,0,0,0,0);
}

std::tuple<uint32_t, uint64_t, uint32_t, uint64_t, uint32_t> USPCUpdate::Query_Search(uint32_t h, uint32_t v) {

	size_t p1 = 0, p2 = 0;
	uint32_t sp_d = UINT32_MAX;
	uint64_t sp_c = 0;

	uint32_t hub_d = UINT32_MAX;
	uint64_t hub_c = 0;
	uint32_t hub_pos = UINT32_MAX;

	while (p1 < cL_[h].size() && p2 < cL_[v].size()) {
		const uint32_t w1 = LEExtractV(cL_[h][p1]);
		const uint32_t w2 = LEExtractV(cL_[v][p2]);

		if (w2 == h) {
			hub_d = LEExtractD(cL_[v][p2]);
			hub_c = LEExtractC(cL_[v][p2]);
			hub_pos = p2;
			return std::make_tuple(sp_d, sp_c, hub_d, hub_c, hub_pos);
		} 
		
		if (rank_[w1] < rank_[w2]) ++p1;
		else if (rank_[w1] > rank_[w2]) ++p2;
		else {
		const uint32_t d = LEExtractD(cL_[h][p1]) +
			LEExtractD(cL_[v][p2]);

		if (d < sp_d) {
			sp_d = d;
			sp_c = static_cast<uint64_t>(LEExtractC(cL_[h][p1])) *
				LEExtractC(cL_[v][p2]);
		} else if (d == sp_d) {
			uint64_t c = static_cast<uint64_t>(LEExtractC(cL_[h][p1])) *
				LEExtractC(cL_[v][p2]);
			sp_c += c;
		}
		++p1; ++p2;
		}
	}

	return std::make_tuple(sp_d, sp_c, UINT32_MAX, 0, p2);
}

uint32_t USPCUpdate::Query_Distance(uint32_t hub, uint32_t v) {
	size_t p1 = 0, p2 = 0;
	uint32_t sp_d = UINT32_MAX;

	while (p1 < cL_[v].size() && p2 < cL_[hub].size()) {
		const uint32_t w1 = LEExtractV(cL_[hub][p1]);
		const uint32_t w2 = LEExtractV(cL_[v][p2]);

		if (rank_[w1] < rank_[w2]) ++p1;
		else if (rank_[w1] > rank_[w2]) ++p2;
		else {
			const uint32_t d = LEExtractD(cL_[hub][p1]) +
				LEExtractD(cL_[v][p2]);

			if (d < sp_d) {
				sp_d = d;
			} 
			++p1; ++p2;
		}
	}

	return sp_d;
}

}	// namespace spc
