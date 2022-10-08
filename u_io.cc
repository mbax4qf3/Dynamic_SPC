#include "u_io.h"

#include <algorithm>
#include <cinttypes>
#include <cstdint>
#include <cstdio>
#include <set>

#include "macros.h"

void GraphRead(const std::string& filename,
               spc::Graph& graph,
               uint32_t& n, uint32_t& m) {
  ASSERT(graph.empty());
  FILE* gfile = fopen(filename.c_str(), "r");
  fscanf(gfile, "%" SCNu32 " %" SCNu32, &n, &m);
  // check the # of vertices
  spc::NormalV(n);
  // construct graph
  graph.resize(n);
  std::set<std::pair<uint32_t, uint32_t>> edges;
  for (uint32_t e = 0; e < m; ++e) {
    uint32_t v1, v2;
    fscanf(gfile, "%" SCNu32 " %" SCNu32, &v1, &v2);
    ASSERT(v1 < n && v2 < n);
    
    // check
    ASSERT(v1 != v2);
    if (v1 > v2) std::swap(v1, v2);
    // ASSERT(0 == edges.count({v1, v2}));
    if (0 == edges.count({v1, v2})) {
      edges.insert({v1, v2});
      graph[v1].push_back(v2);
      graph[v2].push_back(v1);
    }
  }
  for (uint32_t v = 0; v < n; ++v) {
    std::sort(graph[v].begin(), graph[v].end());
  }
  fclose(gfile);
}
