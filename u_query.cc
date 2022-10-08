#include <unistd.h>
#include <chrono>
#include <cinttypes>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <fstream>
#include <queue>
#include <set>
#include <string>
#include <vector>
#include <thread>
#include <omp.h>

#include "progressbar.h"
#include "macros.h"
#include "u_io.h"
#include "u_label.h"
#include "u_spc.h"

int main(int argc, char** argv) {
    std::string lfilename; // label file
    std::string qfilename; // query file
    std::string afilename; // answer file
    std::string gfilename; // graph file
    std::string ufilename; // update edges; graph should be gfile + ufile
    std::string index_Tag; // index flag
    int option = -1;
    while (-1 != (option = getopt(argc, argv, "l:a:q:g:t:u:"))) {
        switch (option) {
            case 'l':
                lfilename = optarg; break;
            case 'q':
                qfilename = optarg; break;
            case 'a':
                afilename = optarg; break;
            case 'g':
                gfilename = optarg; break;
            case 'u':
                ufilename = optarg; break;
            case 't':
                index_Tag = optarg; break;
        }
    }

    printf("label file: %s\n", lfilename.c_str());
    printf("graph file: %s\n", gfilename.c_str());
    printf("update file (if any): %s\n", ufilename.c_str());
    printf("answer file: %s\n", afilename.c_str());

    // read index
    spc::USPCQuery uspc;
    spc::USPCBasic uspb;
    if (index_Tag == "n")
        uspc.IndexRead(lfilename); // cL_ and dL_ are read separatly and then merged, used for querying with an original index
    else
        uspc.IndexRead_UPD(lfilename); // cL_ and dL_ are merged before, used for querying with an updated index

    // read queries
    FILE* file = fopen(qfilename.c_str(), "r");
    uint32_t num_queries = 0;
    fscanf(file, "%" SCNu32, &num_queries);
    std::vector<std::pair<uint32_t, uint32_t>> queries;

    for (uint32_t q = 0; q < num_queries; ++q) {
        uint32_t v1, v2;
        fscanf(file, "%" SCNu32 " %" SCNu32, &v1, &v2);
        queries.push_back(std::make_pair(v1, v2));
    }

    fclose(file);

    // if we need BFS (for correctness proof)
    std::vector<std::pair<uint32_t, uint64_t>> results_bfs;
    if (gfilename != "n") {
        
        printf("BFS Querying:\n");
        uint32_t n, m;
        spc::Graph graph;

        GraphRead(gfilename, graph, n, m);

        // if query under an updated graph, insert or delete edges from ori graph first
        if (ufilename != "n") {
            FILE* file_u = fopen(ufilename.c_str(), "r");

            uint32_t num_update = 0;
            std::vector<std::tuple<uint32_t, uint32_t, char>> upd_edges;

            fscanf(file_u, "%" SCNu32, &num_update);
            
            for (uint32_t u = 0; u < num_update; ++u) {
                uint32_t v1, v2;
                char upd_type;
                fscanf(file_u, "%" SCNu32 " %" SCNu32 " %c", &v1, &v2, &upd_type);
                uspc.UpdateGraph(graph, v1, v2, upd_type);
                upd_edges.push_back(std::make_tuple(v1, v2, upd_type));
            }

            fclose(file_u);
        }

        // BFS answer quering
        std::string bfsafilename = afilename.substr(0,7) + "bibfs_" + afilename.substr(7,afilename.size()-7);
        std::ofstream bfsafile;
        bfsafile.open(bfsafilename.c_str());
        auto btotal = std::chrono::steady_clock::now() - std::chrono::steady_clock::now();

        progressbar bar(queries.size());
        for (int i = 0; i < queries.size(); ++i) {
            auto query = queries[i];
            bar.update();
            const uint32_t v1 = query.first;
            const uint32_t v2 = query.second;
            std::pair<uint32_t, uint64_t> result;

            const auto beg_bfs = std::chrono::steady_clock::now();

            result = uspc.bi_BFS_Count(graph, v1, v2);

            const auto end_bfs = std::chrono::steady_clock::now();
            const auto dif_bfs = end_bfs - beg_bfs;
            
            btotal += dif_bfs;

            bfsafile << v1 << "\t" << v2 << "\t" << result.first << "\t" << result.second << "\t" 
            << std::chrono::duration<double, std::micro>(dif_bfs).count() << "\n";
            results_bfs.push_back(result);
        }
        bfsafile.close();

        printf("\nBFS costs \033[47;31m%f microseconds\033[0m in average\n",
            std::chrono::duration<double, std::micro>(btotal).count() / num_queries);
    }

    // Compute the results by hub labeling
    printf("Hub Labeling Querying:\n");
    std::vector<std::pair<uint32_t, uint64_t>> results;
    std::ofstream afile;
    afile.open(afilename.c_str());
    auto qtotal = std::chrono::steady_clock::now() - std::chrono::steady_clock::now();

    progressbar bar(queries.size());
    for (int i = 0; i < queries.size(); ++i) {
        auto query = queries[i];
        bar.update();
        const uint32_t v1 = query.first;
        const uint32_t v2 = query.second;
        std::pair<uint32_t, uint64_t> result;

        const auto beg = std::chrono::steady_clock::now();

        result = uspc.Count(v1, v2);
        
        const auto end = std::chrono::steady_clock::now();
        const auto dif = end - beg;
        qtotal += dif;

        afile << v1 << "\t" << v2 << "\t" << result.first << "\t" << result.second << "\t" 
        << std::chrono::duration<double, std::micro>(dif).count() << "\n";
        results.push_back(result);
    }

    afile.close();

    ASSERT(results.size() == num_queries);
    printf("\nQuery costs \033[47;31m%f microseconds\033[0m in average\n",
                std::chrono::duration<double, std::micro>(qtotal).count() / num_queries);
}
