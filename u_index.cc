#include <unistd.h>
#include <algorithm>
#include <chrono>
#include <cinttypes>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <queue>
#include <string>
#include <iostream>
#include <fstream>

#include "macros.h"
#include "u_io.h"
#include "u_label.h"
#include "u_spc.h"

bool ReadBool(const char tmp) {
  ASSERT('y' == tmp || 'n' == tmp);
  return 'y' == tmp;
}

int main(int argc, char** argv) {
  // initialize the options
  std::string gfilename;
  std::string lfilename;
  std::string ifilename;
  std::string osname;
  spc::USPCIndex::OrderScheme os = spc::USPCIndex::OrderScheme::kInvalid;

  {
    int option = -1;
    while (-1 != (option = getopt(argc, argv, "g:l:o:f:i:"))) { //s:e:i:
      switch (option) {
        case 'g':
          gfilename = optarg; break;
        case 'l':
          lfilename = optarg; break;
        case 'f': // info
          ifilename = optarg; break;
        case 'o': // ordering
          osname = optarg;
          if ("degree" == osname) {
            os = spc::USPCIndex::OrderScheme::kDegree;
          } else {
            os = spc::USPCIndex::OrderScheme::kInvalid;
          }
          break;
      }
    }
    printf("graph file: %s\n", gfilename.c_str());
    printf("label file: %s\n", lfilename.c_str());
    printf("ordering: %s\n", osname.c_str());
  }

  // read the graph
  uint32_t n, m;
  spc::Graph graph;
  GraphRead(gfilename, graph, n, m);

  const auto beg = std::chrono::steady_clock::now();

  // build index
  spc::USPCIndex spc;
  spc.set_os(os);

  spc.BuildIndex(graph);

  const auto end = std::chrono::steady_clock::now();
  const auto dif = end - beg;
  printf("\nindex construction costs \033[47;31m%f ms\033[0m\n",
         std::chrono::duration<double, std::milli>(dif).count());
  
  auto label_num = spc.IndexWrite(lfilename);


  std::ofstream ifile;
  ifile.open(ifilename.c_str());
  ifile << "Index time: " << std::chrono::duration<double, std::milli>(dif).count() << std::endl 
  << "Index Num: " << label_num << std::endl;
}
