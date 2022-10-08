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
#include <iomanip>
#include <thread>
#include <omp.h>
#include <mutex>
#include <condition_variable>

#include "progressbar.h"
#include "macros.h"
#include "u_io.h"
#include "u_label.h"
#include "u_spc.h"

using namespace std::chrono_literals;

int main(int argc, char** argv) {
    std::string lfilename;
    std::string ufilename; // udate edges
    std::string ifilename; // info file
    std::string newlfilename; // updated index file
    int option = -1;
    while (-1 != (option = getopt(argc, argv, "l:n:u:i:"))) {
        switch (option) {
            case 'l':
            lfilename = optarg; break;
            case 'n':
            newlfilename = optarg; break;
            case 'u':
            ufilename = optarg; break;
            case 'i':
            ifilename = optarg; break;
        }
    }
    
    printf("ori label file: %s\n", lfilename.c_str());
    printf("new label file: %s\n", newlfilename.c_str());
    printf("update file: %s\n", ufilename.c_str());

    // read index
    spc::USPCUpdate uspu;

    // read update edges
    FILE* file_u = fopen(ufilename.c_str(), "r");

    uint32_t num_update = 0;
    // std::vector<std::tuple<uint32_t, uint32_t, char>> upd_edges;

    fscanf(file_u, "%" SCNu32, &num_update);
    
    // Write info title
    std::ofstream ifile;
    ifile.open(ifilename.c_str(), std::ios::app);


    // read index
    uspu.IndexRead(lfilename);

    auto during = std::chrono::steady_clock::now() - std::chrono::steady_clock::now();

    std::cout << "\nStart update: \n";
    std::cout << "There are " << num_update << " updates\n";

    progressbar bar(num_update);

    const auto beg_mp = std::chrono::steady_clock::now();

    for (int i = 0; i < num_update; ++i) {
        
        bar.update();

        uint32_t v1, v2;
        char upd_type;

        fscanf(file_u, "%" SCNu32 " %" SCNu32 " %c", &v1, &v2, &upd_type);


        if (upd_type == 'i') {
            
            const auto beg = std::chrono::steady_clock::now();
            auto inc_info = uspu.Inc_SPC(v1, v2); // inc update
            const auto end = std::chrono::steady_clock::now();
            const auto dif = end - beg;
            
            // write info + update time into info file
            ifile << std::get<0>(inc_info) << std::setw(12) << std::get<1>(inc_info) << std::setw(12) 
            << std::get<2>(inc_info) << std::setw(12) << std::get<3>(inc_info) << std::setw(12) 
            << std::get<4>(inc_info) << std::setw(12) << std::get<5>(inc_info) << std::setw(12) 
            << std::get<6>(inc_info) << std::setw(12) 
            << std::chrono::duration<double, std::milli>(dif).count() << std::setw(12) << "\n";

            during += dif;
                
        } else if (upd_type == 'd'){   
            
            // decremental update
            const auto beg = std::chrono::steady_clock::now();

            auto dec_info = uspu.Dec_SPC(v1, v2); // dec update

            const auto end = std::chrono::steady_clock::now();
            const auto dif = end - beg;

            ifile << v1 << std::setw(12) << v2 << std::setw(12) 
            << std::get<0>(dec_info) <<  std::setw(12) << std::get<1>(dec_info) << std::setw(12) // rank
            << std::get<2>(dec_info) << std::setw(12) << std::get<3>(dec_info) << std::setw(12)  // aff
            << std::get<4>(dec_info) << std::setw(12) << std::get<5>(dec_info) << std::setw(12)  // rec
            << std::get<6>(dec_info) << std::setw(12) << std::get<7>(dec_info) << std::setw(12)  // renew
            << std::get<8>(dec_info) << std::setw(12) << std::get<9>(dec_info) << std::setw(12)  // insert/remove
            << std::chrono::duration<double, std::milli>(dif).count() << "\n"; // time
            during += dif;

        }
    }

    const auto end_mp = std::chrono::steady_clock::now();
    const auto dif_mp = end_mp - beg_mp;
    std::cout << std::endl;
    auto label_num = uspu.IndexWrite(newlfilename);

    // write average update time into info file
    ifile << "Average: " << std::chrono::duration<double, std::milli>(during).count()/num_update << "\n";

    std::cout << "Average: " << std::chrono::duration<double, std::milli>(dif_mp).count()/num_update << "\n";

    ifile.close();

}
