#ifndef U_IO_H_
#define U_IO_H_

#include <string>

#include "u_label.h"

void GraphRead(const std::string& filename, spc::Graph& graph,
               uint32_t& n, uint32_t& m);

#endif
