# Dynamic_SPC

## Efficiently Answering Shortest Path Counting on Dynamic Graphs

We study the maintenance of the 2-hop labeling for shortest path counting. We adopt the SPC-Index constructed by HP-SPC as the underlying index and aim to maintain the index regarding the topological changes to the graph. We propose two dynamic algorithms including incremental algorithm IncSPC and decremental algorithm DecSPC, and also optimization for the isolated vertex. An example graph and its related files are provided, run "*sh dspc_0.sh*" to run the examples.

****
* Author: Qingshuai Feng
* Co-author: You Peng
****

## Folders:
### graph: graph files
|File name|Description|
|---|----|
|0.txt|original graph|
---
### label:	label files
---
### query: query files
|File name|Description|
|---|----|
|0_q.txt|query pairs|
---
###	info: information files
|File name|Description|
|---|----|
|0_ori.txt|index contruction time and # of label entries for graph 0.txt|
|0_inc.txt|update info for incremnetal update|
|0_dec.txt|update info for decremnetal update|
---
### answer: answers files
---
---
## Files:
|File name|Description|
|---|----|
|macros.h|macros operations|
|u_label.h|define labels|
|u_io.cc & u_io.h|read graph|
|u_spc.h & u_spc.cc|all implementations|
|u_index.cc|building index|
|u_query.cc|query|
|u_update.cc|update |
|dspc_0.sh|script for running|
|Makefile|Makefile|

## Execution: (Examples see run_script.sh)
### ./u_index:
|Parameters|Type|Description|
|--|--|---|
|g|string|graph_file|
|l|string|label_file|
|o|degree|ordering|
|f|string|info_file|

### ./u_query:
|Parameters|Type|Description|
|--|--|---|
|l|string|label_file|
|q|string|query_folder|
|a|string|answer_folder|
|g|string|graph_file|
|t|char|index_merge_flag|
|u|string|update_file|

### ./u_update:
|Parameters|Type|Description|
|--|--|---|
|l|string|label_file|
|n|string|updated_label_folder|
|u|string|update_file|
|i|string|info_file|
