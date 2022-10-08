#! /bin/bash -e
make
echo "^^^^Compile completed^^^^"
echo " "

echo "Graph 0 Indexing (Ori)"
./u_index -g graph/0.txt -l label/0_ori -o degree -f info/0_ori.txt
echo ""

echo "Querying (Ori)"
./u_query -l label/0_ori -q query/0_q.txt -a answer/0_ori.txt -g graph/0.txt -t n -u n
echo ""

echo "Inc Updating"
./u_update -l label/0_ori -n label/0_inc -u update/0_inc.txt -i info/0_inc.txt
echo ""

echo "Querying (Inc Updated)"
./u_query -l label/0_inc -q query/0_q.txt -a answer/0_inc.txt -g graph/0.txt -t y -u update/0_inc.txt 
echo ""

echo "Dec Updating"
./u_update -l label/0_ori -n label/0_dec -u update/0_dec.txt -i info/0_dec.txt
echo ""

echo "Querying (Dec Updated)"
./u_query -l label/0_dec -q query/0_q.txt -a answer/0_dec.txt -g graph/0.txt -t y -u update/0_dec.txt
echo ""
