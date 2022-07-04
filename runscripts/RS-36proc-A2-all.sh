#!/bin/bash

pwd
echo "matrix"
pwd
./main ./data/trees/ treesListNoBiggest.txt ./clusters/clusterList_single_small.txt 0 0 1 LMW EX MD

echo "normal random"
 ./main ./data/random_trees/random/ treesList.txt ./clusters/clusterList_single_small.txt 0 0 1 LMW EX MD

pwd
cd data/random_trees
directories="10_children all_large large_makespan_weights 20_children all_small large_node_weights 3_children large_edge_weights"

pwd
 for dir in $directories; do
        cd "$dir"
        treeDir=$(pwd)
        cd ../../..
        echo "!! $treeDir"
        echo "!! $cluster"
        ./main $treeDir /treesList.txt ./clusters/clusterList_single_small.txt 0 0 1 LMW EX MD
        cd data/random_trees
 done


