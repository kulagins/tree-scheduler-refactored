#!/bin/bash

cd real_Trees/random_trees
directories="10_children all_large large_makespan_weights 20_children all_small large_node_weights 3_children large_edge_weights"

for dir in $directories; do
    cd "$dir"
    
    echo "" > treesList.txt
    
    for subdir in */;do
        cd "$subdir"
        for file in *; do
            echo "$subdir$file" >> ../treesList.txt
        done
        cd ..
    done
    cd ..
done

: '
echo "small with extra output"
echo "2-level"
echo "heterogeneous"

pwd
for dir in $directories; do
    cd "$dir"
    treeDir=$(pwd)
    cd ../../..
    echo "$treeDir"
    ./res "$treeDir" treesList.txt 1 0 1 -1 0 clusters/cluster-27proc-1-3.json
    cd real_Trees/random_trees
done
#./main real_Trees/random_trees/ treesList.txt 1 0 1 -1 0 clusters/cluster-27proc-1-3.json
' 
