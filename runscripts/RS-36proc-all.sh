#!/bin/bash

clusters="36proc_1-1,5-3.json  36proc-1-3.json  hom-18proc-1,5.json  hom-18proc-3.json  hom-27proc-1,5.json  hom-27proc-1.json  hom-36proc-1.json  hom-9proc-3.json"
#"36proc_0,5-1-1,5-3.json"

echo "matrix"
pwd
./main ./data/trees/ treesListNoBiggest.txt ./clusters/36-proc/36proc_0,5-1-1,5-3.json 0 0 0   
for cluster in $clusters; do 
    echo "!! $cluster"
    ./main ./data/trees/ treesListNoBiggest.txt ./clusters/36-proc/$cluster 0 1 0      
done


echo "normal random"
 ./main ./data/random_trees/random/ treesList.txt ./clusters/36-proc/36proc_0,5-1-1,5-3.json 0 0 0   
for cluster in $clusters; do 
    echo "!! $cluster"
    ./main ./data/random_trees/random/ treesList.txt ./clusters/36-proc/$cluster 0 1 0   
done


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
        ./main $treeDir /treesList.txt ./clusters/36-proc/36proc_0,5-1-1,5-3.json 0 0 0           
        cd data/random_trees
    done
for cluster in $clusters; do 
    for dir in $directories; do
        cd "$dir"
        treeDir=$(pwd)
        cd ../../..
        echo "!! $treeDir"
        echo "!! $cluster"
        ./main $treeDir /treesList.txt ./clusters/36-proc/$cluster 0 1 0           
        cd data/random_trees
    done
done
