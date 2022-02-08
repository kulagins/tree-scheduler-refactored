#!/bin/bash

cd real_Trees/random_trees
#"10_children all_large large_makespan_weights 20_children all_small large_node_weights 3_children large_edge_weights"
directories="large_makespan_weights"

: 
echo "2-level"
echo "heterogeneous"

pwd
for dir in $directories; do
    cd "$dir"
    treeDir=$(pwd)
    cd ../../..
    echo "!! $treeDir"
    ./main $treeDir /treesList.txt clusters/cluster-27proc-1-3.json 0 0 0   
    cd real_Trees/random_trees
done


echo "many small"
for dir in $directories; do
    cd "$dir"
    treeDir=$(pwd)
    cd ../../..
    echo "!! $treeDir"§
    ./main $treeDir /treesList.txt clusters/cluster-hom-27proc-1.json 0 0 0
   
    cd real_Trees/random_trees
done

echo "few big"
for dir in $directories; do
    cd "$dir"
    treeDir=$(pwd)
    cd ../../..
    echo "!! $treeDir"
    ./main $treeDir /treesList.txt clusters/cluster-hom-18proc-3.json 0 0 0
   
    cd real_Trees/random_trees
done


echo "3-level"
echo "heterogeneous"

pwd
for dir in $directories; do
    cd "$dir"
    treeDir=$(pwd)
    cd ../../..
    echo "!! $treeDir"§
    ./main $treeDir /treesList.txt clusters/cluster_27proc_1-1,5-3.json 0 0 0
    cd real_Trees/random_trees
done
echo "many small"
for dir in $directories; do
    cd "$dir"
    treeDir=$(pwd)
    cd ../../..
    echo "!! $treeDir"§
    ./main $treeDir /treesList.txt clusters/cluster-hom-27proc-1.json 0 0 0
   
    cd real_Trees/random_trees
done
echo "avg avg"
for dir in $directories; do
    cd "$dir"
    treeDir=$(pwd)
    cd ../../..
    echo "!! $treeDir"§
    ./main $treeDir /treesList.txt clusters/cluster-hom-18proc-1,5.json 0 0 0
    cd real_Trees/random_trees
done
echo "few big"
for dir in $directories; do
    cd "$dir"
    treeDir=$(pwd)
    cd ../../..
    echo "!! $treeDir"§
    ./main $treeDir /treesList.txt clusters/cluster-hom-9proc-3.json 0 0 0
    cd real_Trees/random_trees
done

