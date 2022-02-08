 #!/bin/bash
 
clusters="36proc_0,5-1-1,5-3.json  36proc_1-1,5-3.json  36proc-1-3.json  hom-18proc-1,5.json  hom-18proc-3.json  hom-27proc-1,5.json  hom-27proc-1.json  hom-36proc-1.json  hom-9proc-3.json"

echo "matrix"
pwd
for cluster in $clusters; do 
    echo "!! $cluster"
    ./main ./real_Trees/trees/ treesListNoBiggest.txt ./clusters/36-proc/$cluster 0 0 0      
done


echo "normal random"

for cluster in $clusters; do 
    echo "!! $cluster"
    ./main ./real_Trees/random_trees/old_trees/ treesList.txt ./clusters/36-proc/$cluster 0 0 0   
done


pwd
cd real_Trees/random_trees
directories="10_children all_large large_makespan_weights 20_children all_small large_node_weights 3_children large_edge_weights"

pwd
for cluster in $clusters; do 
    for dir in $directories; do
        cd "$dir"
        treeDir=$(pwd)
        cd ../../..
        echo "!! $treeDir"
        echo "!! $cluster"
        ./main $treeDir /treesList.txt ./clusters/36-proc/$cluster 0 0 0           
        cd real_Trees/random_trees
    done
done
