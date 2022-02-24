 #!/bin/bash
#echo "small with extra output"
echo "2-level"
echo "heterogeneous"
./main ./data/random_trees/random/ treesList.txt ./clusters/cluster-27proc-1-3.json 0 0 0
echo "homogeneous many small"
./main ./data/random_trees/random/ treesList.txt ./clusters/cluster-hom-27proc-1.json 0 0 0
echo "homogeneous few big"
./main ./data/random_trees/random/ treesList.txt ./clusters/cluster-hom-18proc-3.json 0 0 0
echo "3-level"
echo "heterogeneous"
./main ./data/random_trees/random/ treesList.txt ./clusters/cluster_27proc_1-1,5-3.json 0 0 0
echo "homogeneous many small"
./main ./data/random_trees/random/ treesList.txt ./clusters/cluster-hom-27proc-1.json 0 0 0
echo "homogeneous average average"
./main ./data/random_trees/random/ treesList.txt ./clusters/cluster-hom-18proc-1,5.json 0 0 0
echo "homogeneous few big"
./main ./data/random_trees/random/ treesList.txt ./clusters/cluster-hom-9proc-3.json 0 0 0

