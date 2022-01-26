 #!/bin/bash
echo "small with extra output"
echo "2-level"
echo "heterogeneous"
./main ./real_Trees/random_trees/old_trees/ treesList.txt 1 0 1 -1 0 ./clusters/cluster-27proc-1-3.json
echo "homogeneous many small"
./main ./real_Trees/random_trees/old_trees/ treesList.txt 0 0 1 -1 1 ./clusters/cluster-hom-27proc-1.json
echo "homogeneous few big"
./main ./real_Trees/random_trees/old_trees/ treesList.txt 0 0 1 -1 3 ./clusters/cluster-hom-18proc-3.json
echo "3-level"
echo "heterogeneous"
./main ./real_Trees/random_trees/old_trees/ treesList.txt 1 0 2 -1 0 ./clusters/cluster_27proc_1-1,5-3.json
echo "homogeneous many small"
./main ./real_Trees/random_trees/old_trees/ treesList.txt 0 0 2 -1 1 ./clusters/cluster-hom-27proc-1.json
echo "homogeneous average average"
./main ./real_Trees/random_trees/old_trees/ treesList.txt 0 0 2 -1 2 ./clusters/cluster-hom-18proc-1,5.json
echo "homogeneous few big"
./main ./real_Trees/random_trees/old_trees/ treesList.txt 0 0 2 -1 3 ./clusters/cluster-hom-9proc-3.json

