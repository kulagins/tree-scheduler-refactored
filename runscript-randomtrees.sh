 #!/bin/bash
echo "small with extra output"
echo "2-level"
echo "heterogeneous"
./main /home/kulaginasv/workspace/MemComRefactored/tree-scheduler-refactored/real_Trees/random_trees/ treesList.txt /home/kulaginasv/workspace/MemComRefactored/tree-scheduler-refactored/clusters/cluster-27proc-1-3.json 0 0 0
echo "homogeneous many small"
./main /home/kulaginasv/workspace/MemComRefactored/tree-scheduler-refactored/real_Trees/random_trees/ treesList.txt /home/kulaginasv/workspace/MemComRefactored/tree-scheduler-refactored/clusters/cluster-hom-27proc-1.json 0 0 0 
echo "homogeneous few big"
./main /home/kulaginasv/workspace/MemComRefactored/tree-scheduler-refactored/real_Trees/random_trees/ treesList.txt /home/kulaginasv/workspace/MemComRefactored/tree-scheduler-refactored/clusters/cluster-hom-18proc-3.json 0 0 0
echo "3-level"
echo "heterogeneous"
./main /home/kulaginasv/workspace/MemComRefactored/tree-scheduler-refactored/real_Trees/random_trees/ treesList.txt /home/kulaginasv/workspace/MemComRefactored/tree-scheduler-refactored/clusters/cluster_27proc_1-1,5-3.json 0 0 0
echo "homogeneous many small"
./main /home/kulaginasv/workspace/MemComRefactored/tree-scheduler-refactored/real_Trees/random_trees/ treesList.txt /home/kulaginasv/workspace/MemComRefactored/tree-scheduler-refactored/clusters/cluster-hom-27proc-1.json 0 0 0
echo "homogeneous average average"
./main /home/kulaginasv/workspace/MemComRefactored/tree-scheduler-refactored/real_Trees/random_trees/ treesList.txt /home/kulaginasv/workspace/MemComRefactored/tree-scheduler-refactored/clusters/cluster-hom-18proc-1,5.json 0 0 0
echo "homogeneous few big"
./main /home/kulaginasv/workspace/MemComRefactored/tree-scheduler-refactored/real_Trees/random_trees/ treesList.txt /home/kulaginasv/workspace/MemComRefactored/tree-scheduler-refactored/clusters/cluster-hom-9proc-3.json 0 0 0

