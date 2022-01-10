 #!/bin/bash
echo "small with extra output"
echo "2-level"
echo "heterogeneous"
./main /home/kulaginasv/workspace/MemComRefactored/tree-scheduler-refactored/real_Trees/trees/ treesListNoBiggest.txt 1 0 1 -1 0 /home/kulaginasv/workspace/MemComRefactored/tree-scheduler-refactored/clusters/cluster-27proc-1-3.json
#echo "homogeneous many small"
#./main /home/kulaginasv/workspace/MemComRefactored/tree-scheduler-refactored/real_Trees/trees/ treesListNoBiggest.txt 0 0 1 -1 1 /home/kulaginasv/workspace/MemComRefactored/tree-scheduler-refactored/clusters/cluster-hom-27proc-1.json
#echo "homogeneous few big"
#./main /home/kulaginasv/workspace/MemComRefactored/tree-scheduler-refactored/real_Trees/trees/ treesListNoBiggest.txt 0 0 1 -1 3 /home/kulaginasv/workspace/MemComRefactored/tree-scheduler-refactored/clusters/cluster-hom-18proc-3.json
echo "3-level"
echo "heterogeneous"
./main /home/kulaginasv/workspace/MemComRefactored/tree-scheduler-refactored/real_Trees/trees/ treesListNoBiggest.txt 1 0 2 -1 0 /home/kulaginasv/workspace/MemComRefactored/tree-scheduler-refactored/clusters/cluster_27proc_1-1,5-3.json
#echo "homogeneous many small"
#./main /home/kulaginasv/workspace/MemComRefactored/tree-scheduler-refactored/real_Trees/trees/ treesListNoBiggest.txt 0 0 2 -1 1 /home/kulaginasv/workspace/MemComRefactored/tree-scheduler-refactored/clusters/cluster-hom-27proc-1.json
#echo "homogeneous average average"
#./main /home/kulaginasv/workspace/MemComRefactored/tree-scheduler-refactored/real_Trees/trees/ treesListNoBiggest.txt 0 0 2 -1 2 /home/kulaginasv/workspace/MemComRefactored/tree-scheduler-refactored/clusters/cluster-hom-18proc-1,5.json
#echo "homogeneous few big"
#./main /home/kulaginasv/workspace/MemComRefactored/tree-scheduler-refactored/real_Trees/trees/ treesListNoBiggest.txt 0 0 2 -1 3 /home/kulaginasv/workspace/MemComRefactored/tree-scheduler-refactored/clusters/cluster-hom-9proc-3.json

