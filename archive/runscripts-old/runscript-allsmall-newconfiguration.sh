 #!/bin/bash
echo "small with extra output"
echo "2-level"
echo "heterogeneous"
./main /home/kulaginasv/workspace/MemComRefactored/tree-scheduler-refactored/real_Trees/trees/ treesListAllSmall.txt 1 0 1 -1 0
echo "homogeneous many small"
./main /home/kulaginasv/workspace/MemComRefactored/tree-scheduler-refactored/real_Trees/trees/ treesListAllSmall.txt 0 0 1 -1 1
echo "homogeneous average average"
./main /home/kulaginasv/workspace/MemComRefactored/tree-scheduler-refactored/real_Trees/trees/ treesListAllSmall.txt 0 0 1 -1 2
echo "homogeneous few big"
./main /home/kulaginasv/workspace/MemComRefactored/tree-scheduler-refactored/real_Trees/trees/ treesListAllSmall.txt 0 0 1 -1 3
echo "3-level"
echo "heterogeneous"
./main /home/kulaginasv/workspace/MemComRefactored/tree-scheduler-refactored/real_Trees/trees/ treesListAllSmall.txt 1 0 2 -1 0
echo "homogeneous many small"
./main /home/kulaginasv/workspace/MemComRefactored/tree-scheduler-refactored/real_Trees/trees/ treesListAllSmall.txt 0 0 2 -1 1
echo "homogeneous average average"
./main /home/kulaginasv/workspace/MemComRefactored/tree-scheduler-refactored/real_Trees/trees/ treesListAllSmall.txt 0 0 2 -1 2
echo "homogeneous few big"
./main /home/kulaginasv/workspace/MemComRefactored/tree-scheduler-refactored/real_Trees/trees/ treesListAllSmall.txt 0 0 2 -1 3

