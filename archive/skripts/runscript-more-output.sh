 #!/bin/bash
echo "small with extra output"
echo "2-level"
echo "heterogeneous"
./main-temp /home/kulaginasv/workspace/MemComRefactored/tree-scheduler-refactored/real_Trees/trees/ treesListAllSmall.txt 1 1 1 -1 0
echo "homogeneous many small"
./main-temp /home/kulaginasv/workspace/MemComRefactored/tree-scheduler-refactored/real_Trees/trees/ treesListAllSmall.txt 0 1 1 -1 1
echo "homogeneous average average"
./main-temp /home/kulaginasv/workspace/MemComRefactored/tree-scheduler-refactored/real_Trees/trees/ treesListAllSmall.txt 0 1 1 -1 2
echo "homogeneous few big"
./main-temp /home/kulaginasv/workspace/MemComRefactored/tree-scheduler-refactored/real_Trees/trees/ treesListAllSmall.txt 0 1 1 -1 3
echo "3-level"
echo "homogeneous many small"
./main-temp /home/kulaginasv/workspace/MemComRefactored/tree-scheduler-refactored/real_Trees/trees/ treesListAllSmall.txt 0 1 2 -1 1
echo "homogeneous average average"
./main-temp /home/kulaginasv/workspace/MemComRefactored/tree-scheduler-refactored/real_Trees/trees/ treesListAllSmall.txt 0 1 2 -1 2
echo "homogeneous few big"
./main-temp /home/kulaginasv/workspace/MemComRefactored/tree-scheduler-refactored/real_Trees/trees/ treesListAllSmall.txt 0 1 2 -1 3
echo "heterogeneous"
./main-temp /home/kulaginasv/workspace/MemComRefactored/tree-scheduler-refactored/real_Trees/trees/ treesListAllSmall.txt 1 1 2 -1 0

