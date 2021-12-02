 #!/bin/bash
echo "heterogeneous"
./main /home/kulaginasv/workspace/MemComRefactored/tree-scheduler-refactored/real_Trees/trees/ treesListSmall.txt 1 0 10 500 0
echo "homogeneous no adaptation"
./main /home/kulaginasv/workspace/MemComRefactored/tree-scheduler-refactored/real_Trees/trees/ treesListSmall.txt 0 0 10 500 0
echo "homogeneous nmany small"
./main /home/kulaginasv/workspace/MemComRefactored/tree-scheduler-refactored/real_Trees/trees/ treesListSmall.txt 0 0 10 500 1
echo "homogeneous average average"
./main /home/kulaginasv/workspace/MemComRefactored/tree-scheduler-refactored/real_Trees/trees/ treesListSmall.txt 0 0 10 500 2
echo "homogeneous few big"
./main /home/kulaginasv/workspace/MemComRefactored/tree-scheduler-refactored/real_Trees/trees/ treesListSmall.txt 0 0 10 500 3



