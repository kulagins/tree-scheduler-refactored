#!/bin/bash

echo "normal random"
./main ./data/random_trees/random/ treesListSmall.txt ./clusters/clusterList_single_small.txt 0 0 1 LMW EX MW
./main ./data/random_trees/random/ treesListSmall.txt ./clusters/clusterList_single_small.txt 0 0 1 LMW EX MD
./main ./data/random_trees/random/ treesListSmall.txt ./clusters/clusterList_single_small.txt 0 0 1 LMW EX CP

./main ./data/random_trees/random/ treesListSmall.txt ./clusters/clusterList_single_small.txt 0 0 1 LMW FFT MW
./main ./data/random_trees/random/ treesListSmall.txt ./clusters/clusterList_single_small.txt 0 0 1 LMW FFT MD
./main ./data/random_trees/random/ treesListSmall.txt ./clusters/clusterList_single_small.txt 0 0 1 LMW FFT CP

./main ./data/random_trees/random/ treesListSmall.txt ./clusters/clusterList_single_small.txt 0 0 1 LMW M MW
./main ./data/random_trees/random/ treesListSmall.txt ./clusters/clusterList_single_small.txt 0 0 1 LMW M MD
./main ./data/random_trees/random/ treesListSmall.txt ./clusters/clusterList_single_small.txt 0 0 1 LMW M CP

./main ./data/random_trees/random/ treesListSmall.txt ./clusters/clusterList_single_small.txt 0 0 1 CP EX MW
./main ./data/random_trees/random/ treesListSmall.txt ./clusters/clusterList_single_small.txt 0 0 1 CP EX MD
./main ./data/random_trees/random/ treesListSmall.txt ./clusters/clusterList_single_small.txt 0 0 1 CP EX CP

./main ./data/random_trees/random/ treesListSmall.txt ./clusters/clusterList_single_small.txt 0 0 1 CP FFT MW
./main ./data/random_trees/random/ treesListSmall.txt ./clusters/clusterList_single_small.txt 0 0 1 CP FFT MD
./main ./data/random_trees/random/ treesListSmall.txt ./clusters/clusterList_single_small.txt 0 0 1 CP FFT CP

./main ./data/random_trees/random/ treesListSmall.txt ./clusters/clusterList_single_small.txt 0 0 1 CP M MW
./main ./data/random_trees/random/ treesListSmall.txt ./clusters/clusterList_single_small.txt 0 0 1 CP M MD
./main ./data/random_trees/random/ treesListSmall.txt ./clusters/clusterList_single_small.txt 0 0 1 CP M CP
