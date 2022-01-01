#ifndef inputenums_h
#define inputenums_h

enum HeterogenityLevels {
    homogeneus, memoryHeteregeneus, heterogeneus
};
enum ClusteringModes {
    treeDependent, staticClustering
};
enum HeterogeneousAdaptationMode {
    noAdaptation, fewBig, average, manySmall
};
#endif