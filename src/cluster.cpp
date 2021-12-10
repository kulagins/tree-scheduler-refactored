#include "../include/cluster.h"
#include "../include/tree.h"
#include <numeric>

using namespace std;

Cluster *Cluster::fixedCluster = NULL;

vector<double> Cluster::build3LevelMemorySizes(double minMem, double maxMem, unsigned int num_processors) {
    // cout << "minimal memory per processor" << minMem << ", maximal " << maxMem << endl;
    double cumulativeMem = 0;
    vector<double> memSizes(num_processors);
    memSizes.resize(num_processors);
    for (int k = 0; k < num_processors / 3; k++) {
        memSizes[k] = minMem;
        cumulativeMem += memSizes[k];
    }
    for (int k = num_processors / 3; k < 2 * num_processors / 3; k++) {
        memSizes[k] = (minMem + maxMem) / 2;
        cumulativeMem += memSizes[k];
    }
    for (int k = 2 * num_processors / 3; k < num_processors; k++) {
        memSizes[k] = maxMem;
        cumulativeMem += memSizes[k];
    }
    // cout << "cumulative mem in system: " << cumulativeMem << " with " << num_processors << " processors." << endl;
    return memSizes;
}

vector<double> Cluster::build3LevelMemorySizes(vector<double> memories, vector<unsigned int> processorGroupSizes) {
    // cout << "minimal memory per processor" << minMem << ", maximal " << maxMem << endl;
    double cumulativeMem = 0;
    double completeProcessorCount = std::accumulate(processorGroupSizes.begin(), processorGroupSizes.end(), 0);
    vector<double> memSizes(completeProcessorCount);
    memSizes.resize(completeProcessorCount);
    for (int i = 0; i < processorGroupSizes.size(); i++) {
        int nextGroupSize = processorGroupSizes.at(i);
        double nextMemSize = memories.at(i);
        for (int i = 0; i < nextGroupSize; i++) {
            memSizes.push_back(nextMemSize);
            cumulativeMem += nextMemSize;
        }
    }
    // cout << "cumulative mem in system: " << cumulativeMem << " with " << num_processors << " processors." << endl;
    return memSizes;
}

vector<double> Cluster::buildHomogeneousMemorySizes(double memSize, unsigned int num_processors) {

    vector<double> memSizes(num_processors);
    memSizes.resize(num_processors);
    for (int k = 0; k < num_processors; k++) {
        memSizes[k] = memSize;
    }

    // cout << "cumulative mem in system: " << num_processors * memSize << " with " << num_processors << " processors."
    //     << endl;
    return memSizes;
}


std::map<int, int> Cluster::buildProcessorSpeeds(int num_processors) {
    std::map<int, int> procSpeeds;
    for (int k = 0; k < num_processors / 3; k++) {
        procSpeeds.insert(pair<int, int>(k, 1));
    }
    for (int k = num_processors / 3; k < 2 * num_processors / 3; k++) {
        procSpeeds.insert(pair<int, int>(k, 2));
    }
    for (int k = 2 * num_processors / 3 + 1; k < num_processors; k++) {
        procSpeeds.insert(pair<int, int>(k, 3));
    }

    return procSpeeds;
}

void Cluster::buildHetStaticClusterWithConfiguration(int configNUmber, double maxMinMem, double maxEdgesToMakespanWeights) {
    vector<unsigned int> processorCounts;
    vector<double> mems;
    double maxMemInCluster = maxMinMem * 1.2;
    int num_processors;
    vector<double> memorySizes;
    switch (configNUmber) {
        case 1:
            //pearlmutter
            // 1536 small nodes memory 256GB
            // 40 middle ndoes 512 GB
            // 4 large nodes 1024 GB memory
            num_processors = 1580;
            processorCounts.insert(processorCounts.end(), {4, 40, 1536});
            maxMemInCluster = maxMinMem * 1.2;
            mems.insert(mems.end(), {maxMemInCluster, maxMemInCluster / 2, maxMemInCluster / 4});
            memorySizes = build3LevelMemorySizes(mems, processorCounts);
            break;
        case 2:
            // juwel
            // 2271 small nodes 96GB
            // 296 middle nodes 192GB
            // 16 big nodes 768GB
            num_processors = 2583;
            processorCounts.insert(processorCounts.end(), {16, 296, 2271});
            vector<double> mems{maxMemInCluster, maxMemInCluster / 4, maxMemInCluster / 8};
            memorySizes =  build3LevelMemorySizes(mems, processorCounts);
            break;
    }

    Cluster *cluster = new Cluster(num_processors, true);
    map<int, int> processor_speeds = Cluster::buildProcessorSpeeds(num_processors);
    cluster->setHomogeneousBandwidth(maxEdgesToMakespanWeights*10);
    Cluster::setFixedCluster(cluster);
    Cluster::getFixedCluster()->setMemorySizes(memorySizes);
}

void Cluster::buildHomStaticClusterWithConfiguration(int configNUmber, double maxMinMem, double maxEdgesToMakespanWeights, HeterogeneousAdaptationMode adaptationMode  ) {
    vector<unsigned int> processorCounts;
    vector<double> mems;
    double maxMemInCluster = maxMinMem * 1.2;
    int num_processors;
    vector<double> memorySizes;
    switch (configNUmber) {
        case 1:
            //pearlmutter
            // 1536 small nodes memory 256GB
            // 40 middle ndoes 512 GB
            // 4 large nodes 1024 GB memory
            switch(adaptationMode){
                case manySmall:
                    num_processors = 1536;
                    memorySizes = buildHomogeneousMemorySizes(maxMemInCluster, num_processors);

                    break;
                case average:
                    num_processors = 44;
                    memorySizes = buildHomogeneousMemorySizes(maxMemInCluster*2, num_processors);
                    break;
                case fewBig:
                    num_processors = 4;
                    memorySizes = buildHomogeneousMemorySizes(maxMemInCluster*4, num_processors);
                    break;
                default: cout<<"adaptation mode "<<adaptationMode<<endl; throw "bad adaptation mode";
            }
            break;
        case 2:
            // juwel
            // 2271 small nodes 96GB
            // 296 middle nodes 192GB
            // 16 big nodes 768GB
            switch(adaptationMode){
                case manySmall:
                    num_processors = 2271;
                    memorySizes = buildHomogeneousMemorySizes(maxMemInCluster, 2271);

                    break;
                case average:
                    num_processors = 312;
                    memorySizes = buildHomogeneousMemorySizes(maxMemInCluster*4, num_processors);
                    break;
                case fewBig:
                    num_processors = 16;
                    memorySizes = buildHomogeneousMemorySizes(maxMemInCluster*8, num_processors);
                    break;
            }
            break;
    }

    Cluster *cluster = new Cluster(num_processors, true);
    map<int, int> processor_speeds = Cluster::buildProcessorSpeeds(num_processors);
    cluster->setHomogeneousBandwidth(maxEdgesToMakespanWeights*10);
    Cluster::setFixedCluster(cluster);
    Cluster::getFixedCluster()->setMemorySizes(memorySizes);
}

void Cluster::buildHomogeneousCluster(double CCR, unsigned int num_processors, Tree *treeobj,
                                      HeterogeneousAdaptationMode mode) {
    // mode 0: build a cluster that uses all nodes with smallest memory
    // mode 1: build a cluster that uses 2/3 of nodes with middle amount of memory
    // mode 2: build a cluster that uses 1/3 of nodes with big memory

    double minMem;
    double maxoutd;
    schedule_traversal *temp_schedule;
    buildTreeDepHomBandwidths(CCR, num_processors, treeobj, minMem, maxoutd, temp_schedule);
    vector<double> memorySizes;
    switch (mode) {
        case noAdaptation:
            memorySizes = Cluster::buildHomogeneousMemorySizes(maxoutd, num_processors);
            break;
        case manySmall:
            memorySizes = Cluster::buildHomogeneousMemorySizes(min(maxoutd, minMem), num_processors);
            break;
        case average:
            memorySizes = Cluster::buildHomogeneousMemorySizes((maxoutd + minMem) / 2, num_processors * 2 / 3);
            break;
        case fewBig:
            memorySizes = Cluster::buildHomogeneousMemorySizes(max(maxoutd, minMem), num_processors / 3);
            break;
        default:
            throw std::runtime_error("Bad mode for creating homogeneous cluster.");

    }

    //Fix, for now we consider the non-homog cluster homogeneuos
    Cluster::getFixedCluster()->setMemorySizes(memorySizes);
    delete temp_schedule;
}

void Cluster::buildMemHetTreeDepCluster(double CCR, unsigned int num_processors, Tree *treeobj) {
    double minMem;
    double maxoutd;
    schedule_traversal *temp_schedule;
    buildTreeDepHomBandwidths(CCR, num_processors, treeobj, minMem, maxoutd, temp_schedule);

    vector<double> memorySizes = Cluster::build3LevelMemorySizes(min(maxoutd, minMem), max(maxoutd, minMem),
                                                                 num_processors);
    Cluster::getFixedCluster()->setMemorySizes(memorySizes);
    delete temp_schedule;
}


void Cluster::buildTreeDepHomBandwidths(double CCR, unsigned int num_processors, Tree *treeobj, double &minMem,
                                        double &maxoutd,
                                        schedule_traversal *&temp_schedule) {
    maxoutd = MaxOutDegree(treeobj, true);
    temp_schedule = new schedule_traversal();
    Cluster *cluster = new Cluster(num_processors, true);
    map<int, int> processor_speeds = Cluster::buildProcessorSpeeds(num_processors);
    cluster->SetBandwidth(CCR, treeobj);
    Cluster::setFixedCluster(cluster);
    Cluster::getFixedCluster()->SetBandwidth(CCR, treeobj);
}

Processor *Cluster::getFirstFreeProcessor() {
    for (vector<Processor *>::iterator iter = this->processors.begin(); iter < this->processors.end(); iter++) {
        if (!(*iter)->isBusy)
            return (*iter);
    }
    throw std::out_of_range("No free processor available anymore!");
}

//todo rewrite with flag
bool Cluster::hasFreeProcessor() {
    for (vector<Processor *>::iterator iter = this->processors.begin(); iter < this->processors.end(); iter++) {
        if (!(*iter)->isBusy)
            return true;
    }
    return false;
}

Processor *Cluster::getFirstFreeProcessorOrSmallest() {
    for (vector<Processor *>::iterator iter = this->processors.begin(); iter < this->processors.end(); iter++) {
        if (!(*iter)->isBusy)
            return (*iter);
    }
    return this->processors.back();
}

Processor *Cluster::getLastProcessor() {
    return this->processors.back();
}

void Processor::assignTask(Task *taskToBeAssigned) {
    this->assignedTask = taskToBeAssigned;
    taskToBeAssigned->setAssignedProcessor(this);
    this->isBusy = true;
}

void Processor::assignTaskId(unsigned int taskToBeAssigned) {
    this->assignedTaskId = taskToBeAssigned;
    this->isBusy = true;
}

int Processor::getAssignedTaskId() const {
    return assignedTaskId;
}

Task *Processor::getAssignedTask() const {
    return assignedTask;
}

void Cluster::assignTasksForIds(Tree *tree) {
    for (Processor *p: getProcessors()) {
        if (p->getAssignedTaskId() != -1 && p->getAssignedTask() == nullptr) {
            p->assignTask(tree->getTask(p->getAssignedTaskId()));
        }
    }

}






