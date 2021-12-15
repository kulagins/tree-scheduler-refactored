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

vector<double> Cluster::buildNLevelMemorySizes(vector<double> memories, vector<unsigned int> processorGroupSizes) {
    // cout << "minimal memory per processor" << minMem << ", maximal " << maxMem << endl;
    double cumulativeMem = 0;
    //double completeProcessorCount = std::accumulate(processorGroupSizes.begin(), processorGroupSizes.end(), 0);
    vector<double> memSizes;
    // memSizes.resize(completeProcessorCount);
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

void Cluster::buildStatic2LevelCluster(double maxMinMem, double maxEdgesToMakespanWeights) {
    vector<unsigned int> processorCounts;
    vector<double> mems;
    double maxMemInCluster = maxMinMem * 1.1;
    int num_processors;
    vector<double> memorySizes;

    num_processors = 50;
    processorCounts.insert(processorCounts.end(), {10, 40});
    mems.insert(mems.end(), {maxMemInCluster*8, maxMemInCluster});
    memorySizes = buildNLevelMemorySizes(mems, processorCounts);


    auto *cluster = new Cluster(num_processors, false);
    map<int, int> processor_speeds = Cluster::buildProcessorSpeeds(num_processors);
    cluster->setHomogeneousBandwidth(maxEdgesToMakespanWeights * 10);
    Cluster::setFixedCluster(cluster);
    Cluster::getFixedCluster()->setMemorySizes(memorySizes);
}

void Cluster::buildStatic3LevelCluster(double maxMinMem, double maxEdgesToMakespanWeights) {
    vector<unsigned int> processorCounts;
    vector<double> mems;
    double maxMemInCluster = maxMinMem * 1.1;
    int num_processors;
    vector<double> memorySizes;

    num_processors = 1475;
    processorCounts.insert(processorCounts.end(), {10, 10, 30});
    mems.insert(mems.end(), {maxMemInCluster*8, maxMemInCluster * 2, maxMemInCluster });
    memorySizes = buildNLevelMemorySizes(mems, processorCounts);


    auto *cluster = new Cluster(num_processors, false);
    map<int, int> processor_speeds = Cluster::buildProcessorSpeeds(num_processors);
    cluster->setHomogeneousBandwidth(maxEdgesToMakespanWeights * 10);
    Cluster::setFixedCluster(cluster);
    Cluster::getFixedCluster()->setMemorySizes(memorySizes);
}

void
Cluster::buildHomStatic2LevelCluster(double maxMinMem, double maxEdgesToMakespanWeights,
                                     HeterogeneousAdaptationMode adaptationMode) {
    vector<unsigned int> processorCounts;
    vector<double> mems;
    double maxMemInCluster = maxMinMem * 1.1;
    int num_processors;
    vector<double> memorySizes;

    switch (adaptationMode) {
        case manySmall:
            cout << "many small" << endl;
            num_processors = 50;
            memorySizes = buildHomogeneousMemorySizes(maxMemInCluster, num_processors);
            break;
        case average:
            throw "No average processors in a 2-step cluster!";
        case fewBig:
            cout << "big" << endl;
            num_processors = 10;
            memorySizes = buildHomogeneousMemorySizes(maxMemInCluster * 8, num_processors);
            break;
        default:
            throw "invalid Clustering-Adaption mode";
            break;
    }


    Cluster *cluster = new Cluster(num_processors, true);
    map<int, int> processor_speeds = Cluster::buildProcessorSpeeds(num_processors);
    cluster->setHomogeneousBandwidth(maxEdgesToMakespanWeights * 10);
    Cluster::setFixedCluster(cluster);
    Cluster::getFixedCluster()->setMemorySizes(memorySizes);
}

void
Cluster::buildHomStatic3LevelCluster(double maxMinMem, double maxEdgesToMakespanWeights,
                                     HeterogeneousAdaptationMode adaptationMode) {
    vector<unsigned int> processorCounts;
    vector<double> mems;
    double maxMemInCluster = maxMinMem * 1.1;
    int num_processors;
    vector<double> memorySizes;

    switch (adaptationMode) {
        case manySmall:
            cout << "many small" << endl;
            num_processors = 50;
            memorySizes = buildHomogeneousMemorySizes(maxMemInCluster, num_processors);
            break;
        case average:
            cout << "avg" << endl;
            num_processors = 20;
            memorySizes = buildHomogeneousMemorySizes(maxMemInCluster * 2, num_processors);
            break;
        case fewBig:
            cout << "big" << endl;
            num_processors = 10;
            memorySizes = buildHomogeneousMemorySizes(maxMemInCluster * 8, num_processors);
            break;
        default:
            throw "invalid Clustering-Adaption mode";
            break;
    }


    Cluster *cluster = new Cluster(num_processors, true);
    map<int, int> processor_speeds = Cluster::buildProcessorSpeeds(num_processors);
    cluster->setHomogeneousBandwidth(maxEdgesToMakespanWeights * 10);
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

void Processor::setAssignedTaskId(int assignedTaskId) {
    Processor::assignedTaskId = assignedTaskId;
}

void Processor::setAssignedTask(Task *assignedTask) {
    Processor::assignedTask = assignedTask;
}

double Processor::getOccupiedMemorySize() const {
    return occupiedMemorySize;
}

void Processor::setOccupiedMemorySize(double occupiedMemorySize) {
    Processor::occupiedMemorySize = occupiedMemorySize;
}

void Cluster::assignTasksForIds(Tree *tree) {
    for (Processor *p: getProcessors()) {
        if (p->getAssignedTaskId() != -1 && p->getAssignedTask() == nullptr) {
            p->assignTask(tree->getTask(p->getAssignedTaskId()));
        }
    }

}

string Cluster::getUsageString() {
    string out = "Used memory per processor \n";
    out += "Processor number; Memory Size; Makespan weight of assigned node; Cost of assigned node; Makespan cost of assigned node; Usage percentage \n";
    int i = 0;
    for (Processor *p: getProcessors()) {
        if (p->getAssignedTask() != nullptr) {
            out += i + " " + to_string(p->getMemorySize()) + " " +
                   to_string(p->getAssignedTask()->getMakespanWeight()) + " "
                   + to_string(p->getAssignedTask()->getCost()) + " " +
                   to_string(p->getAssignedTask()->getMakespanCost()) + " "
                   + to_string(p->getAssignedTask()->getCost() * 100 / p->getMemorySize()) + "% \n";

        } else {
            out += "No task assigned   0%";
        }
    }
    return out;

}

string Cluster::getShortUsageString() {
    double percentageUsedProcesors;
    string out = "";
    int i = 0;
    for (Processor *p: getProcessors()) {
        if (p->getAssignedTask() != nullptr) {
            percentageUsedProcesors++;
            out += i + " " + to_string(p->getMemorySize()) + " " +
                   to_string(p->getAssignedTask()->getMakespanWeight()) + " "
                   + to_string(p->getAssignedTask()->getCost()) + " " +
                   to_string(p->getAssignedTask()->getMakespanCost()) + " "
                   + to_string(p->getAssignedTask()->getCost() * 100 / p->getMemorySize()) + "% \n";

        } else {
            out += "No task assigned   0%";
        }
    }
    percentageUsedProcesors = percentageUsedProcesors * 100 / getProcessors().size();
    out += "percentage used: " + to_string(percentageUsedProcesors);
    return out;

}

void Cluster::clean() {
    for (Processor *p: getProcessors()) {
        p->setAssignedTask(nullptr);
        p->setAssignedTaskId(-1);
        p->isBusy = false;
    }

}





