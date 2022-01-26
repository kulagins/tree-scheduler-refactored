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
    double maxMemInCluster = maxMinMem;
    int num_processors;
    vector<double> memorySizes;

    num_processors = 27;
    processorCounts.insert(processorCounts.end(), {18, 9});
    mems.insert(mems.end(), {maxMemInCluster * 2, maxMemInCluster});
    memorySizes = buildNLevelMemorySizes(mems, processorCounts);

    BuildFixedClusterWithMemories(maxEdgesToMakespanWeights, num_processors, memorySizes);
    Cluster::getFixedCluster()->isMemoryHomogeneous = false;
}


void Cluster::buildStatic3LevelCluster(double maxMinMem, double maxEdgesToMakespanWeights) {
    vector<unsigned int> processorCounts;
    vector<double> mems;
    double maxMemInCluster = maxMinMem;
    int num_processors;
    vector<double> memorySizes;

    num_processors = 27;
    processorCounts.insert(processorCounts.end(), {9, 9, 9});
    mems.insert(mems.end(), {maxMemInCluster * 2, maxMemInCluster * 1.5, maxMemInCluster});
    memorySizes = buildNLevelMemorySizes(mems, processorCounts);


    BuildFixedClusterWithMemories(maxEdgesToMakespanWeights, num_processors, memorySizes);
    Cluster::getFixedCluster()->isMemoryHomogeneous = false;
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


void Cluster::BuildFixedClusterWithMemories(double maxEdgesToMakespanWeights, int num_processors,
                                            vector<double> &memorySizes) {
    Cluster *cluster = new Cluster(num_processors, true);
    map<int, int> processor_speeds = buildProcessorSpeeds(num_processors);
    cluster->setHomogeneousBandwidth(maxEdgesToMakespanWeights);
    setFixedCluster(cluster);
    getFixedCluster()->setMemorySizes(memorySizes);
}

Processor *Cluster::getBiggestFreeProcessor() {
    sort(this->processors.begin(), this->processors.end(),
         [](Processor *lhs, Processor *rhs) { return lhs->getMemorySize() > rhs->getMemorySize(); });
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
    this->assignedTaskId = taskToBeAssigned->getId();
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
        p->setOccupiedMemorySize(0);
    }

}

Processor *Cluster::smallestFreeProcessorFitting(double requiredMem) {
    //TODO method for only free processors
    int min = std::numeric_limits<int>::max();
    Processor *minProc = nullptr;
    for (Processor *proc: (this->getProcessors())) {
        if (proc->getMemorySize() >= requiredMem && !proc->isBusy && min > proc->getMemorySize()) {
            min = proc->getMemorySize();
            minProc = proc;
        }
    }
    return minProc;
}

bool Cluster::areAllUnassigned(Task *currentQNode, const Tree *tree) {
    vector<Task *> *childrenvector = currentQNode->getParent()->getChildren();
    bool mergeThreeNodes = (currentQNode->getChildren()->empty()) & (childrenvector->size() == 2);
    bool areAllUnassigned = true;
    if (mergeThreeNodes) {
        areAllUnassigned = tree->getTask(childrenvector->front()->getOtherSideId())->getAssignedProcessor() == NULL &&
                           tree->getTask(childrenvector->back()->getOtherSideId())->getAssignedProcessor() == NULL &&
                           tree->getTask(currentQNode->getParent()->getOtherSideId())->getAssignedProcessor() == NULL;
    } else {

        areAllUnassigned = tree->getTask(currentQNode->getOtherSideId())->getAssignedProcessor() == NULL &&
                           tree->getTask(currentQNode->getParent()->getOtherSideId())->getAssignedProcessor() == NULL;
    }
    return areAllUnassigned;
}

Processor *Cluster::findSmallestFittingProcessorForMerge(Task *currentQNode, const Tree *tree, double requiredMemory) {
    Processor *optimalProcessor = nullptr;
    double optimalMemorySize = std::numeric_limits<double>::max();
    vector<Task *> *childrenvector = currentQNode->getParent()->getChildren();
    bool mergeThreeNodes = (currentQNode->getChildren()->empty()) & (childrenvector->size() == 2);
    vector<Processor *> eligibleProcessors;
    if (mergeThreeNodes) {
        eligibleProcessors.push_back(tree->getTask(childrenvector->front()->getOtherSideId())->getAssignedProcessor());
        eligibleProcessors.push_back(tree->getTask(childrenvector->back()->getOtherSideId())->getAssignedProcessor());
        eligibleProcessors.push_back(
                tree->getTask(currentQNode->getParent()->getOtherSideId())->getAssignedProcessor());
    } else {
        eligibleProcessors.push_back(tree->getTask(currentQNode->getOtherSideId())->getAssignedProcessor());
        eligibleProcessors.push_back(
                tree->getTask(currentQNode->getParent()->getOtherSideId())->getAssignedProcessor());
    }
    eligibleProcessors.erase(std::remove_if(eligibleProcessors.begin(),
                                            eligibleProcessors.end(),
                                            [](Processor *proc) { return proc == NULL; }),
                             eligibleProcessors.end());
    for (auto eligibleProcessor: eligibleProcessors) {
        if (eligibleProcessor != NULL && eligibleProcessor->getMemorySize() > requiredMemory &&
            eligibleProcessor->getMemorySize() < optimalMemorySize) {
            optimalProcessor = eligibleProcessor;
            optimalMemorySize = eligibleProcessor->getMemorySize();
        }
    }
    return optimalProcessor;
}

