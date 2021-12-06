#include "../include/cluster.h"
#include "../include/tree.h"

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






