#include "../include/cluster.h"

using namespace std;

Cluster *  Cluster::fixedCluster = NULL;

vector<double> Cluster::buildMemorySizes(double maxoutd, double minMem, unsigned int num_processors)
{
    cout << "max deg " << maxoutd << ", MinMem " << minMem << endl;
    double cumulativeMem = 0;
    vector<double> memSizes(num_processors);
    memSizes.resize(num_processors);
    maxoutd = maxoutd / 4;
    //cout << "minProc " << maxoutd << " " << (maxoutd + minMem) / 2 << " " << minMem << endl;
    for (int k = 0; k < num_processors / 3; k++)
    {
        memSizes[k] = maxoutd; 
        cumulativeMem += memSizes[k];
    }
    for (int k = num_processors / 3; k < 2 * num_processors / 3; k++)
    {
        memSizes[k] = (maxoutd + minMem) / 2;
        cumulativeMem += memSizes[k];
    }
    for (int k = 2 * num_processors / 3; k < num_processors; k++)
    {
        memSizes[k] = minMem;
        cumulativeMem += memSizes[k];
    }
    cout << "cumulative mem in system: " << cumulativeMem << endl;
    return memSizes;
}

std::map<int, int> Cluster::buildProcessorSpeeds(int num_processors)
{
    std::map<int, int> procSpeeds;
    for (int k = 0; k < num_processors / 3; k++)
    {
        procSpeeds.insert(pair<int, int>(k, 1));
    }
    for (int k = num_processors / 3; k < 2 * num_processors / 3; k++)
    {
        procSpeeds.insert(pair<int, int>(k, 2));
    }
    for (int k = 2 * num_processors / 3 + 1; k < num_processors; k++)
    {
        procSpeeds.insert(pair<int, int>(k, 3));
    }

    return procSpeeds;
}

Processor* Cluster::getFirstFreeProcessor()

{
    for (vector<Processor*>::iterator iter = this->processors.begin(); iter < this->processors.end(); iter++) {
        if (!(*iter)->isBusy)
            return (*iter);
    }
    throw std::out_of_range("No free processor available anymore!");
}

void Cluster::SetBandwidth(double CCR, unsigned long tree_size, double *ewghts, double *timewghts)
{
    double sum_edges = 0;
    double sum_weights = 0;
    for (unsigned int i = 1; i <= tree_size; ++i)
    {
        sum_edges = sum_edges + ewghts[i];
        sum_weights = sum_weights + timewghts[i];
    }
    if(this->isHomogeneous()){
        this->setHomogeneousBandwidth(sum_edges / (sum_weights * CCR));
    }

}







