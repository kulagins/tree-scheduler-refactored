#include "cluster.h"

using namespace std;

static std::map<int, int> buildProcessorSpeeds(int num_processors)
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






