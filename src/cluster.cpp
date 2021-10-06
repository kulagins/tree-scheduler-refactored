#include "cluster.h"

using namespace std;

Processor *Cluster::getFirstFreeProcessor()
{

    for (vector<Processor *>::iterator iter = this->processors.begin(); iter < processors.end(); iter++)
    {
        if(!(*iter)->isBusy)
        return (*iter);
    }
    return nullptr;
}