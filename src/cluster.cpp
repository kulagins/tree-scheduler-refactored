#include "cluster.h"

using namespace std;

Processor* getFirstFreeProcessor()

{
    for (vector<Processor*>::iterator iter = this->processors.begin(); iter < this->processors.end(); iter++) {
        if (!(*iter)->isBusy)
            return (*iter);
    }
    throw std::out_of_range("No free processor available anymore!");
}
