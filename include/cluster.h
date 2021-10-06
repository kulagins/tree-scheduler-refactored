#include <ostream>
#include <iostream>
#include <stdio.h>
#include <list>

#include <vector>
#include <forward_list>
#include <assert.h>
#include <map>

using namespace std;

class Processor
{
protected:
    double memorySize;
    double processorSpeed;
    int assignedTask;

public:
    bool isBusy;
    Processor()
    {
        this->memorySize = 0;
        this->processorSpeed = 1;
        isBusy = false;
        assignedTask = -1;
    }
    Processor(double memorySize)
    {
        this->memorySize = memorySize;
        this->processorSpeed = 1;
        isBusy = false;
        assignedTask = -1;
    }
    Processor(double memorySize, double processorSpeed)
    {
        this->memorySize = memorySize;
        this->processorSpeed = processorSpeed;
        isBusy = false;
        assignedTask = -1;
    }
    double getMemorySize()
    {
        return memorySize;
    }
    double getProcessorSpeed()
    {
        return processorSpeed;
    }
    void assignTask(int taskId){
        this->assignedTask = taskId;
        this->isBusy = true;
    }
};

class Cluster
{
protected:
    bool isMemoryHomogeneous;
    bool isProcessorHomogeneous;
    bool isBandwidthHomogenenous;
    vector<Processor *> processors;
    vector<vector<double>> bandwidths;

public:
    Cluster()
    {
        this->isMemoryHomogeneous = this->isProcessorHomogeneous = this->isBandwidthHomogenenous = true;

        processors.resize(3, new Processor());
        bandwidths.resize(3);
        for (unsigned long i = 0; i < bandwidths.size(); i++)
        {
            bandwidths.at(i).resize(3, 1);
        }
    }
    Cluster(double clusterSize, bool isMemoryHomogeneous)
    {
        this->isMemoryHomogeneous = isMemoryHomogeneous;
        this->isProcessorHomogeneous = this->isBandwidthHomogenenous = true;

        processors.resize(clusterSize);
        bandwidths.resize(clusterSize);
        for (unsigned long i = 0; i < processors.size(); i++)
        {
            bandwidths.at(i).resize(clusterSize, 1);
        }
    }
    Cluster(vector<double> memorySizes)
    {
        this->isMemoryHomogeneous = false;
        this->isProcessorHomogeneous = this->isBandwidthHomogenenous = true;

        processors.resize(memorySizes.size());
        for (unsigned long i = 0; i < bandwidths.size(); i++)
        {
            //TODO init only upper half
            processors.at(i) = new Processor(memorySizes.at(i));
        }
        bandwidths.resize(memorySizes.size());
        for (unsigned long i = 0; i < bandwidths.size(); i++)
        {
            //TODO init only upper half
            bandwidths.at(i).resize(memorySizes.size(), 1);
        }
    }

    ~Cluster()
    {
        bandwidths.resize(0);

        for (unsigned long i = 0; i < processors.size(); i++)
        {
            delete processors.at(i);
        }
    }

public:
    bool isHomogeneous()
    {
        return isMemoryHomogeneous && isProcessorHomogeneous && isBandwidthHomogenenous;
    }

    vector<Processor *> getProcessors()
    {
        return this->processors;
    }
    double getBandwidthBetween(int firstProcessor, int secondProcessor)
    {
        return this->bandwidths.at(firstProcessor).at(secondProcessor);
    }
    Processor * getFirstFreeProcessor();
};