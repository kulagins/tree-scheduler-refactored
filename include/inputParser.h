#ifndef inputParse_h
#define inputParse_h

#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <json.hpp>
#include "cluster.h"


using json = nlohmann::json;
using namespace std;

class InputParser {
protected:

    string workingDirectory;
    string pathToTree;
    string pathToCluster;
    ClusteringModes clusteringMode;

protected:
    bool buildSmallCluster;
    bool verbose;
public:

    InputParser(int argc, char **argv) {
        //cout << argc << endl;
        //cout << argv << endl;
        if (argc < 7) errorFunction(0);

        this->workingDirectory = argv[1];
        this->pathToTree = this->workingDirectory + argv[2];

        int cluMode = atoi(argv[3]);
        switch (cluMode) {
            case 0:
                this->clusteringMode = treeDependent;
                break;
            case 1:
                this->clusteringMode = staticClustering;
                break;

            default:
                cout<<"clustering mode "<<cluMode<<endl;
                errorFunction(2);
                break;
        }

        this->buildSmallCluster = (bool)atoi((argv[4]));
        this->verbose = (bool)atoi((argv[5]));
        this->pathToCluster = argv[6];
    }

    string getWorkingDirectory() {
        return this->workingDirectory;
    }

    string getPathToTreeList() {
        return this->pathToTree;
    }
    string getPathToCluster(){
        return this->pathToCluster;
    }
    ClusteringModes getClusteringMode() {
        return this->clusteringMode;
    }
    bool getVerbosity(){
        return verbose;
    }
    bool getBuildSmallClusters(){
        return buildSmallCluster;
    }

    void setClusterFromFile(double normedMemory){
        cout <<this->pathToCluster<<endl;
        ifstream inputFile(this->pathToCluster);
        json clusterDescription;
        inputFile >> clusterDescription;

        vector<unsigned int> processorCounts;
        vector<double> mems;
        vector<double> speeds;
        
        for (auto& element : clusterDescription["groups"]) {
            processorCounts.push_back(element["number"]);
            mems.push_back(element["memory"]);
            speeds.push_back(element["speed"]);
        }

        if (this->getClusteringMode() == treeDependent){
            for (auto it = begin(mems); it != end(mems); it++)
            {
                *it = (*it) * normedMemory;
            }
        }

        Cluster *cluster = new Cluster(&processorCounts,&mems,&speeds);
       cout << cluster->getPrettyClusterString()<<endl;
        Cluster::setFixedCluster(cluster);

    }
    void setClusterFromFileWithShrinkingFactor(double normedMemory, double shrinkingFactor){
        cout<<"small cluster"<<endl;
        ifstream inputFile(this->pathToCluster);
        json clusterDescription;
        inputFile >> clusterDescription;

        vector<unsigned int> processorCounts;
        vector<double> mems;
        vector<double> speeds;

        for (auto& element : clusterDescription["groups"]) {
            int numberProcessors = element["number"];
            processorCounts.push_back(ceil(numberProcessors / shrinkingFactor));
            mems.push_back(element["memory"]);
            speeds.push_back(element["speed"]);
        }

        if (this->getClusteringMode() == treeDependent){
            for (auto it = begin(mems); it != end(mems); it++)
            {
                *it = (*it) * normedMemory;
            }
        }

        Cluster *cluster = new Cluster(&processorCounts,&mems,&speeds);
      cout << cluster->getPrettyClusterString()<<endl;
        Cluster::setFixedCluster(cluster);

    }


    void errorFunction(int reason) {
        switch (reason) {
            case 0:
                cout << "Please give the appropriate amount of commant line arguments" << endl;
                exit(-1);
                break;
            case 1:
                cout << "Please give an appropriate value for the Heterogenity Level" << endl;
                exit(-1);
                break;
            case 2:
                cout << "Please give an appropriate value for the Clustering Mode" << endl;
                exit(-1);
                break;
            default:
                exit(-1);
                break;
        }
    }



};
#endif