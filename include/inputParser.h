#ifndef inputParse_h
#define inputParse_h

#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <json.hpp>
#include "inputEnums.h"
#include "cluster.h"


using json = nlohmann::json;
using namespace std;

class InputParser {
protected:

    string workingDirectory;
    string pathToTree;
    HeterogenityLevels heterogenityLevel;
    ClusteringModes clusteringMode;
    HeterogeneousAdaptationMode adaptationMode;
protected:
    double clusteringDependendArg1;
    double clusteringDependendArg2;
public:
    InputParser(int argc, char **argv) {

        //cout << argc << endl;
        //cout << argv << endl;
        if (argc < 7) errorFunction(0);

        this->workingDirectory = argv[1];
        this->pathToTree = this->workingDirectory + argv[2];

        int hetLevel = atoi(argv[3]);
        switch (hetLevel) {
            case 0:
                this->heterogenityLevel = homogeneus;
                break;
            case 1:
                this->heterogenityLevel = memoryHeteregeneus;
                break;
            case 2:
                this->heterogenityLevel = heterogeneus;
                break;

            default:
                errorFunction(1);
                break;
        }

        int cluMode = atoi(argv[4]);
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

        this->clusteringDependendArg1 = atof(argv[5]);
        this->clusteringDependendArg2 = atof(argv[6]);
        if (argv[7]) {
            int adaptMode = atof(argv[7]);
            switch (adaptMode) {
                case 0:
                    this->adaptationMode = noAdaptation;
                    break;
                case 1:
                    this->adaptationMode = manySmall;
                    break;
                case 2:
                    this->adaptationMode = average;
                    break;
                case 3:
                    this->adaptationMode = fewBig;
                    break;

                default:
                    errorFunction(1);
                    break;
            }

        
        }

    }

    string getWorkingDirectory() {
        return this->workingDirectory;
    }

    string getPathToTreeList() {
        return this->pathToTree;
    }

    HeterogenityLevels getHeterogenityLevel() {
        return this->heterogenityLevel;
    }

    ClusteringModes getClusteringMode() {
        return this->clusteringMode;
    }
    HeterogeneousAdaptationMode getAdaptationMode() const {
        return adaptationMode;
    }

    double getCCR() {
        switch (this->clusteringMode) {
            case treeDependent:
                return this->clusteringDependendArg1;
                break;
            default:
                throw "This parameter doesn't exisit within the chosen clustering mode";
                break;
        }
    }

    double getNPR() {
        switch (this->clusteringMode) {
            case treeDependent:
                return this->clusteringDependendArg2;
                break;
            default:
                throw "This parameter doesn't exisit within the chosen clustering mode";
                break;
        }
    }

    int getStaticClusterConfigurationNumber() {
        switch (this->clusteringMode) {
            case staticClustering:
                return this->clusteringDependendArg1;
            default:
                return this->clusteringDependendArg1;
                //throw "This parameter doesn't exisit within the chosen clustering mode";
        }
    }

    void setClusterFromFile(string filePath, double normedMemory){
        ifstream inputFile(filePath);
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