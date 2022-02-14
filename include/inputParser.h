#ifndef inputParse_h
#define inputParse_h

#include <sstream>
#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <json.hpp>
#include <vector>
#include "cluster.h"



using json = nlohmann::json;
using namespace std;

class InputParser {
protected:

    string workingDirectory;
    string pathToTree;
    string clusterListDir;
    ClusteringModes clusteringMode;
    vector<string> *clusterList;
    vector<string>::iterator clusterIterator;

protected:
    bool buildSmallCluster;
    bool verbose;
public:

    InputParser(int argc, char **argv) {
        //cout << argc << endl;
        //cout << argv << endl;
        if (argc < 7) errorFunction(0);

        this->workingDirectory = argv[1];
        cout << workingDirectory<<endl;
        this->pathToTree = this->workingDirectory + argv[2];
        string clusterlistPath = argv[3];
        this->clusterList = initClusterList(clusterlistPath);

        
        resetClusterIterator();

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

        this->buildSmallCluster = (bool)atoi((argv[5]));
        this->verbose = (bool)atoi((argv[6]));
        
    }

    ~InputParser(){
        delete this->clusterList;
    }

    string getWorkingDirectory() {
        return this->workingDirectory;
    }

    string getPathToTreeList() {
        return this->pathToTree;
    }
    string getPathToCluster(){
        return this->clusterListDir+ *(this->clusterIterator);
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
        cout <<this->getPathToCluster()<<endl;
        ifstream inputFile(this->getPathToCluster());
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
       //cout << cluster->getPrettyClusterString()<<endl;
        Cluster::setFixedCluster(cluster);

    }
    void setClusterFromFileWithShrinkingFactor(double normedMemory, double shrinkingFactor){
        cout<<"small cluster"<<endl;
        ifstream inputFile(this->getPathToCluster());
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

    vector<string> * initClusterList(string path){
        vector<string> *list = new vector<string>();
        if(path.substr(path.find_last_of(".") + 1) == "json"){
            
            this->clusterListDir = path.substr(0,path.find_last_of("/")+1); 
            list->push_back(path.substr(path.find_last_of("/")+1));
            return list;
        } 

        ifstream OpenFile(path);
        string line;
        stringstream line_stream;

        while (getline(OpenFile, line)) {
            line_stream.clear();
            line_stream.str(line);

            if (!line_stream.str().empty()) {
                string p;
                line_stream >> p;
                list->push_back(p);
            }

        }
        OpenFile.close();
        this->clusterListDir = path.substr(0,path.find_last_of("/")+1);
        return list;
    }

    bool nextCluster(){
        this->clusterIterator ++;
        return (this->clusterIterator != this->clusterList->end());
    }

    void resetClusterIterator(){
        this->clusterIterator = this->clusterList->begin();
        if (this->clusterIterator == this->clusterList->end()){
            throw "No Cluster Input Given";
        }
    }
};



#endif