from dataclasses import dataclass, field, InitVar
from re import sub
from sys import argv
from typing import Mapping, Union, List,Tuple
from unicodedata import category
from xml.dom import ValidationErr
import csv
from matplotlib import pyplot as plt
import numpy as np
from copy import deepcopy

@dataclass(frozen=True)
class TreeStatistic():
    name: str
    makespan:float

class State():
    def __init__(self,dimensions : List[str]):
        self.statemap : Mapping[str, str] = {}
        for dimension in dimensions:
            self.statemap[dimension] = ""
    
    def updateState (self,dimension: str, state: str):
        try:
            self.statemap[dimension] = state
        except KeyError:
            raise ValidationErr("This Dimension doesn't exist")

    def getFlatRepresentation(self)->str:
        flatlist = [f"{dimension}:{value}" for dimension, value in self.statemap.items()]
        flatlist.sort()
        return ",".join(flatlist)

    def constructFlatRepresentation(attributes:Mapping[str,str]):
        flatlist = [f"{dimension}:{value}" for dimension, value in attributes.items()]
        flatlist.sort()
        return ",".join(flatlist)
    
    def concatFlatRepresentation(state1:str, state2:str)->str:
        string = state1+","+state2
        print(string)
        flatlist = string.split(",")
        print(flatlist)
        flatlist.sort()
        return ",".join(flatlist)

    def reduceBy(dimension:str,states:List[str])->List[List[str]]:
        superStateToIndex:Mapping[str, List[str]] = {}
        print(states)
        for state in states:
            stateList = state.split(",")
            print(stateList)
            
            for curDimension in stateList:
                print(curDimension)
                print(dimension)
                if dimension in curDimension:
                    stateList.remove(curDimension)
                    break
            print("statelist:::::")
            print(stateList)
            superState = ",".join(stateList)
            if superState in superStateToIndex:
                superStateToIndex[superState].append(state)
            else: 
                superStateToIndex[superState] = [state]
                
        print(superStateToIndex)
        return superStateToIndex

    

class OutputAnalyszer():
    def __init__(self,argv):
        if len(argv)!=2:
            raise ValidationErr("Wrong amount of command line arguments")

        self.pathToFile = argv[1]
        self.treeNameToObjectMap: Mapping[str, TreeStatistic] = {}
        self.treesInState: Mapping[str, Mapping[str,TreeStatistic]] = {} # state -> Tree[]

        self.stateKeywords:Mapping[str,List[str]] = {            
            "level" :["2-level\n","3-level\n"],
            "treeCategory" : [],
            "clusteringMode" : ["heterogeneous\n","many small\n","avg avg\n","few big\n"]
        }
        self.state = State([key for key in self.stateKeywords])

    def readNextLine(self,file)->Tuple[str,str]:
        line = file.readline()
        while True:
            if line[0:2] == "&&":
                return ("newTree",line)
            if line[0:2] == "!!":
                return ("stateSwitch treeCategory", line)
            if line in self.stateKeywords["clusteringMode"]:
                return ("stateSwitch clusteringMode", line)
            if line in self.stateKeywords["level"]:
                return ("stateSwitch level", line)
            if not line:
                return ("eof", line)
            line = file.readline()
        

    def adaptState(self,dimension:str, line:str)->None:
        if dimension in ["cluteringMode","level"]:self.state.updateState(dimension,line)
        else:
            category = line.split("/")[-1]
            print(category)
            self.state.updateState(dimension, category)
                
        
    def updateTreeStatistic(self, line:str)->None:
        line = line.split(" ")
        treeName = line[1]
        makespan = float(line[2])
        stateString:str = self.state.getFlatRepresentation()
        newTree = TreeStatistic(treeName,makespan)
        if stateString not in self.treesInState: self.treesInState[stateString]= {treeName:newTree}
        else: self.treesInState[stateString][treeName]= newTree

    def getGeometricMean(self)->np.array:
        n = np.zeros(3)
        prod = np.ones(3)
        for _, treeobj in self.treeNameToObjectMap.items():
            for idx, mode in enumerate(["manySmall","fewBig","averageAverage"]):
                if treeobj.heterogeneousRatio[mode] != "not_computed":
                    n[idx] +=1
                    prod[idx] *=treeobj.heterogeneousRatio[mode]
        n = 1/n
        return np.power(prod,n)

    def getArithmeticMean(self)->np.array:
        n = np.zeros(3)
        sum = np.zeros(3)
        for _, treeobj in self.treeNameToObjectMap.items():
            for idx, mode in enumerate(["manySmall","fewBig","averageAverage"]):
                if treeobj.heterogeneousRatio[mode] != "not_computed":
                    n[idx] +=1
                    sum[idx] +=treeobj.heterogeneousRatio[mode]
        return sum/n   

    def analyze(self):
        with open(self.pathToFile) as file:
            line = ""
            while True:
                action, line = self.readNextLine(file)
                if action == "eof":
                    break
                elif action == "newTree":
                    self.updateTreeStatistic(line) 
                else: #action = stateSwitch
                    dimension = (action.split(" "))[1]
                    self.adaptState(dimension,line)

    def writeCSV(self, filePath = "out.csv"):
        superStates:Mapping[str:List[str]] = State.reduceBy("clusteringMode",self.treesInState.keys()) #returns list of list of states that are equialent except for clustering modes
        with open(filePath, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile, delimiter=',')
            writer.writerow(["Name", "Heterogeneous", "Many_Small", "Few_Big", "Avg_Avg","Ratio_Many_Small","Ratio_Few_Big","Ratio_Avg_Avg"])
            
            for superState in superStates:
                writer.writerow([])
                rows = []
                print(superState)
                rows.append(superState.split(","))
                
                hetTrees = self.treesInState[State.concatFlatRepresentation(superState, "clusteringMode:heterogeneous\n")] 
                
                for treename,tree in hetTrees.items():
                    print("KEYS :::::")
                    print(self.treesInState.keys())
                    print(self.treesInState[State.concatFlatRepresentation(superState, "clusteringMode:many Small\n")].items())
                    try:
                        manySmallMakespan = self.treesInState[State.concatFlatRepresentation(superState, "clusteringMode:many Small\n")][treename].makespan 
                        manySmallRatio = tree.makespan/manySmallMakespan
                    except KeyError:
                        manySmallMakespan = -1
                        manySmallRatio = -1
                    
                    try:
                        fewBigMakespan = self.treesInState[State.concatFlatRepresentation(superState, "clusteringMode:few Big\n")][treename].makespan 
                        fewBigRatio = tree.makespan/fewBigMakespan
                    except KeyError:
                        fewBigRatio = -1
                        fewBigMakespan =-1

                    try:
                        avgAvgMakespan= self.treesInState[State.concatFlatRepresentation(superState, "clusteringMode:averageAverage\n")][treename].makespan 
                        avgAvgRatio = tree.makespan/avgAvgMakespan
                    except:
                        avgAvgMakespan = -1
                        avgAvgRatio = -1
                    row = [treename,tree.makespan,manySmallMakespan,fewBigMakespan,avgAvgMakespan,manySmallRatio,fewBigRatio,avgAvgRatio]
                    rows.append(row)
                writer.writerows(rows)
                

        '''
        rows = []
        for treename, treeobj in  self.treeNameToObjectMap.items():
            manySmallRatio = treeobj.heterogeneousRatio["manySmall"]
            avgAvgRatio = treeobj.heterogeneousRatio["averageAverage"]
            fewBigRatio = treeobj.heterogeneousRatio["fewBig"]
            row = [treename,treeobj.heterogeneousMakespan,treeobj.manySmallMakespan,treeobj.fewBigMakespan,treeobj.averageAverageMakespan,manySmallRatio,fewBigRatio,avgAvgRatio]
            rows.append(row)

        with open(filePath, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile, delimiter=',')
            writer.writerow(["Name", "Heterogeneous", "Many_Small", "Few_Big", "Avg_Avg","Ratio_Many_Small","Ratio_Few_Big","Ratio_Avg_Avg"])
            writer.writerow([])
            writer.writerows(rows)
        
    def plot(self):
        x = np.arange(4)
        y_geom =self.getGeometricMean()
        y_arit = self.getArithmeticMean()
        print(y_geom)

        y_geom = np.concatenate((np.array([1]),y_geom),axis=0)
        y_arit = np.concatenate((np.array([1]),y_arit),axis=0)

        fig, ax = plt.subplots(1,2)
        plt.setp(ax, xticks = x,xticklabels=["Het","Many Small","Few Big","Avg Avg"])
        
        ax[0].bar(x, y_geom,color = ["r"]+["b"]*3)
        ax[0].set_title("Geometric Mean")
        ax[1].bar(x, y_arit,color = ["r"]+["b"]*3)
        ax[1].set_title("Arithmetic Mean")
        fig.suptitle("Ratios in comparison to heterogeneous clustering")
        plt.show()
'''        

o = OutputAnalyszer(argv)
o.analyze()
o.writeCSV()
#print(o.treesInState.keys())
#o.plot()


