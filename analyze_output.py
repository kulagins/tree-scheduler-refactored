from dataclasses import dataclass, field, InitVar
from sys import argv
from typing import Mapping, Union
from xml.dom import ValidationErr
import csv
'''
from matplotlib import pyplot as plt
import numpy as np
'''
@dataclass
class TreeStatistic():
    name: str
    heterogeneousMakespan:float = 0
    manySmallMakespan:float = 0
    fewBigMakespan:float = 0
    averageAverageMakespan:float = 0
    heterogeneousRatio: Mapping[str, Union[float,str]] = None

    def __post_init__(self):
        self.heterogeneousRatio:Mapping[str, float] ={} 

    def setHeterogeneous(self, value:float)->None:
        self.heterogeneousMakespan = value
    
    def setManySmall(self,value:float)->None:
        self.manySmallMakespan = value
        self.heterogeneousRatio["manySmall"] = self.heterogeneousMakespan/value if value != -1 else "not_computed"

    def setfewBig(self,value:float)->None:
        self.fewBigMakespan = value
        self.heterogeneousRatio["fewBig"] = self.heterogeneousMakespan/value if value != -1 else "not_computed"

    def setAverageAverage(self,value:float)->None:
        self.averageAverageMakespan = value
        self.heterogeneousRatio["averageAverage"] = self.heterogeneousMakespan/value if value != -1 else "not_computed"


class OutputAnalyszer():
    def __init__(self,argv):
        if len(argv)!=2:
            raise ValidationErr("Wrong amount of command line arguments")
        self.pathToFile = argv[1]
        self.treeNameToObjectMap: Mapping[str, TreeStatistic] = {}
        self.modeSwitchKeywords = ["heterogeneous\n","many small\n","avg avg\n","few big\n",""]

    def readNextLine(self,file):
        line = file.readline()
        while line[0:2] != "&&" and line not in self.modeSwitchKeywords:
            line = file.readline()
        return line

    def adaptState(self,line)->int:
        for idx, keywrd in enumerate(self.modeSwitchKeywords):
            if keywrd == line:
                return idx
        
    def updateTreeStatistic(self,state:int, line:str)->None:
        line = line.split(" ")
        treeName = line[1]
        if state == 0:
            treeStatistic = TreeStatistic(treeName, heterogeneousMakespan=float(line[2]))
            self.treeNameToObjectMap[treeName] = treeStatistic
        else:
            try:
                treeStatistic = self.treeNameToObjectMap[treeName]
                makespanValue = float(line[2])
                if state == 1: treeStatistic.setManySmall(makespanValue)
                elif state == 2: treeStatistic.setfewBig(makespanValue)
                elif state == 3: treeStatistic.setAverageAverage(makespanValue)
                else: raise ValidationErr("state not supported")

            except KeyError:
                raise ValidationErr(f"The tree {treeName} doesn't have all entries.")


    def analyze(self):
        state:int = 0
        with open(self.pathToFile) as file:
            line = "string"
            while line:
                line = self.readNextLine(file)
                if line in self.modeSwitchKeywords:
                    state = self.adaptState(line)
                else: self.updateTreeStatistic(state,line)

    def writeCSV(self, filePath = "out.csv"):
        
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
'''
    def plot(self):
        x = np.arange(4)
'''
        

o = OutputAnalyszer(argv)
o.analyze()
o.writeCSV()


