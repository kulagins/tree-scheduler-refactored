from dataclasses import dataclass
from sys import argv
from typing import Mapping
from xml.dom import ValidationErr


@dataclass
class TreeStatistic():
    name: str
    heterogeneousMakespan:float = 0
    manySmallMakespan:float = 0
    fewBigMakespan:float = 0
    averageAverageMakespan:float = 0

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
            print(line[0:2])
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
                if state == 1: treeStatistic.manySmallMakespan = makespanValue
                elif state == 2: treeStatistic.fewBigMakespan = makespanValue
                elif state == 3: treeStatistic.averageAverageMakespan = makespanValue
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

o = OutputAnalyszer(argv)
o.analyze()
for element in o.treeNameToObjectMap:
    print(o.treeNameToObjectMap[element])



