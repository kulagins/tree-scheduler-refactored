from random import Random, randint, shuffle, uniform
import os
from typing import Tuple
class Node:
    def __init__(self,id, nodeWeight = 0, edgeWeight = 0, makespanWeight = 0) -> None:
        self.id = id
        self.nodeWeight = nodeWeight
        self.edgeWeight = edgeWeight
        self.makespanWeight = makespanWeight
        self.neighbours = []
        self.parent = None
    
    def addNeighbour(self, neighbour):
        self.neighbours.append(neighbour)


"""
Dataclass to store the values that describe the randomness
numberNodes: number of nodes in the Tree
avgNodeBorders: between which borders should the average nodeweight lie
avgEdgeBorders: between which borders should the average edgeweight lie
avgMakespanBorders: between which borders should the average edgeweight lie
nodeSpread: all the nodweigths will lie in an area of +- nodeSpread around the average nodeweight
makespanSpread: all the makespaneigths will lie in an area of +- makespanSpread around the average makespanweight
avgChildren: average number of children for each node
childrenSpread: the number of children for each node will lie in an area of +- childrenSpread around the average avgChildren
"""
class RandomParameters:
    def __init__(self, numberNodes : int, avgNodeBorders : Tuple[int,int], avgMakespanBorders : Tuple[int,int], avgEdgeBorders:Tuple[float,float], nodeSpread : int, makespanSpread : int, avgChildren :int = None, childrenSpread :int = None) -> None:
        self.numberNodes = numberNodes
        self.avgNodeBorders = avgNodeBorders
        self.avgEdgeBorders = avgEdgeBorders
        self.avgMakespanBorders = avgMakespanBorders
        self.nodeSpread = nodeSpread
        self.makespanSpread = makespanSpread
        self.avgChildren = avgChildren
        self.childrenSpread = childrenSpread



class RandomTreeGenerator:
    def __init__(self, randomParameters) -> None:
        self.randomParameters = randomParameters

    def initializeTreeEdges(self,prufer):
        numberVertices = len(prufer) + 2	
        graph = self.initializeRandomNodes()

        occurancesInPrufer = [0] * numberVertices
        for node in prufer:
            occurancesInPrufer[node-1]+=1
        
        for node_p in prufer:
            for node_b in range(1,numberVertices+1):
                if occurancesInPrufer[node_b-1] == 0:
                    occurancesInPrufer[node_p-1] -= 1
                    occurancesInPrufer[node_b-1] = -1
                    graph[node_p-1].addNeighbour(graph[node_b-1])
                    graph[node_b-1].addNeighbour(graph[node_p-1])
                    break
        
        firstFound = False
        for index, value in enumerate(occurancesInPrufer):
            if value > -1:
                if firstFound:
                    node_p = index+1
                    firstFound = True
                else:
                    node_b = index+1
                    graph[node_p-1].addNeighbour(graph[node_b-1])
                    graph[node_b-1].addNeighbour(graph[node_p-1])
        return graph

    def createUniformRandomPrufer(self):
        n = self.randomParameters.numberNodes
        prufer = []
        length = n - 2
        for i in range(length):
            prufer.append(randint(0,length+1)+1)
        return prufer
    
    def createParameterizedChildAmtPrufer(self):
        n = self.randomParameters.numberNodes
        prufer = []
        freeSlots = n - 2
        index = 1
        while freeSlots > 0:
            numberOccurences = self.randomParameters.avgChildren + randint(-self.randomParameters.childrenSpread,self.randomParameters.childrenSpread)
            numberOccurences = min(numberOccurences,freeSlots)
            prufer.extend([index]*numberOccurences)
            index +=1
            freeSlots-= numberOccurences
        shuffle(prufer)
        return prufer
        

    def createRandomTree(self):
        if not self.randomParameters.avgChildren:
            prufer = self.createUniformRandomPrufer()
        else:
            prufer = self.createParameterizedChildAmtPrufer()
        
        graph =  self.initializeTreeEdges(prufer)
        self.topologicalOrder(graph)
        return graph

    def topologicalOrder(self,graph):
        orderedIds = []
        numberNodes = len(graph)
        nodeVisited = {node:False for node in graph}
        queue = [graph[randint(0,len(graph)-1)]]
        while queue:
            node = queue.pop(0)
            nodeVisited[node] = True
            orderedIds.append(node.id)
            for child in node.neighbours:
                if not nodeVisited[child]:
                    queue.append(child)
        
        for newId,oldId in zip([node.id for node in graph],orderedIds):
            graph[oldId-1].id = numberNodes - newId+1
        cnt = 0
        for node in graph:
            for child in node.neighbours:                
                if child.id<node.id:
                    child.parent = node

        
    def initializeRandomNodes(self):
        params = self.randomParameters
        nodeWeightSeed = randint(params.avgNodeBorders[0],params.avgNodeBorders[1])
        makespanWeightSeed = randint(params.avgMakespanBorders[0],params.avgMakespanBorders[1])
        nodes = []
        for i in range (1,params.numberNodes+1):
            nodeWeight = nodeWeightSeed + randint(-params.nodeSpread,params.nodeSpread)
            edgeWeight = uniform(params.avgEdgeBorders[0],params.avgEdgeBorders[1])
            makespanWeight = makespanWeightSeed + randint(-params.makespanSpread,params.makespanSpread)
            node = Node(i,nodeWeight,edgeWeight,makespanWeight)
            nodes.append(node)
        return nodes
    
    def prettyString(self,graph):
        indexMap = {}
        for node in graph:
            if node.parent:
                parent = node.parent.id
            else:
                parent = 0
            indexMap[node.id] = f"{node.id} {parent} {node.nodeWeight} {node.edgeWeight} {node.makespanWeight} \n"
        
        returnString = ""
        for i in range(1, len(graph)+1):
            returnString+=indexMap[i]
        return returnString

# the treelengths list should always be sorted as the file naming depends on that for now
treelengths = [2000]*5+[4000]*5+[10000]*5+[20000]*5+[30000]*5+[50000]*5
#treelengths = [2000]*5
parameterList = [
    #RandomParameters(10,(11,200),(1000,5000),(0.01,0.9),10,1000),
    
	RandomParameters(length,(11,200),(1000,5000),(0.01,0.9),10,500),
	RandomParameters(length,(11,200),(1000,5000),(0.01,0.9),10,500, 3,1),
	RandomParameters(length,(11,200),(1000,5000),(0.01,0.9),10,500, 20,4),
	RandomParameters(length,(11,200),(1000,5000),(0.01,0.9),10,500, 10,2),
	RandomParameters(length,(1100,20000),(1000,5000),(0.01,0.9),1000,500),
	RandomParameters(length,(11,200),(1000,5000),(1.0,90.0),10,500),
	RandomParameters(length,(11,200),(100000,500000),(0.01,0.9),10,50000),
	RandomParameters(length,(1100,20000),(100000,500000),(1.0,90.0),1000,50000),
	RandomParameters(length,(1,20),(100,500),(0.001,0.09),1,5)
]
description = ["old_trees", "3-children", "20-children", "10-children", "large_nw", "large_ew", "large_msw", "all_large", "all_small"]

index = 1
for i,length in enumerate(treelengths):
    index = (1 if i == 0 or treelengths[i-1]!=length else index+ 1)  
    for j, params in enumerate(parameterList):   
        try: os.makedirs(f"{description[j]}/size_{length}")
        except FileExistsError:
            pass
        try:
            file=open(f"{description[j]}/size_{length}/random_tree_{index}","x")
            params.numberNodes = length
            r = RandomTreeGenerator(params)
            g = r.createRandomTree()
            file.write(r.prettyString(g))
            print(f"generated {description[j]}/size_{length}/random_tree_{index}") 
        except FileExistsError:
            print(f"skipped {description[j]}/size_{length}/random_tree_{index}") 
            
"""

# Parameters for different types of tree-generation:

old Trees
RandomParameters(length,(11,200),(1000,5000),(0.01,0.9),10,500),

3-children:
RandomParameters(length,(11,200),(1000,5000),(0.01,0.9),10,500, 3,1),

20-children:
RandomParameters(length,(11,200),(1000,5000),(0.01,0.9),10,500, 20,4),

10-children:
RandomParameters(length,(11,200),(1000,5000),(0.01,0.9),10,500, 10,2),

"large_node_weights",
RandomParameters(length,(1100,20000),(1000,5000),(0.01,0.9),1000,500),

"large_edge_weights",
RandomParameters(length,(11,200),(1000,5000),(1.0,90.0),10,500),

"large_makespan_weights",
RandomParameters(length,(11,200),(100000,500000),(0.01,0.9),10,50000),

"all_large",
RandomParameters(length,(1100,20000),(100000,500000),(1.0,90.0),1000,50000),

"all_small"
RandomParameters(length,(1,20),(100,500),(0.001,0.09),1,5)
"""
