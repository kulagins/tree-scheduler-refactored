from random import Random, randint, uniform
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


class RandomTreeGenerator:

    def initializeTreeEdges(self,prufer):
        numberVertices = len(prufer) + 2	
        graph = self.initializeRandomNodes(numberVertices)

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


        '''
        for node in graph:
            for child in node.neighbours:
                if node.id < child.id:
                    print(f"({node.id},{child.id})", sep="", end = "")
        '''
        return graph


    def createRandomTree(self,n):
        prufer = []
        length = n - 2
        for i in range(length):
            prufer.append(randint(0,length+1)+1)

        graph =  self.initializeTreeEdges(prufer)
        self.topologicalOrder(graph)
        return graph

    def topologicalOrder(self,graph):
        '''
        print()                
        for node in graph:
            print(f"{node.id},", sep="", end = "")
        print()
        '''

        orderedIds = []
        numberNodes = len(graph)
        nodeVisited = {node:False for node in graph}
        queue = [graph[0]]
        while queue:
            node = queue.pop(0)
            nodeVisited[node] = True
            orderedIds.append(node.id)
            for child in node.neighbours:
                if not nodeVisited[child]:
                    queue.append(child)
        
        for newId,oldId in zip([node.id for node in graph],orderedIds):
            graph[oldId-1].id = numberNodes - newId+1
        '''
        for node in orderedIds:
            print(f"{node},", sep="", end = "")
        print()
        '''
        cnt = 0
        for node in graph:
            for child in node.neighbours:                
                if child.id<node.id:
                    child.parent = node
        '''
        for node in graph:
            for child in node.neighbours:
                if node.id > child.id:
                    print(f"({node.id},{child.id})", sep="", end = "")
        print()
        '''

        
    def initializeRandomNodes(self,numberNodes):
        nodeWeightSeed = randint(11,200)
        makespanWeightSeed = randint(1000,5000)
        nodes = []
        for i in range (1,numberNodes+1):
            nodeWeight = nodeWeightSeed + randint(-10,10)
            edgeWeight = uniform(0.01,0.9)
            makespanWeight = makespanWeightSeed + randint(-100,100)
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



treelengths = [15]+ [2000]*4 + [4000]*4 + [10000]*2 + [20000]
r = RandomTreeGenerator()
for index,length in enumerate(treelengths):
    g = r.createRandomTree(length)
    with open(f"random_tree_{index}_{length}","w") as file:
        file.write(r.prettyString(g))

