from cgitb import text
import simexpal
import numpy as np
from matplotlib import pyplot as plt, ticker
from sys import argv
import csv
from typing import List, Set, Tuple
class OutputObject:
    def __init__(self):
        self.tree_category = ""
        self.tree_names = []
        self.cluster_names = []
    
    def init_makespan_matrix(self,num_trees, num_clusters):
        self.makespan_matrix = np.zeros((num_trees,num_clusters), dtype=np.float64)
        self.num_clusters = num_clusters
        self.num_trees = num_trees

    def init_ratio_matrix(self):
        self.ratio_matrix = np.zeros((self.num_trees, self.num_clusters))
        first_col = self.makespan_matrix[:,0]
        for i in range (0,len(self.makespan_matrix[1]) ):
            self.ratio_matrix[:,i] = np.divide(first_col, self.makespan_matrix[:,i])

        self.geom_means = np.power(np.prod(self.ratio_matrix,0),(1/self.num_trees))
        self.arith_means = np.mean(self.ratio_matrix,0)

    def plot(self):
        x = np.arange(self.num_clusters)

        fig, ax = plt.subplots(1,2)
        plt.setp(ax, xticks = x,xticklabels=self.cluster_names)
        
        ax[0].bar(x, self.geom_means,color = ["r"]+["b"]*(self.num_clusters-1))
        ax[0].set_title("Geometric Mean")
        ax[1].bar(x, self.arith_means,color = ["r"]+["b"]*(self.num_clusters-1))
        ax[1].set_title("Arithmetic Mean")
        fig.suptitle("Ratios in comparison to first cluster")
        plt.show()
    
    def writeCSV(self, filePath = "out.csv"):
        
        rows = [[self.tree_names[i]]+ self.makespan_matrix[i,:].tolist() + self.ratio_matrix[i,:].tolist() for i in range(self.num_trees) ]
        rows += [["geometric_mean"]+[""]*self.num_clusters+self.geom_means.tolist()]
        with open(filePath, 'a+', newline='') as csvfile:
            writer = csv.writer(csvfile, delimiter=',')
            writer.writerow([self.tree_category])
            writer.writerow(["tree"]+self.cluster_names+["ratio "+elem for elem in self.cluster_names])
            writer.writerow([])
            writer.writerows(rows)
            writer.writerow("")

def getAllCategories(results:List[OutputObject])->List[str]:
    categories = []
    for result in results:
        if not result.tree_category in categories:
            categories.append(result.tree_category)
    return categories

def getAllClusters(results:List[OutputObject])->List[str]:
    clusters = []
    for result in results:
        for clusterName in result.cluster_names[1:]:
            if clusterName not in clusters:
                clusters.append(clusterName)
    return clusters

#returns a tuple of the cluster that it comparing to and the diffference in mean
def findCategoryClusterMean(results:List[OutputObject], category:str, clusterName:str)->Tuple[str, float]:
    for result in results:
        if result.tree_category == category:
            if clusterName in result.cluster_names:
                idx = result.cluster_names.index(clusterName)
                return (result.cluster_names[0],result.geom_means[idx])
    return ("",0)

def plot(results:List[OutputObject])->None:
    categories = getAllCategories(results)
    #categories.sort()
    clusters = getAllClusters(results)
    #clusters.sort()

    fig, ax = plt.subplots(1,len(clusters))
    max_dist = 0
    min_dist = 0
    colourmap = plt.cm.get_cmap('hsv', len(clusters))
    colours = []
    for i in range(len(clusters)):
        colours.append(colourmap(i))
    for index,cluster in enumerate(clusters):
        print(cluster)
        x = np.arange(len(categories))
          
        #ax[index].tick_params(bottom=False)
        geomMeans = []
        for category in categories:
            (comparingCluster, mean ) = findCategoryClusterMean(results, category,cluster)
            geomMeans.append(mean)

        geomMeans = np.array(geomMeans)-1
        print(geomMeans)
        max_dist = max(max_dist,np.max(geomMeans))
        min_dist = min(min_dist,np.min(geomMeans))
        print(max_dist)
        ax[index].bar(x, geomMeans,color = colours)
        ax[index].set_title(f"{cluster} in comparison to \n {comparingCluster}")
        #ax[index].xaxis.set_major_locator(ticker.NullLocator())
        plt.setp(ax[index], xticks = x,xticklabels=categories)
        

    distance = 1.1*max(abs(max_dist),abs(min_dist))
    for i in range(len(clusters)):
        ax[i].set_ylim(ymin=-distance, ymax=distance)
    fig.suptitle("Geometric Means for different categories and clusters")
    
    plt.show()

def writeCSV(results:List[OutputObject], filePath = "out.csv"):
    with open(filePath, 'w') as f:
        f.write("")
    print(results)
    for result in results:
        result.writeCSV(filePath)
        

def parse(file) -> OutputObject:
    lines = file.readlines()
    output_information = OutputObject()
    for index, line in enumerate(lines):
        if index == 0:
            output_information.tree_category = line.split("/")[-2]
            continue
        line = line.split("\t")
        if index == 1:
            cluster_names = line[1:-1]
            output_information.cluster_names = [long_name.split("/")[-1] for long_name in cluster_names]
            output_information.init_makespan_matrix(len(lines)-2, len(line)-2)
            continue
        output_information.tree_names.append(line[0])
        for j in range(1,len(line)-1):
            output_information.makespan_matrix[index-2][j-1] = float(line[j])
    output_information.init_ratio_matrix()
    return output_information
if len(argv) >1:
    csv_file = argv[1]
else:
    csv_file = "out.csv"

'''
input_file = argv[1]
with open(input_file) as f:
    o = parse(f)
o.plot()
o.writeCSV(csv_file)
'''
cfg = simexpal.config_for_dir()
results :List[OutputObject] = []
for successful_run in cfg.collect_successful_results():
	with successful_run.open_output_file() as f:
		results.append(parse(f))


for result in results:
    print(result.tree_category)

writeCSV(results, csv_file)
plot(results)
'''
results = []
for successful_run in cfg.collect_successful_results(parse_fn=parse):
    print
    with successful_run.open_output_file() as f:
        parse(f) 
    '''