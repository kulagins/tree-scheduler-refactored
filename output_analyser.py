import simexpal
import numpy as np
from matplotlib import pyplot as plt
from sys import argv
import csv

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
        rows += [["geometric mean"]+[""]*self.num_clusters+self.geom_means.tolist()]
        with open(filePath, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile, delimiter=',')
            writer.writerow(["tree"]+self.cluster_names+["ratio "+elem for elem in self.cluster_names])
            writer.writerow([])
            writer.writerows(rows)


def parse(file):
    lines = file.readlines()
    output_information = OutputObject()
    for index, line in enumerate(lines):
        if index == 0:
            output_information.tree_category = line
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

input_file = argv[1]
if len(argv) >2:
    csv_file = argv[2]
else:
    csv_file = "out.csv"

with open(input_file) as f:
    o = parse(f)
o.plot()
o.writeCSV(csv_file)

"""
cfg = simexpal.config_for_dir()
results = []
for successful_run in cfg.collect_successful_results():
    
    with successful_run.open_output_file() as f:
        parse(f) """
