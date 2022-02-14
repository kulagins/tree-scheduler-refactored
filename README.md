# tree-scheduler

Modification of MemComJournal for heterogeneous clusters

## Getting Started
### How to compile
In the root directory, do `make allmain`.
### Input files
#### Tree Description
Input files are in the directory `./real_trees`. `trees` contain the trees generated from sparse matrices, while `random_trees` contain various random trees generated automatically.

All trees have the same format. Each line is dedicated to a single unique task.
Each line contains `id parent_id node_weight makespan_weight edge_weight`

`id` is the number of the current task, `parent_id` is the id of the (unique) parent of the current task. 
`node_weight` is the amount of memory required to execute this task on a processor,
`edge_weight` is the amount of data needed to be transferred from parent to current task in order to start the execution on the current node.
`makespan_weight` represents the runtime length of the execution of the task.

All values are numbers without a measurement unit.

When you create your own trees to run the program with, please be aware of the following assumptions that are made on the input file:
1. Every task has to be declared _after_ its children in the input file. This particularly means that the root-node of the tree has to be declared at last.
2. The ID of a node has to be higher than the ID of its children.
3. The children of a node have to habe ascending IDs without any gaps.

#### Cluster Description
The `cluster-file` is a JSON-File that specifies the cluster that the programm should run at. The It specifies a list of groups of processors, each group consists of n equivalent processors running at the same speed and with the same memory.

```code
{
    "groups":
    [
        {
            "number": unsigned int
            "memory": double
            "speed": double
        }, ...
    ]
}
```
Depending on the `clustering-mode`, the `memory` and `speed` values are interpreted differently: with static clustering, the numbers are given as absolute values. If the clustering mode is dependent on the tree, then the speed value stays an absolute number and the memory value acts as a scalar changing the processors' memory relative to the maximum memory usage of a single task in the tree.

The `cluster-file-path` can be called in two diffrent variants: either this is the path to a single JSON-File which contains the description for a cluster, or it is the path to a txt-File which then contains the paths to multiple Clusters, seperated by newlines. The base-path of each of this cluster is the directory where the txt-File ist stored.

### How to run
The executable requires the following input parameters:

- `clustering-mode` : How many processors exist, how much memory does a processor have. For more information, please refer to the Cluster-Description, explained above.
  - `clustering-mode = 0`: the clustering-mode is dependent on the tree instance.
  - `clustering-mode = 1`: the clustering-mode is static.
- `build-small-clusters`: has yet to be implemented
- `verbosity`: Should all of the debugging output be printed when the program runs


The basic call is then 
```Shell
./main trees-directory trees-list cluster-file-path clustering-mode build-small-clusters verbosity
```

### How to run tests and add new ones
To compile all existing tests, do ```make alltests``` in root.

The following tests have been written already and can be executed in the `test/` directory:
- `cluster_unittest`
- `task_tree_unittest`

Add new tests in the directory `test`. Any tests added in the existing file `test.cpp` will be executed without further changes. If you add new test files, add them in the Makefile as new targets after line `100`.

### How to run experiments using simexpal
1. Install simexpal using `pip3 install simexpal`
2. run experiments using `simex experiments run --all`

You can delete the results using `simex experiments purge --all` and show the program output using `simex experiments print --all`.

#### How to add new experiments to simexpal

Depending on what new parameter you want to add, you have to change a different part of the `experiments.yml`:
|Parameter|Descritption|
|---|---|
|`treesList`| under the `experiments`-key, add a new entry with a useful name and change the path to the treesList.|
|`clusterList` |make sure that the new `clusterList` is located in `./clusters/` and add it's name to the `instances`-`items`-key|


## Frequently Asked Questions and Common Errors

A collection of common errors for compiling and running the project that are easy to fix.

#### `heuristics.a: No such file or directory` when calling `make all` in the root-directory

Check if the `lib`-Directory exists in your repository, so the `heuristics.a` file can be written into it.

#### In VS Code, when running the build task: "Cannot build and debug because the active file is not a C or C++ source file", or a similar notice that the program doesn't exist on x84-Architecture

Check if you are running the build task while having the `main.cpp` open.

#### Wie sind Bäume in der Datei geordnet und wie werden sie in der `Tree`-Klasse umgesetzt?
Der Root hat in der Datei den größten id. Er steht auch immer an der letzten Stelle, aber eigentlich sollen wir uns darauf nicht verlassen. Er soll den id 1 bekommen und an der Stelle 0 stehen.
Der Knoten mit dem zweitgrößten id ist sein erstes Kind. Er soll den id 2 bekommen und an der Stelle 1 stehen.
Weiter können Kinder von Root oder vom ersten Kind stehen, aber sie sind weiterhin absteigend durchnummeriert und sollen umgedreht dann die ids 3,4 etc bekommen.
Der allerletzte Blatt ist solcher, dessen id 1 im Dokument ist. Wir sollten wiederum uns eigentlich nicht darauf verlassen, dass er an der 1 Zeile steht.
Siehe auch Issue #18
