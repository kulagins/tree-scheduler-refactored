# tree-scheduler

Modification of MemComJournal for heterogeneous clusters

## Getting Started
### How to compile
1. In the root directory, do `make all` to compile files for the include
3. Move to the `./examples/` directory and compile the files there via `make all`. After the succesful compilation, there will be an executable file `./res` in the root directory

### How to run
The executable requires the following input parameters:

-  `heterogeneity`:
   -  `heterogeneity = 0`: both, memory and computing power are homogeneous
   -  `heterogeneity = 1`: the memory is heterogeneous, but computing power is still homogeneous
   -  `heterogeneity = 2`: both, memory and computing power are heterogeneous
- `clustering-mode` : How many processors exist, how much memory does a processor have
  - `clustering-mode = 0`: the clustering-mode is dependent on the tree instance, we will, just like in the MemComJournal, give the CCR, NPR-values as an input
  - `clustering-mode = 1`: the clustering-mode is static, we will give the number of processors and their memory as input

The basic call is then 
```Shell
./res trees-directory trees-list heterogeneity clustering-mode cluster_arg_1 cluster_arg_2
```
Where `cluster_arg_1` and `cluster_arg_2` are the parameters that we give for the `clustering-mode`.

Plugging in a fixed value for the `clustering-mode` then yields:

```Shell
./res trees-directory trees-list heterogeneity 0 CCR NPR 
```
if we use the tree-dependent clustering, or:
```Shell
./res trees-directory trees-list heterogeneity 1 number_processors processor_memory
```
if we use the static clustering.


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
