# tree-scheduler

Modification of MemComJournal for heterogeneous clusters


## Frequently Asked Questions and Common Errors

A collection of common errors for compiling and running the project that are easy to fix.

#### `heuristics.a: No such file or directory` when calling `make all` in the root-directory

Check if the `lib`-Directory exists in your repository, so the `heuristics.a` file can be written into it.

#### In VS Code, when running the build task: "Cannot build and debug because the active file is not a C or C++ source file", or a similar notice that the program doesn't exist on x84-Architecture

Check if you are running the build task while having the `main.cpp` open.
