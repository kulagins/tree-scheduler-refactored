# tree-scheduler

Modification of MemComJournal for heterogeneous clusters


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
