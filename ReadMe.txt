DESCRIPTION
------------
The files provided in this project contain the python source code
referent to the algorithms described in the paper "MultiAspect Graphs:
Algebraic representation and algorithms", available at xxx.

Folder “src” contains the source code
Folder “MAG Edge Lists” contains examples of edge lists that can be
used to generate Multi-Aspect Graphs.

In particular, the file “Algorithms.py” contains the source code
referent to the implementation of all algorithms shown in the paper.
File “SupportAlgorithms.py” contains supplementary algorithms that are
not presented at the paper, but are useful for the creation and
manipulation of Multi-Aspect Graphs (MAGs). It includes the function
“EdgesToMAG”, which receives the name of a file containing a list of
MAG edges and returns a MAG constructed from these edges.


REQUIREMENTS
------------
To run the algorithms presented, the following modules are required:
 * NumPy (http://www.numpy.org)
 * SciPy (http://www.scipy.org)