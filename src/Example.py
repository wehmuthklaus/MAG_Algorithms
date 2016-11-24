import SupportAlgorithms as sup
import Algorithms as alg

H = sup.EdgesToMAG("../MAG_Edge_Lists/MAG_EX1_Edges.txt") 	# Builds MAG H from an edge list
J,T = alg.AdjMatrix(H)						# Computes Adjacency matrix (J) and companion tuple (T) for MAG H
V,D,P = alg.BFS(J,T, (2,1,1))					# Runs BFS starting from composite vertex (2,1,1)

print(V)							# Prints reached composite vertices 
print(D)							# Prints distances
print(P)							# Prints predecessors
