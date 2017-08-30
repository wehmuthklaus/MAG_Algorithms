import scipy.sparse as sps
import numpy as np

def Betweenness(J_H, T_H):
   n = 1
   for i in range(len(T_H)):
      if T_H[i] > 0:
         n = n * T_H[i]				# n = |V(H)|
   Cb = [0] * n					# Betweenness centrality list
   for s in range(n):				# for each s in V(H)
      if J_H.getrow(s).nnz > 0:			# if s has at least one neighbor - many composite vertices don't
         S = []					# S is an empty stack
         P = [[] for x in range(n)]		# P has an empty list for each w in V(H)
         sigma = [0] * n			# list for counting shortest paths
         sigma[s] = 1				# by definition, each vertex has 1 shortest path to itself
         d = [-1] * n				# list of distances
         d[s] = 0				# distance from a vertex to itself is 0
         Q = [s]				# BFS queue containing vertex s
         while len(Q) > 0:			# Do BFS while vertex queue is not empty
            v = Q.pop(0)			# v <- first element of Q
            S.append(v)				# Pushes v into stack S
            v_neighb = J_H.getrow(v).indices	# neighbors of v
            for w in v_neighb:			# for each neighbor w of v do
               if d[w] < 0:			# if w found for the first time
                  Q.append(w)			# enqueue w in Q
                  d[w] = d[v]+1			# register distance from s to w
               if d[w] == d[v]+1:		# if there is a shortest path to w via v
                  sigma[w] = sigma[w]+sigma[v]	# counts shortest paths
                  P[w].append(v)		# v is predecessor of w
         delta = [0] * n			# list for computing dependencies
         while len(S) > 0:			# for each vertex stacked in S
            w = S.pop()				# S returns vertices in order of non-increasing distance from s
            for v in P[w]:			# for each predecessor of w
               delta[v] = delta[v]+((float(sigma[v])/float(sigma[w])) * (1 + delta[w]))	# calc of dependencies
            if w != s:
               Cb[w] = Cb[w] + delta[w]		# Calculates betwwennes of w
   return Cb


def SubBetweenness(J_H, T_H, zeta):
   Mz = SubDetMatrix(T_H, zeta)
   n = 1
   nS = 1
   for i in range(len(T_H)):
      n = n * T_H[i]				# n = |V(H)|
      if zeta[i] != 0:
         nS = nS * T_H[i]			# nS = |Vz(H)|
   Tz = SubCompTuple(T_H,zeta)
   Cb = [0] * nS				# Betwwenness centrality list
   for s in range(nS):				# for each s in Vz(H)
      color = [0] * n				# Color of vertices
      colorS = [0] * nS				# Color of sub-determined vertices
      S = []					# S is an empty stack
      P = [[] for x in range(nS)]		# P has an empty list for each w in Vz(H)
      sigma = [0] * nS				# list for counting shortest paths
      sigma[s] = 1				# by definition, each vertex has 1 shortest path to itself
      d = [-1] * nS				# list of distances
      d[s] = 0					# distance from a vertex to itself is 0
      Q = []					# Sub-determined BFS queue
      for idx in Mz.getrow(s).indices:
         if J_H.getrow(idx).nnz > 0:
            color[idx] = 1			# vertex idx is black
            Q.append(idx)			# Vertices with sub-determination s and not zero degree
      while len(Q) > 0:				# Do sub-determined BFS while queue not empty
         v = Q.pop(0)				# v <- first element of Q
         vn = InvD(v+1, T_H)			# vn is v in numeric tuple form
         vS = D(vn, Tz)-1			# vS is v sub-determined by Tz
         if colorS[vS] == 0:			# vS is white
            colorS[vS] = 1			# Turn vS black
            S.append(vS)			# pushes vS into stack S
         v_neighb = J_H.getrow(v).indices	# neighbors of v
         for w in v_neighb:			# for each neighbor w of v do
            if color[w] == 0:			# if w found for first time
               color[w] = 1			# turn w black
               Q.append(w)			# insert w in BFS queue
            wn = InvD(w+1, T_H)			# wn is w in numeric tuple form
            wS = D(wn, Tz)-1			# wS is w sub-determined by Tz
            if d[wS] < 0:			# if wS is found for the first time
               d[wS] = d[vS]+1			# determine distance from s to wS
            if d[wS] == d[vS]+1:		# if there is a shortest path to wS via vS
               sigma[wS] = sigma[wS]+sigma[vS]	# counts shortest paths
               P[wS].append(vS)			# vS is predecessor of wS
      delta = [0] * nS				# list for computing dependencies
      while len(S) > 0:				# for each vertex stacked in S
         w = S.pop()				# S returns vertices in order of non-increasing distance from s
         for v in P[w]:				# for each predecessor of w
            delta[v] = delta[v]+((float(sigma[v])/float(sigma[w])) * (1 + delta[w]))	# calc of dependencies
         if w != s:
            Cb[w] = Cb[w] + delta[w]		# Calculates betwwennes of w
   return Cb

def SubDetMatrix(T_H, zeta):
   """ Returns the Sub-determination matrix of a MAG, given a sub-determination zeta
   :input: Companion tuple T_H of MAG H and sub-determination tuple zeta
   :output: Sub-determination matrix of H under zeta 
   The Sub-determination matrix is in CSR form -- see scipy.sparse
   Usage: Considering zeta = (1,1,0)
   SubDetMatrix(T, (1,1,0))
   SubDetMatrix(T,zeta)
   """
   Tz = SubCompTuple(T_H,zeta)
   n = 1
   m = 1
   for i in range(len(T_H)):
      n = n * T_H[i]			# n = |V(H)|
      if zeta[i] != 0:
         m = m * T_H[i]			# m = |Vz(H)|
   Mz = sps.dok_matrix((m,n), dtype = np.int64)
   for j in range(n):
      u = InvD(j+1, T_H)		# j+1 because python range starts at 0
      i = D(u, Tz)
      Mz[i-1,j] = 1			# i-1 because python matrix index starts as 0
   return Mz.tocsr()

def SubCompTuple(tau_H, zeta):
   """ Returns the sub-determined companion tuple of a given companion tuple
   :input: companion tuple tau_H and sub-determination tuple zeta
   :output: Sub-determined companion tuple Tz
   Usage: Considering zeta = (1,1,0)
   SubCompTuple(T, (1,1,0))	or	SubCompTuple(T, zeta)
   """
   p = len(tau_H)
   Tz = []
   for i in range(p):
      Tz.append(tau_H[i] * zeta[i])
   return tuple(Tz) 

def D(v, tau_H):		
   """ Returns the numerical representation of a given composite vertex v
   :input: composite vertex v and MAG's companion tuple tau_H
   :output: Numerical representation of v
   The composite vertex v has to be in the form of a numerical tuple, which can be obtained from idx or idxV functions
   The numerical representation returned is an integer number ranging from 1 to |V(H)|
   This number corresponds to the row or column of v in the MAG's adjacency matrix - considering that the first element is numbered 1
   The idxV function can be found in SupportAlgorithms.py
   Usage: Considering v = (2,1,1)
   D((2,1,1), T)	or	D(v,T)
   """
   T = tau_H
   if len(tau_H) > len(v):
      T = [x for x in tau_H if x > 0]	# makes len(T) == len(v)
   p = len(T)
   d = 0
   w = 1
   for i in range(p):
      if T[i] != 0:
         d = d + ((v[i]-1) * w) # v[i]-1 because python tuples start at index 0
         w = w * T[i]
   return d+1			# d+1 because python tuples start at index 0

def InvD(d, tau_H):
   """ Returns numerical tuple of the composite vertex
   :input: d - numerical representation of v and companion tuple of the MAG
   :output: composite vertex as numerical tuple
   d is the integer corresponding to the composite vertex
   This function is the inverse of function D
   Usage: Considering d = 2
   InvD(2,T)	or	InvD(d,T)
   """
   p = len(tau_H)
   w = [1]
   j = 0
   for i in range(p):
      if tau_H[i] != 0:
         w.append(w[j] * tau_H[i])
         j = j+1
   v = []
   p = len(w)-1
   for i in range(p):
      v.append(((d-1) % w[i+1])/w[i] + 1) #d-1 and +1 because python tuples start at index 0
   return tuple(v)

