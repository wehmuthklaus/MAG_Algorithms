__author__ = 'Klaus Wehmuth'
# MAG Algorithms
# 
# The algorithms present on this file were written to exemplify the algorithms showed in the paper "MultiAspect Graphs: Algebraic representation and algorithms"
# All algorithms are were written to best resemble the algorithms shown on the paper. 
# The performance of the algorithms was not taken into account.
#
# The algorithms presented on the paper and also in this file are in part based on well known graph algorithms, which can be found the literature, as for instance in "Introduction to Algorithms" by Cormen et. al, McGraw-Hill Higher Education 2dn edition, 2001
#
# In all input and output of algorithms presented on this file, tuples, lists and matrices are considered to start at index 1
# This is not the behavior of python, so it is necessary to add and subtract "1" at some points of the algorithms
# The algorithms were written this way to ease the visual verification of results in small examples
# For real applications, it is more convenient to use idexing starting from 0
#
# The usage examples consider the following:
# H is the MAG T example of the paper
# J is the adjacency matrix of H
# T is the Companion tuple of H

import scipy.sparse as sps
import numpy as np
import sys as sys
import SupportAlgorithms as spa

sys.setrecursionlimit(500000)	# increased recursion limit for DFS algorithms

def CompTuple(A_H):
   """ Returns the companion tuple of a MAG
   :input: The Aspects list of the MAG
   :output: Companion tuple of H
   Usage:
   CompTuple(H[0])
   """
   p = len(A_H)
   T = []
   for i in range(p):
      T.append(len(A_H[i]))
   return tuple(T)

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

def CompIdx(H, zeta=None):
   """ Returns the compact indices for a given MAG H.
       The compact form of a MAG's adjacency matrix is the matrix of the MAG without the trivial components,
       i.e. the adjacency matrix containig only vertices with degree > 0.
       The compact index toComp is used to convert the vertex numerical representation to compact form
       The compact index toFull is used to convert from compact form to the full representation form
       IF a sub-determination (zeta) is provided, the sub-determined indices are returned.
       :input the MAG H and optionally a sub-determination zeta
       :output the compact indices toFull and toComp
       Usage: CompIdx(H)   or   CompIdx(H, zeta)
   """ 
   A = aspectsToDictList(H[0])		# Aspects as a list of dictionaries 
   E = H[1]
   T = CompTuple(A)
   if zeta != None:
      T = SubCompTuple(T, zeta)
   vts = set()
   for e in E:
      en = idxDict(e,A)			# edge as numerical tuple
      u = D(pi_o(en), T)
      v = D(pi_d(en), T)
      vts.add(u-1)
      vts.add(v-1)
   toFull =  sorted(list(vts))
   toComp = dict()
   for i in range(len(vts)):
      toComp[toFull[i]] = i
   return toFull, toComp

def D(v, tau_H, toComp=None):		
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
         d = d + ((v[i]-1) * w) 	# v[i]-1 because python tuples start at index 0
         w = w * T[i]
   if toComp == None:
      return d+1			# d+1 because python tuples start at index 0
   return toComp[d]+1

def InvD(dn, tau_H, toFull=None):
   """ Returns numerical tuple of the composite vertex
   :input: d - numerical representation of v and companion tuple of the MAG
   :output: composite vertex as numerical tuple
   d is the integer corresponding to the composite vertex
   This function is the inverse of function D
   Usage: Considering d = 2
   InvD(2,T)	or	InvD(d,T)
   """
   if toFull == None:
      d = dn
   else:
      d = toFull[dn-1]+1
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
 
def RMatrix(H):
   A = aspectsToDictList(H[0])		# Aspects as a list of dictionaries 
   E = H[1]
   T = CompTuple(A)
   n = 1
   for v in T:
      n = n*v
   vts = set()
   for e in E:
      en = idxDict(e,A)			# edge as numerical tuple
      u = D(pi_o(en), T)
      v = D(pi_d(en), T)
      vts.add(u-1)
      vts.add(v-1)
   vts = sorted(list(vts))
   m = len(vts)
   R_H = sps.dok_matrix((n,m), dtype = np.int64)
   toComp = dict()
   for i in range(len(vts)):
      R_H[vts[i],i] = 1
      toComp[vts[i]] = i
   return R_H.tocsr()

def AdjMatrix(H):
   """ Returns the Compact Adjacency matrix and companion tuple of a MAG
   :input: MAG H = (A,E)
   :output: Adjacency matrix of H, Companion tuple of H, toCompV, toFull and toComp
   The Adjacency matrix is in CSR form -- see scipy.sparse
   Usage:
   AdjMatrix(H)
   """
   A = aspectsToDictList(H[0])		# Aspects as a list of dictionaries 
   E = H[1]				# Edge set
   n = 1
   tf, tc = CompIdx(H)
   n = len(tf)				# Number of composite nodes with edges
   T = CompTuple(A)
   J_H = sps.dok_matrix((n,n), dtype = np.int64)
   for e in E:
      en = idxDict(e,A)			# edge as numerical tuple
      u = D(pi_o(en), T, tc)
      v = D(pi_d(en), T, tc)
      J_H[u-1,v-1] = 1			# u-1 and v-1 because python matrix index starts at 0
   return J_H.tocsr(), T, tf, tc

def SubDetMatrix(T_H, zeta, tf, tzc):
   """ Returns the Sub-determination matrix of a MAG, given a sub-determination zeta
   :input: MAG H and sub-determination tuple zeta
   :output: Sub-determination matrix of H under zeta 
   The Sub-determination matrix is in CSR form -- see scipy.sparse
   Usage: Considering zeta = (1,1,0)
   SubDetMatrix(T, (1,1,0))
   SubDetMatrix(T,zeta)
   """
   Tz = SubCompTuple(T_H,zeta)
   n = len(tf)			# n = number of connected vertices
   m = len(tzc)			# m = number of connected sub-determined vertices
   Mz = sps.lil_matrix((m,n), dtype = np.int64)
   for j in range(n):
      u = InvD(j+1, T_H, tf)	# j+1 because python range starts at 0
      i = D(u, Tz, tzc)
      Mz[i-1,j] = 1		# i-1 because python matrix index starts as 0
   return Mz.tocsr()

def Degree(H):
   """ Returns the indegree and outdegree of all composite vertices on a MAG H
   :input: MAG H = (A,E)
   :output: indegree and outdegree of all composite vertices
   The indegree and outdegree are returned as two tuples
   Each tuple has |V(H)| entries
   Usage:
   Degree(H)
   """
   A = aspectsToDictList(H[0])		# Aspects as a list of dictionaries
   E = H[1]				# Edge set
   n = 1
   for i in range(len(A)):
      n = n * len(A[i])			# n = |V(H)|
   T = CompTuple(A)
   indegree = [0] * n			# indegree <- list of n zeros
   outdegree = [0] * n
   for e in E:
      en = idxDict(e,A)			# edge as numerical tuple
      o = D(pi_o(en),T) - 1		# -1 because python list starts with index 0
      d = D(pi_d(en),T) - 1		# -1 because python list starts with index 0
      indegree[d] = indegree[d] + 1
      outdegree[o] = outdegree[o] + 1
   return tuple(indegree), tuple(outdegree)

def SubDetDegree(H, zeta):
   """ Returns the indegree and outdegree of all sub-determined composite vertices on a MAG H under a given sub-determination zeta
   :input: MAG H = (A,E) and sub-determination tuple zeta
   :output: indegree and outdegree of all composite vertices
   The indegree and outdegree are returned as two tuples
   Each tuple has |Vz(H)| entries - the number of sub-determined composite vertices under zeta
   Usage: Considering zeta = (1,1,0)
   SubDetDegree(H, (1,1,0))
   SubDetDegree(H,zeta)
   """
   A = aspectsToDictList(H[0])		# Aspects as a list of dictionaries
   E = H[1]				# Edge set
   n = 1
   for i in range(len(A)):
      if zeta[i] != 0:
         n = n * len(A[i])		# n = |Vz(H)|
   T = CompTuple(A)
   Tz = SubCompTuple(T,zeta)
   indegree = [0] * n			# indegree <- list of n zeros
   outdegree = [0] * n
   for e in E:
      en = idxDict(e,A)			# edge as numerical tuple
      o = D(pi_o(en),Tz) - 1		# -1 because python list starts with index 0
      d = D(pi_d(en),Tz) - 1		# -1 because python list starts with index 0
      indegree[d] = indegree[d] + 1
      outdegree[o] = outdegree[o] + 1
   return tuple(indegree), tuple(outdegree)

def SubDetDegreeSepLoops(H, zeta):
   """ Returns the indegree and outdegree of all composite vertices on a MAG H under sub-determination zeta, with self-loops in separeta tuples
   :input: MAG H = (A,E) and sub-determination tuple zeta
   :output: indegree and outdegree and self loops of all composite vertices
   The indegree and outdegree are returned as two tuples
   Each tuple has |V(H)| entries
   Usage: Considering zeta = (1,1,0)
   SubDetDegreeSepLoops(H, (1,1,0))
   SubDetDegreeSepLoops(H,zeta)
   """
   A = aspectsToDictList(H[0])		# Aspects as a list of dictionaries
   E = H[1]				# Edge set
   n = 1
   for i in range(len(A)):
      if zeta[i] != 0:
         n = n * len(A[i])		# n = |Vz(H)|
   T = CompTuple(A)
   Tz = SubCompTuple(T,zeta)
   indegree = [0] * n			# indegree <- list o n zeros
   outdegree = [0] * n
   sdegree = [0] * n
   for e in E:
      en = idxDict(e,A)			# edge as numerical tuple
      o = D(pi_o(en),Tz) - 1		# -1 because python list starts with index 0
      d = D(pi_d(en),Tz) - 1		# -1 because python list starts with index 0
      if d != o:
         indegree[d] = indegree[d] + 1
         outdegree[o] = outdegree[o] + 1
      else:
         sdegree[d] = sdegree[d] + 1
   return tuple(indegree), tuple(outdegree), tuple(sdegree)

def BFS(J_H, T_H, s, toComp):
   """ Executes a sub-determined Breadth First Search on a given MAG starting from a given composite vertex, under a sub-determination
   Returns a tuple indicating the sub-determined composite vertices reached, a tuple with the distances from the origin vertex, and a tuple with the predecessor of each reached sub-determined composite vertex
   :input: Adjacency matrix J_H, Companion tuple T_H, sub-determination tuple zeta, and start composite vertex s (in numeric tuple format)
   :output: reached sub-determined vertices, distances, and predecessors in three distinct tuples
   Each tuple has |Vz(H)| entries - the number of sub-determined composite vertices under zeta
   The tuple of reached vertices has 1 if the corresponding composite vertex was reached by the BFS and 0 if it was not reached
   The tuple of distances has -1 for not reached composite vertices, and number of hops from s for reached composite vertices
   The tuple of predecessors has -1 for composite vertices which do not have a predecessor on the BFS tree, and the numeric representation of the predecessor if there is one
   Usage: Considering s = (2,1,1) and zeta = (1,1,0)
   BFS_Sub(J,T,(1,1,0),(2,1,1),tf,tc,tzc)
   BFS_Sub(J,T,zeta,s,tf,tc,tzc)
   tf and tc are generated by the AdjMatrix function
   tzc is generated by the CompIdx function
   """
   n = len(toComp)			# n = nro vertices with edges
   vertices = [0] * n			# vertices <- list of n zeros
   distance = [-1] * n			# -1 => infinity
   pred = [-1] * n			# -1 => Nil
   color = [0] * n			# all vertices are white			   
   frontier = []			# empty queue
   vertices[D(s,T_H, toComp)-1] = 1	# python list starts at index 0
   distance[D(s,T_H, toComp)-1] = 0	# python list starts at index 0
   frontier.append(D(s,T_H, toComp)-1)	# python list starts at index 0
   while len(frontier) > 0:
      u = frontier.pop(0)		# pop dequeues frontier
      for v in J_H.getrow(u).indices:	# out neighbors - J_H is a CSR matrix
         if color[v] == 0:		# if v is white
            color[v] = 1		# make v gray
            vertices[v] = 1
            distance[v] = distance[u]+1
            pred[v] = u+1		# +1 because python lists start at 0
            frontier.append(v)
      color[u] = 2			# make u black
   return tuple(vertices), tuple(distance), tuple(pred)

def BFS_Sub(J_H,T_H,zeta,s,tf,tc,tzc):
   """ Executes a sub-determined Breadth First Search on a given MAG starting from a given composite vertex, under a sub-determination
   Returns a tuple indicating the sub-determined composite vertices reached, a tuple with the distances from the origin vertex, and a tuple with the predecessor of each reached sub-determined composite vertex
   :input: Adjacency matrix J_H, Companion tuple T_H, sub-determination tuple zeta, and start composite vertex s (in numeric tuple format)
   :output: reached sub-determined vertices, distances, and predecessors in three distinct tuples
   Each tuple has |Vz(H)| entries - the number of sub-determined composite vertices under zeta
   The tuple of reached vertices has 1 if the corresponding composite vertex was reached by the BFS and 0 if it was not reached
   The tuple of distances has -1 for not reached composite vertices, and number of hops from s for reached composite vertices
   The tuple of predecessors has -1 for composite vertices which do not have a predecessor on the BFS tree, and the numeric representation of the predecessor if there is one
   Usage: Considering s = (2,1,1) and zeta = (1,1,0)
   BFS_Sub(J,T,(1,1,0),(2,1,1))
   BFS_Sub(J,T,zeta,s)
   """
   n = len(tf)					# n = number of composite vertices with degree > 0
   nS = len(tzc)				# nS = number of sub-determined composite vertices with degree > 0
   Tz = SubCompTuple(T_H,zeta)
   vertices = [0] * nS				# vertices <- list of nS zeros
   distance = [-1] * nS				# -1 => infinity
   pred = [-1] * nS				# -1 => Nil
   colorS = [0] * nS				# all sub-determined composite vertices are white
   color = [0] * n				# all composite vertices are white
   frontier = []				# empty queue
   for v in range(n):				# for every v in V(H)
      if J_H.getrow(v).nnz > 0:			# if v has out edges
         v = InvD(v, T_H, tf)			# v in numerico tuple form
         if D(v,Tz,tzc) == D(s,Tz,tzc):		# where D(v,Tz) = D(s,Tz)
            frontier.append(D(v,T_H,tc)-1)	# -1 because python list starts at 0
            color[D(v,T_H,tc)-1] = 1		# make v gray
   vertices[D(s,Tz,tzc)-1] = 1			# -1 because python list starts at 0
   distance[D(s,Tz,tzc)-1] = 0
   while len(frontier) > 0:
      u = frontier.pop(0)			# pop dequeues frontier
      for v in J_H.getrow(u).indices:		# out neighbors - J_H is a CSR matrix
         if color[v] == 0:			# if v is white
            color[v] = 1			# make v gray
            frontier.append(v)
            vn = InvD(v+1, T_H, tf)		# vn is v in numeric tuple form
            vS = D(vn, Tz, tzc)-1		# vS is v sub-determined by Tz
            if colorS[vS] == 0:			# if Vs is white
               colorS[vS] = 1			# make vS gray
               vertices[vS] = 1
               un = InvD(u+1, T_H, tf)		# un is u in tuple form
               uS = D(un, Tz, tzc)-1		# uS is u sub-determined by Tz
               if distance[vS] == -1:		# not reached yet
                  distance[vS] = distance[uS] + 1
                  pred[vS] = uS + 1
      color[u] = 2				# make u black
   return tuple(vertices), tuple(distance), tuple(pred)

def DFS(J_H, T_H, toFull):
   """ Executes a Depth first search on a given MAG
   Returns a tuple indicating the discovery time of each composite vertex, a tuple indicating the finish time of each composite vertex and a tuple indicating the predecessor of each composite vertex
   :input: Adjacency matrix J_H, and Companion tuple T_H
   :output: dicovery times, finish times, and predecessors in three distinct tuples
   Each tuple has |V(H)| entries
   Usage:
   DFS(J,T)
   """
   n = len(toFull)			# n = number of composite vertices with degree > 0
   color = [0] * n			# color <- list of n zeros => all composite vertices are white
   discTime = [-1] * n
   finTime = [-1] * n
   pred = [-1] * n
   time = [0]				# time in mutable form, for passing as reference
   for u in range(n):
      if color[u] == 0:
         DFS_Visit(u, time, color, discTime, finTime, pred, J_H)
   return tuple(discTime), tuple(finTime), tuple(pred)

def DFS_Visit(u, time, color, discTime, finTime, pred, J_H):
   """ Helper function to recursively determine each DFS tree 
       Should not be used directly
   """
   color[u] = 1				# make u gray
   discTime[u] = time[0]
   time[0] = time[0]+1
   for v in J_H.getrow(u).indices:	# out neighbors - J_H is a CSR matrix
      if color[v] == 0:
         pred[v] = u+1			# +1 because python list starts at 0
         DFS_Visit(v, time, color, discTime, finTime, pred, J_H)
   color[u] = 2				# make u black
   finTime[u] = time[0]
   time[0] = time[0]+1
   return

def DFS_Sub(J_H, T_H, zeta, tf, tzf, tc, tzc):
   """ Executes a sub-determined Depth first search on a given MAG
   Returns a tuple indicating the discovery time of each sub-determined composite vertex, a tuple indicating the finish time of each sub-determined composite vertex and a tuple indicating the predecessor of each sub-determined composite vertex
   :input: Adjacency matrix J_H, Companion tuple T_H, and sub-determination tuple zeta
   :output: dicovery times, finish times, and predecessors in three distinct tuples
   Each tuple has |Vz(H)| entries - the numaber of sub-determined composite vertices under zeta
   Usage: Considering zeta = (1,1,0)
   DFS_Sub(J,T,(1,1,0))
   DFS_Sub(J,T,zeta)
   """
   Tz =  SubCompTuple(T_H, zeta)
   Mz = SubDetMatrix(T_H, zeta, tf, tzc)
   Jz = Mz * J_H * Mz.transpose()
   n = len(tzf)				# n = number of sub-determined composite vertices with degree > 0
   color = [0] * n			# color <- list of n zeros => all sub-determined composite vertices are white
   discTime = [-1] * n
   finTime = [-1] * n
   pred = [-1] * n
   time = [0]				# time in mutable form, for passing as reference
   for u in range(n):
      if color[u] == 0:
         un = InvD(u+1, Tz, tzf)	# +1 beacause python list starts at 0
         vertices,d,p = BFS_Sub(J_H,T_H,zeta,un,tf,tc,tzc)
         DFS_Sub_Visit(u, vertices, time, color, discTime, finTime, pred, Jz)
   return tuple(discTime), tuple(finTime), tuple(pred)

def DFS_Sub_Visit(u, vertices, time, color, discTime, finTime, pred, Jz):
   """ Helper function to recursively determine each DFS tree 
       Should not be used directly
   """
   color[u] = 1				# make u gray
   discTime[u] = time[0]
   time[0] = time[0]+1
   for v in Jz.getrow(u).indices:	# out neighbors - Jz is a CSR matrix
      if color[v] == 0 and vertices[v] != 0:
         pred[v] = u+1
         DFS_Sub_Visit(v, vertices, time, color, discTime, finTime, pred, Jz)
   color[u] = 2				# make u black
   finTime[u] = time[0]
   time[0] = time[0]+1
   return

def strongConnComps(J_H, T_H, toFull):
   n = len(toFull)			# n = number of composite vertices with degree > 0
   color = [0] * n			# color <- list of n zeros 
   discTime = [-1] * n
   finTime = [-1] * n
   pred = [-1] * n
   time = [0]				# time in mutable form, for passing as reference
   for u in range(n):
      if color[u] == 0:
         DFS_Visit(u, time, color, discTime, finTime, pred, J_H)

   ft = []
   for i in range(n):
      ft.append((finTime[i],i))
   ft.sort(reverse=True)
   JT = J_H.transpose()

   color = [0] * n			# color <- list of n zeros 
   discTime = [-1] * n
   finTime = [-1] * n
   pred = [-1] * n
   time = [0]				# time in mutable form, for passing as reference
   for i in range(n):
      u = ft[i][1]
      if color[u] == 0:
         DFS_Visit(u, time, color, discTime, finTime, pred, JT)

   C = spa.getDFSTreesVts(discTime, finTime)
   return C

def boolTupleToCV(bt, A_H, T, toFull):
   cvs = []
   for i in range(len(bt)):
      if bt[i] != 0:
         v = toFull[i]+1
         cvs.append(toCV(v, A_H, T))
   return cvs

def vTupleToCV(vt, A_H, T, toFull):
   cvs = []
   for i in range(len(vt)):
      if vt[i] >= 0:
         v = toFull[vt[i]-1]+1
         cvs.append(toCV(v, A_H, T))
      else:
         cvs.append('-')
   return cvs

def toCV(v, A_H, T_H, toFull=None):
   return invIdxV(InvD(v, T_H, toFull), A_H, T_H)

def invIdxV(v, A_H, T_H):
   """ Transforms numerical composite vertex into symbolic composite vertex
   :input composite vertex v in its numerical form, the MAG's Aspect list A_H and companion tuple T_H
   :output composite vertex on its symbolic tuple form
   This function just take the symbolic name associated to each index and constructs the symbolic tuple
   It is the inverse of idxV
   """ 
   p = len(v)
   vi = []
   j = 0
   for i in range(p):
      while T_H[j] == 0:
         j = j+1
      n = v[i]
      s = sorted(list(A_H[j]))		# Aspect as a sorted list
      vi.append(s[n-1])
      j = j+1
   return tuple(vi)

def idx(e, A_H):
   """ Transforms an edge in symbolic form to its corresponding numerical form
   :input edge e in symbolic tuple form and the MAG's Aspect list A_H
   :output edge e in numeric tuple form
   This function just indexes each symbolic aspect element to its corresponding integer index
   Usage: Considering e = ('Loc3', 'Bus', 'T1', 'Loc2', 'Bus', 'T2') 
   idx(('Loc3', 'Bus', 'T1', 'Loc2', 'Bus', 'T2'), H[0])
   idx(e,H[0])
   """
   p = len(e) / 2
   ei = []
   for i in range(len(e)):
      a = e[i]			#Aspect element in symbolic form
      s = sorted(list(A_H[i%p]))#Aspect as a sorted list
      for j in range(len(A_H[i%p])):
         if a == s[j]:
            ei.append(j+1)
            break
   return tuple(ei)

def idxV(v, A_H):
   """ Transforms a composite vertex in symbolic form to its corresponding numerical form
   :input composite vertex v in symbolic tuple form and the MAG's Aspect list A_H
   :output composite vertex v in numeric tuple form
   This function just indexes each symbolic aspect element to its corresponding integer index
   """
   p = len(v)
   vi = []
   for i in range(p):
      a = v[i]
      s = sorted(list(A_H[i]))		# Aspect as a sorted list
      for j in range(len(A_H[i])):
         if a == s[j]:
            vi.append(j+1)
            break
   return tuple(vi)

def pi_o(e):
   """ Returns the origin composite vertex of a given edge e
   :input edge e
   :output origin composite vertex of e
   Usage: Considering e = ('Loc3', 'Bus', 'T1', 'Loc2', 'Bus', 'T2')
   pi_o(('Loc3', 'Bus', 'T1', 'Loc2', 'Bus', 'T2'))
   pi_o(e)
   """
   p = len(e) / 2
   return e[0:p]

def pi_d(e):
   """ Returns the destination composite vertex of a given edge e
   :input edge e
   :output destination composite vertex of e
   Usage: Considering e = ('Loc3', 'Bus', 'T1', 'Loc2', 'Bus', 'T2')
   pi_d(('Loc3', 'Bus', 'T1', 'Loc2', 'Bus', 'T2'))
   pi_d(e)
    """
   p = len(e) / 2
   return e[p: 2*p]

def aspectsToDictList(A_H):
   A = []
   for i in range(len(A_H)):		# for each aspect
      d = dict()			# new dictionary
      s = sorted(list(A_H[i]))		# Aspect as a sorted list
      for j in range(len(s)):		# for each aspect element
         d[s[j]] = j+1
      A.append(d)
   return A

def idxDict(e, D_A):
   p = len(e) / 2			# p <- nro of aspects
   ei = []
   for i in range(len(e)):		# for each edge element
      a = e[i]				# a <- edge element
      n = D_A[i%p][a]			# numeric index of element a
      ei.append(n)
   return tuple(ei)

