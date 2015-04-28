__author__ = 'Klaus Wehmuth'
# Support MAG Algorithms
#
# The algorithms present on this file are not explicitly shown in the paper "MultiAspect Graphs: Algebraic representation and algorithms".
# However, they perform useful support functions for MAG manipulation.
#
# In all input and output of algorithms presented on this file, tuples, lists and matrices are considered to start at index 1
# This is not the behavior of python, so it is necessary to add and subtract "1" at some points of the algorithms
# The algorithms were written this way to ease the visual verification of results in small examples
# For real applications, it is more convenient to use idexing starting from 0
#

import itertools as its
import scipy.sparse as sps
import numpy as np
import Algorithms as alg

def EdgesToMAG(edgeFile):
   """ Builds a MAG from its edge list
   :input string with the name of the file containing the MAG's edge list
   :output MAG H = (A,E) constructed from the list
   """
   f = open(edgeFile, "rb")
   nroAspects = -1
   l = f.readline()
   while len(l) > 0:
      i = l.find("(")
      j = l.find(")")
      if i >= 0 and j >= 0:
         ed = l[i+1:j].replace(" ","").split(",")
         if (len(ed) % 2) != 0:
            print "invalid edge"
            f.close()
            return None
         if nroAspects == -1:
            nroAspects = len(ed) / 2
            E = []
            As = []
            for n in range(nroAspects):
               As.append([])
         if len(ed) != (2 * nroAspects):
            print "Edge " + repr(ed) + " has invalid nro of aspects"
            f.close()
            return None
         E.append(tuple(ed))
         for n in range(2 * nroAspects):
            As[n % nroAspects].append(ed[n])
      l = f.readline()
      A = []
   f.close()
   for n in range(nroAspects):
      A.append(set(sorted(As[n])))
   return A,set(E)

def inducedSubMAG(H, Ai):
   """ Returns the subMAG of MAG H induced by the aspect sublist Ai
   :input MAG H and aspects sublist Ai
   :output The subMAG induced by Ai
   """
   A = H[0]					# Original MAG Aspect list
   E = H[1]					# Original MAG Edge set
   T = alg.CompTuple(A)
   n = 1
   for i in range(len(A)):			# n is the number of composite vertices on the original MAG
      n = n * len(A[i])
   cv = [0] * n					# list of n = |V(H)| integers
   for v in its.product(*Ai):			# every composite vertex from aspetcs in Ai
      cv[alg.D(idxV(v,A), T)-1] = 1
   Ei = set()
   for e in E:
      o = alg.D(idxV(alg.pi_o(e),A),T)-1	# origin cv numeric representation
      d = alg.D(idxV(alg.pi_d(e),A),T)-1	# destination cv numeric representation
      if cv[o] == 1 and cv[d] == 1:
         Ei.add(e)
   return Ai,Ei

def cvToAdjMatrixRow(v, A_H):
   """ Returns the matrix line corresponding to a composite vertex
   :input composite vertex v in its symbolic tuple form, and MAG's Aspect list A_H
   :output numeric representation of the composite vertex
   The numeric representation of a composite vertex directly indicates its corresponding row and column on the adjacency matrix
   """
   T_H = alg.CompTuple(A_H)
   return alg.D(idxV(v, A_H), T_H)

def adjMatrixRowToCv(n, A_H):
   """ Returns the composite vertex corresponding to a matrix line
   :input numeric representation of the composite vertex n, and MAG's Aspect list A_H
   :output composite vertex in its symbolic tuple form
   The numeric representation of a composite vertex directly indicates its corresponding row and column on the adjacency matrix
   This function is the inverse of cvToAdjMatrixRow
   """
   T_H = alg.CompTuple(A_H)
   return invIdxV(alg.InvD(n, T_H), A_H)

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

def invIdxV(v, A_H):
   """ Transforms numerical composite vertex into symbolic composite vertex
   :input composite vertex v in its numerical form and the MAG's Aspect list A_H
   :output composite vertex on its symbolic tuple form
   This function just take the symbolic name associated to each index and constructs the symbolic tuple
   It is the inverse of idxV
   """ 
   p = len(v)
   vi = []
   for i in range(p):
      n = v[i]
      s = sorted(list(A_H[i]))		# Aspect as a sorted list
      vi.append(s[n-1])
   return tuple(vi)

def edgeStringToTuple(es):
   es = es.replace('(', '').replace(')', '')
   es = es.replace(' ','')
   return tuple(es.split(','))

