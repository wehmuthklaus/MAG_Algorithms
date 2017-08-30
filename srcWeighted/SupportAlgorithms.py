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
   W = dict()
   Eid = dict()
   while len(l) > 0:
      i = l.find("(")	# get edges
      j = l.find(")")
      if i >= 0 and j >= 0:
         ed = l[i+1:j].replace(" ","").split(",")
         if (len(ed) % 2) != 0:
            print "invalid edge " + repr(ed)
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

         i = l.find("%") # get edge ids
         j = i + l[i+1: -1].find("%") + 1
         Eid[tuple(ed)] = l[i+1: j]

         h = l.find("< ") # get edge weights
         k = l.find(" >")
         wg = 1
         if h >= 0 and k > h:
            wg = int(l[h+1:k])
         W[tuple(ed)] = wg
         for n in range(2 * nroAspects):
            As[n % nroAspects].append(ed[n])
      l = f.readline()
      A = []
   f.close()
   for n in range(nroAspects):
      A.append(set(sorted(As[n])))
   return (A,set(E)), W, Eid

def subDeterminedMAG(H, zeta):
   """Returns a sub-determined MAG, i.e. a MAG with 
   aggregated aspects according to the sub-determination 
   tuple zeta. The tuple zeta has the same number of
   elements as the number of aspects on the original MAG H, 
   and for each aspect carries a value 1 or 0, so that only 
   the aspects marked with 1 will remain on the resulted MAG. 
   The resulting sub-determined MAG has the edges of the original 
   MAG projected over the reduced aspect structure given by the 
   sub-determination.
   input: MAG H and sub-determination tuple zeta
   output: Sub-determined MAG Hz
   """ 
   E = H[1]
   Ez = set()
   asps = list(zeta) + list(zeta)
   n = len(asps)
   for e in E:
      ez = []
      for i in range(n):
         if asps[i] != 0:
            ez.append(e[i])
      if alg.pi_o(ez) != alg.pi_d(ez): 
         Ez.add(tuple(ez))
   A = buildAspectListFromEdges(Ez)
   return A,Ez

def inducedSubMAG(H,W,Eid, Ai, reduceAspects=True):
   """ Returns the subMAG of MAG H induced by the aspect sublist Ai
   :input MAG H and aspects sublist Ai
   :output The subMAG induced by Ai
   """
   E = H[1]
   Ei = set()
   Wi = dict()
   Ein = dict()
   for e in E:
      if checkAspects(e,Ai):
         print e
         Ei.add(e)
         Wi[e] = W[e]
         Ein[e] = Eid[e]
   if not reduceAspects:
      return (Ai,Ei),Wi, Ein
   Ar = buildAspectListFromEdges(Ei)
   return (Ar,Ei),Wi, Ein

def checkAspects(e,A):
   le = len(e)
   la = len(A)
   if le != 2*la:				# Wrong sizes
      return False
   for i in range(le):
      if e[i] not in A[i%la]:
         return False				# Edge Aspect not in A
   return True

def buildAspectListFromEdges(E):
   Ar = []
   n = len(list(E)[0])/2			# number of aspects in Ai
   for i in range(n):				# build induced aspects
      A = {e[i] for e in E}
      A = A.union({e[i+n] for e in E})
      Ar.append(A)
   return Ar

def showAdjacency(AdjList, v, AH):
   vn = cvToAdjMatrixRow(v, AH) - 1
   adj = AdjList[vn]
   al = []
   for u in adj:
      al.append(adjMatrixRowToCv(u+1, AH))
   return al

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

def primeiroCV(edg, aero, cam): #ok
   l = sorted([int(e[2]) for e in edg if e[0] == aero and e[1] == cam and e[4] == cam])
   if len(l) == 0:
      return None
   tm = "%05d" % l[0]
   return (aero, cam, tm)


