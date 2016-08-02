"""" 
BipartiteGraph Utilities File specifically for LDPC-Channel coding

Contains Bit Node, Check Node, Graph and Belief Propogation classes

"""""

import numpy
import math

def extrSum(vec):
	sum=0
	op=range(len(vec))
	for i in vec:
		sum+=i
	for i in range(0,len(vec)):
		op[i]=sum-vec[i]
	return op

def extrProd(vec):
	op=[None]*len(vec)
	fprod=op
	bprod=op
	fprod[0]=1
	bprod[len(vec)-1]=1
	
	for i in range(1,len(vec)):
		fprod[i]=fprod[i-1]*vec[i-1]
	for i in range(len(vec)-1,0,-1):
		bprod[i-1]=bprod[i]*vec[i]
	for i in range(0,len(vec)):
		op[i]=fprod[i]*bprod[i]
	return op

def SoftThresh(vec,thresh):
	for i in range(0,len(vec)):
		if vec[i]>thresh:
			vec[i]=thresh
		elif vec[i]<-thresh:
			vec[i]=-thresh
		
	return vec


class Edge:
	"""Bit Node to Check Node, Edge class"""

	def __init__(self, n1, n2,name):
		"""Constructor. Takes bit and check node IDs as arguments"""
		self.ID=name
		self.node1=n1 #Variable Node ID
		self.node2=n2 #Check Node ID
		self.msg1=0
		self.msg2=0

	def getNodes(self):
		"""Returns a list containing the bit and check nodes for this edge"""
		return [self.node1, self.node2]

	def hasNodes(self, n1, n2):
		"""Takes two node IDs. Returns true if the IDs match the two nodes of this edge in that order."""
		if(self.node1==n1 and self.node2==n2):
			return True
			
		return False

	def getID(self):
		return self.ID


class BitNode:
	""" Variable node class."""

	def __init__(self, name):
		"""Constructor. Takes a node ID"""
		self.ID=name
		self.neighbors=[]
		self.edges=[]
		self.degree= 0
		self.chan=0

	def addNeighbor(self, nb,ed):
		"""Adds a list of neighbors to the current list. Takes a list of node IDs"""
		if((nb not in self.neighbors)):# & (ed not in self.edges)):
			self.neighbors.append(nb)
			self.edges.append(ed)
			self.degree+=1

	
	def addNeighbors(self, nbs,eds):
		"""Adds a list of neighbors to the current list. Takes lists of node and edge IDs"""
		if (type(nbs) is list and type(eds) is list):
			assert len(nbs)==len(eds)
		for i in range(0, len(nbs)):
			if(nbs[i] not in self.neighbors & eds[i] not in self.edges):
				self.neighbors.append(nbs[i])
				self.edges.append(eds[i])
				self.degree+=1

	def replaceNeighbors(self, nbs,eds):
		"""Replaces the current list of neighbors. Takes a list of node IDs"""
		if (type(nbs) is list and type(eds) is list):
			assert size(nbs)==size(eds)
		self.neighbors=nbs
		self.edges=eds
		self.degree=len(nbs)

	def removeNeighbor(self, nb,ed):
		"""Removes the specified node ID from the list of neighbors, if it exists"""
		if(nb in self.neighbors):
			self.neighbors.remove(nb)
			self.edges.remove(ed)
		assert(len(self.neighbors)==len(self.edges))
		self.degree=len(self.neighbors)

	def getID(self):
		"""Returns node ID"""
		return self.ID

	def getNeighbors(self):
		"""Returns list of neighbor IDs"""
		return self.neighbors
	
	def getEdges(self):
		"""Returns list of edge IDs"""
		return self.edges


class ChkNode:
	""" Check node class."""

	def __init__(self, name):
		"""Constructor. Takes a node ID"""
		self.ID=name
		self.neighbors=[]
		self.edges=[]
		self.degree= 0
		

	def addNeighbor(self, nb,ed):
		"""Adds a list of neighbors to the current list. Takes a list of node IDs"""
		if((nb not in self.neighbors)):# & (ed not in self.edges)):
			self.neighbors.append(nb)
			self.edges.append(ed)
			self.degree+=1

	
	def addNeighbors(self, nbs,eds):
		"""Adds a list of neighbors to the current list. Takes a list of node IDs"""
		if (type(nbs) is list and type(eds) is list):
			assert len(nbs)==len(eds)
		for i in range(0, len(nbs)):
			if((not nbs[i] in self.neighbors) & (eds[i] not in self.edges)):
				self.neighbors.append(nbs[i])
				self.edges.append(eds[i])
				self.degree+=1

	def replaceNeighbors(self, nbs,eds):
		"""Replaces the current list of neighbors. Takes a list of node IDs"""
		if (type(nbs) is list and type(eds) is list):
			assert len(nbs)==len(eds)	
		
		self.neighbors=nbs
		self.edges=eds
		self.degree=len(nbs)

	def removeNeighbor(self, nb,ed):
		"""Removes the specified node ID from the list of neighbors, if it exists"""
		if(nb in self.neighbors):
			self.neighbors.remove(nb)
			self.edges.remove(ed)
		assert(len(self.neighbors)==len(self.edges))
		self.degree=len(self.neighbors)
			

	def getID(self):
		"""Returns node ID"""
		return self.ID

	def getNeighbors(self):
		"""Returns list of neighbor IDs"""
		return self.neighbors
		
	def getEdges(self):
		"""Returns list of edge IDs"""
		return self.edges



class BipartGraph(object):
	"""Bipartite graph class.
	"""

	def __init__(self, eds):
		"""Constructor. Takes a list of edge primitives, which is a list of two node IDs.
		Iterates through the edges, creates nodes for unique node IDs, and adds all edges and nodes.
		"""
		self.size = 0
		self.BitNodes = []
		self.ChkNodes = []
		self.edges = []
		for i in range(0, len(eds)):
			self.addEdge(eds[i])
			self.size=len(self.edges)
		
	def addEdge(self, edgep):
		"""Adds an edge into the graph, and updates neighbors & degrees of relevant nodes.
		Takes an edge primitive, a list of two node IDs
		"""
		if(not self.containsEdge(edgep)):
			no1 = self.getBitNode(edgep[0])
			no2 = self.getChkNode(edgep[1])
			newEdge =self.getEdge(edgep[0], edgep[1])
			no1.addNeighbor(no2.getID(), newEdge)
			no2.addNeighbor(no1.getID(), newEdge)
		
	
	def removeEdge(self, edgep):
		"""Removes an edge from the graph, and updates the neighbors of relevant nodes
		Takes an edge primitive.
		"""

		for i in range(0, len(self.edges)):
			if(self.edges[i].hasNodes(edgep[0], edgep[1])):
				edID=self.edges[i].getID()
				self.edges.pop(i)
				break
		no1 = self.getBitNode(edgep[0])
		no2 = self.getChkNode(edgep[1])
		no1.removeNeighbor(edgep[1],edID)
		no2.removeNeighbor(edgep[0],edID)
		
	def getBitNode(self, name):
		"""Checks if a given Bit Node ID exists in the graph. If not, it creates and adds a Bit Node for the given ID. Returns the Bit Node"""
		for i in range(0, len(self.BitNodes)):
			if(self.BitNodes[i].getID()==name):
				return self.BitNodes[i]
		newNode = BitNode(name)
		self.BitNodes.append(newNode)
		return self.BitNodes[len(self.BitNodes)-1]

	def getChkNode(self, name):
		"""Checks if a given Chk Node ID exists in the graph. If not, it creates and adds a Chk Node for the given ID. Returns the Chk Node"""
		for i in range(0, len(self.ChkNodes)):
			if(self.ChkNodes[i].getID()==name):
				return self.ChkNodes[i]
		newNode = ChkNode(name)
		self.ChkNodes.append(newNode)
		return self.ChkNodes[len(self.ChkNodes)-1]

	def getEdge(self,n1,n2):	
		"""Checks if a given Edge ID exists in the graph. If not, it creates and adds a Bit Node for the given ID. Returns the Bit Node"""
		newEdge = Edge(n1,n2,len(self.edges))
		self.edges.append(newEdge)
		return self.edges[len(self.edges)-1]

	

	def containsEdge(self, edgep):
		"""Checks for an edge in the graph. Takes an edge primitive, which is a list of two node IDs. First ID is bit node, second ID is of Check node"""
		for e in self.edges:
			if(e.hasNodes(edgep[0], edgep[1])):
				return True

		return False


	
	def printGraph(self):
		"""Prints graphs nodes, neighbors, and edges. Should use ___str___()"""

		print("\n\t ---EDGES---")
		print("\n Bit Node \t Check Node \t Edge \n")
		for i in range(0, len(self.edges)):
			currentNodes = self.edges[i].getNodes()
			print(str(currentNodes[0])+" \t \t   "+str(currentNodes[1])+" \t\t  "+str(self.edges[i].getID()))

	def __str__(self):
		"""Prints graphs nodes, neighbors, and edges"""
		fullstring = "\n\t Bit Nodes: "
		for i in self.BitNodes:
		    fullstring +=str(i.getId()) + "  "
		    
		fullstring+= "\n\t Check Nodes: "
		for i in self.ChkNodes:
		    fullstring +=str(i.getId()) + "  "

		print("\n\t ---EDGES---")
		print("\n Bit Nodes  Check Nodes \n")
		for i in range(0, len(self.edges)):
			currentNodes = self.edges[i].getNodes()
			print(str(currentNodes[0])+" ---- "+str(currentNodes[1]))
		
		fullstring+="\n\tEDGES---EDGES---\n Bit Nodes  Check Nodes \n"
		for i in range(0, len(self.edges)):
			currentNodes = self.edges[i].getNodes()
			fullstring += str(currentNodes[0])+" ---- "+str(currentNodes[1])+"\n"

		return fullstring


class BP(BipartGraph):
	"""Belief Propogation decoder on a Bipartite graph class.
	"""
	def __init__(self,eds,ch):
		"""Constructor. Takes a graph and channel output vector as arguments"""
		BipartGraph.__init__(self,eds)
		for i in range(0,len(self.BitNodes)):
		    self.BitNodes[i].chan=ch[i]
		for e in self.edges:
			e.msg1=0
			e.msg2=0
		

	def BitOp(self):
		"""Computes the BP operation at a Bit Node. Takes the Edge ID as argument"""
		THRESH=20
		for i in self.BitNodes:
			msgvec=[]
			for j in i.edges:
				msgvec.append(j.msg2)
			msgvec.append(i.chan)
			extr=extrSum(msgvec)
			extr=SoftThresh(extr,THRESH)
			for j in range(0,len(i.edges)):
				i.edges[j].msg1=extr[j]
			
		
	def ChkOp(self):
		for i in self.ChkNodes:
			msgvec=[]
			for j in i.edges:
				msgvec.append(math.tanh(0.5*j.msg1))
			extr=extrProd(msgvec)
			for j in range(0,len(i.edges)):
				i.edges[j].msg2=2.0*math.atanh(extr[j])
			
	def BPIterations(self,MaxIter):
		for i in range(0,MaxIter):
			self.BitOp()
			self.ChkOp()
			
	def HardDec(self):
		op=[]
		for i in self.BitNodes:
			sum=i.chan
			for j in i.edges:
				sum+=j.msg2
				
			op.append(numpy.sign(sum))
		return op
		
		
		

		
	

		
	