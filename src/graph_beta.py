#!/usr/bin/python
# -*- coding: utf-8 -*-

from enum import Enum, IntEnum, unique
import re
from os.path import dirname, abspath
#from utils import *

'''
Define a enumous class node type
try:
    @unique
    class NodeType(Enum):
        metabolite = 0
        reaction = 1
except ValueError as e:
    print(e)
'''

class Node(object):
    """
    A class of nodes in the metabolic network.
    Specifically, every metabolite and reaction are considered as nodes.
    
    Arguments & Attributes:
        node_type: 0-metabolite, 1-reaction
        name:
        id:
        reverse:
    """
    def __init__(self, name, id, node_type, reverse):
        """
        Initialize a single node.
        """
        if str(id)[0] in 'mM0':
            self.type = '0'
        elif str(id)[0] in 'rR1':
            self.type = '1'
        else:
            raise ValueError()
        if isinstance(name, str) and isinstance(id, str): # and isinstance(reverse, str):
            self.name = name
            self.id = id
            self.reverse = reverse
        else:
            raise TypeError()

    def __str__(self):
        """
        This function returns a string representation for the node.
        node:[vertex, vertex]
        """
        return str(self.name)+','+str(self.id)+','+str(self.type)+','+str(self.reverse)


class SubGraph(object):
    """
    A class of graph.

    Attributes:
        net_name :
        node_list: A list of valid nodes.
                   If no list is given, the graph is instantiated empty.
        graph    : A dictionary with a name of the node as a key, and a list
                   of nodes which are connected to the key node as the value.

    Methods:
        addNode: 
        addEdge:
        addRxn :
    """
    def __init__(self, network_name):
        """
        Initialize function.
        """
        self.net_name = network_name
        self.node_list = []
        self.graph = {}
    
    def addNode(self, node):
        """
        This function adds a node to the graph, Probabily not use.
        """
        if isinstance(node, Node):
            self.node_list.append(node)
        else:
            raise TypeError()

    def addRxn(self, itermlist):
        i = 1
        while itermlist[i] != '-':
            self.addEdge(itermlist[i], itermlist[0])
            i += 1
        i += 1
        while i < len(itermlist):
            self.addEdge(itermlist[0], itermlist[i])
            i += 1

    def addEdge(self, node1, node2):
        if node1 in self.graph:
            self.graph[node1].append(node2)
        else:
            self.graph[node1] = [node2]

    def isEdge(self, node1, node2):
        if isinstance(node1, Node) and isinstance(node2, Node):
            if node2 in self.graph[node1]:
                return True
            else:
                return False
        else:
            raise TypeError()
    def length(self):
        return len(self.node_list)

class Graph(object):
    """
    A class of graph, which was a list of subgraphs.

    Attributes:
        subGraphs: A list of SubGraph object

    Method:
        read : Read a .met file, which is calling when initializing.
        num_of_subGraphs : Returns the number of subgraphs.
    """
    def __init__(self, filename):
        path = dirname(dirname(abspath(__file__))) + '/data/ALLSPECIES_tsv/' + filename
        self.filename = filename
        self.subGraphs = self.__read(path)

    def __read(self,filename,separator = '\t'):
        """
        Read a network file, and return a list of subgraphs
        
        list_of_subgraph = []
        with open(filename, 'r') as file:
            
            ### Handling the first subgraph
            while line = file.readline():
                line = line.strip() 
                
                ###Judge if the line is empty or annotation
                
                if not len(line) or line.startswith('//'):
                    network_name = line[2:].strip()
                    break
            # network_name = file.readline().strip()[2:]
            sub = SubGraph(network_name)
        """
        list_of_subgraph = []
        with open(filename, 'r') as file:
            network_name = file.readline().strip()[2:]
            sub = SubGraph(network_name)
            for line in file:
                if line.strip().startswith('##'):
                    """
                    Handling a new subgraph
                    """
                    list_of_subgraph.append(sub)
                    sub = SubGraph(line.strip()[2:])
                elif line.strip().startswith('#'):
                    """
                    Handling a reaction entry
                    """
                    itemrlist = line[1:].strip().split(separator)
                    sub.addRxn(itemrlist)
                else:
                    """
                    Handling a node entry
                    """
                    itermlist = line.strip().split(separator)
                    if len(itermlist) is 3:
                        node_id,name, reverse = itermlist
                        node_type = 1
                    elif len(itermlist) is 2:
                        node_id, name = itermlist
                        node_type,reverse = 1,0
                    node = Node(name, node_id, node_type, reverse)
                    sub.addNode(node)
            list_of_subgraph.append(sub)
        return list_of_subgraph

    def num_of_subGraphs(self):
        return len(self.subGraphs)
    def __str__(self):
        s = ''
        for i in self.subGraphs:
            s1 = ''
            for j in i.node_list:
                s1 += str(j)+'\n'
            s += s1
        return s
    '''
    def graph2dot(self):
        graph_list = []
        for sub_graph in self.subGraphs:
            sub = []
            sub.append(sub_graph.net_name)
            node_list = []
            for i in sub_graph.node_list:
                s = ''
                s = s + i.name + '\t' + i.id + '\t'
                if i.node_type is NodeType(0):
                    s += '0'
                else:
                    s += '1'
                node_list.append(s)
            sub.append(node_list)
    '''


class Align_node(object):
    """
    A class of Node in the AlignGraph

    Attributes:
        node1: node in the first graph
        node2: node in the second graph
        weight: the node similarity
    """
    def __init__(self, node1, node2, weight):
        if isinstance(node1, Node) and isinstance(node2, Node):
            self.node1 = node1
            self.node2 = node2
            self.weight = weight
        else:
            raise TypeError()

Weight = [1,0.5,0.2]
Threshod = 0.25

class Align_sub_graph(object):
    """
    A class of Align subgraph, which is the element of aligngraph.

    Attributes:
        node_list: a list of Align_Node
        align_graph: a dictionary represents the adjacency list of the graph
        graphs: the original graph

    Method:
        __add_node:
        __add_edge:
    """
    def __init__(self, sub_graph1, sub_graph2, sim):
        self.node_list = []
        self.align_graph = {}
        self.sub_graph1 = sub_graph1
        self.sub_graph2 = sub_graph2
        self.align_sub_graph1 = sub_graph1
        self.align_sub_graph1.node_list = []
        self.align_sub_graph1.graph = {}
        for node1 in sub_graph1.node_list:
            for node2 in sub_graph2.node_list:
                self.__add_node(node1, node2, sim)
        self.__add_edge()

    def __add_node(self, node1, node2, sim):
        if sim(node1, node2) >= Threshod:
            node = Align_node(node1, node2, sim(node1, node2))
            if node not in self.node_list:
                self.node_list.append(node)
                self.align_sub_graph1.node_list.append(node1)
            else:
                pass
        else:
            pass

    def __add_edge(self):
        if len(self.node_list) <= 1:
            return
        else:
            for i in xrange(len(self.node_list)-1):
                for j in xrange(i+1,len(self.node_list)):
                    if self.sub_graph1.isEdge(self.node_list[i].node1,self.node_list[j].node1) and self.sub_graph2.isEdge(self.node_list[i].node2,self.node_list[j].node2):
                        self.align_sub_graph1.addEdge(self.node_list[i].node1,self.node_list[j].node1)
                        if self.node_list[i] in self.align_graph.keys():
                            self.align_graph[self.node_list[i]].append(self.node_list[j])
                            self.align_graph[self.node_list[j]].append(self.node_list[i])
                        else:
                            self.align_graph[self.node_list[i]] = []
                            self.align_graph[self.node_list[i]].append(self.node_list[j])
                            self.align_graph[self.node_list[j]].append(self.node_list[i])
                            '''
                    elif self.sub_graph1.isEdge(self.node_list[i][0],self.node_list[j][0])\
                    or self.sub_graph2.isEdge(self.node_list[i][1],self.node_list[j][1]):
                        if self.align_graph[self.node_list[i]]:
                            self.align_graph[self.node_list[i]][self.node_list[j]] = Weight[1]
                            self.align_graph[self.node_list[j]][self.node_list[i]] = Weight[1]
                        else:
                            self.align_graph[self.node_list[i]] = {}
                            self.align_graph[self.node_list[i]][self.node_list[j]] = Weight[1]
                            self.align_graph[self.node_list[j]][self.node_list[i]] = Weight[1]
                    else:
                        if self.align_graph[self.node_list[i]]:
                            self.align_graph[self.node_list[i]][self.node_list[j]] = Weight[2]
                            self.align_graph[self.node_list[j]][self.node_list[i]] = Weight[2]
                        else:
                            self.align_graph[self.node_list[i]] = {}
                            self.align_graph[self.node_list[i]][self.node_list[j]] = Weight[2]
                            self.align_graph[self.node_list[j]][self.node_list[i]] = Weight[2]
                            '''

    def __cal_ELI(self):
        pass
        

class AlignmentGraph(object):
    
    ### A class for alignment graph
    
    def __init__(self, graph1, graph2, sim):
        self.align_node_list = []
        self.align_sub_graph = []
        self.graph1 = graph1
        self.graph2 = graph2

        for sub in graph1:
            self.align_sub_graph.append(Align_sub_graph(sub, graph2[0],sim).align_sub_graph1)        
        l = 0
        m = None
        for sub in align_sub_graph:
            if sub.length() > l:
                m = sub
                l = m.length()
        if m.length() > sqrt(graph2.length()):
            b = List2dot()
            b.list2dot([a.subGraphs[0]], a.filename)






