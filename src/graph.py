
from enum import Enum, IntEnum, unique
import re


'''
Define a enumous class node type
'''
try:
    @unique
    class NodeType(Enum):
        metabolite = 0
        reaction = 1
except ValueError as e:
    print(e)

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
        try:
            self.type = NodeType(node_type)
        except ValueError as e:
            if node_type[0] in 'mM0':
                self.type = NodeType(0)
            elif node_type[0] in 'rR1':
                self.type = NodeType(1)
            else:
                raise ValueError()
        if isinstance(name, str) and isinstance(id, str) and isinstance(reverse, str):
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
        self.subGraphs = self.__read(filename)

    def __read(self,filename):
        """
        Read a network file, and return a list of subgraphs
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
                    itemrlist = line.strip().split(',')
                    sub.addRxn(itemrlist)
                else:
                    """
                    Handling a node entry
                    """
                    itermlist = line.strip().split(',')
                    name, node_id, node_type, reverse = itermlist
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
