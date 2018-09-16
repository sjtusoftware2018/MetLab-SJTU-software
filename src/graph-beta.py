#!/usr/bin/python
# -*- coding: utf-8 -*-

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
                   If no list is given, the graph is initialized empty.
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

class Align_node(object):
    """
    """
    def __init__(self, node1, node2, weight):
        if isinstance(node1, Node) and isinstance(node2, Node):
            self.node1 = node1
            self.node2 = node2
            self.weight = weight
        else:
            raise TypeError()

Weight = [1,0.5,0.2]

class Align_sub_graph(object):
    """
    """
    def __init__(self, sub_graph1, sub_graph2, sim):
        self.node_list = []
        self.align_graph = {}
        self.graphs = [sub_graph1, sub_graph2]
        for node1 in sub_graph1.node_list:
            for node2 in sub_graph2.node_list:
                self.__add_node(node1, node2, sim)
        self.__add_edge()

    def __add_node(self, node1, node2, sim):
        node = Align_node(node1, node2, sim(node1, node2))
        if node not in self.node_list:
            self.node_list.append(node)
        else:
            pass

    def __add_edge(self):
        for i in xrange(len(self.node_list)-1):
            for j in xrange(i+1,len(self.node_list)):
                if self.sub_graph1.isEdge(self.node_list[i][0],self.node_list[j][0])\
                and self.sub_graph2.isEdge(self.node_list[i][1],self.node_list[j][1]):
                    if self.align_graph[self.node_list[i]]:
                        self.align_graph[self.node_list[i]][self.node_list[j]] = Weight[0]
                        self.align_graph[self.node_list[j]][self.node_list[i]] = Weight[0]
                    else:
                        self.align_graph[self.node_list[i]] = {}
                        self.align_graph[self.node_list[i]][self.node_list[j]] = Weight[0]
                        self.align_graph[self.node_list[j]][self.node_list[i]] = Weight[0]
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

    def __cal_ELI(self):
        


class AlignmentGraph(object):
    """A class for alignment graph
    """
    def __init__(self, graph1, graph2, sim):
        self.align_node_list = []
        self.sub_graph = []

        for sub in graph1:
            Align_sub_graph(sub, graph2[0], sim)
        for 

    def __sim():

    def __construct_align_graph():

    def __optimize():

    def get_alignment_result():

    def _show_graph():

    def _present_graph():


def _show_graph(objs, edge_func, swap_source_target,
                max_depth=3, extra_ignore=(), filter=None, too_many=10,
                highlight=None, filename=None, extra_info=None,
                refcounts=False, shortnames=True, output=None,
                cull_func=None):
    if not _isinstance(objs, (list, tuple)):
        objs = [objs]

    is_interactive = False
    if filename and output:
        raise ValueError('Cannot specify both output and filename.')
    elif output:
        f = output
    elif filename and filename.endswith('.dot'):
        f = codecs.open(filename, 'w', encoding='utf-8')
        dot_filename = filename
    elif IS_INTERACTIVE:
        is_interactive = True
        f = StringIO()
    else:
        fd, dot_filename = tempfile.mkstemp(prefix='objgraph-',
                                            suffix='.dot', text=True)
        f = os.fdopen(fd, "w")
        if getattr(f, 'encoding', None):
            # Python 3 will wrap the file in the user's preferred encoding
            # Re-wrap it for utf-8
            import io
            f = io.TextIOWrapper(f.detach(), 'utf-8')
    f.write('digraph ObjectGraph {\n'
            '  node[shape=box, style=filled, fillcolor=white];\n')
    queue = []
    depth = {}
    ignore = set(extra_ignore)
    ignore.add(id(objs))
    ignore.add(id(extra_ignore))
    ignore.add(id(queue))
    ignore.add(id(depth))
    ignore.add(id(ignore))
    ignore.add(id(sys._getframe()))   # this function
    ignore.add(id(sys._getframe().f_locals))
    ignore.add(id(sys._getframe(1)))  # show_refs/show_backrefs
    ignore.add(id(sys._getframe(1).f_locals))
    for obj in objs:
        f.write('  %s[fontcolor=red];\n' % (_obj_node_id(obj)))
        depth[id(obj)] = 0
        queue.append(obj)
        del obj
    gc.collect()
    nodes = 0
    while queue:
        nodes += 1
        # The names "source" and "target" are reversed here because
        # originally there was just show_backrefs() and we were
        # traversing the reference graph backwards.
        target = queue.pop(0)
        tdepth = depth[id(target)]
        f.write('  %s[label="%s"];\n' % (_obj_node_id(target),
                                         _obj_label(target, extra_info,
                                                    refcounts, shortnames)))
        h, s, v = _gradient((0, 0, 1), (0, 0, .3), tdepth, max_depth)
        if inspect.ismodule(target):
            h = .3
            s = 1
        if highlight and highlight(target):
            h = .6
            s = .6
            v = 0.5 + v * 0.5
        f.write('  %s[fillcolor="%g,%g,%g"];\n'
                % (_obj_node_id(target), h, s, v))
        if v < 0.5:
            f.write('  %s[fontcolor=white];\n' % (_obj_node_id(target)))
        if hasattr(getattr(target, '__class__', None), '__del__'):
            f.write('  %s->%s_has_a_del[color=red,style=dotted,'
                    'len=0.25,weight=10];\n' % (_obj_node_id(target),
                                                _obj_node_id(target)))
            f.write('  %s_has_a_del[label="__del__",shape=doublecircle,'
                    'height=0.25,color=red,fillcolor="0,.5,1",fontsize=6];\n'
                    % (_obj_node_id(target)))
        if tdepth >= max_depth:
            continue
        if cull_func is not None and cull_func(target):
            continue
        neighbours = edge_func(target)
        ignore.add(id(neighbours))
        n = 0
        skipped = 0
        for source in neighbours:
            if id(source) in ignore:
                continue
            if filter and not filter(source):
                continue
            if n >= too_many:
                skipped += 1
                continue
            if swap_source_target:
                srcnode, tgtnode = target, source
            else:
                srcnode, tgtnode = source, target
            elabel = _edge_label(srcnode, tgtnode, shortnames)
            f.write('  %s -> %s%s;\n' % (_obj_node_id(srcnode),
                                         _obj_node_id(tgtnode), elabel))
            if id(source) not in depth:
                depth[id(source)] = tdepth + 1
                queue.append(source)
            n += 1
            del source
        del neighbours
        if skipped > 0:
            h, s, v = _gradient((0, 1, 1), (0, 1, .3), tdepth + 1, max_depth)
            if swap_source_target:
                label = "%d more references" % skipped
                edge = "%s->too_many_%s" % (_obj_node_id(target),
                                            _obj_node_id(target))
            else:
                label = "%d more backreferences" % skipped
                edge = "too_many_%s->%s" % (_obj_node_id(target),
                                            _obj_node_id(target))
            f.write('  %s[color=red,style=dotted,len=0.25,weight=10];\n'
                    % edge)
            f.write('  too_many_%s[label="%s",shape=box,height=0.25,'
                    'color=red,fillcolor="%g,%g,%g",fontsize=6];\n'
                    % (_obj_node_id(target), label, h, s, v))
            f.write('  too_many_%s[fontcolor=white];\n'
                    % (_obj_node_id(target)))
    f.write("}\n")

    if output:
        return

    if is_interactive:
        return graphviz.Source(f.getvalue())
    else:
        # The file should only be closed if this function was in charge of
        # opening the file.
        f.close()
        print("Graph written to %s (%d nodes)" % (dot_filename, nodes))
        _present_graph(dot_filename, filename)


def _present_graph(dot_filename, filename=None):
    """Present a .dot file to the user in the requested fashion.
    If ``filename`` is provided, runs ``dot`` to convert the .dot file
    into the desired format, determined by the filename extension.
    If ``filename`` is not provided, tries to launch ``xdot``, a
    graphical .dot file viewer.  If ``xdot`` is not present on the system,
    converts the graph to a PNG.
    """
    if filename == dot_filename:
        # nothing to do, the user asked for a .dot file and got it
        return
    if not filename and _program_in_path('xdot'):
        print("Spawning graph viewer (xdot)")
        subprocess.Popen(['xdot', dot_filename], close_fds=True)
    elif _program_in_path('dot'):
        if not filename:
            print("Graph viewer (xdot) not found, generating a png instead")
            filename = dot_filename[:-4] + '.png'
        stem, ext = os.path.splitext(filename)
        cmd = ['dot', '-T' + ext[1:], '-o' + filename, dot_filename]
        dot = subprocess.Popen(cmd, close_fds=False)
        dot.wait()
        if dot.returncode != 0:
            # XXX: shouldn't this go to stderr or a log?
            print('dot failed (exit code %d) while executing "%s"'
                  % (dot.returncode, ' '.join(cmd)))
        else:
            print("Image generated as %s" % filename)
    else:
        if not filename:
            print("Graph viewer (xdot) and image renderer (dot) not found,"
                  " not doing anything else")
        else:
print("Image renderer (dot) not found, not doing anything else")