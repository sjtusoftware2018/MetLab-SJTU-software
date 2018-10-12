#!/usr/bin/python2
# -*- encoding:utf-8 -*-
#

from utils import *
from graph_beta import *
from list2dot1 import *


a = Graph('iYL1228.met')
b = Graph('test.met')
# for i in a.subGraphs[0].node_list:
#    print(i.name)
c = AlignmentGraph(a,b,SIM.sim())
d = List2dot()
d.list2dot(c.align_sub_graph)