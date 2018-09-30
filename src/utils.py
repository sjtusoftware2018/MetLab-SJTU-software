#!/usr/bin/python
# -*- coding: utf-8 -*-

###############################
########## node similarity ####
###############################

from rdkit import Chem,DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols
from difflib import SequenceMatcher
from os.path import dirname, abspath
import re

f1=open(dirname(dirname(abspath(__file__))) + '/data/node_list/allmetabolites.txt','r')
set_metabolites=f1.read()
f2=open(dirname(dirname(abspath(__file__))) + '/data/node_list/allreactions.txt','r')
set_reactions=f2.read()

class SIM():
    #Get the ratio of two nodes to be compared.
    # threshold : 0.2
    # mode : metabolite/reaction
    # node1,node2 : The nodes whose similarity will be compared to calculate the ratio.
    # ratio : The result of comparison of two nodes.
    def __init__(self,node1,node2):
        self.threshold = 0.2
        self.mode=''
        ratio=self.sim(node1,node2)
        # print 'The ratio of two',self.mode,'is:',ratio
        print("The ratio of two %s is %f" % (self.mode, ratio))
        self.mode=''

    #Find the Smiles or ECnumber in two files
    # re_node : regular expression of node in f1,f2
    # name : An array of the matching results.
    def get_name(self,node):
        _node = node + '\t(.+?)\r\n'
        re_node = re.compile(_node)
        self.mode='metabolites'
        name = re.findall(re_node,set_metabolites)
        if name == []:
            self.mode='reactions'
            name= re.findall(re_node,set_reactions)
        return name[0]
    
    #Calculate the similarity of two nodes
    # name1,name2 : The corresponding Smiles or EC number of two nodes.
    def sim(self,node1,node2):
        name1 = self.get_name(node1)
        name2 = self.get_name(node2)
        #compare the similarity of metabolites
        #Calculate the ratio of two Smiles.
        if self.mode == 'metabolites':
            try:
                ms = [Chem.MolFromSmiles(name1),Chem.MolFromSmiles(name2)]
                fps = [FingerprintMols.FingerprintMol(x) for x in ms]
                ratio = DataStructs.FingerprintSimilarity(fps[0],fps[1])
                return ratio
            except:
                ratio = self.sdiff(node1,node2)
                return ratio
        #compare the similarity of enzymes
        #Calculate the ratio of two EC number
        else: 
            ECnumber1 = name1.split('\t')
            ECnumber2 = name2.split('\t')
            if ECnumber1 != 'no EC number' and ECnumber2 != 'no EC number':
                ratio=0
                for i in ECnumber1:
                    for j in ECnumber2:
                        ei = i.split('.')
                        ej = j.split('.')
                        temp_ratio = 0
                        if ei[0]==ej[0]:
                            temp_ratio = 0.25
                            if ei[1]==ej[1]:
                                temp_ratio = 0.5
                                if ei[2]==ej[2]:
                                    temp_ratio = 0.75
                                    if ei[3]==ej[3]:
                                        temp_ratio = 1
                        if temp_ratio > ratio:
                            ratio = temp_ratio#choose the max ratio of multiple ECnumber matchs
                return ratio
            else:#If one enzyme has no EC number, we choose to compare their name strings.
                ratio = self.sdiff(node1,node2)
                return ratio
        
    #compare the name string and get the match ratio
    def sdiff(self,s1,s2):
        seq = SequenceMatcher(None,s1,s2)
        ratio = seq.ratio()#0<ratio<1
        if ratio >=0.5:
            return self.threshold
        else:   
            return 0 
    
if __name__ == '__main__':
    node1 = input("please input the name of node1:")
    node2 = input("please input the name of node2:")
    ratio = SIM(node1,node2)

###################################
########## output dot file ########
###################################

class List2dot(object):

    #创建一个类，将列表转化为.dot文件

    number=1  #定义一个类的静态变量，方便之后给节点命名

    id_dict={}  #定义一个字典，记录node_id与number的对应关系
    

    def create_point(self,f,name,node_type):

        #将创建一个点的语句写入.dot文件
        #f是打开的文件，name是节点在图中的名字
        #node_type是节点类型，为'0'表示是代谢物，为'1'表示是反应,为'0'时用默认形状，为'1'时用矩形
        if node_type=='0':
            f.write('\t\t'+str(List2dot.number)+'[label='+'"'+name+'"'+'];\n')
        else:
            f.write('\t\t'+str(List2dot.number)+'[label='+'"'+name+'"'+',shape=rectangle];\n')
        List2dot.number=List2dot.number+1
        

    def create_edge(self,f,id1,id2):

        #将创建一个边的语句写入.dot语句
        #id1,id2是两个节点的变量名
        
        f.write('\t\t'+id1+'->'+id2+';\n')
        

    def list2subgraph(self,f,list1):

        #将给出的列表转换为dot文件中子图的代码

        id_dict={}  #定义一个字典，记录node_id与number的对应关系
    
        f.write('\t'+'subgraph cluster_'+list1[0]+'{\n'+'\t\t'+'label='+'"'+list1[0]+'"'+';\n')#创建一个子图

        for node in list1[1]:
            #创建所有的节点
            one_node=node.split(',')#将节点的字符串分割
            id_dict[one_node[1]]=str(List2dot.number)
            self.create_point(f,one_node[0],one_node[2])

        for key in list1[2]:#遍历字典
            for i in list1[2][key]:
                self.create_edge(f,id_dict[key],id_dict[i])

        f.write('\t'+'}\n')
        

    def list2dot(self,graph_list,graph_name='graph_name',file_name='graph_name.dot'):

        #这里的graph_list是列表的列表，即其中的每一项都是一个列表代表一个子图
        #graph_name即整个图的名字，file_name即生成的文件的名字，这两变量可以设置缺省值
        
        f=open(file_name,'w')#创建一个.dot文件
        f.write('digraph graph_name{\n\tlabel='+'"'+graph_name+'"'+';\n')
        for listi in graph_list:
            self.list2subgraph(f,listi)
        f.write('}\n')
        f.close()
        List2dot.number=1


#以下是一个简单的例子用来验证上面的代码    
L=List2dot()
graph_list=[['test_graph1',['n1,1,1','n2,2,0','n3,3,1','n4,4,0','n5,5,0'], \
          {'1':['5'],'2':['1'],'3':['1','2'],'4':[],'5':['3','4']}],['test_graph2',['n6,1,0','n7,2,0','n8,3,0','n9,4,0','n10,5,0'], \
          {'1':['5'],'2':['1'],'3':['1','2'],'4':[],'5':['3','4']}]]
L.list2dot(graph_list,'test_graph','graph_name.dot')
    
    


