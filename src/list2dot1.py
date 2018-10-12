#coding: UTF-8

class List2dot(object):

    #创建一个类，将列表转化为.dot文件

    number=1  #定义一个类的静态变量，方便之后给节点命名
    

    def create_point(self,f,name,node_type):

        #将创建一个点的语句写入.dot文件
        #f是打开的文件，name是节点在图中的名字
        #node_type是节点类型，为'0'表示是代谢物，为'1'表示是反应,为'0'时用默认形状，为'1'时用矩形
        if node_type is '0':
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
    
        f.write('\t'+'subgraph cluster_'+filter(str.isalnum, list1.net_name)+'{\n'+'\t\t'+'label='+'"'+list1.net_name+'"'+';\n')#创建一个子图

        for node in list1.node_list:
            #创建所有的节点
            
            id_dict[node.id]=str(List2dot.number)
            self.create_point(f,node.id,node.type)
        '''
        for key in list1[2]:#遍历字典
            for i in list1[2][key]:
                self.create_edge(f,id_dict[key],id_dict[i])
        '''

        for k,v in list1.graph.items():
            for i in v:
                self.create_edge(f,id_dict[k],id_dict[i])

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

