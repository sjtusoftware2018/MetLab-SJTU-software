
class Graph_to_dot(object):
	s=0
	def make_dotfile(self,Graph_x):
		s=s+1
		f=open('dot'+str(s)+'.dot','w')

		f.write('digraph graph1{\n')

		j=0


		for subG in Graph_x.subGraphs:
			j=j+1
			f.write('	subgraph cluster_subgraph'+str(j)+'{\n')

			for node_i in subG.nodelist:

				f.write('		'+node_i.name+';\n')

			for key in subG.graph.keys():
				for value in subG.graph[key]:
					f.write('		'+key.name+'->'+value.name+';\n')

			f.write('	}\n')

		f.write('}\n')
		f.close()

		return 'dot'+str(s)+'.dot'




