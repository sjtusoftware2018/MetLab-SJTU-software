from rdkit import Chem,DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols
from difflib import SequenceMatcher
import re
# f1 : The data file of metabolites
# f2 : The data file of reactions
f1=open('allmetabolites.txt','r')
set_metabolites=f1.read()
f2=open('allreactions.txt','r')
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
		print 'The ratio of two',self.mode,'is:',ratio
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
	node1 = raw_input("please input the name of node1:")
	node2 = raw_input("please input the name of node2:")
	ratio = SIM(node1,node2)
