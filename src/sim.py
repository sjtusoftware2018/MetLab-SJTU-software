from rdkit import Chem,DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols
from difflib import SequenceMatcher
import re
f1=open('../data/node_list/allmetabolites.txt','r')
set_metabolites=f1.read()
f2=open('../data/node_list/allreactions.txt','r')
set_reactions=f2.read()

threshold = 0.2
mode=''
def get_name(node):#find the Smiles or ECnumber in two files
	_node = node + '\t(.+?)\r\n'
	re_node = re.compile(_node)
	global mode
	mode='metabolite'
	name = re.findall(re_node,set_metabolites)
	if name == []:
		mode='reaction'
		name= re.findall(re_node,set_reactions)
	return name[0]

def sim(node1,node2):
	name1 = get_name(node1)
	name2 = get_name(node2)
	if mode == 'metabolite':#compare the similarity of metabolites
		try:
			ms = [Chem.MolFromSmiles(name1),Chem.MolFromSmiles(name2)]
			fps = [FingerprintMols.FingerprintMol(x) for x in ms]
			ratio = DataStructs.FingerprintSimilarity(fps[0],fps[1])
			return ratio
		except:
			ratio = sdiff(node1,node2)
			return ratio
	elif mode == 'reaction': #compare the similarity of enzymes
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
		else:
			ratio = sdiff(node1,node2)
			return ratio
		
def sdiff(s1,s2):#compare the string and get the match ratio
	seq = SequenceMatcher(None,s1,s2)
	ratio = seq.ratio()#0<ratio<1
	if ratio >= 0.5:
		return threshold
	else:
		return 0