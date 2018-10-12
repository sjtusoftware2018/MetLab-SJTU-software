#!/usr/bin/python

from rdkit import Chem,DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols
from difflib import SequenceMatcher
import sys


if __name__ == '__main__':
    if len(sys.argv) is not 3:
        print "Error!"

    else:
        try:
            # print sys.argv
            r = DataStructs.FingerprintSimilarity(FingerprintMols.FingerprintMol(Chem.MolFromSmiles(sys.argv[1])), FingerprintMols.FingerprintMol(Chem.MolFromSmiles(sys.argv[2])))
            print(r)
            # return r
        except:
            # return 0
            print "Error"