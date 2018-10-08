#!/usr/bin/python2
# -*- encoding:utf-8 -*-
#

from utils import *
from graph import *
from optparse import OptionParser

if __name__ == '__main__':
    parser = OptionParser(description = 'Metabolic Network Alignment')
    parser.add_option('-m', '--model', default = 'all', help = 'The whole genetic metabolic network name (default: all)')
    parser.add_option('-n', '--number', default = 1, help = 'The number of alignment results to be output (default: 1)')
    #parser.add_option()
    (opts, args) = parser.parse_args()

    if opts.m == 'all':
        
    elif opts.m in metabolic_network:
        

    else:
        print "Wrong model name!"
        exit(1)

metabolic_network = ['iAF692.met',
                     'iAF987.met',
                     'iAF1260.met',
                     'iCHOv1.met',
                     'iECO111_1330.met',
                     'iECUMN_1333.met',
                     'iHN637.met',
                     'iIT341.met',
                     'iJB785.met',
                     'iJN678.met',
                     'iJN746.met',
                     'iLB1027_lipid.met',
                     'iLJ478.met',
                     'iMM904.met',
                     'iMM1415.met',
                     'iNF517.met',
                     'iNJ661.met',
                     'iPC815.met',
                     'iRC1080.met',
                     'iSB619.met',
                     'iSbBS512_1146.met',
                     'iSBO_1134.met',
                     'iSDY_1059.met',
                     'iSF_1195.met',
                     'iSSON_1240.met',
                     'iYL1228.met',
                     'iYO844.met',
                     'STM_v1_0.met']


