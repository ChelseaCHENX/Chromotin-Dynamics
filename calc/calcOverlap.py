import os
os.chdir('/home/fangyuan/projects/hic/codes')

from utils import *
from prody import *
from hicparsers import parseBulkMap

# calc overlap from GNM
# calc gnm
chroms = ['17']
resol = '50kb'
binsize = resol2bin(resol)

for x in ['02','03','04','05','06','51']:
    for chrom in chroms:

        ## load GNM results 
        gnm = loadModel(f'../data/GNM/mcf7_{sp}.chr{chrom}.{resol}.sqrtvc.gnm.npz')
        C = calcCrossCorr(gnm)

        cov = calcCovariance 
        corr = calcCrossCorr(gnm[:500])
        sqfs = calcSqFlucts(gnm[:500])

        