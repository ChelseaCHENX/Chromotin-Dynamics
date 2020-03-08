import os
os.chdir('/home/fangyuan/hic/codes')

from utils import *
from prody import * # hic conda env: prody uptodate from git: v1.11 
from hicparsers import parseBulkMap

# calc overlap from GNM
# calc gnm
meta = 'tamR'
resol = '50kb'
binsize = resol2bin(resol)
chroms = [str(i) for i in range(1,23)] + ['X']
# chroms = ['16']

# for sp in ['02','03','04','05','06','51']:
for f in os.listdir(f'../data/contact/{meta}'):
    sp = f.strip('.txt')
    for chrom in chroms:
        if not os.path.isfile(f'../data/GNM/{meta}_{sp}.chr{chrom}.{resol}.vc.gnm.npz'):
            fpath = f'../data/contact/{meta}/{f}' # from unbalanced cool files
            try:
                M, bin = parseBulkMap(fpath, chroms=chrom, bin=binsize)

                hic = HiC('population', M, bin) # prody.chromatin.hic
                hic.normalize(method=VCnorm)

                hic.masked = True
                
                gnm = hic.calcGNM(zeros=True)  
                gnm.masked = False
                gnm.fixTail(M.shape[0]) # add zeros at tail to make length == chrom length
                saveModel(gnm, f'../data/GNM/{meta}_{sp}.chr{chrom}.{resol}.vc.gnm.npz', matrices=True)
                print(f'Done for {sp} chr{chrom}')
            except:
                print(f'Failed for {sp} chr{chrom}')
