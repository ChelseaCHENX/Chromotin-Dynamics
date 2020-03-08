import os
os.chdir('/home/fangyuan/projects/hic/codes')
from prody import *
from hicparsers import parseBulkMap
import matplotlib.pyplot as plt
import numpy as np
import pandas
import cooler
import os 

# after balancing 
for f in os.listdir('../data/cool_unbalanced'):
    if f.endswith('.cool'):
        sp = f.split('.')[0]
        fpath = f'../data/cool_unbalanced/{f}'

        c = cooler.Cooler(fpath)
        df = c.pixels(join=True)[:]
        df.to_csv(f'../data/contact/{sp}.txt', sep='\t', index=False, header=False)

'''
## process coontact map with iterative balancing method
f = 'SRR6493702.GRCh38.mapq_30.50000.cool'
fpath = f'../data/cool/{f}'
c = cooler.Cooler(fpath)
c.bins()[:5]

chr_sel = '17'
ar1 = c.matrix(sparse=True, balance=False).fetch(chr_sel).toarray()
ar2 = c.matrix(sparse=True, balance=True).fetch(chr_sel).toarray()[::-1]


fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(111)
im = ax.matshow(np.log10(ar1), cmap='YlOrRd')
fig.colorbar(im)
plt.title('Raw contact map')

fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(111)
im = ax.matshow(np.log10(ar2), cmap='YlOrRd')
fig.colorbar(im)
plt.title('Balanced (VC) contact map')

### process GNM
chroms = ['17']
binsize = 50000
sample = 'SRR6493702'
display = True 
matched = False

chrom = ['17']
fname = f'../data/contact/{sample}.txt'
M, bin = parseBulkMap(fname, chroms=chrom, bin=binsize)

hic = HiC('population', M, bin)
hic.normalize(method=SQRTVCnorm)

if display:
    from matplotlib.pyplot import *
    ion()
    hic.view()

hic.masked = True
gnm = hic.calcGNM(zeros=True)  
calcCovariance 
cov = calcCrossCorr(gnm[:500])
sqfs = calcSqFlucts(gnm[:500])

fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(111)
im = ax.matshow(np.log10(cov), cmap='YlOrRd')
fig.colorbar(im)
plt.title('GNM covariance map')

###
calcSpectralOverlap(gnm[:500], gnm[:100], weighted=True, turbo=False)
'''