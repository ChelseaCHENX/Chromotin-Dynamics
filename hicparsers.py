from numpy import *
import h5py
from collections import defaultdict
from collections.abc import MutableMapping

class contactdict(MutableMapping):
    """A dictionary that applies an arbitrary key-altering
       function before accessing the keys"""

    def __init__(self, *args, **kwargs):
        self.store = dict()
        self.update(dict(*args, **kwargs))  # use the free update to set keys

    def __getitem__(self, key):
        if key in self.store:
            return self.store[key]
        elif key[::-1] in self.store:
            key = key[::-1]
            data = [[y, x, v] for x, y, v in self.store[key]] 
            return data

        return self.store[key]

    def __setitem__(self, key, value):
        self.store[self._transkey(key)] = value

    def __delitem__(self, key):
        del self.store[self._transkey(key)]

    def __iter__(self):
        return iter(self.store)

    def __contains__(self, key):
        return self._transkey(key) in self.store

    def __len__(self):
        return len(self.store)

    def _transkey(self, key):
        if key[::-1] in self.store:
            return key[::-1]
        return key

    def __repr__(self):
        return repr(self.store)

    def add(self, key, contact):
        if key not in self.store:
            if key[::-1] in self.store:
                key = key[::-1]
                contact = [contact[1], contact[0], contact[-1]]
        
        if key in self.store:
            self.store[key].append(contact)
        else:
            self.store[key] = [contact]


def parseBulkMap(filename, chroms=None, bin=None):
    data_dir = contactdict()
    lengths = defaultdict(int)
    
    if chroms is not None:
        if isinstance(chroms, str):
            chroms = [chroms]
    else:
        chroms = [str(i+1) for i in range(19)] + ['X']

    n_chroms = len(chroms)

    with open(filename, 'r') as f:
        for line in f:
            txt = line.strip()
            if txt.startswith('chrom1'):
                continue

            info = txt.split('\t')
            chr_A, start_A, end_A = info[:3]
            chr_B, start_B, end_B = info[3:6]
            num_obs = info[6]

            chA = chr_A # chA = chr_A[3:]
            chB = chr_B # chB = chr_B[3:]

            if not chA in chroms or not chB in chroms:
                continue

            if bin is None:
                bin = int(end_A) - int(start_A) + 1
            else:
                bin = int(bin)

            posA = (int(start_A) - 1) // bin
            posB = (int(start_B) - 1) // bin
            value = int(num_obs)

            if posA + 1 > lengths[chA]:
                lengths[chA] = posA + 1

            if posB + 1 > lengths[chB]:
                lengths[chB] = posB + 1

            chpair = (chA, chB)
            data_dir.add(chpair, [posA, posB, value])
            if chA == chB and posA != posB:
                data_dir.add(chpair, [posB, posA, value])

    
    #total_length = sum(list(lengths.values()))
    #labels = [''] * total_length

    map_mat = []
    for i in range(n_chroms):
        map_row = [None] * n_chroms
        map_mat.append(map_row)

    for i, chA in enumerate(chroms):
        for j, chB in enumerate(chroms):
            if i < j:
                continue
            lA = lengths[chA]
            lB = lengths[chB]

            M = zeros((lA, lB))

            data = data_dir[chA, chB]
            for x, y, v in data:
                # try:
                M[x, y] = v
                # except IndexError:
                #     print(M.shape)

            map_mat[i][j] = M
            if i != j:
                map_mat[j][i] = M.T

    rows = []
    for map_row in map_mat:
        R = hstack(map_row)
        rows.append(R)
    
    M = vstack(rows)

    return M, bin

def parseSCMap(filename, chroms=None, bin=None, contact_value=False):
    if chroms is not None:
        if isinstance(chroms, str):
            chroms = [chroms]
    else:
        chroms = [str(i+1) for i in range(19)] + ['X']

    col = 3 if contact_value else 2

    n_chroms = len(chroms)
    map_mat = []
    for i in range(n_chroms):
        map_row = [None] * n_chroms
        map_mat.append(map_row)

    with h5py.File(filename, 'r') as data_dir:
        if bin is None:
            bin = data_dir['bin_size'][()]
        bin = int(bin)

        L = []
        C = []
        last = 0
        for chrom in chroms:
            pos = data_dir['/bp_pos/%s'%chrom][()]
            coord = unique((pos - 1) // bin)
            l = max(coord) + 1

            coord += last
            last = max(coord) + 1

            C.extend(coord)
            L.append(l)

        for i, chri in enumerate(chroms):
            for j, chrj in enumerate(chroms):
                M = zeros((L[i], L[j]))
                map_mat[i][j] = M

        for i, chri in enumerate(chroms):
            for j, chrj in enumerate(chroms):
                ctt_dir = '/expr_contacts/%s/%s'%(chri, chrj)
                if ctt_dir in data_dir:
                    M = map_mat[i][j]
                    raw = data_dir[ctt_dir][()]
                    X = (raw[:, 0] - 1) // bin
                    Y = (raw[:, 1] - 1) // bin
                    V = raw[:, col]

                    for x, y, v in zip(X, Y, V):
                        M[x, y] += v

                    if i == j:
                        D = diag(diag(M))
                        M = M + M.T - D
                    
                    map_mat[i][j] = M
                    if i != j:
                        map_mat[j][i] = M.T
    
    rows = []
    for map_row in map_mat:
        R = hstack(map_row)
        rows.append(R)
    
    M = vstack(rows)

    return M, C, bin

def parseSCDistMap(filename, chroms=None, bin=None):
    if chroms is not None:
        if isinstance(chroms, str):
            chroms = [chroms]
    else:
        chroms = [str(i+1) for i in range(19)] + ['X']

    n_chroms = len(chroms)
    map_mat = []
    for i in range(n_chroms):
        map_row = [None] * n_chroms
        map_mat.append(map_row)

    with h5py.File(filename, 'r') as data_dir:
        if bin is None:
            bin = data_dir['bin_size'][()]
        bin = int(bin)

        L = []
        last = 0
        for chrom in chroms:
            pos = data_dir['/bp_pos/%s'%chrom][()]
            l = unique((pos - 1) // bin) + last
            last = max(l) + 1
            L.extend(l)

        for i, chri in enumerate(chroms):
            for j, chrj in enumerate(chroms):
                ctt_dir = '/dists/%s/%s'%(chri, chrj)
                if ctt_dir in data_dir:
                    M = data_dir[ctt_dir][()]
                    
                    map_mat[i][j] = M
                    if i != j:
                        map_mat[j][i] = M.T
    
    rows = []
    for map_row in map_mat:
        R = hstack(map_row)
        rows.append(R)
    
    M = vstack(rows)

    return M, L, bin

def parseSCContactMap(filename, chrom=None, n_loci=None, k=2.30624):
    from distance import dist2contact

    D, C, bin = parseSCDistMap(filename, chroms=chrom)

    # convert distance to contacts
    n = len(D)
    I = eye(n)

    D += I

    M = dist2contact(D, k)

    cap = dist2contact(1.15, k)
    M[M > cap] = cap
    M = floor(M)

    # add margins
    if n_loci is None:
        n_loci = max(C) + 1
    M_ = zeros((n_loci, n_loci))
    M_[ix_(C, C)] = M

    return M_, bin

def getSCLengths(filename, chroms=None, bin=None, mapped=False):
    if chroms is not None:
        if isinstance(chroms, str):
            chroms = [chroms]
    else:
        chroms = [str(i+1) for i in range(19)] + ['X']

    n_chroms = len(chroms)
    map_mat = []
    for i in range(n_chroms):
        map_row = [None] * n_chroms
        map_mat.append(map_row)

    with h5py.File(filename, 'r') as data_dir:
        if bin is None:
            bin = data_dir['bin_size'][()]
        bin = int(bin)

        L = []
        for chrom in chroms:
            pos = data_dir['/bp_pos/%s'%chrom][()]
            if mapped:
                l = (max(pos) - min(pos)) // bin
            else:
                l = max(pos) // bin
            l += 1
            L.append(l)

    return L

def parseContactFile_GSE48262(filename, bin, chroms=None):
    data_dir = contactdict()
    lengths = defaultdict(int)
    bin = int(bin)
    
    if chroms is not None:
        if isinstance(chroms, str):
            chroms = [chroms]
    else:
        chroms = [str(i+1) for i in range(19)] + ['X']

    n_chroms = len(chroms)

    skip = True
    with open(filename, 'r') as f:
        for line in f:
            if skip:
                skip = False
                continue

            txt = line.strip()
            info = txt.split('\t')
            chr_A, pos_A = info[:2]
            chr_B, pos_B = info[2:4]

            chA = chr_A
            chB = chr_B

            if not chA in chroms or not chB in chroms:
                continue

            x = int(pos_A) // bin
            y = int(pos_B) // bin
            value = 1

            if x + 1 > lengths[chA]:
                lengths[chA] = x + 1

            if y + 1 > lengths[chB]:
                lengths[chB] = y + 1

            chpair = (chA, chB)
            data_dir.add(chpair, [x, y, value])
            if chA == chB and x != y:
                data_dir.add(chpair, [y, x, value])

    
    #total_length = sum(list(lengths.values()))
    #labels = [''] * total_length

    map_mat = []
    for i in range(n_chroms):
        map_row = [None] * n_chroms
        map_mat.append(map_row)

    for i, chA in enumerate(chroms):
        for j, chB in enumerate(chroms):
            if i < j:
                continue
            lA = lengths[chA]
            lB = lengths[chB]

            M = zeros((lA, lB))

            data = data_dir[chA, chB]
            for x, y, v in data:
                # try:
                M[x, y] = v
                # except IndexError:
                #     print(M.shape)

            map_mat[i][j] = M
            if i != j:
                map_mat[j][i] = M.T

    rows = []
    for map_row in map_mat:
        R = hstack(map_row)
        rows.append(R)
    
    M = vstack(rows)

    return M

if __name__ == '__main__':
    from matplotlib.pyplot import *
    from prody import *

    ion()

    # filename = '../../Data/mouse/GSM2123564_Haploid_mESC_population_hic.txt'
    # #data = loadtxt(filename, dtype=object, delimiter='\t')
    # M0, bin = parseBulkMap(filename, chroms=['17'], bin=100e3)

    # figure()
    # showMatrix(M0, percentile=1)

    # filename = '../../Data/mouse/Cell-1.hdf5'
    # M, bin = parseSCMap(filename, chroms=['17'], contact_value=True)

    # figure()
    # showMatrix(M)

    # filename = '../../Data/mouse/GSE48262_cell-1.txt'
    # M = parseContactFile_GSE48262(filename, 50e3, chroms=['17'])

    # figure()
    # showMatrix(M)

    # L = getSCLengths('../../Data/mouse/Cell-1.hdf5')
