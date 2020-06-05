

import numpy as np
curr_int = np.int16

def align_to_seed(seqs_seed, lea, lepmin, lepmax, seqs_new, weights=None):
    import matlab.engine
    eng = matlab.engine.start_matlab()
    eng.addpath(rootf + '/Align_utils')
    if weights is None:
        alignc_new = eng.balign_peps_seed(seqs_seed, lea, lepmin, lepmax, seqs_new)
    else:
        alignc_new = eng.balign_peps_seed(seqs_seed, lea, lepmin, lepmax, seqs_new, weights) # lea -> final length, lepmin to lepmax -> length included in the profile itself
    eng.quit()
    msa_al = np.asarray(alignc_new).astype(curr_int) - 1 # needed for different indexing convention in matlab
    return msa_al


import sys
sys.path.append(rootf + '/PGM3/source/')
sys.path.append(rootf + '/PGM3/utilities/')
sys.path.append(rootf + '/Align_utils/')

import argparse
parser = argparse.ArgumentParser() 
parser.add_argument('-sseed', nargs = '*', type = str, required=True)
parser.add_argument('-sseqs', nargs = '*', type = str, required=True)
parser.add_argument('-SA', nargs = '*', type = int, required=True)
parser.add_argument('-SAmin', nargs = '*', type = int, required=True)
parser.add_argument('-SAmax', nargs = '*', type = int, required=True)
parser.add_argument('-yw', nargs = '*', type = int, required = True)
args = parser.parse_args()

seqs_seed_name = args.sseed[0]
seqs_new_name = args.sseqs[0]

seqs_seed = []
with open(seqs_seed_name) as f:
    for line in f:
        linesplit = line.strip().split('\t')
        seqs_seed.append(linesplit[0])

seqs_new = []
with open(seqs_new_name) as f:
    for line in f:
        linesplit = line.strip().split('\t')
        seqs_new.append(linesplit[0])

SA = args.SA[0]
SAmin = args.SAmin[0]
SAmax = args.SAmax[0]
yw = args.yw[0]

if yw == 1:
    seqs_new_n = align_to_seed(seqs_seed, SA, SAmin,SAmax, seqs_new, yw)
if yw == 0:
    seqs_new_n = align_to_seed(seqs_seed, SA, SAmin,SAmax, seqs_new) 

np.savetxt(rootf + '/Align_utils/aligned_temp.txt', seqs_new_n, fmt='%i')
