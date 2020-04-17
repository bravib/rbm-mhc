"""

"""

import pandas as pd
import numpy as np
import re


aa = ['A', 'C', 'G', 'T']
aadict = {aa[k]:k for k in range(len(aa))}
for k, key in enumerate(['a', 'c', 'g', 't']):
    aadict[key] = aadict[aa[k]]

curr_int = np.int16
curr_float = np.float32


# gives the same color to Watson-Crick canonical base pairs (A-T in red, C-G in blue), maybe one can do better
def aa_color(letter):
    if letter in ['C', 'G']:
        return 'blue'
    elif letter in  ['A', 'T']:
        return 'red'
    else:
        return 'black'


def load_DNA_MSA(filename, with_labels=False, with_counts=False, remove_insertions = True, drop_duplicates=True):
    count = 0
    current_seq = ''
    all_seqs = []
    if with_labels:
        all_labels = []
    if with_counts:
        all_counts = []
    with open(filename, 'r') as f:
        for line in f:
            if line[0] == '>':
                all_seqs.append(current_seq)
                current_seq = ''
                if with_labels:
                    all_labels.append(line[1:].replace('\n','').replace('\r',''))
                if with_counts:
                    all_counts.append(int(re.split("-", line)[-1]))
            else:
                current_seq += line.replace('\n', '').replace('\r', '')
                count += 1
        all_seqs.append(current_seq)
        all_seqs = np.array(list(map(lambda x: [aadict[y] for y in x], all_seqs[1:])), dtype=curr_int, order="c")
    if remove_insertions:
        all_seqs = np.asarray(all_seqs[:, ((all_seqs == -1).max(0) == False) ],dtype=curr_int,order='c')
    if drop_duplicates:
        all_seqs = pd.DataFrame(all_seqs).drop_duplicates()
        if with_labels:
            all_labels = np.array(all_labels)[all_seqs.index]
        all_seqs = np.array(all_seqs)
    if with_labels:
        return all_seqs, np.array(all_labels)
    elif with_counts:
        return all_seqs, np.array(all_counts)
    else:
        return all_seqs

    
    
def write_FASTA(filename,all_data,all_labels=None):
    sequences = num2seq(all_data)
    if all_labels is None:
        all_labels = ['S%s'%k for k in range(len(sequences))]
    with open(filename,'w') as fil:
        for seq, label in zip(sequences,all_labels):
            fil.write('>%s\n'%label)
            fil.write('%s\n'%seq)
    return 'done'


def seq2num(string):
    if type(string) == str:
        return np.array([aadict[x] for x in string])[np.newaxis,:]
    elif type(string) ==list:
        return np.array([[aadict[x] for x in string_] for string_ in string])


def num2seq(num):
    return [''.join([aa[x] for x in num_seq]) for num_seq in num]


def distance(MSA,verbose=False):
    B = MSA.shape[0]
    N = MSA.shape[1]
    distance = np.zeros([B,B])
    for b in range(B):
        if verbose:
            if b%1000 ==0:
                print(b)
        distance[b] =  ((MSA[b] != MSA).mean(1))
        distance[b,b] = 2.
    return distance

def count_neighbours(MSA,threshold = 0.1): # Compute reweighting
    B = MSA.shape[0]
    N = MSA.shape[1]
    num_neighbours = np.zeros(B)
    for b in range(B):
        if b%1000 ==0:
            print(b)
        num_neighbours[b] =  ((MSA[b] != MSA).mean(1) < threshold).sum()
    return np.asarray(num_neighbours,dtype=curr_int)



def couplings_to_contacts(couplings,with_gaps=True): # From N x N x n_c x n_c to Average Product Correction.
    if with_gaps:
        W = np.sqrt( (couplings**2).sum(-1).sum(-1) )
    else:
        W = np.sqrt( (couplings[:,:,:-1,:-1]**2).sum(-1).sum(-1) )
    tmp = W
    tmp2 = tmp.sum(1)
    F = tmp - tmp2[np.newaxis,:] * tmp2[:,np.newaxis]/tmp2.sum()
    return F


def compare_contact_maps(F,contact_map, distant = False,return_list=False,filter_alignment =None):

    F2 = F.copy()
    contact_map = contact_map.copy()
    if filter_alignment is not None:
        F2 = F2[filter_alignment,:][:,filter_alignment]
        contact_map = contact_map[filter_alignment,:][:,filter_alignment]

    n_sites = F2.shape[0]
    if distant:
        for i in range(n_sites):
            for j in range(n_sites):
                if np.abs(i-j)<5:
                    F2[i,j] = -1
                    contact_map[i,j] = 0

    n_contacts = int(contact_map.sum()/2)
    contact_map = contact_map.flatten()
    for i in range(n_sites):
        for j in range(i,n_sites):
            F2[i,j]= -1
    l = np.argsort(F2.flatten())[::-1]
    corrects = np.zeros(n_contacts)
    if return_list:
        I,J = np.unravel_index(  l, (n_sites,n_sites) )
        liste_pairs = zip(I,J)

    for i in range(n_contacts):
        if contact_map[l[i]]:
            corrects[i]=1

    if return_list:
        is_correct = corrects.copy()


    corrects = np.cumsum(np.array(corrects))/ np.arange(1,n_contacts+1)

    if return_list:
        return liste_pairs, is_correct, corrects
    else:
        return corrects
