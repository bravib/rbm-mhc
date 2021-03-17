

import argparse
parser = argparse.ArgumentParser() 

parser.add_argument('-hla', nargs = '+', type = str, required= False, help='HLA-I alleles') # For custom dataset, if not given, only scoring by RBM is performed; For IEDB dataset, if not given, it is set to 'Haplotype1'
parser.add_argument('-i', type = str, required=False, help='Name of input file (with peptide sequences)') # if not given, IEDB sequences are used
parser.add_argument('-o', type = str, required=False, help='Name of output folder') # Required for custom datasets. If not given, it is set to 'IEDB_out_Haplotype1' - see data folders available at https://github.com/bravib/rbm-mhc
parser.add_argument('-nameo', type = str, required=False, help='String to appear in file names')
parser.add_argument('-rl', nargs = '+', type = int, required=False, help='Range of peptide lengths') # if not given, it is set to 8,9,10,11
parser.add_argument('-mt', type = int, required=False, default = 1, help='0 disables training of the method and reads an existing model')
parser.add_argument('-deg', type = float, required=False, default = 1, help='Fraction of data for training RBM')
parser.add_argument('-hu', type = int, required=False,  default = 10, help='Number of RBM hidden units')
parser.add_argument('-l12', type = float, required=False, default = 0.001, help='L^1_2 RBM regularization')
parser.add_argument('-niter', type = int, required=False, default = 100, help='Number of iterations for training the RBM')
parser.add_argument('-niterd', type = int, required=False, default = 1000, help='Number of iterations for training the HLA-I classifier')
parser.add_argument('-perc', type = float, required=False, default = 0.1, help='Fraction of labelled data for training the HLA-I classifier')
parser.add_argument('-al', type = int, required=False, default = 1, help='Disable re-iteration of alignment after classification') # Default: 1 re-iterated alignment when multiple peptide lengths
parser.add_argument('-score', type = str, required=False, help='Name of file of peptides to score')
parser.add_argument('-ba', type = int, default = 0, required=False, help='Choose data among positive binding assays instead of MS')

parser.add_argument('-rwnms', type = int, default = 0, required=False, help='Re-weight to correct for differences between MS and non-MS aa frequency')
parser.add_argument('-rwhp', type = int,  default = 0, required=False, help='Re-weight to correct for differences between MS and Human Proteome aa frequency')
parser.add_argument('-gen', type = int,  default = 0, required=False, help='Generate # synthetic peptides of given HLA specificity')
parser.add_argument('-fig', type = int,  default = 0, required=False, help='Print motif figures')
args = parser.parse_args()

if args.mt not in [0,1]:
   print('Binary option -mt require 0 or 1')
   exit()

if args.al not in [0,1]:
   print('Binary option -al require 0 or 1')
   exit()

if args.rwnms not in [0,1]:
   print('Binary option -al require 0 or 1')
   exit()

if args.rwhp not in [0,1]:
   print('Binary option -al require 0 or 1')
   exit()

if args.rwhp == args.rwnms == 1:
   print('Cannot enable re-weighting by both non-MS and Human Proteome aa frequency')
   exit()

if args.o is None and args.i is not None:
   print('Need an output folder')
   exit()

if args.perc == 0:
   print('Percentage of labelled data should be > 0')
   exit()

maketraining = args.mt 
Nrun = 2
if args.al == 0:
    Nrun = 1

makereweighting = 0
if args.rwnms or args.rwhp:
    makereweighting = 1

if args.score is None:
    NAranking = 0
else:
    NAranking = 1

makegeneration = args.gen # generate new peptides with a given HLA-I specificity
makefig = args.fig


if args.rl is None:
    range_len = [8, 9, 10, 11] 
else:
    range_len = list(args.rl)

LL=len(range_len)
if LL == 1:
   CC = 20
else:
   CC = 21

# set SA, SAmin etc. with the lengths provided
if set(range_len) == set([8, 9, 10, 11]) or set(range_len) == set([8, 9, 10]) or set(range_len) == set([8, 9]) or set(range_len) == set([9]):
    SA = 9
    SAmin = 8
    SAmax = 11
elif set(range_len) == set([9, 10, 11]) or set(range_len) == set([9, 10]):
    SA = 9
    SAmin = 9
    SAmax = 11    
else: 
   print('peptide lengths not recognized - see options admitted')
   exit()

import sys,os,pickle
sys.path.append(rootf + '/PGM/source/')
sys.path.append(rootf + '/PGM/utilities/')
sys.path.append(rootf + '/Align_utils/')
from common_imports import set_num_threads
set_num_threads(1) # Set the number of cores. Must be executed before importing numpy&numba.
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import rbm,utilities
import Proteins_utils, RBM_utils, utilities,sequence_logo,plots_utils
import csv
import importlib

from scipy.stats import percentileofscore
import subprocess


#RBM and Decoder parameters  
l12 = args.l12
hu = args.hu
n_iter = args.niter # Number of epochs
n_epochs = args.niterd

# Use of data 
deg  = args.deg # degradation of training set
perc = args.perc
yesreal = 0
if args.i is None:
    yesreal = 0
else:
    yesreal = 1

## Example haplotypes - give ''synthetic individual'' samples
hap1 = ['HLA-A*31:01', 'HLA-B*27:05', 'HLA-C*02:02']
hap2 = ['HLA-A*24:02', 'HLA-B*39:01', 'HLA-C*12:03']
hap3 = ['HLA-A*02:03', 'HLA-B*27:07', 'HLA-C*01:02']
hap4 = ['HLA-A*68:02', 'HLA-B*15:01', 'HLA-C*03:03']
hap5 = ['HLA-A*11:01', 'HLA-B*35:01', 'HLA-C*03:04'] 
hap6 = ['HLA-A*02:01', 'HLA-B*27:02', 'HLA-C*14:02'] 
hap7 = ['HLA-A*29:02', 'HLA-B*44:03', 'HLA-C*16:01']
hap8 = ['HLA-A*02:07', 'HLA-B*40:01', 'HLA-C*01:02']
hap9 = ['HLA-A*02:01', 'HLA-B*40:02', 'HLA-C*15:02'] 
hap10 = ['HLA-A*24:02', 'HLA-B*54:01', 'HLA-C*01:02'] 
hap11 = ['HLA-A*01:01', 'HLA-B*18:01', 'HLA-C*05:01']
hap12 = ['HLA-A*68:02', 'HLA-B*57:01', 'HLA-C*06:02']
hap13 = ['HLA-A*02:04', 'HLA-B*45:01', 'HLA-C*06:02'] 
hap14 = ['HLA-A*29:02', 'HLA-B*58:01', 'HLA-C*08:02']
hap15 = ['HLA-A*01:01', 'HLA-B*35:01', 'HLA-C*04:01']
hap16 = ['HLA-A*32:01', 'HLA-B*27:02', 'HLA-C*02:02']


Haplotype1 = ['HLA-A*01:01', 'HLA-A*03:01','HLA-B*07:02', 'HLA-B*08:01', 'HLA-C*07:02', 'HLA-C*07:01']
Haplotype2 = hap1 + hap2
Haplotype3 = hap3 + hap4
Haplotype4 = hap5 + hap6
Haplotype5 = hap7 + hap8
Haplotype6 = hap9 + hap10
Haplotype7 = hap11 + hap12
Haplotype8 = hap13 + hap14
Haplotype9 = hap15 + hap16
Haplotype10 = ['HLA-A*02:01', 'HLA-A*03:01','HLA-B*51:01', 'HLA-B*08:01', 'HLA-C*01:02', 'HLA-C*07:01']


if args.hla is None:
    if yesreal == 0:
        set_hla = Haplotype1
        str_hla = 'Haplotype1'
        list_hla_l = list(np.unique(set_hla)) # if 2 same alleles, count as 1
        nc = len(list_hla_l)
    if yesreal == 1:
        nc = 1
else:
    set_hla = args.hla
    str_hla = 'custom_haplotype'
    if len(list(np.unique(set_hla))) != len(list(set_hla)):
        list_hla_l = list(np.unique(set_hla)) # if 2 same alleles, count as 1
    else:
        list_hla_l = list(set_hla)
    nc = len(list_hla_l)
    for y in range(len(set_hla)):
        if 'HLA-A' not in set_hla[y] and 'HLA-B' not in set_hla[y] and 'HLA-C' not in set_hla[y]:
            print('HLA allele format not admitted')
            exit()

# Plots material
#matplotlib.use('Agg')
from matplotlib.colors import LogNorm
import matplotlib as mpl
#from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import patches
from pandas.plotting import table
mpl.rcParams['font.family'] = ['Garuda']
mpl.rcParams['font.serif'] = ['Garuda-Oblique']

# Plots options
off = 0.15
s1 = 14
s2 = 16
s3 = 2
s4 =12
colors = ['red','orange', 'deepskyblue', 'blue', 'springgreen', 'forestgreen']

# definitions of functions 

def perf_measure(y_actual, y_hat): 
    tp = 0
    fp = 0
    tn = 0
    fn = 0

    for i in range(len(y_hat)): 
        if y_actual[i]==y_hat[i]==1:
           tp += 1
        if y_hat[i]==1 and y_actual[i]!=y_hat[i]:
           fp += 1
        if y_actual[i]==y_hat[i]==0:
           tn += 1
        if y_hat[i]==0 and y_actual[i]!=y_hat[i]:
           fn += 1
        
    den = float(tp+fp+tn+fn)
    if den !=0:
        acc = float(tp+tn)/float(tp+fp+tn+fn)
    else:
        acc=0
    den = float(tp+fp)
    if den !=0:   
        prec= float(tp)/float(tp+fp)
    else:
        prec=0
    den = float(tn+fp)
    if den !=0:
        spec= float(tn)/float(tn+fp)
    else:
        spec=0
    den = float(tp+fn)   
    if den !=0:
        sens= float(tp)/float(tp+fn)
    else:
        sens=0

    return(acc, prec, spec, sens)


def convert_number(seqs): # convert to numbers already aligned seqs
    aa = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V',  'W', 'Y','-']
    aadict = {aa[k]: k for k in range(len(aa))} 
    
    msa_num = np.array(list(map(lambda x: [aadict[y] for y in x], seqs[0:])), dtype=curr_int, order="c")
    
    return msa_num

def convert_letter(seqs_n): # convert to numbers already aligned seqs
    aa = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V',  'W', 'Y','-']
    aadictinv = {k: aa[k] for k in range(len(aa))} 
    seqs=[]
    if type(seqs_n[0]) == curr_int:
        seqs.append(''.join([aadictinv[e] for e in seqs_n]))
    else:
        for t in range(len(seqs_n)):
            seqs.append(''.join([aadictinv[e] for e in seqs_n[t]]))
    return seqs

def overlap_seqs(list1,list2):
    overlap=[]
    for i in range(len(list1)):
        if list1[i] in list2:
            overlap.append(list1[i])
    return overlap

def flatten_list(listoflist):
    listoflist_fl = [];
    for l in range(len(listoflist)):
        for u in range(len(listoflist[l])):
            listoflist_fl.append(listoflist[l][u])
    return listoflist_fl


import random 
def Rand(start, end, num): 
    res = []  
    for j in range(num*100): 
        r = random.randint(start, end)
        if r not in res:
            res.append(r)  
        if len(res) == num:
            break
    return res

def partition(list_in, n):
    return [list_in[i::n] for i in range(n)]

curr_float = np.float32
curr_int = np.int16
def average_n(config,q):
    B = config.shape[0]
    N = config.shape[1]
    out = np.zeros((N,q) ,dtype = curr_float)
    for b in list(range(B)):
        for n in list(range(N)):
            out[n, config[b,n]] +=1
    out/=B
    return out


# keras and sklearn
import keras
from keras.models import Sequential
from keras.layers import Dense, Dropout
from keras import regularizers
import matlab.engine
import sklearn
from sklearn import metrics
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.cluster import AgglomerativeClustering
from sklearn.mixture import GaussianMixture
cov_type = ['spherical', 'diag', 'tied', 'full']

# tools for aminoacids
aa = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V',  'W', 'Y','-']
aadict_rev = {k:aa[k] for k in range(len(aa))}
aadict = {aa[k]: k for k in range(len(aa))} 

pseudocount = 0.000001

# import IEDB and read relevant columns
cond = yesreal == 1 and nc!=1
if yesreal == 0 or makereweighting or cond:
    filename_lab = rootf + '/mhc_ligand_full.csv'
    iedb = pd.read_csv(filename_lab, sep=',')
    head = iedb.columns # select relevant indices
    index_type = 10
    index_antigen = 11
    index_quality = 83
    index_species = 39
    hla = list(iedb['MHC'].values)
    antigen = list(iedb['Epitope.2'].values)
    quality = list(iedb[head[index_quality]].values)

    # define a given set of HLA types and of quality of binding assessment
    quality_pos = ['Positive-High', 'Positive', 'Positive-Intermediate', 'Positive-Low']

    # conditions of selection of MS data
    condition0 = iedb[head[index_type]] == 'Linear peptide'
    condition1 = iedb[head[index_species]].isin(['human (Homo sapiens)', 'Homo sapiens', 'Homo sapiens (human)']) 
    condition2 = iedb['MHC.3'] == 'I'
    condition3 = iedb[head[index_quality]].isin(quality_pos) 
    condition5 = iedb['Assay.1'].isin(['cellular MHC/mass spectrometry', 'mass spectrometry', 'secreted MHC/mass spectrometry'])
    condition_mon = iedb['MHC.2'] == 'Single allele present'
    condition_aps = iedb['MHC.2'] == 'Allele specific purification'

    # conditions of selection of binding assay data
    list_meth_ba = ['binding assay', 'cellular MHC/T cell inhibition','cellular MHC/competitive/fluorescence', 'cellular MHC/competitive/radioactivity','cellular MHC/direct/fluorescence', 'lysate MHC/direct/radioactivity','purified MHC/competitive/fluorescence','purified MHC/competitive/radioactivity', 'purified MHC/direct/fluorescence', 'purified MHC/direct/radioactivity']
    condition_ba1 = iedb['Assay.3'] == 'nM'
    condition_ba2 = iedb['Assay.1'].isin(list_meth_ba)
    c_ex11 = iedb['Assay.1'] == 'purified MHC/direct/radioactivity'
    c_ex12 = iedb['Assay.2'] == 'dissociation constant KD'
    c_ex21 = iedb['Assay.1'] == 'purified MHC/direct/fluorescence'
    c_ex22 = iedb['Assay.2'] == 'half maximal effective concentration (EC50)'
    c_ex31 = iedb['Assay.1'] == 'cellular MHC/direct/fluorescence'
    c_ex32 = iedb['Assay.2'] == 'half maximal effective concentration (EC50)'

    nminseq = 50 # Num min seqs before removing allele
    nminseqMS = 300 # Num min seqs before switching to all MS 
    nminseqAL = 300 # Num min seqs before switching to all techs
    NSmin = 10

    if yesreal:
        nminseq = 5 
        nminseqMS = 30
        nminseqAL = 30

    # Estimate Frequencies in different sets for reweighting factor: background frequency is obtained following Bassani-Sternberg et al. PLoS Comput. Biol. 2017

    if makereweighting:
        condition4 = iedb['MHC'].isin(['HLA-A*02:01', 'HLA-A*02:05', 'HLA-A*02:06', 'HLA-A*02:20', 'HLA-A*25:01', 'HLA-A*26:01','HLA-A*29:02', 'HLA-B*08:01', 'HLA-B*14:01', 'HLA-B*14:02', 'HLA-C*03:03', 'HLA-C*07:04'])
        # exclude alleles that have specificity at P4-P7
        # select ms
        iedb_selms = iedb[condition0 & condition1 & condition2 & condition3 & condition5 & ~condition4]
        # select nms
        iedb_selnms = iedb[condition0 & condition1 & condition2 & condition3 & ~condition5]
        # Find the HLAs in common - to guarantee that one is correcting for the bias in the technique not in the composition
        hla_corr = overlap_seqs(list(np.unique(iedb_selms['MHC'].values)), list(np.unique(iedb_selnms['MHC'].values)))
        # narrow them down to the hlas that have at least 100 seqs in each subset - MS and NON-MS
        hla_corr_100ms = []
        for u in range(len(hla_corr)):
            iedb_sel0 = iedb_selms[iedb_selms['MHC'] == hla_corr[u]]
            if len((iedb_sel0.values)[:,index_antigen]) >= 100:
                hla_corr_100ms.append(hla_corr[u])

        hla_corr_100nms = []
        for u in range(len(hla_corr)):
            iedb_sel0 = iedb_selnms[iedb_selnms['MHC'] == hla_corr[u]]
            if len((iedb_sel0.values)[:,index_antigen]) >= 100:
                hla_corr_100nms.append(hla_corr[u])

        # finally take overlap
        hla_corr_fin = overlap_seqs(hla_corr_100ms, hla_corr_100nms)
        iedb_ms = iedb_selms[iedb_selms['MHC'].isin(hla_corr_fin)]
        iedb_nms = iedb_selnms[iedb_selnms['MHC'].isin(hla_corr_fin)]

        seqs_ms = (iedb_ms.values)[:, index_antigen]
        seqs_msL = [a for a in seqs_ms if len(a)==9 and 'X' not in a and 'Z' not in a] 

        # 9-mers only, to have clear identification of non-anchor sites
        list_counts_aa = []
        for a in aa:
            countaa = 0
            totL = 0
            for ss in range(len(seqs_msL)): # Limit to non-anchor positions
                s = ''.join(list(seqs_msL[ss][3:7]))
                totL += len(s)
                if a in s:
                    countaa += s.count(a)
            list_counts_aa.append(countaa)
        aafreq_ms = [float(list_counts_aa[t])/float(totL) for t in range(len(aa)-1)]
        aafreq_ms.append(0.05) # add a uniform-distribution value for the gap

        if args.rwnms:
        # Find frequencies of non-MS data
            seqs_nms = (iedb_nms.values)[:, index_antigen]
            seqs_nmsL = [a for a in seqs_nms if len(a) == 9 and 'X' not in a and 'Z' not in a] 

            list_counts_aa = []
            for a in aa:
                countaa = 0
                totL = 0
                for ss in range(len(seqs_nmsL)): # Limit estimation to non-anchor positions
                    s = ''.join(list(seqs_nmsL[ss][3:7]))
                    totL += len(s)
                    if a in s:
                        countaa += s.count(a)
                list_counts_aa.append(countaa)

            aafreq_nms = [float(list_counts_aa[t])/(totL) for t in range(len(aa)-1)]

        if args.rwhp:
            aafreq_nms = [0.06973193229784201, 0.022227300211810606, 0.04776998270837931, 0.07077023506787744, 0.035831258011111006, 0.06569891873175196, 0.026058299821958418, 0.04291157865732822, 0.056774727153814875, 0.09904393617954255, 0.02195800091948729, 0.03549834564548375, 0.06341291522288485, 0.0477736944581562, 0.05668205186949188, 0.0842078354011136, 0.054770500734393386, 0.059924778626702334, 0.012472466417845177, 0.026161320621617003]

        aafreq_nms.append(0.05)

# Parameters re-alignment
mex = 1 #exp HMM score for weighting alignment score
yw = 1 #re-weights seqs in building profile
dth = 0 # threshold in decoder probability for re-weighting step

# Re-adjust options
if LL == 1:
    Nrun = 1    
    print('No need for re-alignment with 1 length')

if CC==20:
    Na = CC
if CC ==21:
    Na = CC-1

if makereweighting:
    strw = 'YesAlpha'
else:
    strw = 'NoAlpha' 
    
if maketraining == 0:
    Nrun = 1

if args.o is None:
    out_fold = rootf + '/data/IEDB_out_'+ str_hla # name folder  which allotype
else:
    out_fold = rootf +'/' + args.o

if os.path.exists(out_fold) is False:
    os.mkdir(out_fold)

if args.nameo is None:
   out_par = 'YesReal_' + str(yesreal) + '_hu_' + str(hu) + '_l12' + str(l12) + '_AL' + str(np.mean(range_len)) + '_SA' + str(SA) + '_RW' + strw + '_TR' + str(deg)
else:
   out_par = args.nameo


## name of a log file ##
log_file = out_fold + '/Log_' + out_par + '.txt' 


if nc == 1 and yesreal == 1:
    filename = out_fold + '/' + args.i + '.txt'
    iedb_data_fl = []
    with open(filename) as f:
        for line in f:
            linesplit = line.strip().split('\t')
            nogap1=linesplit[0].replace(' ','')
            if len(nogap1) in range_len:
                iedb_data_fl.append(nogap1)
                

else:
    f = open(log_file,'a+')
    f.write('Data retrieval from IEDB' + '\n') 
    f.close()
    
    iedb_data_T = [] # _T temporary variables
    iedb_cat_T = []
    iedb_qua_T = []
    size_nc =[]
    for lcou in range(len(list(np.unique(set_hla)))):
        if lcou ==0:
            l = lcou
        condition4 = iedb['MHC'] == list_hla_l[l]
        if args.ba == 1:

            iedb_hla_pos = iedb[~(c_ex11 & c_ex12) & ~(c_ex21 & c_ex22) & ~(c_ex31 & c_ex32) & condition0 & condition1 & condition2 & condition4 & condition_ba1 & condition_ba2]
            size_nc.append(len(iedb_hla_pos))
            an_chla = (iedb_hla_pos.values)[:,index_antigen]
            ba_chla = (iedb_hla_pos.values)[:,85]
            f = open(log_file,'a+')
            f.write('Total data ' + list_hla_l[l] + ' = ' + str(len(an_chla)) + '\n')
            f.close()
            if len(an_chla) < nminseqAL: # at least 100, otherwise appeal to other techniques
               iedb_hla_pos = iedb[condition0 & condition1 & condition2 & condition3 & condition4]
               an_chla = (iedb_hla_pos.values)[:,index_antigen]
               ba_chla = (iedb_hla_pos.values)[:,85]
               f = open(log_file,'a+')
               f.write('in MS less ' + str(nminseqAL) + ', increase with all techniques ' + str(len(an_chla))+ '\n')
               f.close()
        else:
            iedb_hla_pos = iedb[condition0 & condition1 & condition2 & condition3 & condition4 & condition5 & condition_mon]
          #  iedb_hla_pos = iedb[condition0 & condition1 & condition2 & condition3 & condition4 & condition5]
            size_nc.append(len(iedb_hla_pos))
            an_chla = (iedb_hla_pos.values)[:,index_antigen]
            f = open(log_file,'a+')
            f.write('Total data ' + list_hla_l[l] + ' = ' + str(len(an_chla)) + '\n')
            f.close()
            if len(an_chla) < nminseqMS:
                iedb_hla_pos0 = iedb[condition0 & condition1 & condition2 & condition3 & condition4 & condition5]
                tech1 = iedb_hla_pos0['MHC.2'] == 'Allele specific purification'
                tech2 = iedb_hla_pos0['MHC.2'] == 'Single allele present'
                iedb_hla_pos = iedb_hla_pos0[tech1 | tech2]
                an_chla = (iedb_hla_pos.values)[:,index_antigen]
                f = open(log_file,'a+')
                f.write('in monoallelic less ' + str(nminseqMS) + ', increase with all MS-allele specific purification ' + str(len(an_chla)) + '\n')
                f.close()
                if len(an_chla) < nminseqMS:
                    iedb_hla_pos = iedb[condition0 & condition1 & condition2 & condition3 & condition4 & condition5]
                    an_chla = (iedb_hla_pos.values)[:,index_antigen]
                    f = open(log_file,'a+')
                    f.write('in monoallelic less ' + str(nminseqMS) + ', increase with all MS ' + str(len(an_chla)) + '\n')
                    f.close()
                    if len(an_chla) < nminseqAL: # at least 100, otherwise appeal to other techniques
                        iedb_hla_pos = iedb[condition0 & condition1 & condition2 & condition3 & condition4]
                        an_chla = (iedb_hla_pos.values)[:,index_antigen]
                        f = open(log_file,'a+')
                        f.write('in MS less ' + str(nminseqAL) + ', increase with all techniques ' + str(len(an_chla))+ '\n')
                        f.close()

        for si in range(len(an_chla)):
            seq = an_chla[si]
            if '+' in seq:
                ii=seq.index('+')
                an_chla[si] =seq[0:ii-1]

        an_chla9 = [a for a in an_chla if 'X' not in a and 'Z' not in a and len(a) in range_len]
        if args.ba == 1:
            an_chla9 = [an_chla[a] for a in range(len(an_chla)) if 'X' not in an_chla[a] and 'Z' not in an_chla[a] and len(an_chla[a]) in range_len and float(ba_chla[a]) < 500]
        iedb_data_T.append(list(an_chla9))
        iedb_cat_T.append(list(map(int,list(l*np.ones(len(an_chla9))))))
        indices = []
        for u in range(len(an_chla9)):
            indices.append(list(an_chla).index(an_chla9[u]))
        iedb_qua_T.append(list((iedb_hla_pos.values)[indices,index_quality]))

        if len(an_chla9) < nminseq or len([a for a in an_chla9 if len(a) == SA]) < 1:
            print('overall less ' + str(nminseq) + ', remove allele ' + list_hla_l[l] + ' from deconvolution')
            f = open(log_file,'a+')
            f.write('overall less ' + str(nminseq) + ', remove allele ' + list_hla_l[l] + ' from deconvolution \n')
            f.close()
            list_hla_l.remove(list_hla_l[l])
            nc = len(list_hla_l)
            if len(iedb_qua_T) > l:
                iedb_qua_T.remove(iedb_qua_T[l])
            if len(iedb_data_T) > l:
                iedb_data_T.remove(iedb_data_T[l])
            if len(iedb_cat_T) > l:
                iedb_cat_T.remove(iedb_cat_T[l])
        else:
            l+=1

    # clean sample from overlap between classes to get uniquely labelled peps

    iedb_data = []
    iedb_cat = []
    iedb_qua = []
    iedb_ind = []
    iedb_data_val = []
    iedb_cat_val = []
    iedb_qua_val = []
    iedb_ind_V = []

    pwm_refA =[] # Pwm built only on all data

    count=0
    count_val=0
    for l in range(len(iedb_data_T)):
        listd=[]
        listc=[]
        listq=[]
        listi=[]

        listd_val=[]
        listc_val=[]
        listq_val=[]
        listi_val=[]

        lo = list(range(0, len(iedb_data_T)))
        lo.remove(l)
        listo = [iedb_data_T[t] for t in lo]
        list_other = flatten_list(listo)
        list_u = Rand(0,(len(iedb_data_T[l])-1), int(round(deg*len(iedb_data_T[l]))))
        for u in list_u:
            if iedb_data_T[l][u] not in listd and iedb_data_T[l][u] not in list_other:
           # if iedb_data_T[l][u] not in listd:
                listd.append(iedb_data_T[l][u])
                listc.append(iedb_cat_T[l][u])
                listq.append(iedb_qua_T[l][u])
                listi.append(count)
                count+=1

        for u in range(len(iedb_data_T[l])):
            if u not in list_u and iedb_data_T[l][u] not in listd and iedb_data_T[l][u] not in list_other:
                listd_val.append(iedb_data_T[l][u])
                listc_val.append(iedb_cat_T[l][u])
                listq_val.append(iedb_qua_T[l][u])
                listi_val.append(count_val)
                count_val+=1

        iedb_data.append(listd)
        f = open(log_file,'a+')
        f.write('Final data ' + list_hla_l[l] + ' = ' + str(len(listd)) + ' for average length ' + str(np.mean(range_len)) + '\n')
        f.close()
        if makefig: 
            if SA not in [len(s) for s in listd]:
                aux_str = ''.join(['A' for i in range(SA)])
                T=len(listd)
                listd.append(aux_str) 

                if LL == 1:
                    temp_nn = convert_number(listd)
                else: 
                    name_mat = rootf + '/Align_utils/align_seqpy.py'
                    name_seqs = rootf + '/Align_utils/seqs_str.txt'
                    with open(name_seqs, 'w') as out_f:
                        for u in range(len(listd)):
                            out_f.write(listd[u] + '\n')

                    subprocess.call('python3 -W ignore ' + name_mat +' -ss ' + name_seqs + ' -SA ' + str(SA) + ' -SAmin ' + str(SAmin) + ' -SAmax ' + str(SAmax), shell=True)
                    name_al = rootf + '/Align_utils/aligned_temp.txt'
                    temp_nn = np.loadtxt(name_al)
                    subprocess.call('rm ' + name_seqs, shell=True)
                    subprocess.call('rm ' + name_al, shell=True)

                an_chla9_n = np.copy(temp_nn[0:T])
            else:
                if LL == 1:
                    an_chla9_n = convert_number(listd)
                else: 
                    name_mat = rootf + '/Align_utils/align_seqpy.py'
                    name_seqs = rootf + '/Align_utils/seqs_str.txt'
                    with open(name_seqs, 'w') as out_f:
                        for u in range(len(listd)):
                            out_f.write(listd[u] + '\n')
                    subprocess.call('python3 -W ignore ' + name_mat +' -ss ' + name_seqs + ' -SA ' + str(SA) + ' -SAmin ' + str(SAmin) + ' -SAmax ' + str(SAmax), shell=True)
                    name_al = rootf + '/Align_utils/aligned_temp.txt'
                    an_chla9_n = np.loadtxt(name_al)
                    subprocess.call('rm ' + name_seqs, shell=True)
                    subprocess.call('rm ' + name_al, shell=True)

            pwm = average_n(an_chla9_n.astype(np.int16), CC)
            pwm_refA.append(pwm) # alignment allotype by allotype

            stringt  = list_hla_l[l] + ' # seqs ' + str(len(iedb_data[l]))
            fig = sequence_logo.Sequence_logo(pwm, figsize=(15,3), ylabel= 'bits', title = stringt, show=True, ticks_every=5, ticks_labels_size=20, title_size=24);

            pwm_ng = average_n(an_chla9_n.astype(np.int16), CC)
            fig = sequence_logo.Sequence_logo(pwm_ng, figsize=(15,3), ylabel= 'bits', title = stringt, show=True, ticks_every=5, ticks_labels_size=20, title_size=24);
            pathfig = out_fold + '/IEDB_allotype_' + str(l) + '_AL' + str(np.mean(range_len)) + '_SA' + str(SA) + '.pdf'
            fig.savefig(pathfig)
        plt.close('all')

        iedb_cat.append(listc)
        iedb_qua.append(listq)
        iedb_ind.append(listi)
        iedb_data_val.append(listd_val)
        iedb_cat_val.append(listc_val)
        iedb_qua_val.append(listq_val)
        iedb_ind_V.append(listi_val)

    pwm_random = float(1)/float(CC)*np.ones((SA,CC))

    iedb_data_fl = flatten_list(iedb_data)
    iedb_cat_fl = flatten_list(iedb_cat)
    iedb_qua_fl = flatten_list(iedb_qua)

    if deg<1:
        iedb_data_fl_val = flatten_list(iedb_data_val)
        iedb_cat_fl_val = flatten_list(iedb_cat_val)
        iedb_qua_fl_val = flatten_list(iedb_qua_val)
        with open(out_fold + '/validation.txt', 'w') as out_f:
            for u in range(len(iedb_data_fl_val)):
                out_f.write(iedb_data_fl_val[u] + '\n')

    # select 10% of data, those with high confidence for the supervised part
    if yesreal:
        filename = out_fold + '/' + args.i + '.txt'
        ms_seqs = []
        with open(filename) as f:
            for line in f:
                linesplit = line.strip().split('\t')
                nogap1=linesplit[0].replace(' ','')
                if len(nogap1) in range_len:
                    ms_seqs.append(nogap1)

        iedb_ind_tr=[]
        for l in range(len(iedb_ind)):
            listi = []
            ncli = float(1)/nc
            fullperc= (float(ncli))*perc
            limit = round(len(ms_seqs)*fullperc) # amount if labelled data == perc*size dataset
            if len(iedb_ind[l])>limit:
                for p in range(len(iedb_ind[l])):
                    cond_qua = (iedb_qua_fl[iedb_ind[l][p]] == quality_pos[0] or iedb_qua_fl[iedb_ind[l][p]] == quality_pos[1])
                    if iedb_data_fl[iedb_ind[l][p]] in ms_seqs and cond_qua:
                        listi.append(iedb_ind[l][p])
                    if len(list(np.unique(listi))) > limit:
                            break
                if len(list(np.unique(listi))) < limit:
                    for t in range(10000):
                        r = np.copy(random.randint(0,len(iedb_ind[l])-1))
                        if iedb_qua_fl[iedb_ind[l][r]] == quality_pos[0] or iedb_qua_fl[iedb_ind[l][r]] == quality_pos[1]:
                            listi.append(iedb_ind[l][r])
                        if len(list(np.unique(listi))) > limit:
                            break
            else:
                limit = len(iedb_ind[l])-2
                for p in range(len(iedb_ind[l])):
                    cond_qua = (iedb_qua_fl[iedb_ind[l][p]] == quality_pos[0] or iedb_qua_fl[iedb_ind[l][p]] == quality_pos[1])
                    if iedb_data_fl[iedb_ind[l][p]] in ms_seqs and cond_qua:
                        listi.append(iedb_ind[l][p])
                    if len(list(np.unique(listi))) > limit:
                            break
                if len(list(np.unique(listi))) < limit:
                    for t in range(10000):
                        r = np.copy(random.randint(0,len(iedb_ind[l])-1))
                        listi.append(iedb_ind[l][r])
                        if len(list(np.unique(listi))) > limit:
                            break
            iedb_ind_tr.append(list(np.unique(listi)))
                     
    else: 
        iedb_ind_tr=[]
        for l in range(len(iedb_ind)):
            listi = []
            ncli = float(1)/nc
            fullperc = (float(ncli))*perc
            limit = round(len(iedb_data_fl)*fullperc)
            #limit = round(len(iedb_ind[l])*perc)
            if limit > NSmin:
                for t in range(10000):
                    r = np.copy(random.randint(0,len(iedb_ind[l])-1))
                    if iedb_qua_fl[iedb_ind[l][r]] == quality_pos[0] or iedb_qua_fl[iedb_ind[l][r]] == quality_pos[1]:
                        listi.append(iedb_ind[l][r])
                    if len(list(np.unique(listi))) > limit:
                        break
                iedb_ind_tr.append(list(np.unique(listi)))
            else:
                f = open(log_file,'a+')
                f.write('Exception case')
                f.close()
                limit = len(iedb_ind[l])-2
                for t in range(10000):
                    r = np.copy(random.randint(0,len(iedb_ind[l])-1))
                    listi.append(iedb_ind[l][r])
                    if len(list(np.unique(listi))) > limit:
                        break
                iedb_ind_tr.append(list(np.unique(listi)))
        
    iedb_ind_val = [[a for a in iedb_ind[l] if a not in iedb_ind_tr[l]] for l in range(len(iedb_ind))]

    if yesreal:     

        full_data = []
        full_ind_tr = []
        full_ind_val = []
        full_cat = []
        len_tr = 0
        for c in range(nc):
              full_data.append([iedb_data_fl[iedb_ind_tr[c][u]] for u in range(len(iedb_ind_tr[c]))])
              full_cat.append([iedb_cat_fl[iedb_ind_tr[c][u]] for u in range(len(iedb_ind_tr[c]))])                    
              
        for c in range(nc):    
              full_ind_tr.append([u for u in range(len_tr, len_tr + len(iedb_ind_tr[c]))])
              len_tr += len(iedb_ind_tr[c]) 
        ms_seqadd=[ms_seqs[ns] for ns in range(len(ms_seqs)) if ms_seqs[ns] not in flatten_list(full_data)]
        full_data.append(ms_seqadd)

        len_val = len(flatten_list(full_data)) - len_tr
        full_cat.append([0]*len_val)
        full_ind_val = partition(list(range(len_tr, len_tr + len_val)), nc)
        iedb_data_fl = flatten_list(full_data)
        iedb_cat_fl = flatten_list(full_cat)
        iedb_ind_tr = list(full_ind_tr)
        iedb_ind_val = list(full_ind_val)
        
        f = open(log_file,'a+')
        for cl in range(len(iedb_ind_tr)):
            f.write('Amount of labelled peptides for ' + list_hla_l[cl] + ' = ' + str(len(iedb_ind_tr[cl])) + '\n')
        f.write('Length full sample ' + str(len(iedb_data_fl)) + '\n')
        f.close()

        iedb_ind=[]
        for h in range(nc):
            iedb_ind.append(iedb_ind_tr[h] + iedb_ind_val[h])

    # Build PWMs from 10% of data only
    pwm_refD = [] # Pwm built only on 10% of data
    for cl in range(nc):
        temp_msa = [iedb_data_fl[i] for i in iedb_ind_tr[cl]]
        if SA not in [len(s) for s in temp_msa]:
              aux_str = ''.join(['A' for i in range(SA)])
              T = len(temp_msa)
              temp_msa.append(aux_str)
              if LL == 1:
                temp_nn = convert_number(temp_msa)
              else: 
                name_mat = rootf + '/Align_utils/align_seqpy.py'
                name_seqs = rootf + '/Align_utils/seqs_str.txt'
                with open(name_seqs, 'w') as out_f:
                    for u in range(len(temp_msa)):
                        out_f.write(temp_msa[u] + '\n')    
                                          
                subprocess.call('python3 -W ignore ' + name_mat +' -ss ' + name_seqs + ' -SA ' + str(SA) + ' -SAmin ' + str(SAmin) + ' -SAmax ' + str(SAmax), shell=True)
                name_al = rootf + '/Align_utils/aligned_temp.txt'
                temp_nn = np.loadtxt(name_al)
                subprocess.call('rm ' + name_seqs, shell=True)
                subprocess.call('rm ' + name_al, shell=True)
              seqs_n = np.copy(temp_nn[0:T])
        else:
            if LL == 1:
                seqs_n = convert_number(temp_msa)
            else: 
                name_mat = rootf + '/Align_utils/align_seqpy.py'
                name_seqs = rootf + '/Align_utils/seqs_str.txt'
                with open(name_seqs, 'w') as out_f:
                    for u in range(len(temp_msa)):
                        out_f.write(temp_msa[u] + '\n') 
                subprocess.call('python3 -W ignore ' + name_mat +' -ss ' + name_seqs + ' -SA ' + str(SA) + ' -SAmin ' + str(SAmin) + ' -SAmax ' + str(SAmax), shell=True)
                name_al = rootf + '/Align_utils/aligned_temp.txt'
                seqs_n = np.loadtxt(name_al)
                subprocess.call('rm ' + name_seqs, shell=True)
                subprocess.call('rm ' + name_al, shell=True)
        pwm = average_n(seqs_n.astype(np.int16), CC)
        pwm_refD.append(pwm) 

    # Set the PWM to use as reference for the labelling of clusters
    pwm_ref = list(pwm_refD)
    
    '''
    if yesreal:
        for cl in range(nc):
            name_fi = out_fold + '/labelled_seqs_'+ out_par + '_' + str(cl)
            with open(name_fi  + '.txt', 'w') as out_f:
                for r in range(len(iedb_ind_tr[cl])):
                    out_f.write(iedb_data_fl[iedb_ind_tr[cl][r]] + '\n')
    '''
   

pwm_random = float(1)/float(CC)*np.ones((SA,CC))

# Train a RBM: datatset
temp_msa = list(iedb_data_fl)

# alignment of the dataset:
if SA not in [len(s) for s in temp_msa]:
      aux_str = ''.join(['A' for i in range(SA)])
      T = len(temp_msa)
      temp_msa.append(aux_str)
      if LL == 1:
            temp_nn = convert_number(temp_msa)
      else: 
            name_mat = rootf + '/Align_utils/align_seqpy.py'
            name_seqs = rootf + '/Align_utils/seqs_str.txt'
            with open(name_seqs, 'w') as out_f:
              for u in range(len(temp_msa)):
                  out_f.write(temp_msa[u] + '\n')
            subprocess.call('python3 -W ignore ' + name_mat +' -ss ' + name_seqs + ' -SA ' + str(SA) + ' -SAmin ' + str(SAmin) + ' -SAmax ' + str(SAmax), shell=True)
            name_al = rootf + '/Align_utils/aligned_temp.txt'
            temp_nn = np.loadtxt(name_al)
            subprocess.call('rm ' + name_seqs, shell=True)
            subprocess.call('rm ' + name_al, shell=True)
      iedb_data_fl_nall = np.copy(temp_nn[0:T])
else:
      if LL == 1:
            iedb_data_fl_nall = convert_number(iedb_data_fl)
      else: 
            name_mat = rootf + '/Align_utils/align_seqpy.py'
            name_seqs = rootf + '/Align_utils/seqs_str.txt'
            with open(name_seqs, 'w') as out_f:
             for u in range(len(iedb_data_fl)):
                 out_f.write(iedb_data_fl[u] + '\n')
            subprocess.call('python3 -W ignore ' + name_mat +' -ss ' + name_seqs + ' -SA ' + str(SA) + ' -SAmin ' + str(SAmin) + ' -SAmax ' + str(SAmax), shell=True)
            name_al = rootf + '/Align_utils/aligned_temp.txt'
            iedb_data_fl_nall = np.loadtxt(name_al)
            subprocess.call('rm ' + name_seqs, shell=True)
            subprocess.call('rm ' + name_al, shell=True)

iedb_data_fl_ntemp = np.copy(iedb_data_fl_nall).astype(np.int16) # temporary alignment for iterative training

for runs in range(Nrun):

    B = len(iedb_data_fl_ntemp)
    all_data = np.copy(iedb_data_fl_ntemp)

    if makereweighting:
        pc_ms = 0.0015 # pseudo-count adapted on cysteine frequency to avoid too large weights \alpha
        if args.rwhp:
            pc_nms = pc
        if args.rwnms:
            pc_nms = 0.0015
        alphas_rbm=[]
        for m in range(len(all_data)):
            alpha = 1
            for t in range(len(all_data[m])):
                alpha *= float(aafreq_nms[all_data[m][t]]+ pc_nms)/float(aafreq_ms[all_data[m][t]] + pc_ms)
                alphas_rbm.append(alpha) 

        all_weights = np.copy(alphas_rbm)
    else:
        all_weights =np.ones(len(all_data))

    n_v = SA # Number of visible units; = # sites in alignment.
    n_h = hu # Number of hidden units.
    visible = 'Potts' # Nature of visible units potential. Here, Potts states...
    n_cv = CC # With n_cv = 21 colors (all possible amino acids + gap)
    name_r = out_fold + '/model_' + out_par + '.data'

    hidden = 'dReLU' # Nature of hidden units potential. Here, dReLU potential. 'Gaussian'
    seed = 0 # Random seed (optional)
    all_data = all_data.astype(np.int16)
    seed = utilities.check_random_state(0) # shuffle data

    batch_size = 100 # Size of mini-batches (and number of Markov chains used). Default: 100. 
    learning_rate = 0.05 # Initial learning rate (default: 0.1)
    decay_after = 0.5 # Decay learning rate after 50% of iterations (default: 0.5)
    l1b = l12 # L1b regularization. Default : 0.
    N_MC = 15 # Number of Monte Carlo steps between each update
    l2f = 1/len(all_data)

    if maketraining:
        RBM = rbm.RBM(visible = visible, hidden = hidden,n_v = n_v,n_h = n_h, n_cv = n_cv, random_state = seed)
        RBM.fit(all_data.astype(np.int16), weights = all_weights.astype('float32'), batch_size = batch_size, n_iter = n_iter, l1b = l1b, l2_fields = l2f, N_MC = N_MC, decay_after = decay_after, verbose = 1)
        RBM_utils.saveRBM(name_r,RBM)
    else:
        RBM = RBM_utils.loadRBM(name_r)

    if nc == 1: 
        if deg<1:
            name_wes = out_fold + '/validation.txt'
    
            if os.path.exists(name_wes) is False:
                print('File of sequences to score not found, thus exiting ')
                exit()
            wespepss=[]
            with open(name_wes) as f:
                for line in f:
                    linesplit = line.strip().split('\n')
                    wespepss.append(linesplit[0])

            out_fold_NA = out_fold + '/scoring'
            if os.path.exists(out_fold_NA) is False:
                os.mkdir(out_fold_NA)

            wespepswt = [a for a in wespepss if len(a) in range_len]
            if LL == 1:
                wespeps_n = convert_number(wespepswt)
            else:
                name_mat = rootf + '/Align_utils/align_to_seedpy.py'
                name_seqs = rootf + '/Align_utils/seqs_str.txt'
                with open(name_seqs, 'w') as out_f:
                    for u in range(len(wespepswt)):
                        out_f.write(wespepswt[u] + '\n')

                name_seed = rootf + '/Align_utils/seqs_seed.txt'
                with open(name_seed, 'w') as out_f:
                    for u in range(len(iedb_data_fl)):
                        out_f.write(iedb_data_fl[u] + '\n')

                subprocess.call('python3 -W ignore ' + name_mat + ' -sseed ' + name_seed + ' -sseqs ' + name_seqs +' -SA ' + str(SA) + ' -SAmin ' + str(SAmin) + ' -SAmax ' + str(SAmax) + ' -yw 0 ', shell = True)
                name_al = rootf + '/Align_utils/aligned_temp.txt'
                wespeps_n = np.loadtxt(name_al)
                subprocess.call('rm ' + name_al, shell = True)
                subprocess.call('rm ' + name_seed, shell = True)
                subprocess.call('rm ' + name_seqs, shell = True)

            likwes = RBM.likelihood(wespeps_n.astype(np.int16))
            wespeps_g = convert_letter(wespeps_n.astype(np.int16))
            lenwes = [len(a) for a in wespepswt]

            NAlist = zip(wespepswt, wespeps_g, likwes, lenwes, list(np.repeat(list_hla_l[0],len(likwes))) , likwes)
            NAstr = 'Sequences' + '\t'+ 'Seqs core' + '\t' + 'RBM Log-Lik.' + '\t' + 'Length' +   '\t' + 'Pred. HLA class' + '\t'  + 'Total RBM-MHC score'

            with open(out_fold_NA + '/validation_scored_' +  out_par + '.txt', 'w') as out_f:
                out_f.write(NAstr + '\n')
                writer = csv.writer(out_f, delimiter='\t')
                writer.writerows(NAlist)

        if NAranking:

            name_wes = out_fold + '/' + args.score + '.txt'        
            if os.path.exists(name_wes) is False:
                print('File of sequences to score not found, thus exiting ')
                exit()
            wespepss=[]
            with open(name_wes) as f:
                for line in f:
                    linesplit = line.strip().split('\n')
                    wespepss.append(linesplit[0])

            out_fold_NA = out_fold + '/scoring'
            if os.path.exists(out_fold_NA) is False:
                os.mkdir(out_fold_NA)

            wespepswt = [a for a in wespepss if len(a) in range_len]
            if LL == 1:
                wespeps_n = convert_number(wespepswt)
            else:
                name_mat = rootf + '/Align_utils/align_to_seedpy.py'
                name_seqs = rootf + '/Align_utils/seqs_str.txt'
                with open(name_seqs, 'w') as out_f:
                    for u in range(len(wespepswt)):
                        out_f.write(wespepswt[u] + '\n')

                name_seed = rootf + '/Align_utils/seqs_seed.txt'
                with open(name_seed, 'w') as out_f:
                    for u in range(len(iedb_data_fl)):
                        out_f.write(iedb_data_fl[u] + '\n')

                subprocess.call('python3 -W ignore ' + name_mat + ' -sseed ' + name_seed + ' -sseqs ' + name_seqs +' -SA ' + str(SA) + ' -SAmin ' + str(SAmin) + ' -SAmax ' + str(SAmax) + ' -yw 0 ', shell = True)
                name_al = rootf + '/Align_utils/aligned_temp.txt'
                wespeps_n = np.loadtxt(name_al)
                subprocess.call('rm ' + name_al, shell = True)
                subprocess.call('rm ' + name_seed, shell = True)
                subprocess.call('rm ' + name_seqs, shell = True)

            likwes = RBM.likelihood(wespeps_n.astype(np.int16))
            wespeps_g = convert_letter(wespeps_n.astype(np.int16))
            lenwes = [len(a) for a in wespepswt]

            NAlist = zip(wespepswt, wespeps_g, likwes, lenwes)
            NAstr = 'Sequences' + '\t'+ 'Seqs core' + '\t' + 'RBM Log-Lik.' + '\t' + 'Length'

            with open(out_fold_NA + '/' + args.score + '_scored_' +  out_par + '.txt', 'w') as out_f:
                out_f.write(NAstr + '\n')
                writer = csv.writer(out_f, delimiter='\t')
                writer.writerows(NAlist)

        if makegeneration == 1:
            print('Sequence HLA-specific generation is implemented for multiallelic samples')
        quit()


    list_categories = np.zeros((len(iedb_data_fl),len(list_hla_l)))
    for m in range(len(iedb_data_fl)):
        list_categories[m, iedb_cat_fl[m]] = 1

    I = RBM.vlayer.compute_output(iedb_data_fl_ntemp.astype(np.int16), RBM.weights)
    if runs == 0:
        Inall = np.copy(I)

    datade = np.vstack((I[iedb_ind_tr[0]], I[iedb_ind_tr[1]]))
    for l in list(range(2,len(iedb_ind_tr))):    
        datade = np.vstack((datade, I[iedb_ind_tr[l]]))

    labelsde = np.vstack((list_categories[iedb_ind_tr[0]], list_categories[iedb_ind_tr[1]]))
    for l in list(range(2,len(iedb_ind_tr))):
        labelsde = np.vstack((labelsde, list_categories[iedb_ind_tr[l]]))

    seed = utilities.check_random_state(0) # shuffle data
    permutation = np.argsort(seed.rand(labelsde.shape[0]))
    datade = datade[permutation] # Shuffle data.
    labelsde = labelsde[permutation] # Shuffle data.

    model = Sequential([Dense(len(iedb_ind), input_shape=(hu,), activation='softmax')])
    model.compile(optimizer='Adam', loss='categorical_crossentropy', metrics=['categorical_accuracy'])

    fac=1.5
    val_data = datade[int(round(len(datade)/fac)):]
    val_label = labelsde[int(round(len(labelsde)/fac)):]
    train_data = datade[:int(round(len(datade)/fac))]
    train_label = labelsde[:int(round(len(labelsde)/fac))]
    name_w = out_fold + '/weights_' + out_par + '_perc' + str(perc) + '.h5'
    if maketraining:
        best = keras.callbacks.ModelCheckpoint(name_w, monitor='val_categorical_accuracy', verbose=0, save_best_only=True, save_weights_only = True, mode='auto', period=1)
        model.fit(train_data, train_label, epochs =  n_epochs, batch_size=64, shuffle=True, callbacks = [best], validation_data=(val_data, val_label))
        
    model.load_weights(name_w)

    # Put forward new labels
    labels_new = [] # look at the new labels after the classification
    values_new = []
    pred_new = model.predict(I)
    for t in range(len(pred_new)):
        labels_new.append(np.argmax(pred_new[t])) 
        values_new.append(np.max(pred_new[t])) 

    labels_ssrbm = list(labels_new)
    for i in range(len(labels_ssrbm)): 
        if i in flatten_list(iedb_ind_tr):
             labels_ssrbm[i] = iedb_cat_fl[i]

    
    if not(yesreal):

      # Precision on validation set 
        mat_auc_list = []
        mat_aucR_list = []
        mat_acc_list = []
        mat_sens_list = []
        mat_spec_list = []
        mat_prec_list = []
        for l1 in range(nc):
            listtemp = list(range(nc))
            listtemp.remove(l1)
            for l2 in listtemp:   
                
                labels = np.hstack((np.ones((len(iedb_ind_val[l1]))), np.zeros((len(iedb_ind_val[l2]))) )) 
                pred1 = (np.array(labels_new)[iedb_ind_val[l1]] == l1).astype(np.int16)
                pred2 = (np.array(labels_new)[iedb_ind_val[l2]] == l1).astype(np.int16)
                scores = np.hstack((pred1, pred2))
                fpr, tpr, thresholds = metrics.roc_curve(labels, scores)
                mat_auc_list.append(metrics.auc(fpr, tpr))    
                (acc,prec,spec,sens) = perf_measure(labels, scores) 
                mat_acc_list.append(acc)
                mat_sens_list.append(sens)
                mat_spec_list.append(spec)
                mat_prec_list.append(prec)

                pred1 = model.predict(I[iedb_ind_val[l1]])
                pred2 = model.predict(I[iedb_ind_val[l2]])
                scores = np.hstack((pred1[:,l1], pred2[:,l1]))
                fpr, tpr, thresholds = metrics.roc_curve(labels, scores)
                mat_aucR_list.append(metrics.auc(fpr, tpr))

        aucv_rbm = np.mean(mat_auc_list)    
        aucRv_rbm = np.mean(mat_aucR_list)  
        accv_rbm = np.mean(mat_acc_list)
        precv_rbm =np.mean(mat_prec_list)
        sensv_rbm =np.mean(mat_sens_list) # this is what I want to increase
        specv_rbm =np.mean(mat_spec_list)

        f = open(log_file,'a+')
        perf = ' AUC ' + str(aucv_rbm)+ ' AUC Rank ' +str(aucRv_rbm) + ' Accuracy ' +str(accv_rbm) + ' Precision ' + str(precv_rbm) + ' Sensitivity ' + str(sensv_rbm)  + ' Specificity ' + str(specv_rbm) + '\n'
        f.write('Unlabelled only, iteration ' + str(runs) + ' RBM \n') 
        f.write(perf)
        f.close()

        # Precision on full dataset #
        mat_aucR = np.zeros((nc,nc))
        mat_auc_list=[]
        mat_aucR_list =[]
        mat_acc_list = []
        mat_sens_list = []
        mat_spec_list = []
        mat_prec_list = []
        for l1 in range(nc):
            listtemp = list(range(nc))
            listtemp.remove(l1)
            for l2 in listtemp:   
                labels = np.hstack((np.ones((len(iedb_ind[l1]))), np.zeros((len(iedb_ind[l2]))) )) 
                pred1 = np.ones((len(iedb_ind_tr[l1])))
                pred1 = np.hstack((pred1,(np.array(labels_new)[iedb_ind_val[l1]] == l1).astype(np.int16)))
                pred2 = np.zeros((len(iedb_ind_tr[l2])))
                pred2 = np.hstack((pred2, (np.array(labels_new)[iedb_ind_val[l2]] == l1).astype(np.int16) ))
                scores = np.hstack((pred1, pred2))
                fpr, tpr, thresholds = metrics.roc_curve(labels, scores)
                mat_auc_list.append(metrics.auc(fpr, tpr))    
                (acc,prec,spec,sens) = perf_measure(labels,scores) 
                mat_acc_list.append(acc)
                mat_sens_list.append(sens)
                mat_spec_list.append(spec)
                mat_prec_list.append(prec)

                pred1 = model.predict(I[iedb_ind_val[l1]])
                pred1 = np.vstack((pred1,list_categories[iedb_ind_tr[l1]]))
                pred2 = model.predict(I[iedb_ind_val[l2]])
                pred2 = np.vstack((pred2,list_categories[iedb_ind_tr[l2]]))
                scores = np.hstack((pred1[:,l1], pred2[:,l1]))
                fpr, tpr, thresholds = metrics.roc_curve(labels, scores)
                mat_aucR_list.append(metrics.auc(fpr, tpr))
                mat_aucR[l1,l2] = metrics.auc(fpr, tpr) 

        auc_rbm = np.mean(mat_auc_list)    
        aucR_rbm = np.mean(mat_aucR_list)  
        acc_rbm =np.mean(mat_acc_list)
        prec_rbm =np.mean(mat_prec_list)
        sens_rbm =np.mean(mat_sens_list)
        spec_rbm =np.mean(mat_spec_list)
        f = open(log_file,'a+')
        perf = ' AUC ' + str(auc_rbm)+ ' AUC Rank ' + str(aucR_rbm) + ' Accuracy ' +str(acc_rbm) + ' Precision ' + str(prec_rbm) + ' Sensitivy ' + str(sens_rbm)  + ' Specificity ' + str(spec_rbm) + '\n'
        f.write('Full set, iteration ' + str(runs) + ' RBM \n')
        f.write(perf)
        f.close()

      # re-align and re-run     
        if Nrun > 1 and nc > 1:
                
            seqs_seed_all = []
            seqs_new_all_n = []
            scores_new_all = []
            labels_newv = [labels_new[ff] for ff in flatten_list(iedb_ind_val)]
            values_newv = [values_new[ff] for ff in flatten_list(iedb_ind_val)]
            seqs_new = [iedb_data_fl[ff] for ff in flatten_list(iedb_ind_val)]
            for cl in range(nc):
                cc = np.where(np.array(labels_newv) == cl)[0] # For each class, align them to the corresponding seed
                seqs_val_new = []
                weights_val_new = []
                for u in range(len(cc)):
                    if values_newv[cc[u]] > dth:
                        seqs_val_new.append(iedb_data_fl[flatten_list(iedb_ind_val)[cc[u]]])
                        weights_val_new.append(values_newv[cc[u]])          
                seqs_seedt = [iedb_data_fl[u] for u in iedb_ind_tr[cl]]
                seqs_seed = list(seqs_seedt + seqs_val_new)
                weightsw = list([1]*len(iedb_ind_tr[cl]) + weights_val_new)
                namew = rootf + '/Align_utils/weights.txt'
                
                with open(namew, 'w') as out_f:
                    for r in range(len(weightsw)):
                        out_f.write(str(weightsw[r]) + '\n')
                if SA not in [len(s) for s in seqs_seedt]:
                    aux_str = ''.join(['A' for i in range(SA)])
                    T = len(seqs_seed)
                    seqs_seed.append(aux_str)
                    if LL == 1:
                        temp_nn = convert_number(seqs_seedt)
                    else: 
                        name_mat = rootf + '/Align_utils/align_seqpy.py'
                        name_seqs = rootf + '/Align_utils/seqs_str.txt'
                        with open(name_seqs, 'w') as out_f:
                            for u in range(len(seqs_seedt)):
                                out_f.write(seqs_seedt[u] + '\n')  
                        subprocess.call('python3 -W ignore ' + name_mat +' -ss ' + name_seqs + ' -SA ' + str(SA) + ' -SAmin ' + str(SAmin) + ' -SAmax ' + str(SAmax), shell=True)
                        name_al = rootf + '/Align_utils/aligned_temp.txt'
                        temp_nn = np.loadtxt(name_al)
                        subprocess.call('rm ' + name_seqs, shell=True)
                        subprocess.call('rm ' + name_al, shell = True)
                    seqs_seed_n = np.copy(temp_nn[0:T])
                    seqs_seedt = seqs_seedt[0:T]
                else:
                    if LL == 1:
                        seqs_seed_n = convert_number(seqs_seedt)
                    else: 
                        name_mat = rootf + '/Align_utils/align_seqpy.py'
                        name_seqs = rootf + '/Align_utils/seqs_str.txt'
                        with open(name_seqs, 'w') as out_f:
                            for u in range(len(seqs_seedt)):
                                out_f.write(seqs_seedt[u] + '\n')
                        subprocess.call('python3 -W ignore ' + name_mat + ' -ss ' + name_seqs + ' -SA ' + str(SA) + ' -SAmin ' + str(SAmin) + ' -SAmax ' + str(SAmax), shell=True)
                        name_al = rootf + '/Align_utils/aligned_temp.txt'
                        seqs_seed_n = np.loadtxt(name_al)
                        subprocess.call('rm ' + name_seqs, shell = True)
                        subprocess.call('rm ' + name_al, shell = True)
                
                name_mat = rootf + '/Align_utils/align_to_seedpy.py'
                name_seed = rootf + '/Align_utils/seqs_seed.txt'
                with open(name_seed, 'w') as out_f:
                    for u in range(len(seqs_seed)):
                        out_f.write(seqs_seed[u] + '\n')
                name_seqs = rootf + '/Align_utils/seqs_str.txt'
                with open(name_seqs, 'w') as out_f:
                    for u in range(len(seqs_new)):
                        out_f.write(seqs_new[u] + '\n')                      
                subprocess.call('python3 -W ignore ' + name_mat + ' -sseed ' + name_seed + ' -sseqs ' + name_seqs + ' -SA ' + str(SA) + ' -SAmin ' + str(SAmin) + ' -SAmax ' + str(SAmax) + ' -yw ' + str(yw), shell = True)
                name_al = rootf + '/Align_utils/aligned_temp.txt'
                seqs_new_n = np.loadtxt(name_al)
                subprocess.call('rm ' + name_seed, shell = True)
                subprocess.call('rm ' + name_seqs, shell = True)
                subprocess.call('rm ' + name_al, shell = True)   

                scores_new = list(np.loadtxt(rootf + '/Align_utils/scores.txt'))
                seqs_seed_all.append(seqs_seedt)
                if cl==0:
                    seqs_seed_all_n = np.copy(seqs_seed_n)
                else:
                    seqs_seed_all_n = np.vstack((seqs_seed_all_n, seqs_seed_n))

                seqs_new_all_n.append(seqs_new_n)
                scores_new_all.append(scores_new)

            seqs_seed_all = flatten_list(seqs_seed_all)

            iedb_data_fl_repl = []
            scores_repl = []
            for cl in range(nc):
                iedb_data_fl_al_new = np.zeros((len(iedb_data_fl), SA))
                scores_temp = np.zeros((len(iedb_data_fl)))
                for g in range(len(seqs_seed_all)):
                    ii = iedb_data_fl.index(seqs_seed_all[g])
                    iedb_data_fl_al_new[ii,:] = np.copy(seqs_seed_all_n[g])
                    scores_temp[ii] = 1
                for h in range(len(seqs_new)):
                    jj = iedb_data_fl.index(seqs_new[h])
                    iedb_data_fl_al_new[jj,:] = np.copy(seqs_new_all_n[cl][h])
                    scores_temp[jj] = np.copy(scores_new_all[cl][h])
                if cl==0:
                    iedb_data_fl_ntemp = np.copy(iedb_data_fl_al_new.astype(np.int16)) # nc replicates of the system
                else:
                    iedb_data_fl_ntemp = np.vstack((iedb_data_fl_ntemp, iedb_data_fl_al_new.astype(np.int16)))
                iedb_data_fl_repl.append(iedb_data_fl_al_new.astype(np.int16))
                scores_repl.append(scores_temp)

            opt_score = []
            for u in range(len(iedb_data_fl)):
                scorem = [scores_repl[cl][u] for cl in range(nc)]
                opt_score.append(scorem)


            iedb_data_fl_ntemp = np.copy(iedb_data_fl_repl[0])
            for i in range(len(flatten_list(iedb_ind_val))):
                ii = flatten_list(iedb_ind_val)[i]
                if len((iedb_data_fl)) != 0:
                    vec_scores = [np.exp(opt_score[ii][cl])/sum(np.exp(opt_score[ii])) for cl in range(nc)]
                    vec_scores = [(opt_score[ii][cl]**mex) for cl in range(nc)]
                    maxvs = np.argmax(vec_scores)
                    iedb_data_fl_ntemp[ii,:] = np.copy(iedb_data_fl_repl[maxvs][ii, :].astype(np.int16))


    iedb_data_fl_nf = np.copy(iedb_data_fl_ntemp)     

    if not(yesreal):
        iedb_data_fl_g = convert_letter(iedb_data_fl_nf)
        OLlist = zip(iedb_data_fl, iedb_data_fl_g, iedb_cat_fl, labels_ssrbm)
        OLstr = 'Data' + '\t'+ 'Data aligned' + '\t' + 'IEDB HLA' + '\t'+ 'RBM-MHC HLA'
        with open(out_fold + '/Table_classification_' +  out_par + '.txt', 'w') as out_f:
            out_f.write(OLstr + '\n')
            writer = csv.writer(out_f, delimiter='\t')
            writer.writerows(OLlist)

    if yesreal:
        iedb_ind_val0 = list(iedb_ind_val)
        iedb_ind_val=[]
        labels_newv = [labels_new[ff] for ff in flatten_list(iedb_ind_val0)]
        for cl in range(nc):
            cc = np.where(np.array(labels_newv) == cl)[0]
            iedb_ind_val.append([flatten_list(iedb_ind_val0)[u] for u in cc])

        iedb_ind=[]
        for h in range(nc):
            iedb_ind.append(iedb_ind_tr[h] + iedb_ind_val[h])
       
        labels_rbm = []
        labels_tr = [0]*len(iedb_data_fl)
        for i in range(len(iedb_data_fl)):
            for cl in range(nc):
                if i in iedb_ind[cl]:
                    labels_rbm.append(cl)
                    if i in iedb_ind_tr[cl]:
                        labels_tr[i] = 1
 
        iedb_data_fl_g = convert_letter(iedb_data_fl_nf)
        CLlist = zip(iedb_data_fl,iedb_data_fl_g, labels_tr, labels_rbm)
        CLstr = 'Data' + '\t'+ 'Data Core' + '\t'+ 'Training' + '\t' + 'RBM-MHC label' 
        fold_real = out_fold + '/' + args.i + '_classification_' +  out_par + '.txt'
        with open(fold_real, 'w') as out_f:
            out_f.write(CLstr + '\n')
            writer = csv.writer(out_f, delimiter='\t')
            writer.writerows(CLlist)

if deg<1:

    ## read wes ##
    name_wes = out_fold + '/validation.txt'
    
    if os.path.exists(name_wes) is False:
        print('File of sequences to score not found, thus exiting ')
        exit()
    wespepss=[]
    with open(name_wes) as f:
        for line in f:
            linesplit = line.strip().split('\n')
            wespepss.append(linesplit[0])

    ## Create folders for output ##
    out_fold_NA = out_fold + '/scoring'
    if os.path.exists(out_fold_NA) is False:
        os.mkdir(out_fold_NA)

    wespepswt = [a for a in wespepss if len(a) in range_len]
    
    if LL == 1:
        wespeps_n = convert_number(wespepswt)
    else:
        name_mat = rootf + '/Align_utils/align_to_seedpy.py'
        name_seqs = rootf + '/Align_utils/seqs_str.txt'
        with open(name_seqs, 'w') as out_f:
            for u in range(len(wespepswt)):
                out_f.write(wespepswt[u] + '\n')

        name_seed = rootf + '/Align_utils/seqs_seed.txt'
        with open(name_seed, 'w') as out_f:
             for u in range(len(iedb_data_fl)):
                 out_f.write(iedb_data_fl[u] + '\n')

        subprocess.call('python3 -W ignore ' + name_mat + ' -sseed ' + name_seed + ' -sseqs ' + name_seqs +' -SA ' + str(SA) + ' -SAmin ' + str(SAmin) + ' -SAmax ' + str(SAmax) + ' -yw 0 ', shell = True)
        name_al = rootf + '/Align_utils/aligned_temp.txt'
        wespeps_n = np.loadtxt(name_al)       
        subprocess.call('rm ' + name_al, shell = True)
        subprocess.call('rm ' + name_seed, shell = True)
        subprocess.call('rm ' + name_seqs, shell = True)
    
    Iwes = RBM.vlayer.compute_output(wespeps_n.astype(np.int16), RBM.weights)

    # Calculate RBM prediction            
    predwes = model.predict(Iwes)
    predwes_max = []
    predwes_labels = []
    for y in range(len(predwes)):
        predwes_max.append(np.max(predwes[y]))
        predwes_labels.append(np.argmax(predwes[y]))

    ent = []
    for t in range(len(predwes)):
       ent.append(-sum(predwes[t]*np.log(predwes[t]+0.000001)))

    # Re-align based on putative class
    wespeps_real=[]
    wespeps_n_real=[]
    wespeps_g_real=[]
    likwes_real=[]
    #scowes_real=[]
    #scowesHMM_real=[]
    lenwes_real=[]
    pred_real=[]
    predl_real=[]
    score_ssrbm=[]
    labels_newv = [labels_new[ff] for ff in flatten_list(iedb_ind_val)]
    values_newv = [values_new[ff] for ff in flatten_list(iedb_ind_val)]
    for cl in range(nc):
        cc = np.where(np.array(labels_newv) == cl)[0]
        ccw = np.where(np.array(predwes_labels) == cl)[0] # For each class, align them to the corresponding seed
        if len(ccw)==1:
            tt = np.copy(int(ccw))
            ccw = np.array([tt,tt])
        if len(cc) and len(ccw):
            wespeps0 = [wespepswt[u] for u in ccw]
            pred0 = [predwes[u][cl] for u in ccw]  
            predl = [predwes_labels[u] for u in ccw]
            if LL == 1:
                wespeps_n = convert_number(wespeps0)
            else:
                seqs_seedt = [iedb_data_fl[iedb_ind_tr[cl][u]] for u in range(len(iedb_ind_tr[cl]))]
                seqs_new =  [iedb_data_fl[flatten_list(iedb_ind_val)[u]] for u in cc]
                pn = [values_newv[u] for u in cc]
                seqs_seed = list(seqs_seedt + seqs_new)
                name_mat = rootf + '/Align_utils/align_to_seedpy.py'
                name_seqs = rootf + '/Align_utils/seqs_str.txt'
                with open(name_seqs, 'w') as out_f:
                    for u in range(len(wespeps0)):
                        out_f.write(wespeps0[u] + '\n')
                
                name_seed = rootf + '/Align_utils/seqs_seed.txt'
                with open(name_seed, 'w') as out_f:
                    for u in range(len(seqs_seed)):
                        out_f.write(seqs_seed[u] + '\n')
                subprocess.call('python3 -W ignore ' + name_mat + ' -sseed ' + name_seed + ' -sseqs ' + name_seqs +' -SA ' + str(SA) + ' -SAmin ' + str(SAmin) + ' -SAmax ' + str(SAmax) + ' -yw 0', shell = True)
                name_al = rootf + '/Align_utils/aligned_temp.txt'
                wespeps_n = np.loadtxt(name_al)
                subprocess.call('rm ' + name_al, shell = True)
                subprocess.c#all('rm ' + name_seqs, shell = True)
                subprocess.call('rm ' + name_seed, shell = True)
                
            wespeps_g = convert_letter(wespeps_n.astype(np.int16))
            #scowes = list(np.loadtxt(rootf + '/Align_utils/scores_cons.txt'))
            #scowes2 = list(np.loadtxt(rootf + '/Align_utils/scores.txt'))
            likwes = RBM.likelihood(wespeps_n.astype(np.int16))
            lenwes = [len(a) for a in wespeps0] 
            wespeps_real.append(wespeps0)
            wespeps_g_real.append(wespeps_g)
            wespeps_n_real.append(wespeps_n)
            likwes_real.append(likwes)
            #scowes_real.append(scowes)
            #scowesHMM_real.append(scowes2)
            lenwes_real.append(lenwes)
            pred_real.append(pred0)
            predl_real.append(predl)

    wespeps_g=[]
    wespeps_n=[]
    likwes=[]
    #scowes=[]
    #scowesHMM=[]
    lenwes=[]
    pred0=[]
    predl=[]
    wespeps_fl = flatten_list(wespeps_real)

    for r in range(len(wespepswt)):
        ii = wespeps_fl.index(wespepswt[r])
        wespeps_g.append(flatten_list(wespeps_g_real)[ii])
        wespeps_n.append(flatten_list(wespeps_n_real)[ii])
        likwes.append(flatten_list(likwes_real)[ii])
        #scowes.append(flatten_list(scowes_real)[ii])
        #scowesHMM.append(flatten_list(scowesHMM_real)[ii])
        lenwes.append(flatten_list(lenwes_real)[ii])
        pred0.append(flatten_list(pred_real)[ii])
        predl.append(flatten_list(predl_real)[ii])


    # print tables for NA
    NAlist = zip(wespepswt, wespeps_g, likwes, lenwes, pred0, predl, likwes - np.array(ent))

    NAstr = 'Sequences' + '\t'+ 'Seqs core' + '\t' + 'RBM Log-Lik.' + '\t' + 'Length' + '\t' + 'Pred. classifier' + '\t'  + 'Pred. HLA class' + '\t'  + 'Total RBM-MHC score'

    with open(out_fold_NA + '/validation_scored_' +  out_par + '.txt', 'w') as out_f:
        out_f.write(NAstr + '\n')
        writer = csv.writer(out_f, delimiter='\t')
        writer.writerows(NAlist)

if NAranking:

    ## read wes ##
    name_wes = out_fold + '/' + args.score + '.txt'
    
    if os.path.exists(name_wes) is False:
        print('File of sequences to score not found, thus exiting ')
        exit()

    wespepss=[]
    with open(name_wes) as f:
        for line in f:
            linesplit = line.strip().split('\n')
            wespepss.append(linesplit[0])

    ## Create folders for output ##
    out_fold_NA = out_fold + '/scoring'
    if os.path.exists(out_fold_NA) is False:
        os.mkdir(out_fold_NA)

    wespepswt = [a for a in wespepss if len(a) in range_len]
    if LL == 1:
        wespeps_n = convert_number(wespepswt)
    else:
        name_mat = rootf + '/Align_utils/align_to_seedpy.py'
        name_seqs = rootf + '/Align_utils/seqs_str.txt'
        with open(name_seqs, 'w') as out_f:
            for u in range(len(wespepswt)):
                out_f.write(wespepswt[u] + '\n')

        name_seed = rootf + '/Align_utils/seqs_seed.txt'
        with open(name_seed, 'w') as out_f:
             for u in range(len(iedb_data_fl)):
                 out_f.write(iedb_data_fl[u] + '\n')

        subprocess.call('python3 -W ignore ' + name_mat + ' -sseed ' + name_seed + ' -sseqs ' + name_seqs +' -SA ' + str(SA) + ' -SAmin ' + str(SAmin) + ' -SAmax ' + str(SAmax) + ' -yw 0 ', shell = True)
        name_al = rootf + '/Align_utils/aligned_temp.txt'
        wespeps_n = np.loadtxt(name_al)
        subprocess.call('rm ' + name_al, shell = True)
        subprocess.call('rm ' + name_seed, shell = True)
        subprocess.call('rm ' + name_seqs, shell = True)
    
    Iwes = RBM.vlayer.compute_output(wespeps_n.astype(np.int16), RBM.weights)

    # Calculate RBM prediction            
    predwes = model.predict(Iwes)
    predwes_max = []
    predwes_labels = []
    for y in range(len(predwes)):
        predwes_max.append(np.max(predwes[y]))
        predwes_labels.append(np.argmax(predwes[y]))

    ent = []
    for t in range(len(predwes)):
       ent.append(-sum(predwes[t]*np.log(predwes[t]+0.000001)))

    # Re-align based on putative class
    wespeps_real=[]
    wespeps_n_real=[]
    wespeps_g_real=[]
    likwes_real=[]
    #scowes_real=[]
    #scowesHMM_real=[]
    lenwes_real=[]
    pred_real=[]
    predl_real=[]
    score_ssrbm=[]
    labels_newv = [labels_new[ff] for ff in flatten_list(iedb_ind_val)]
    values_newv = [values_new[ff] for ff in flatten_list(iedb_ind_val)]
    for cl in range(nc):
        cc = np.where(np.array(labels_newv) == cl)[0]
        ccw = np.where(np.array(predwes_labels) == cl)[0] # For each class, align them to the corresponding seed
        if len(ccw)==1:
            tt = np.copy(int(ccw))
            ccw = np.array([tt,tt])
        if len(cc) and len(ccw):
            wespeps0 = [wespepswt[u] for u in ccw]
            pred0 = [predwes[u][cl] for u in ccw]  
            predl = [predwes_labels[u] for u in ccw] 
            if LL == 1:
                wespeps_n = convert_number(wespeps0)
            else:           
                seqs_seedt = [iedb_data_fl[iedb_ind_tr[cl][u]] for u in range(len(iedb_ind_tr[cl]))]
                seqs_new =  [iedb_data_fl[flatten_list(iedb_ind_val)[u]] for u in cc]
                pn = [values_newv[u] for u in cc]
                seqs_seed = list(seqs_seedt + seqs_new)

                name_mat = rootf + '/Align_utils/align_to_seedpy.py'
                name_seqs = rootf + '/Align_utils/seqs_str.txt'
                with open(name_seqs, 'w') as out_f:
                    for u in range(len(wespeps0)):
                        out_f.write(wespeps0[u] + '\n')
                
                name_seed = rootf + '/Align_utils/seqs_seed.txt'
                with open(name_seed, 'w') as out_f:
                    for u in range(len(seqs_seed)):
                        out_f.write(seqs_seed[u] + '\n')
                subprocess.call('python3 -W ignore ' + name_mat + ' -sseed ' + name_seed + ' -sseqs ' + name_seqs +' -SA ' + str(SA) + ' -SAmin ' + str(SAmin) + ' -SAmax ' + str(SAmax) + ' -yw 0', shell = True)
                name_al = rootf + '/Align_utils/aligned_temp.txt'
                wespeps_n = np.loadtxt(name_al)
                subprocess.call('rm ' + name_al, shell = True)
                subprocess.call('rm ' + name_seqs, shell = True)
                subprocess.call('rm ' + name_seed, shell = True)
            wespeps_g = convert_letter(wespeps_n.astype(np.int16))
            #scowes = list(np.loadtxt(rootf + '/Align_utils/scores_cons.txt'))
            #scowes2 = list(np.loadtxt(rootf + '/Align_utils/scores.txt'))
            likwes = RBM.likelihood(wespeps_n.astype(np.int16))
            lenwes = [len(a) for a in wespeps0] 
            wespeps_real.append(wespeps0)
            wespeps_g_real.append(wespeps_g)
            wespeps_n_real.append(wespeps_n)
            likwes_real.append(likwes)
            #scowes_real.append(scowes)
            #scowesHMM_real.append(scowes2)
            lenwes_real.append(lenwes)
            pred_real.append(pred0)
            predl_real.append(predl)

    wespeps_g=[]
    wespeps_n=[]
    likwes=[]
    #scowes=[]
    #scowesHMM=[]
    lenwes=[]
    pred0=[]
    predl=[]
    wespeps_fl = flatten_list(wespeps_real)

    for r in range(len(wespepswt)):
        ii = wespeps_fl.index(wespepswt[r])
        wespeps_g.append(flatten_list(wespeps_g_real)[ii])
        wespeps_n.append(flatten_list(wespeps_n_real)[ii])
        likwes.append(flatten_list(likwes_real)[ii])
        #scowes.append(flatten_list(scowes_real)[ii])
        #scowesHMM.append(flatten_list(scowesHMM_real)[ii])
        lenwes.append(flatten_list(lenwes_real)[ii])
        pred0.append(flatten_list(pred_real)[ii])
        predl.append(list_hla_l[flatten_list(predl_real)[ii]])


    # print tables for NA
    NAlist = zip(wespepswt, wespeps_g, likwes, lenwes, pred0, predl, likwes - np.array(ent))

    NAstr = 'Sequences' + '\t'+ 'Seqs core' + '\t' + 'RBM Log-Lik.' + '\t' + 'Length' + '\t' + 'Pred. classifier' + '\t'  + 'Pred. HLA class' + '\t'  + 'Total RBM-MHC score'

    with open(out_fold_NA + '/' + args.score + '_scored_' +  out_par + '.txt', 'w') as out_f:
        out_f.write(NAstr + '\n')
        writer = csv.writer(out_f, delimiter='\t')
        writer.writerows(NAlist)


if makegeneration: 
### GENERATE NEW SEQUENCES with MC exploration of space hidden inputs ###

    ## Create folders for output ##
    out_fold_gen = out_fold + '/Generative'
    if os.path.exists(out_fold_gen) is False:
        os.mkdir(out_fold_gen)       

    # MC exploration parameters
    nsteps = 1000 # total steps
    nsample = 100 # sample every 100 steps
    beta = 50

    # gen_data parameters
    N_sequences = 100 # Total number of sequences per chain
    Nchains = 100 # Number of Markov chains
    Nthermalize = 500 # Number of Gibbs sampling steps to perform before the first sample of a chain.
    Nstep = 10 # Number of Gibbs sampling steps between each sample of a chain

    # class by class generation
    I = RBM.vlayer.compute_output(iedb_data_fl_nf.astype(np.int16),RBM.weights)

    ind_synth = []
    cat_synth = []
    for cl in range(nc):
         pred = model.predict(I[iedb_ind[cl]])
         ind_opt = np.argmax([pred[l][cl] for l in range(len(pred))])
         opt_I = np.copy(I[iedb_ind[cl]][ind_opt])
         opt_h = RBM.hlayer.mean_from_inputs(opt_I)
         cond = [(l, opt_h[0][l]) for l in range(hu-1)]
         condRBM1 = RBM_utils.conditioned_RBM(RBM,cond)
         samplev,sampleh = condRBM1.gen_data(Nchains = Nchains, Lchains = int(N_sequences/Nstep), Nthermalize=Nthermalize)

            # Initialization
         I_old = np.copy(opt_I)
         energy_old = -np.log(pred[ind_opt][cl])
         I_total = np.copy(I_old)

            # MC chain
         liste=[]
         listp=[]

         for nt in range(1,nsteps):
            dI = np.zeros(n_h)
            for i in range(len(dI)):
                dI[i] = np.random.normal(0,0.01)
            I_new = np.copy(I_old + dI)
            energy_new = -np.log(model.predict(np.vstack((I_new,I_new)))[0][cl])  
            de = energy_new - energy_old 
            p = np.exp(-beta*de)
            if np.random.uniform(0,1) <= p:
                I_old = np.copy(I_new)
                energy_old = energy_new
            liste.append(energy_old) # should get smaller
            listp.append(model.predict(np.vstack((I_old,I_old)))[0][cl]) # should stay always above a threshold which is set by beta
            I_total = np.vstack((I_total, I_old))  

            if not(nt%nsample):
                hm = RBM.hlayer.mean_from_inputs(I_old.astype(curr_float))
                cond = [(l, hm[0][l]) for l in range(hu-1)]
                condRBM1 = RBM_utils.conditioned_RBM(RBM,cond)
                datav1, datah1 = condRBM1.gen_data(Nchains = Nchains, Lchains = int(N_sequences/Nstep), Nthermalize=Nthermalize)
                samplev = np.vstack((samplev,datav1))

         if cl==0:
             sample_synth = np.vstack((samplev))
             ind_synth.append(list(range(len(samplev))))
             cat_synth.append([cl]*len(samplev))
         else:
             start = len(sample_synth)
             sample_synth = np.vstack((sample_synth, samplev))
             ind_synth.append(list(range(start, len(sample_synth))))
             cat_synth.append([cl]*len(samplev))

        ## print table ##
         sample_synth_g = convert_letter(sample_synth)
         OLlist = zip(sample_synth_g, flatten_list(cat_synth))
         OLstr = 'Synthetic data' + '\t'+ 'HLA'
         with open(out_fold_gen + '/Synthetic_sample_' +  out_par + '.txt', 'w') as out_f:
             out_f.write(OLstr + '\n')
             writer = csv.writer(out_f, delimiter='\t')
             writer.writerows(OLlist)
            

    # Motifs of generated sequences 

    for cl in range(nc):
        pwm = average_n(sample_synth[ind_synth[cl]].astype(np.int16), CC)
        stringt  = list_hla_l[cl] + ' # seqs ' + str(len(sample_synth[ind_synth[cl]]))
        fig = sequence_logo.Sequence_logo(pwm, figsize=(15,3), ylabel= 'bits', title = stringt, show=True, ticks_every=5, ticks_labels_size=20, title_size=24);
        if makefig: 
            pathfig = out_fold_gen + '/allotype_'+ str(cl) + '_' + out_par + '.pdf'
            fig.savefig(pathfig)
        plt.close()
        
