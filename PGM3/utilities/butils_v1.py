import utilities as utilities
import rbm
import numpy as np

curr_int = np.int16
curr_float = np.float32
rootf = '/home/barbara/Barbara_Bravi'

# Check quality of statistics
def calculate_error(RBM, data_tr, N_sequences = 800000, Nstep = 10, background = None):
    N=RBM.n_v
    q=RBM.n_cv

    # Check how moments are reproduced
    # Means
    mudata = RBM.mu_data # empirical averages
    #datav, datah = RBM.gen_data(Nchains = int(100), Lchains = int(N_sequences/100), Nthermalize=int(500), background= background)
    datav, datah = RBM.gen_data(Nchains = int(100), Lchains = int(N_sequences/100), Nthermalize=int(500))
    mugen = utilities.average(datav, c = q, weights=None)

    # Correlations
    covgen = utilities.average_product(datav,datav,c1=q,c2=q) - mugen[:,np.newaxis,:,np.newaxis] * mugen[np.newaxis,:,np.newaxis,:]
    covdata = utilities.average_product(data_tr,data_tr,c1=q,c2=q) - mudata[:,np.newaxis,:,np.newaxis] * mudata[np.newaxis,:,np.newaxis,:]
    fdata= utilities.average_product(data_tr,data_tr,c1=q,c2=q)

    #put to zero the diagonal elements of the covariance
    for i in range(N):
        covdata[i,i,:,:] = np.zeros((q,q))
        fdata[i,i,:,:] = np.zeros((q,q))
        covgen[i,i,:,:] = np.zeros((q,q))

    M=len(data_tr)
    maxp=float(1)/float(M)
    pp = 1
    ps = 0.00001 # pseudocount for fully conserved sites

    # error on frequency
    pp=1
    errm=0
    neffm=0
    for i in range(N):
        for a in range(q):
            neffm += 1
            if mudata[i,a] < maxp:
                errm += np.power((mugen[i,a] - mudata[i,a]),2)/(float(1-maxp)*float(maxp))
            else:
                if mudata[i,a] != 1.0:
                    errm += np.power((mugen[i,a] - mudata[i,a]),2)/(float(1-mudata[i,a])*float(mudata[i,a]))  
                else:
                    errm += np.power((mugen[i,a] - mudata[i,a]),2)/(float(1-mudata[i,a]-ps)*float(mudata[i,a])) 

    errmt = np.sqrt(float(1)/(float(neffm)*float(maxp))*float(errm))
    # rigourously, there would be also the regularization term in the difference errm!

    # error on correlations
    errc=0
    neffc=0
    for i in range(N):
        for j in range(i+1,N):
            for a in range(q):
                for b in range(a+1,q):
                    neffc+=1
                    if covdata[i,j,a,b] < maxp:
                        den = np.power(np.sqrt(float(1-maxp)*float(maxp)) + mudata[i,a] * np.sqrt(mudata[j,b] * (1 - mudata[j,b])) + mudata[j,b] * np.sqrt(mudata[i,a] * (1 - mudata[i,a])),2) 
                        errc += np.power((covgen[i,j,a,b] - covdata[i,j,a,b]),2)/float(den)
                    else:
                        den = np.power(np.sqrt(float(1-fdata[i,j,a,b])*float(fdata[i,j,a,b])) + mudata[i,a] * np.sqrt(mudata[j,b] * (1 - mudata[j,b])) + mudata[j,b] * np.sqrt(mudata[i,a] * (1 - mudata[i,a])),2) 
                        errc += np.power((covgen[i,j,a,b] - covdata[i,j,a,b]),2)/float(den)

    errct = np.sqrt(float(1)/(float(neffc)*float(maxp))*float(errc))
    print('calculate error finish')
    return (errmt, errct)

def convert_number(seqs): # convert to numbers already aligned seqs
    aa = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V',  'W', 'Y','-']
    aadict = {aa[k]: k for k in range(len(aa))} 
    
    msa_num = np.array(list(map(lambda x: [aadict[y] for y in x], seqs[0:])), dtype=curr_int, order="c") ### Here change ####
    
    return msa_num

def convert_letter(seqs_n): # convert to numbers already aligned seqs
    aa = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V',  'W', 'Y','-']
    aadictinv = {k: aa[k] for k in range(len(aa))} 
    
    seqs=[]
    for t in range(len(seqs_n)):
        seqs.append(''.join([aadictinv[e] for e in seqs_n[t]]))
    
    return seqs

def align_seq(seqs,lea, lepmin, lepmax, weights=None):
    import matlab.engine
    eng = matlab.engine.start_matlab()
    eng.addpath(rootf + '/rbm/PGM3/utilities')
    if weights is None:
       alignc = eng.balign_peps(seqs, lea, lepmin, lepmax)
    else:
       alignc = eng.balign_peps(seqs, lea, lepmin, lepmax, weights)
    eng.quit()
    msa_al = np.asarray(alignc).astype(curr_int) - 1 # needed for different indexing convention in matlab
    return msa_al

def align_seqTCR(seqs, lea):
    import matlab.engine
    eng = matlab.engine.start_matlab()
    eng.addpath(rootf + '/rbm/PGM3/utilities')
    alignc = eng.balign(seqs, lea-1)
    eng.quit()
    msa_al = np.asarray(alignc).astype(curr_int) - 1 # needed for different indexing convention in matlab
    return msa_al

def align_to_seed(seqs_seed, lea, lepmin, lepmax, seqs_new, weights=None):
    import matlab.engine
    eng = matlab.engine.start_matlab()
    eng.addpath(rootf + '/rbm/PGM3/utilities')
    if weights is None:
        alignc_new = eng.balign_peps_seed(seqs_seed, lea, lepmin, lepmax, seqs_new)
    else:
        alignc_new = eng.balign_peps_seed(seqs_seed, lea, lepmin, lepmax, seqs_new, weights) # lea -> final length, lepmin to lepmax -> length included in the profile itself
    eng.quit()
    msa_al = np.asarray(alignc_new).astype(int) - 1 # needed for different indexing convention in matlab
    return msa_al

def align_ranto_seed(seqs_seed, lea, lepmin, lepmax, seqs_new, weights=None):
    import matlab.engine
    eng = matlab.engine.start_matlab()
    eng.addpath(rootf + '/rbm/PGM3/utilities')
    if weights is None:
        alignc_new = eng.balign_ranpeps_seed(seqs_seed, lea, lepmin, lepmax, seqs_new)
    else:
        alignc_new = eng.balign_ranpeps_seed(seqs_seed, lea, lepmin, lepmax, seqs_new, weights) # lea -> final length, lepmin to lepmax -> length included in the profile itself
    eng.quit()
    msa_al = np.asarray(alignc_new).astype(curr_int) - 1 # needed for different indexing convention in matlab
    return msa_al

def count_neighbours(MSA,threshold = 0.1): # Compute reweighting
    B = MSA.shape[0]
    N = MSA.shape[1]
    num_neighbours = np.zeros(B)
    for b in range(B):
        if b%1000 ==0:
            print(b)
        num_neighbours[b] =  ((MSA[b] != MSA).mean(1) < threshold).sum()
    return num_neighbours

def float_to_str(f):
    float_string = repr(f)
    if 'e' in float_string:  # detect scientific notation
        digits, exp = float_string.split('e')
        digits = digits.replace('.', '').replace('-', '')
        exp = int(exp)
        zero_padding = '0' * (abs(int(exp)) - 1)  # minus 1 for decimal point in the sci notation
        sign = '-' if f < 0 else ''
        if exp > 0:
            float_string = '{}{}{}.0'.format(sign, digits, zero_padding)
        else:
            float_string = '{}0.{}{}'.format(sign, zero_padding, digits)
    return float_string

