## ssRBM Python3 script ##

Python3 script to implement the ssRBM (semi-supervised RBM) architecture described in Bravi et al. 2020: a method to predict antigen presentation and for HLA-I motif deconvolution 

## Copyright ##
Copyright 2020 - by Barbara Bravi (bbravi.bb@gmail.com)
All rights reserved
     
Permission is granted for anyone to copy, use, or modify this
software for any uncommercial purposes, provided this copyright 
notice is retained, and note is made of any changes that have 
been made. This software is distributed without any warranty, 
express or implied. In no event shall the author or contributors be 
liable for any damage arising out of the use of this software.

The publication of research using this software, modified or not, must include 
appropriate citations to:

## Download and Install the packages ## 

ssRBM is written in Python version 3.6.9

Packages required: biopython, numba, keras (with Theano or tensorFlow backend), scikit-learn, along with standard packages (numpy, cython, matplotlib).
The alignment routines require matlab and matlab engine API for Python https://fr.mathworks.com/help/matlab/matlab_external/install-the-matlab-engine-for-python.html


The ssRBM script , the folder Align_utils and the package PGM3 (RBM implementation in Python3 by Jerome Tubiana, see ...) should be saved in the same folder. 
The path to this folder should be  specified inside setup.py (assigned to NAME_FOLDER). Run this script to set the right path:

python3 setup.py


Create or download in the same folder the file containing antigen sequences annotated by their HLA-specificity to guide motif deconvolution in the supervised learning step. The current implementation searches annotated ligands in the file mhc_ligand_full.csv from the Immune Epitope Database downloadable here:
http://www.iedb.org/database_export_v3.php


## Run the script ##

To run the script, an example command line is: 
python3 ssRBM.py -hla 'HLA-A\*01:01' 'HLA-A*02:01' 'HLA-B*15:01' 'HLA-B*27:02' 'HLA-C*08:02' 'HLA-C*16:01' -rl 9 10 -i 'sample_file' -o 'output_folder' -nameo 'string_output'

This command line reads the peptides of length 9-10 residues from the file NAME_FOLDER/output_folder/sample_file.txt, trains a RBM on them and deconvolves the peptides specifically binding to HLA-A*01:01, HLA-A*02:01, HLA-B*15:01, HLA-B*27:02, HLA-C*08:02, HLA-C*16:01 (HLA alleles expressed in the sample known from HLA typing). The deconvolution is guided by an amount of labelled peptides for these specificities, equal to 0.1 of the sample size, extracted from IEDB. The output (trained RBM, trained classifier, table of peptides with assigned specificity) is saved in NAME_FOLDER/output_folder. 


Options: 

-hla : list of HLA-I alleles characterizing the sample to analyze or, if a sample is not provided, for which data on IEDB should be searched; if for a given sample 
this option is not provided, only a RBM is trained and RBM presentation scores are assigned to the peptides in the file provided by the -score option

-i : Name of the input file with peptide sequences to analyze, saved in format .txt inside the main folder (NAME_FOLDER)

-o : Name of output folder that will be created inside the main folder (NAME_FOLDER)

-nameo: String that will be used in the name of files produced

-rl : Range of peptide lengths in number of residues to consider (default choice is the range 8-11)

-mt : When set to 0, it disables the training of the method and reads an existing model (default = 1)

-deg : Fraction of sample data retained for the RBM training (default = 1); if deg is smaller than 1, the training is done on one part of the dataset and the rest, left for validation, is scored and printed to the file validation.txt  

-hu : Number of RBM hidden units (default = 10)

-l12 : L^1_2 RBM regularization (default = 0.001)

-niter : Number of iterations for the RBM training (default = 100)

-niterd : Number of iterations for training the HLA-I classifier (default = 1000)

-perc : Equivalent fraction of the sample size to be labelled for training the HLA-I classifier (default = 0.1)

-al : When set to 0, it disables the re-iteration of the alignment after classification (default = 1)

-rwnms : When set to 0, it implements the re-weighting scheme to correct for differences between amino acid frequencies of antigens detected by Mass Spectrometry and by other techniques (default = 0)

-rwhp :  When set to 0, it implements the re-weighting scheme to correct for differences between amino acid frequencies of antigens detected by Mass Spectrometry and the human proteome (default = 0)

-score : Name of the file inside the output folder with the list of peptides to score (as peptides containing missense mutations from Whole Exome Sequencing of a tumour sample); the list of peptides and corresponding scores is saved in a subfolder 'scoring' (default = 0)

-gen : It generates 10000 synthetic peptides of each HLA specificity (default = 0) and saves them into the subfolder 'Generative'

-fig : It prints motif figures (default = 0)


## Data retrieval from Immune Epitope Database ###

To guide motif deconvolution in samples of interest for the user, ssRBM routine seeks in IEDB (vita 2019) a given amount of 
peptides labelled with the HLA-specificities expected in the sample. The search first targets monoallelic source data as described above. If less than a given amount of sequences (set by default to 300) are retrieved, the search is first extended to ``Allele Specific Purification'' MS-data, next, if the latter are not available, to all MS-data, next to data obtained by all techniques (thus e.g.\ by binding affinity assays). Regardless of the technique, this labelled set of sequences is chosen among ligands annotated as ``positive'', ``positive-high'' to the HLA-I alleles under consideration.

