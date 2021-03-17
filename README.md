## RBM-MHC Python3 script ##

Python3 script to implement RBM-MHC, a Restricted Boltzmann Machine (RBM)-based method to predict antigen presentation and for HLA-I motif reconstruction. 

## Copyright ##
Copyright 2020 - by Barbara Bravi (bbravi.bb@gmail.com)
All rights reserved
     
Permission is granted for anyone to copy, use, or modify this
software for any uncommercial purposes, provided this copyright 
notice is retained, and note is made of any changes that have 
been made. This software is distributed without any warranty, 
express or implied. In no event shall the author or contributors be 
liable for any damage arising out of the use of this software.

The publication of research using this software, modified or not, must include appropriate citations to: Bravi et al. Cell Systems 2021 https://www.sciencedirect.com/science/article/pii/S2405471220304567

## Download and Install the packages ## 

RBM-MHC is written in Python version 3.6.9

Packages required: biopython, numba, keras (with Theano or tensorFlow backend), scikit-learn, along with standard packages (numpy, cython, matplotlib). The alignment routines require matlab and matlab engine API for Python https://fr.mathworks.com/help/matlab/matlab_external/install-the-matlab-engine-for-python.html

RBM-MHC has a dependency on the RBM implementation Copyright 2020 by Jerome Tubiana, imported as a submodule and downloadable from https://github.com/jertubiana/PGM

The RBM-MHC.py script, the folders Align_utils and PGM should be saved in the same folder. The path to this folder should be specified inside setup.py (assigned to NAME_FOLDER). Run this script to set the right path:

python3 setup.py

Antigen sequences annotated by their HLA-binding specificity are used to guide motif reconstruction. The current implementation searches annotated ligands in the file mhc_ligand_full.csv from the Immune Epitope Database downloadable here:
http://www.iedb.org/database_export_v3.php

## Run the script ##
We assume to have a custom dataset, e.g. a dataset of unannotated peptides from an elution experiment saved in sample_file.txt, and the HLA-I known to be expressed in the sample, e.g. from HLA-typing, are HLA-A\*01:01, HLA-A\*03:01, HLA-B\*07:02, HLA-B\*08:01, HLA-C\*07:01, HLA-C\*07:02. Running the RBM-MHC script allows to build a presentation model able to: 
- assign an HLA-I type to each peptide in sample_file.txt;
- assign presentation scores to a custom list of peptides, e.g. all peptides harbouring cancer-specific mutations, to understand which ones have high likelihood of being presented (see option -score).

To run the script, an example command line is: 

python3 RBM-MHC.py -hla 'HLA-A\*01:01' 'HLA-A\*03:01' 'HLA-B\*07:02' 'HLA-B\*08:01' 'HLA-C\*07:01' 'HLA-C\*07:02' -rl 9 -i 'sample_file' -o 'output_folder' -nameo 'string_output' -score 'peptides-to-score'

This command line reads the peptides of length 9 residues from the file NAME_FOLDER/output_folder/sample_file.txt, trains RBM-MHC on them and predicts the peptides specifically binding to the 6 HLA-I provided. The HLA assignment is guided by an amount of labelled peptides for these specificities, equal to 0.1 of the sample size, extracted from IEDB. The model is used to assign probabilistic scores of presentation to the peptides in the file NAME_FOLDER/output_folder/peptides-to-score.txt. The output (trained RBM and HLA-I classifier, table of peptides with assigned HLA-binding specificity, scored peptides) is saved in NAME_FOLDER/output_folder. 

The 'output_folder' provided as example contains a python notebook that explains how to read and analyze the results from the command above. When considering only fixed-length peptides (as here 9 residues) no alignment, hence no Matlab routine, needs to be called. 
 

## Options ##

-hla : list of HLA-I alleles characterizing the sample to analyze or, if a sample is not provided, for which data on IEDB should be searched; if for a given sample this option is not provided, only a RBM is trained and RBM presentation scores are assigned to the peptides in the file provided by the -score option

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

-ba : When set to 1, the search in IEDB for peptides of known HLA type is performed among data from binding assays (default = 0)

-rwnms : When set to 1, it implements the re-weighting scheme to correct for differences between amino acid frequencies of antigens detected by Mass Spectrometry and by other techniques (default = 0)

-rwhp :  When set to 1, it implements the re-weighting scheme to correct for differences between amino acid frequencies of antigens detected by Mass Spectrometry and the human proteome (default = 0)

-score : Name of the file inside the output folder with the list of peptides to score (as peptides containing missense mutations from Whole Exome Sequencing of a tumour sample); the list of peptides and corresponding scores is saved in a subfolder 'scoring' (default = 0)

-gen : It generates 10000 synthetic peptides of each HLA specificity (default = 0) and saves them into the subfolder 'Generative'

-fig : It prints motif figures (default = 0)

## Data retrieval from Immune Epitope Database ###

To guide motif reconstruction in samples of interest for the user, RBM-MHC uses an amount of peptides equal to (percentage provided by -perc)x(sample size) labelled with their HLA association. The current routine seeks an overlap between the sample and IEDB in such a way that these peptides can be assigned a known HLA preference. If there is little overlap, it adds labelled peptides from IEDB. The search first targets monoallelic source data as described in the paper. If less than a given amount of sequences (set by default to 300) are retrieved, the search is first extended to 'Allele Specific Purification' mass spectrometry (MS) data, next, if the latter are not available, to all MS-data, next to data obtained by all techniques (thus e.g. by binding affinity assays). Regardless of the technique, this labelled set of sequences is chosen among ligands annotated as 'positive', 'positive-high' to the HLA-I alleles under consideration. By the option -ba the search for labelled peptides can be instead performed in IEDB data from binding assays. 

## Trouble-shooting
version of h5py should be < 3 

