# AACBias

This repository contains code for "Improves estimation of maximal microbial growth rate by additionally considering amino acid composition bias". 

Scripts for calculating CUBias and AACBias features from genomic sequences. 
multiple linear regression of those features to the originating species' Doubling Time,and compared prediction model between CUBias and CUBias+AACBias.

Demostration
1. calculate Amino Acid Composition Bias(AACBias);  
  (1) Move into the amino acid compositon bias(AminoAcidBias) directory  
    cd AminoAcidBias  
  (2) Measuring AACBias    
    python3  AACB.py -ribo 'ribosomal gene fileName' -genome 'genome fileName'  
    
2. calculate Codon Usage Bias(CUBias) using python3;  
  (1) Move into the CUBbias-python_script directory  
    cd CUBbias-python_script  
  (2) Measuring CUBias   
    python3 calculate_GrowthRate_features.py -ribo_dir 'Directory path of your bacterial ribosomal gene file(s)'  
                                              -genome_dir 'Directory path of your bacterial genome file(s)'  
                                              -ribo_suffix 'Ribosome gene filename suffix'  
                                              -genome_duffix 'genome filename suffix'  
                                              
3. View the results of model(CUBias VS CUBias+AACBias) comparison in jupyter  
  (1) Move into the data directory  
     cd data  
     open LinearRegression_CUBiasModel_VS_CUBias_AACBiasModel.ipynb  
  
