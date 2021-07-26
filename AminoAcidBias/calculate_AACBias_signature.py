#usage: python3  calculate_AACBias_signature.py > score.txt
from Bio import SeqIO
from AACB import cal_AA_Bias

import numpy as np
from pprint import pprint
import os
import re
import sys


def calculate_index(path_to_genome,path_to_ribosome,genome_suffix='_cds_from_genomic.fna',ribosomal_suffix='_ribosomal_cds.fna'):
	genome_file_list = os.listdir(path_to_genome)
	ribosomal_file_list = os.listdir(path_to_ribosome)

	print('species_name\tAACB_HE\tAACB_consistency')
	for genome_file in  genome_file_list:
		species_name = genome_file.replace(genome_suffix,'')
		genome_file = path_to_genome + genome_file
		ribosomal_file = path_to_ribosome + species_name +ribosomal_suffix 

		AACB_HE = np.median(cal_AA_Bias(genome_file,ribosomal_file))
		AACB_consistency = np.mean(cal_AA_Bias(ribosomal_file,ribosomal_file))

		print(species_name+'\t'+str(AACB_HE)+'\t'+str(AACB_consistency))		

if __name__ == '__main__':

	path_to_ribosome = "The dir of your genome file(s)" #example:"D:\\dir1\\dir2\\genome_files\\"
	path_to_genome = "The dir of your ribosome file(s)" #example:"D:\\dir1\\dir3\\ribosome_files\\"
	genome_suffix = 'genome file suffix' #your genome file suffix,example: "_cds_from_genomic.fna"
	ribosomal_suffix = 'ribosome file suffix' #your ribosome gene file suffix,example: "_ribosomal_cds.fna"

	calculate_index(path_to_genome,path_to_ribosome,genome_suffix='_cds_from_genomic.fna',ribosomal_suffix = '_ribosomal_cds.fna')