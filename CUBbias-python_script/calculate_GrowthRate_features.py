#calculate growth rate features
from Bio import SeqIO
from MILC import cal_MILC
from CPB import cal_CPB
from AACB import cal_AA_Bias

import numpy as np
import os
import re


def calculate_index(path_to_genome,path_to_ribosome,genome_suffix,ribosomal_suffix):
	"""
	calculate growth rate features
	Args:

		path_to_genome: The dir path of your bacterial genome file(s); example: "D:\\dir1\\dir2\\genome_files\\"
		path_to_ribosome: The dir path of your bacterial ribosomal protein file(s); expample: example:"D:\\dir1\\dir3\\ribosome_files\\"
		genome_suffix: your genome file suffix; example: "_cds_from_genomic.fna"
		ribosomal_suffix: #your ribosome gene file suffix; example: "_ribosomal_cds.fna"
	
	Returns:
		#species_name	CUB_HE	CUB_consistency	AAI_HE	AAI_consist; feature values
	"""
	genome_file_list = os.listdir(path_to_genome)
	ribosomal_file_list = os.listdir(path_to_ribosome)

	print('species_name\tCUB_HE\tCUB_consistency\tCPB\tAAI_HE\tAAI_consist')
	for genome_file in  genome_file_list:
		species_name = genome_file.replace(genome_suffix,'')
		genome_file = path_to_genome + genome_file
		ribosomal_file = path_to_ribosome + species_name +ribosomal_suffix 
		#genome sequence as a reference, calculate CUB_HE
		CUB_HE = np.median(cal_MILC(genome_file,ribosomal_file,cutoff=240, rm_start_stop_codons=True))
		#ribosome gene sequence as a reference, calculate CUB_consistency
		CUB_consistency = np.mean(cal_MILC(ribosomal_file,ribosomal_file,cutoff=0, rm_start_stop_codons=True))
		#calulate codon pair bias(CPB)
		CPB = cal_CPB(genome_file,cutoff=240,rm_start_stop_codons=True)
		#calulate amino acid compositon bias
		AAI_HE = np.median(cal_AA_Bias(genome_file,ribosomal_file))
		AAI_consist = np.mean(cal_AA_Bias(ribosomal_file,ribosomal_file))

		print(species_name+'\t'+str(CUB_HE)+'\t'+str(CUB_consistency)+'\t'+str(CPB)+'\t'+str(AAI_HE)+'\t'+str(AAI_consist))


if __name__ == '__main__':

	#path_to_ribosome = "The dir path of your bacterial genome file(s)"
	#path_to_genome = "The dir path of your bacterial ribosomal protein file(s)"
	#genome_suffix = 'your genome file suffix'
	#ribosomal_suffix = 'your ribosome gene file suffix'

	LenArgv = len(sys.argv)

	for i in range(LenArgv):
		if sys.argv[i] == '-ribo_dir':
			path_to_ribosome = sys.argv[i+1]
		if sys.argv[i] == '-genome_dir':
			path_to_genome = sys.argv[i+1]

		if sys.argv[i] == '-ribo_suffix':
			ribosomal_suffix = sys.argv[i+1]
		if sys.argv[i] == '-genome_suffix':
			genome_suffix = sys.argv[i+1]

	calculate_index(path_to_genome,path_to_ribosome,genome_suffix,ribosomal_suffix)
