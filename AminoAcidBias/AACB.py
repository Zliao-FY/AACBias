'''
function: calculate amino acid composition bias
inputs: two fasta files; the first one is a reference DNA fasta file(genome fasta file), second one is a cds fasta file(ribosome gene fasta file)
'''
import sys
from pprint import pprint
from Bio import SeqIO
from Bio.Seq import Seq
from collections import Counter
import numpy as np

nucl = ['A','T','C','G','a','t','c','g']
#signature_AA have three amino acids: E,V,Q
signature_AA = ['E','V','Q']


def check_cds_Right_cutoff(cds,cutoff=240):
	cds_len = len(cds)
	if cds_len < cutoff:
		return False
	if cds_len % 3 != 0:
		return False
	nucl_cds = set(str(cds))
	for n in nucl_cds:
		if n not in nucl:
			return False
	return True

def AAI_reference_sequence(referenceFile):
	#reference sequence
	reference_protein = ''
	for record in SeqIO.parse(referenceFile,'fasta'):
		if not check_cds_Right_cutoff(record.seq):
			continue
		reference_protein += cds_AA_sequence(str(record.seq))

	return str(reference_protein)

def cds_AA_sequence(cdss):
	#delete stop codon
	cds_protein = Seq(cdss).translate(table='1',stop_symbol='*')
	if cds_protein.endswith('*'):
		cds_protein = cds_protein[:-1]
	return str(cds_protein)

def sequence_AA_table(sequence_aa):
	#the frequency of amino acid  in a sequence
	seq_aa_table = Counter(sequence_aa)
	sum_aa = sum(seq_aa_table.values())
	return {aa:aa_value/sum_aa for aa,aa_value in seq_aa_table.items()}

def AA_Bias(ref_AA_table,cds_seq):
	cds_len = len(cds_seq)
	cds_AA_table = sequence_AA_table(cds_seq)

	cds_signature_AA_table = {aa:aa_value for aa, aa_value in cds_AA_table.items() if aa in signature_AA}

	M_AA = sum([aa_freq* np.log(aa_freq/ref_AA_table[aa]) for aa,aa_freq in cds_signature_AA_table.items()])

	return M_AA

def cal_AA_Bias(ref_file,cds_or_ref_file):
	ref_seq = AAI_reference_sequence(ref_file)
	ref_AA_table = sequence_AA_table(ref_seq)

	AA_Bias_value = []

	for record in SeqIO.parse(cds_or_ref_file,'fasta'):
		if not check_cds_Right_cutoff(record.seq):
			continue
		cds_seq = cds_AA_sequence(str(record.seq))
		aai_distance = AA_Bias(ref_AA_table, cds_seq)

		AA_Bias_value.append(aai_distance)
	
	return AA_Bias_value


if __name__ == '__main__':
	LenArgv = len(sys.argv)

	for i in range(LenArgv):
		if sys.argv[i] == '-ribo':
			ribo_file = sys.argv[i+1]
		if sys.argv[i] == '-genome':
			genome_file = sys.argv[i+1]

	AACB_HE = np.median(cal_AA_Bias(genome_file, ribo_file))
	AACB_consistency = np.mean(cal_AA_Bias(ribo_file, ribo_file))

	print('AACB_HE = ' , AACB_HE)
	print('AACB_consistency = ' , AACB_consistency)
