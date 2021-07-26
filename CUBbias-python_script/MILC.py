#usage: python3 predict_doublingTime5 reference_file.fasta cds_file.fasta  > milc.score
#function: This module contains code measuring CUB independent of length and composition(MILC)
#reference file: genome file; cds_file: ribosome file; cutoff=240
#20 standard amino acid
import sys
import re
from collections import Counter
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Data import CodonTable
from Bio.Seq import Seq
from pprint import pprint
from math import log
import warnings

#biopython stardard codonTable
xcodontable = CodonTable.unambiguous_dna_by_id[11]
#codon:aa dictionary
ycodontable = xcodontable.forward_table
#amino acid list
zcodontable = [ycodontable[i] for i in ycodontable]
#codon list except three stop codons
qcodontable = [i for i in ycodontable ]
stop_codons = xcodontable.stop_codons

nucl = ['A','T','C','G','a','t','c','g']
def check_cds_Right_cutoff(cds,cutoff=240):
	cds = str(cds)
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

def aa_codon_table(y_codontable):
	aa_table = {}
	for c,aa in y_codontable.items():
		if aa in aa_table:
			aa_table[aa].append(c)
		else:
			aa_table[aa] = [c]
	return aa_table
aa_codontable = aa_codon_table(ycodontable)

def get_codon_count(cdss):
	#cdss is a string of DNA sequence(cds)
	#return format: {'codon1':count1,'codon2':count2...}
	len_cdss = len(cdss)
	return Counter([cdss[i:i+3] for i in range(0,len_cdss,3)])

def get_codon_freq_inAA(cds_codon_count,cds_used_aa):
	'''
	ie: Pro -> [CCT,CCC,CCA,CCG]
		f(cct) = Occt / (Occt + Occc + Occa +Occg) = Occt/Opro
	return {'CCT relative frequence' : 0.21, 
			etc}

	'''
	codon_freq = {}
	for aa in cds_used_aa:
		for codon in aa_codontable[aa]:
				if codon in cds_codon_count:
					codon_freq[codon] = cds_codon_count[codon] / cds_used_aa[aa]
				else:
					codon_freq[codon] = 0
	return codon_freq


def get_ref_codon_freq_inAA(refFile):
	'''
		calculate the frequence of codon in the synonymous codons list in a reference sequence
	
	'''
	ref_seq = ''.join([str(record.seq) for record in SeqIO.parse(refFile,'fasta') if check_cds_Right_cutoff(record.seq,cutoff=0)])
	ref_codon_count = get_codon_count(ref_seq)
	ref_aa = str(Seq(ref_seq).translate(stop_symbol=''))
	ref_used_aa = Counter(ref_aa)
	return get_codon_freq_inAA(ref_codon_count,ref_used_aa)

def MILC(ref_codon_freq,cdss,rm_start_stop_codons=True):
	if rm_start_stop_codons:
		cdss = cdss[3:-3]

	cdss_len = len(cdss)
	cdss_codon_count = get_codon_count(cdss)

	cdss_aa = Seq(cdss).translate(stop_symbol='')
	cdss_used_aa = Counter(cdss_aa)
	cdss_codon_freq = get_codon_freq_inAA(cdss_codon_count,cdss_used_aa)

	def Ma(aa):
		#calculate Ma 
		m_a = 0
		if aa in cdss_aa:
			for codon in aa_codontable[aa]:
				if codon in cdss_codon_count:
					m_a += 2 * cdss_codon_count[codon] * log(cdss_codon_freq[codon]/ref_codon_freq[codon])
		return m_a
	
	M = sum([Ma(aa) for aa in cdss_used_aa])
	ra = sum([len(aa_codontable[aa])-1 for aa in cdss_used_aa]) #R stop codon 3-1
	L = float(cdss_len/3)-1
	C = ra / L - 0.5
	milc = M/L -C
	return milc


def cal_MILC(reference_file, cds_File, cutoff=240,rm_start_stop_codons=True):
	ref_codon_freq = get_ref_codon_freq_inAA(reference_file)
	milc_list = []
	for record in SeqIO.parse(cds_File,'fasta'):
		if not check_cds_Right_cutoff(record.seq,cutoff=cutoff):
			continue

		#print(record.id+'=',end='\t')
		milc_value = MILC(ref_codon_freq,str(record.seq).upper(),rm_start_stop_codons=rm_start_stop_codons)
		milc_list.append(milc_value)
	return milc_list



if __name__ == '__main__':
	#reference DNA sequence file,expample:genome_file
	reference_file = sys.argv[1]
	#sequence file to be calculated
	cds_file = sys.argv[2]
	#ref_codon_frequence = get_ref_codon_freq_inAA(reference_file)
	print(cal_MILC(reference_file,cds_file))