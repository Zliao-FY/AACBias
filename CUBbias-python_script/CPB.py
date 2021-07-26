#codon pair bias 
#delete start codon and stop codon
#delete special rare amino acid
# 20 standard amino acids

import sys
import re
from collections import Counter
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Data import CodonTable
from Bio.Seq import Seq
from pprint import pprint
from math import log
import pandas as pd
import numpy as np
import traceback

#pd.set_option('display.max_rows',500)
#pd.set_option('display.max_columns',500)
#pd.set_option('display.width',1000)


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
def check_cds_Right_cutoff(cds,cutoff=0):
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

def get_codon(cdss):
	#cdss is coding sequence
	#format: {'codon1':count1,'codon2':count2...}
	#get rid of stop codons
	len_cdss = len(cdss)
	return [cdss[i:i+3] for i in range(0,len_cdss,3)]

def get_AA(cdss):
	if type(cdss) == str:
		cdss = Seq(cdss)
	return list(str(cdss.translate(table='1',stop_symbol='*')))

def CPB(cds_codon_total,cds_aa_total,codon_pair_total,aa_pair_total):

	cds_codon_dict = Counter(cds_codon_total)
	cds_aa_dict = Counter(cds_aa_total)
	aa_pair_dict = Counter(aa_pair_total)

	C12_c = [1 for i in range(0,len(codon_pair_total))]
	codonPT = pd.DataFrame({'codon_pair':codon_pair_total,
							'C12_c':C12_c
							})
	codonPT = codonPT.groupby(by='codon_pair').count()
	C1 = []
	C2 = []
	A1 = []
	A2 = []
	A12 = []

	for c in codonPT.index:
		pair = c.split('_')

		C1.append(pair[0])
		A1.append(str(Seq(pair[0]).translate(table='1')))
		C2.append(pair[1])
		A2.append(str(Seq(pair[1]).translate(table='1')))

	C1_c = [cds_codon_dict[c] for c in C1]
	C2_c = [cds_codon_dict[c] for c in C2]
	A1_c = [cds_aa_dict[a] for a in A1]
	A2_c = [cds_aa_dict[a] for a in A2]

	codonPT['C1'] = C1
	codonPT['C2'] = C2
	codonPT['A1'] = A1
	codonPT['A2'] = A2
	codonPT['A12'] = codonPT['A1']+'_'+codonPT['A2']

	codonPT['C1_c'] = C1_c
	codonPT['C2_c'] = C2_c
	codonPT['A1_c'] = A1_c
	codonPT['A2_c'] = A2_c
	codonPT['A12_c'] = [aa_pair_dict[a_a] for a_a in codonPT['A12']]
	
	def fifter_A12(x):
		if '*' in x:
			return False
		else:
			return True
	codonPT = codonPT[codonPT['A12'].apply(fifter_A12)]

	codonPT['CPS'] = np.log(codonPT['C12_c'] / ((codonPT['C1_c']*codonPT['C2_c']) / ((codonPT['A1_c']*codonPT['A2_c'])) * codonPT['A12_c']))

	cpb = np.sum(codonPT['CPS']) / (len(codonPT)-1)
	codonPT.to_csv('sta-CPB-python.txt',sep='\t')
	return cpb



def cal_CPB(cds_file,cutoff=240,rm_start_stop_codons=True):
	cds_codon_total = []
	cds_aa_total = []
	codon_pair_total = []
	aa_pair_total = []
	for record in SeqIO.parse(cds_file,'fasta'):

		if not check_cds_Right_cutoff(record.seq,cutoff=cutoff):
			continue
			
		if rm_start_stop_codons:
			cds_seq = str(record.seq)[3:-3]
		else:
			cds_seq = str(record.seq)

		cds_codon = get_codon(str(cds_seq))
		cds_aa = get_AA(str(cds_seq))
		cds_codon_total.extend(cds_codon)
		cds_aa_total.extend(cds_aa)

		codon1 = cds_codon[:-1]
		codon2 = cds_codon[1:]
		aa1 = cds_aa[:-1]
		aa2 = cds_aa[1:]
		codon_pair = [codon1[i]+'_'+codon2[i] for i in range(0,len(codon1))]
		aa_pair = [aa1[i]+'_'+aa2[i] for i in range(0,len(aa1))]

		codon_pair_total.extend(codon_pair)
		aa_pair_total.extend(aa_pair)
	return CPB(cds_codon_total,cds_aa_total,codon_pair_total,aa_pair_total)




if __name__ == '__main__':
	cds_file_file = sys.argv[1]
	cpb = cal_CPB(cds_file_file)
	print(cpb)

