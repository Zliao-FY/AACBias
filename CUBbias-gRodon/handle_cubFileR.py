import sys
import re

def handle_CUB_Rfiles(cub_file):
	flag = 0
	file_suffix = '_cds_from_genomic.fna'
	with open(cub_file) as fileobj:
		for line in fileobj:
			line = line.strip()
			if flag == 0:
				if re.match(r'file_name',line):
					headline = line.replace('species_predict.','')
					print('#'+headline)
			elif re.match(r'file_name',line):
				continue
			else:
				dataline = line.replace(file_suffix,'')
				print(dataline)
			flag += 1


if __name__ == '__main__':
	cub_file = sys.argv[1]
	handle_CUB_Rfiles(cub_file)
	