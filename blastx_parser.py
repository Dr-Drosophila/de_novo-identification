"""
================================= FILE INFORMATION ========================================
File name:			blastx_parser.py
Run using:			python3 ../blastx_parser.py -s <sorted_file> -u <uclust_output> -b <bedtools_input>
Input(s):			$uclust_output
Output(s):			$bedtools_input
Description:
	Takes in uncommented BLASTX data and compares to UCLUST data.
	Merges the names of the query and subject so can use [bedtools merge] to merge overlapping alignments.

Author: 	Chinmay Rele
Date:		2019/01/25

===================================== IMPORTS =============================================
"""

import sys # to import a file from the commmand line
import getopt # TO PASS IN INPUT AND OUTPUT PARAMETERS
from Bio import SeqIO # Need Seq.IO for parsing and making FASTA file

"""
========= FUNCTIONS ===============================================================
"""

if __name__ == "__main__":
	try:
		opts, args = getopt.getopt( sys.argv[1:],"hs:u:b:",["sorted_blast=","uclust_output=","bedtools_input="] )
	except getopt.GetoptError:
		print( 'python3 blastx_parser.py -s <blastx_file> -u <uclust_output> -b <bedtools_input>' )
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			print( 'python3 blastx_parser.py -s <blastx_file> -u <uclust_output> -b <bedtools_input>' )
			sys.exit()
		elif opt in ("-s", "--sorted_blast"):
			# ENCAPSULATE THE SORTED BLAST FILE
			blastx_file = arg
		elif opt in ("-u", "--uclust_output"):
			# ENCAPSULATE THE UCLUSTED DATA
			uclust = arg
		elif opt in ("-b", "--bedtools_input"):
			# ENCAPSULATE THE OUTPUT THAT HAS TO PUT INTO BEDTOOLS
			bedtools_merge = arg
	print( 'no main SORTED BLAST is:', blastx_file )
	print( 'no main UCLUST OUTPUT is:', uclust )
	print( 'no main BEDTOOLS INPUT is:', bedtools_merge )
	print(  )

"""
======================================== CODE =============================================
"""

# OPENING THE DATA
data = open(blastx_file, 'r')
data_len = len(open(blastx_file, 'r').readlines())

# INITIALIZING EMPTY LISTS TO HOLD DATA
reads = []
likely_genes = []

# ENTERING DATA FROM BLAST INTO A LIST
for i in range(data_len):
	line = data.readline()
	read = line[0:-1].split('\t')
	reads.append(read)

# needed = [["query_id+subjectid", "alignment length", "q. start", "q. end", "query length"]]
needed = []

for item in reads:
	temp = []
	# following indeces based on GFF format.
	name = item[0] + "_|_|_|_|_" + item[1]
	q_start = item[2]
	q_end = item[3]
	s_start = item[4]
	s_end = item[5]
	q_len = item[6]
	temp.append(name)
	if float(q_start) > float(q_end):
		temp.append(q_end)
		temp.append(q_start)
	elif float(q_start) < float(q_end):
		temp.append(q_start)
		temp.append(q_end)
	temp.append(q_len)
	needed.append(temp)

# IF THE SEQUENCE IS LIKELY A GENE, IT WILL NOT BE COPIED INTO THE FINAL FILE
# OPENS A FILE READY TO WRITE
with open(bedtools_merge, 'w') as f:
	for item in needed:
		print('\t'.join(item), file=f)

"""
---------------------------
	END OF FILE
---------------------------
"""
