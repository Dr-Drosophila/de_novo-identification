"""
========= FILE INFORMATION =========================================================
File name:			assembly_summary.py
Run using:			python3 assembly_summary.py -a <assembly> -s <species_short> -o <output_file>
Input(s):			-a <assembly>
Output(s):			-s <species_short>
					-o <output_file>
Description:
	finds the summary of the assembly:
		Total size
		Total contigs
		average contig size
		N_50
	NOTE: Need to echo -e "species\tassembly_size\tcontig_count\tavg_contig_size\tN_50" to output file before starting running.

Author: 	Chinmay Rele
Date:		2019/03/14
Version:	0.0.001

========= IMPORTS =================================================================
"""

import sys # to import a file from the commmand line
import getopt # TO PASS IN INPUT AND OUTPUT PARAMETERS
from Bio import SeqIO # Need Seq.IO for parsing and making FASTA file
import statistics # to calculate statistics such as mean

"""
========= FUNCTIONS ===============================================================
"""

if __name__ == "__main__":
	string = "python3 assembly_summary.py \n\t-a <assembly> \n\t-s <species_short> \n\t-o <output_file>"
	try:
		opts, args = getopt.getopt( sys.argv[1:],"ha:s:o:",["assembly=","species_short=","output_file="] )
	except getopt.GetoptError:
		print( string )
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			print( string )
			sys.exit()
		elif opt in ("-a", "--assembly"):
			assembly = arg
		elif opt in ("-s", "--species_short"):
			species_short = arg
		elif opt in ("-o", "--output_file"):
			output_file = arg
	print( 'assembly is:\t\t', assembly )
	print( 'species_short is:\t', species_short )
	print( 'output_file is:\t\t', output_file )
	print(  )

"""
========= CODE ====================================================================
"""


# Total size
assembly_size = 0

# Total contigs
contig_count = 0

# average contig size
contigs_len = []

# N_50

for record in SeqIO.parse( assembly, 'fasta' ):
	assembly_size += len( record.seq )
	contig_count += 1
	contigs_len.append( len(record.seq) )

# average contig size
avg_contig_size = statistics.mean( contigs_len )

contigs_len.sort(reverse=True)

n_50 = None
n_50_stop = 0
for i in range( len(contigs_len) ):
	if n_50_stop > assembly_size/2:
		break
	elif n_50_stop <= assembly_size/2:
		n_50_stop += contigs_len[i]
		n_50 = contigs_len[i]

final_output = [ species_short, str(assembly_size), str(contig_count), str(avg_contig_size), str(n_50) ]
with open( output_file, "a" ) as out_file:
	print( "\t".join(final_output), file=out_file )
