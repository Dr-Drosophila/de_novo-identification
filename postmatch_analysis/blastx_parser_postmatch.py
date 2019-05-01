"""
========= FILE INFORMATION =========================================================
File name:			blastx_parser_postmatch.py
Run using:			python3 blastx_parser_postmatch.py <blastx_out<>
Input(s):			<blastx_out>
Output(s):			-
Description:
	Calculates the coverage values of blastx file

Author: 	Chinmay Rele
Date:		2019/04/18
Version:	0.0.1

========= IMPORTS =================================================================
"""

from tqdm import tqdm as pbar				# to print progress bar on terminal
import time									# debug
from pprint import pprint 					# to print items to screen enatly
import sys									# input data
import operator								# to sort

"""
========= FUNCTIONS ===============================================================
"""

# sorting table by index
def sort_table(table, col=0):
    return sorted(table, key=operator.itemgetter(col))

"""
========= CODE ====================================================================
"""

# column assignments
# 	actual assignments are base-1
# 	Python requeres base-0
qid = 0
sid = 1
qstart = 2
qend = 3
sstart = 4
send = 5
qlen = 6

# alignments
# 		format = [ qid, sid, alignment_score ]
alignments = []

blastx_out = sys.argv[-1]

with open( blastx_out, "r" ) as blast_file:
	for line in pbar( blast_file, desc=blastx_out ):
		if "# " not in line:
			line = line.split()
			alignment = [ line[qid], line[sid], (float(line[send])-float(line[sstart]))/float(line[qlen]) ]
			alignments.append( alignment )

sorted_data = []

for row in sort_table( alignments, 2):
	item = [ row[0], row[1], str(row[2]) ]
	sorted_data.append( item )

# long alignments
long_align = []

# print items with large alignments
for item in pbar( alignments, desc="alignments" ):
	if item[2] >= 0.5:
		long_align.append( item )

with open( "postmatch_alignment_scores.tab", "w" ) as outfile:
	for item in pbar( sorted_data, desc="writing to ouput" ):
		print( "\t".join( item ), file=outfile )
