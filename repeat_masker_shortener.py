"""
================================= FILE INFORMATION ========================================
File name:			repeat_masker_shortener.py
Run using:			python3 ../repeat_masker_shortener.py $repeatmasker_gff $repeatmasker_shortened_gff
Input(s):			rmask_no_genes_raw( Raw file that contains RepeatMasker data )
					-n <repeatmasker_out> repeatmasker.out file; to find anotations
Output(s):			rmask_no_genes_short( Shortened version of RepeatMasker data { seq_id+"_|_|_"+annotation_id, start_alignment, end_alignment } )
					pkl file of unknown dictionary
Description:
	Combines data from RepeatMasker into a single sequenceid so we can run mergeBed on it.

Author: 	Chinmay Rele
Date:		2019/01/25


# # GFF FORMAT
# 00	seqname - name of the chromosome or scaffold; chromosome names can be given with or without the 'chr' prefix. Important note: the seqname must be one used within Ensembl, i.e. a standard chromosome name or an Ensembl identifier such as a scaffold ID, without any additional content such as species or assembly. See the example GFF output below.
# 01	source - name of the program that generated this feature, or the data source (database or project name)
# 02	feature - feature type name, e.g. Gene, Variation, Similarity
# 03	start - Start position of the feature, with sequence numbering starting at 1.
# 04	end - End position of the feature, with sequence numbering starting at 1.
# 05	score - A floating point value.
# 06	strand - defined as + (forward) or - (reverse).
# 07	frame - One of '0', '1' or '2'. '0' indicates that the first base of the feature is the first base of a codon, '1' that the second base is the first base of a codon, and so on..
# 08	attribute - A semicolon-separated list of tag-value pairs, providing additional information about each feature.


===================================== IMPORTS =============================================
"""

import sys # to import a file from the commmand line
import getopt # TO PASS IN INPUT AND OUTPUT PARAMETERS
from Bio import SeqIO # Need Seq.IO for parsing and making FASTA file
import collections # TEMP TO SEE IF THERE ARE REPEATS IN ANY LISTS
import re # USED TO SPLIT STRINGS TO SEPARATE THE NAME AND ALIGNMENT OF THE REPEATMASKER LIBRARY ANNOTATION.
import pickle # to export the dictionary to a file
from pprint import pprint
from tqdm import tqdm as pbar

"""
========= FUNCTIONS ===============================================================
"""

if __name__ == "__main__":
	# print( "from main" )
	help_str = "python3 ../repeat_masker_shortener.py -i <rmask_no_genes_raw> -o <rmask_no_genes_short>"
	try:
		opts, args = getopt.getopt( sys.argv[1:],"hi:o:",["help", "input=", "repeatmasker_out=" ] )
	except getopt.GetoptError:
		print( help_str )
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			print( help_str )
			sys.exit(2)
		elif opt in ("-i", "--input"):
			rmask_no_genes_raw = arg
		elif opt in ("-o", "--output"):
			rmask_no_genes_short = arg

	print( 'no main rmask_no_genes_raw is:', rmask_no_genes_raw )
	print( 'no main rmask_no_genes_short is:', rmask_no_genes_short )
	print()


"""
======================================== CODE =============================================
"""

# create a dictionary that has annotated repeats
annotations = []
with open( rmask_no_genes_raw, "r" ) as in_file:
	for line in in_file:
		if "SW" not in line and "score" not in line and len(line) >= 14:
			line = line.strip().split()
			temp = [ line[4], line[5], line[6], line[10], line[9] ]
			alignment = []
			rename = temp[0] + "_|_|_" + temp[3] + "::" + temp[4]
			alignment.append( rename )
			alignment.append( temp[1] )
			alignment.append( temp[2] )
			annotations.append(alignment)

# print the annotations to the output files

with open( rmask_no_genes_short, "w" ) as outfile:
	for item in pbar( annotations, desc="writing annotations" ):
		print( "\t".join( item ), file=outfile )

sys.exit( "test" )
# ==============================================================================================
# CODE THE BELOW NEATLY; ALSO comment
# ==============================================================================================

raw_data = open(rmask_no_genes_raw, "r")
len_raw_data = len( open(rmask_no_genes_raw, "r").readlines() )




dict_annotations = {}
for item in annotations:
	if item[0] not in list( dict_annotations.keys() ):
		dict_annotations[ item[0] ] = [ item[1] ]
	elif item[0] in list( dict_annotations.keys() ):
		if item[1] not in dict_annotations[ item[0] ]:
			dict_annotations[ item[0] ].append( item[1] )

[ print(key) for key in list(dict_annotations.keys()) ]



# REPEATMASKER GFF FORMAT INDEXES
indx_seq_name = 0
indx_source = 1
indx_feature = 2
indx_start = 3
indx_end = 4
indx_score = 5
indx_strand = 6
indx_frame = 7
indx_attribute = 8

# TEMP LIST TO HOLD RAW DATA IN.
# 	WILL BE COMPUTED ON LATER
data = []

for i in range( len_raw_data ):
	line = raw_data.readline()[:-1]
	if "rnd-" in line:
		array = line.split("\t")
		temp = []
		temp.append( array[indx_seq_name] )
		temp.append( array[indx_source] )
		temp.append( array[indx_feature] )
		temp.append( int(array[indx_start]) )
		temp.append( int(array[indx_end]) )
		temp.append( float(array[indx_score]) )
		temp.append( array[indx_strand] )
		temp.append( array[indx_frame] )
		temp.append( array[indx_attribute] )
		data.append( temp )

for item in data:
	temp = item[0].split("#")
	TE_class = item[8].split("\"")[1][6:]
	# print( TE_class )
	printi = False
	if temp[1] == "Unknown":
		print(item)
		printi = True
		for key, value in dict_annotations.items():
			if TE_class in value:
				temp[1] = key
				print(key)
	item[0] = "#".join(temp)
	print(item)

# SHORTENED VERSION OF DATA THAT HAS THE SEQUENCE NAME, THE ALIGNMENT LENGTH AND THE REPEATMASKER LIBRARY ANNOTATION
shortened_data = []

for item in data:
	temp = []
	s = item[8]
	result = re.search('\"(.*)\"', s).group(1)
	temp.append( item[0] + "_|_|_" + result )
	temp.append( str(item[3]) )
	temp.append( str(item[4]) )
	shortened_data.append( temp )

# # TEST TO SEE WHAT HAS BEEN ADDED TO shortened_data
# print( "shortened_data[0]: ", shortened_data[0] )
# print( "shortened_data[1]: ", shortened_data[1] )
# print( "shortened_data[2]: ", shortened_data[2] )

with open(rmask_no_genes_short, "w") as f:
	for item in shortened_data:
		string_to_add = "\t".join(item)
		print( string_to_add, file=f )


"""
---------------------------
	END OF FILE
---------------------------
"""
