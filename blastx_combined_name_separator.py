"""
================================= FILE INFORMATION ========================================
File name:			blastx_combined_name_separator.py
Run using:			python3 ../blastx_combined_name_separator.py $cutoff repeat_gene_alignments.bed $uclust_output $uclust_no_genes $uclust_genes
Input(s):			cutoff, shortened_merge_bed, uclust
Output(s):			uclust_genes, uclust_no_genes
Description:
	# 1. GOES THROUGH EACH OF THE ITEMS
	# 2. SEPARATES THE NAME OF THE REPEAT AND THE NAME OF THE GENE
	# 3. IF THERE ARE MULTIPLE PAIRS, THEN ADDS THE ALIGHNMENT SCORES AND RETURNS ONLY 1 PAIR
	# 4. RETURNS SEQUENCE OF REPEATS THAT ARE GENES AND NOT GENES (LIKELY)
	# 5. CONTAINS THE CUTOFF VALUE AS A VARIABLE. CAN BE CHANGED AT A LATER DATE IF NEED BE.

Author: 	Chinmay Rele
Date:		2019/01/25

===================================== IMPORTS =============================================
"""

import sys # to import a file from the commmand line
import getopt # TO PASS IN INPUT AND OUTPUT PARAMETERS
from Bio import SeqIO # Need Seq.IO for parsing and making FASTA file
import collections # TEMP TO SEE IF THERE ARE REPEATS IN ANY LISTS
from pprint import pprint

"""
========= FUNCTIONS ===============================================================
"""

## get information from bash
if __name__ == "__main__":
	help_string = "python3 ../blastx_combined_name_separator.py \n\t -c <cutoff> \n\t -s <shortened_mergeBed> \n\t -u <uclust> \n\t -p <uclust_no_genes> \n\t -g <uclust_genes>"
	try:
		opts, args = getopt.getopt( sys.argv[1:],"hc:s:u:p:g:",["help", "cutoff=", "shortened_mergeBed=", "uclust=", "uclust_no_genes=", "uclust_genes="] )
	except getopt.GetoptError:
		print( help_string )
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			print( help_string )
			sys.exit()
		elif opt in ("-i", "--intput_file"):
			input_file = arg
		elif opt in ("-s", "--shortened_mergeBed"):
			# THE FILE THAT HAS HAS BEEN RUN THROUGH BEDTOOLS MERGE. HAS ALIGNMENT SCORES FOR REPEATS.
			shortened_mergeBed = arg
		elif opt in ("-c", "--cutoff"):
			# THE CUTOFF VALUE:
			cutoff = arg
			cutoff = cutoff.strip()
			cutoff = float(cutoff)
		elif opt in ("-u", "--uclust"):
			# THE RAW UCLUST DATA THAT HAS ALL THE REPEATS IN IT. MAY HAVE GENES ALSO
			uclust = arg
		elif opt in ("-p", "--uclust_no_genes"):
			# PARSED UCLUST DATA THAT DOES NOT HAVE ANY GENES IN IT.
			uclust_no_genes = arg
		elif opt in ("-g", "--uclust_genes"):
			# PARSED UCLUST DATA THAT DOES NOT HAVE ANY GENES IN IT.
			uclust_genes = arg

	print( 'cutoff is:\t\t', str(cutoff) )
	print( 'shortened_mergeBed is:\t', shortened_mergeBed )
	print( 'uclust is:\t\t', uclust )
	print( 'uclust_no_genes is:\t', uclust_no_genes )
	print( 'uclust_genes is:\t', uclust_genes )
	print(  )

"""
======================================== CODE =============================================
"""


# create a list of annotations
gene_annotations = {}
annotations = []
with open( shortened_mergeBed, "r" ) as merged:
	for line in merged:
		line = line.strip().split()
		name = line[0].split( "_|_|_|_|_" )[0]
		annotations.append( name )
		alignment = float( line[1] )
		# print( name, alignment )
		if name not in gene_annotations:
			gene_annotations[ name ] = [ alignment, line[0].split( "_|_|_|_|_" )[1] ]
		if name in gene_annotations:
			if gene_annotations[ name ][0] < alignment:
				gene_annotations[ name ] = [ alignment, line[0].split( "_|_|_|_|_" )[1] ]
			elif gene_annotations[ name ][0] >= alignment:
				continue

# create set of annotations to loop through to see > cutoff overlap.
annotations = list(set( annotations ))

for rep in annotations:
	if gene_annotations[ rep ][0] < cutoff:
		del gene_annotations[ rep ]
	# print( gene_annotations[ rep ] )

with open( uclust_no_genes, "w" ) as no_genes, open( uclust_genes, "w" ) as genes:
	for record in SeqIO.parse(uclust, "fasta"):
		if record.id in gene_annotations:
			record.id += " ~= " + gene_annotations[ record.id ][1]
			SeqIO.write( record, genes, "fasta")
		elif record.id not in gene_annotations:
			SeqIO.write( record, no_genes, "fasta")


print( "\n\n\tTEST\n\n" )
sys.exit()


# ==============================================================================================
# CODE THE BELOW NEATLY; ALSO comment
# ==============================================================================================



# OPENING THE DATA
data = open(shortened_mergeBed, 'r')
data_len = len(open(shortened_mergeBed, 'r').readlines())
# print("len(shortened_mergeBed) = " + str(data_len) )

# INITIALIZING EMPTY LISTS TO HOLD DATA
raw_lines = []

# looping throuch each datapoint, and combining the names of the RepeatModeler identifier and the name of the identified repeat.
for i in range(data_len):
	line = data.readline()
	line = line.replace("_|_|_|_|_", " ")
	read = line[0:-1].split('\t')
	temp = read[0].split(' ')
	raw_lines.append(temp)


for item in raw_lines:
	item[2] = float(item[2])

# print( "raw_lines[1:3]: ", raw_lines[1:3] )

# for item in raw_lines:
# 	print(item)
# print()

# print( "len(raw_lines): ", len(raw_lines) )

# if two items are the same, get rid of one of them
for item1 in raw_lines:
	for item2 in raw_lines:
		if item1[0] == item2[0] and item1[1] == item2[1] and raw_lines.index(item1) != raw_lines.index(item2) :
			# print(item1, "\t", item2)
			item1[2] += item2[2]
			raw_lines.remove(item2)

# print( "len(raw_lines): ", len(raw_lines) )

lines = raw_lines

genes = []
for item in lines:
	temp = []
	if item[2] > cutoff:
		# print( "cutoff, item[2]: ", cutoff, item[2] )
		temp.append( item[0] )
		temp.append( item[1] )
		temp.append( item[2] )
		genes.append( temp )

pinko = set( [ item[0] for item in genes ] )
print( "pinko: ", pinko )
print( "len(pinko): ", len(pinko) )

print( "genes[0:3]: ", genes[0:3] )
print( "len(genes): ", len(genes) )

names = []

# remove duplicates
for i in range(1000):
	for item1 in genes:
		for item2 in genes:
			if item1[0] == item2[0] and genes.index(item1) != genes.index(item2):
				if item1[2] >= item2[2]:
					genes.remove(item2)
				elif item1[2] < item1[2]:
					genes.remove(item1)

# print()
print( "len(genes): ", len(genes) )

# if score is less that cutoff, then add to genes, else add to no_genes
not_genes = []
for item in lines:
	if item[2] <= cutoff:
		not_genes.append(item[0])

print( "len(not_genes): ", len(list(set(not_genes))) )
print( "not_genes: ", list(set(not_genes)) )

# find the length of the uclust files (the number of sequences)
uclust_len = 0
for record in SeqIO.parse(uclust, "fasta"):
	uclust_len+=1

print( "uclust_len: " + str(uclust_len) )

# IF THE SEQUENCE IS LIKELY A GENE, IT WILL NOT BE COPIED INTO THE FINAL FILE
# OPENS A FILE READY TO WRITE
adding_to_no_genes = []
no_genes_count = 0
with open(uclust_no_genes, 'w') as f:
	for record in SeqIO.parse(uclust, "fasta"):
		# CHECKS IF THE SEQUENCE ID IS NOT A likely_gene:
		if record.id in not_genes:
			# IF IT IS NOT, THEN WRITES THE SEQUENCE TO THE FILE
			print( record.format("fasta"), file=f )
			adding_to_no_genes.append(record.id)
			no_genes_count+=1

print(no_genes_count)

adding_to_genes = []
genes_count = 0
with open(uclust_genes, 'w') as f:
	for record in SeqIO.parse(uclust, "fasta"):
		# CHECKS IF THE SEQUENCE ID IS A likely_gene:
		for pair in genes:
			if record.id == pair[0]:
				record.id = record.id + " ~= " + pair[1]
				# print( record.id )
				print(record.format("fasta"), file=f)
				adding_to_no_genes.append(record.id)
				genes_count += 1

# print( "adding_to_genes: ", adding_to_genes )
# print( "genes_count: ", str(genes_count) )
# print( "adding_to_no_genes: ", adding_to_no_genes )
# print( "no_genes_count: ", str(no_genes_count) )


"""
---------------------------
	END OF FILE
---------------------------
"""
