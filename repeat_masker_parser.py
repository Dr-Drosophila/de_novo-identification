"""
================================= FILE INFORMATION ========================================
File name:			repeat_masker_parser.py
Run using:			python3 ../repeat_masker_parser.py $repeat_masker_shortened_merged_bed $threshold $uclust_no_genes $uclust_postmatch $uclust_summary
Input(s):			rmask_no_genes_raw( Raw file that contains RepeatMasker data )
Output(s):			uclust_post_match( FASTA file that has aligned names next to RepeatModeler names ); uclust_summary( Summary file of aligned sequences (parsable) )
					<chimera_list>
Description:
	READS IN THE RepeatMasker OUTPUT TO BE ABLE TO ALIGN THEM AND SEE SEQUENCE ALIGNMENT SIMILARITIES TO THE RepeatMasker Library.

Author: 	Chinmay Rele
Date:		2019/01/25

"""

#
#    ,,
#    db                                                    mm
#                                                          MM
#  `7MM  `7MMpMMMb.pMMMb.  `7MMpdMAo.  ,pW"Wq.  `7Mb,od8 mmMMmm  ,pP"Ybd
#    MM    MM    MM    MM    MM   `Wb 6W'   `Wb   MM' "'   MM    8I   `"
#    MM    MM    MM    MM    MM    M8 8M     M8   MM       MM    `YMMMa.
#    MM    MM    MM    MM    MM   ,AP YA.   ,A9   MM       MM    L.   I8
#  .JMML..JMML  JMML  JMML.  MMbmmd'   `Ybmd9'  .JMML.     `Mbmo M9mmmP'
#                            MM
#                          .JMML.


import sys # to import a file from the commmand line
from os import system as bash				# to run bash from python
from Bio import SeqIO # Need Seq.IO for parsing and making FASTA file
import getopt # TO PASS IN INPUT AND OUTPUT PARAMETERS
import collections # TEMP TO SEE IF THERE ARE REPEATS IN ANY LISTS
import re # USED TO SPLIT STRINGS TO SEPARATE THE NAME AND ALIGNMENT OF THE REPEATMASKER LIBRARY ANNOTATION.
from pprint import pprint
from tqdm import tqdm as pbar



#
#      ,...                                           ,,
#    .d' ""                                   mm      db
#    dM`                                      MM
#   mMMmm  `7MM  `7MM  `7MMpMMMb.   ,p6"bo  mmMMmm  `7MM   ,pW"Wq.  `7MMpMMMb.  ,pP"Ybd
#    MM      MM    MM    MM    MM  6M'  OO    MM      MM  6W'   `Wb   MM    MM  8I   `"
#    MM      MM    MM    MM    MM  8M         MM      MM  8M     M8   MM    MM  `YMMMa.
#    MM      MM    MM    MM    MM  YM.    ,   MM      MM  YA.   ,A9   MM    MM  L.   I8
#  .JMML.    `Mbod"YML..JMML  JMML. YMbmd'    `Mbmo .JMML. `Ybmd9'  .JMML  JMML.M9mmmP'
#
#

if __name__ == "__main__":
	help_string = "python3 ../repeat_masker_parser.py \n\t -r <rmask_no_genes_short[in]> \n\t -t <threshold[num]> \n\t -u <uclust_no_genes[in]> \n\t -p <uclust_postmatch[out]> \n\t -s <summary_file[out]>"
	try:
		opts, args = getopt.getopt( sys.argv[1:],"hr:t:u:p:s:",["help", "rmask_no_genes_short=", "threshold=", "uclust_no_genes=", "uclust_postmatch=", "summary_file=" ] )
	except getopt.GetoptError:
		print( help_string )
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			print( help_string )
			sys.exit()
		elif opt in ("-r", "--rmask_no_genes_short"):
			# INPUT FILE THAT HAS ONLY THE NAME OF THE UCLUST DATA, THE START AND END ALIGNMENT OF THE UCLUST DATA NAD THE NAME OF THE LIBRARY ANNOTATION THAT HAS BEEN MERGED
			rmask_no_genes_short = arg
		elif opt in ("-t", "--threshold"):
			# THRESHOLD FOR ALIGNMENT MATCH
			threshold = arg
			threshold = threshold.strip()
			threshold = float(threshold)
		elif opt in ("-u", "--uclust_no_genes"):
			# PARSED UCLUST DATA THAT DOES NOT HAVE ANY GENES IN IT.
			uclust_no_genes = arg
		elif opt in ("-p", "--uclust_postmatch"):
			# UCLUST DATA THAT HAS BEEN SCREENED AND REPEATS IDENTIFIED.
			#  	REPEATS THAT REPEATMODELER MISSED HAVE NOW BEENIDENTIFIED.
			uclust_postmatch = arg
		elif opt in ("-s", "--summary_file"):
			# SUMMARY FILE
			# 	HAS THE SUMMARY OF ALL IDENTIFIED AND UNIDENTIFIED UCLUST SEQUENCES
			# 	UCLUST NAME, REPEATMASKER NAME, ALIGNMENT SCORE.
			summary_file = arg
	print()
	print( 'rmask_no_genes_short is:\t', rmask_no_genes_short )
	print( 'threshold is:\t\t\t', str(threshold) )
	print( 'uclust_no_genes is:\t\t', uclust_no_genes )
	print( 'uclust_postmatch is:\t\t', uclust_postmatch )
	print( 'summary_file is:\t\t', summary_file )

	print( "Hello world" )
	print(  )


#
#                            ,,
#                          `7MM
#                            MM
#   ,p6"bo   ,pW"Wq.    ,M""bMM   .gP"Ya
#  6M'  OO  6W'   `Wb ,AP    MM  ,M'   Yb
#  8M       8M     M8 8MI    MM  8M/////
#  YM.    , YA.   ,A9 `Mb    MM  YM.    ,
#   YMbmd'   `Ybmd9'   `Wbmd"MML. `Mbmmd'
#
#


# open the repeatmasker data, and find the sequence, what it aligned to, and the length of the alignment
rmask_data = []
with open( rmask_no_genes_short, "r" ) as rmask, open( "rmod_tomerge.bed", "w" ) as rmod:
	for line in rmask:
		line = line.strip().split()
		line[1] = int( line[1] )
		line[2] = int( line[2] )
		temp = [ line[0], line[2]-line[1] ]
		extra_temp = [ line[0].split( "_|_|_" )[0], str(line[1]), str(line[2]), line[0].split( "_|_|_" )[1] ]
		print( "\t".join(extra_temp), file=rmod )
		rmask_data.append( temp )
rmask_data.sort()

print( "len(rmask_data):", len(rmask_data) )

# merge the alignments only; distinct the aligned items
bash( "cat rmod_tomerge.bed | sort -k1,1 -k2,2n | bedtools merge -c 4 -o distinct > rmod_merged.bed" )

# "add" those alignments up and create a dictionary of ITEMS
putative_chimeras = {}
with open( "rmod_merged.bed", "r" ) as rmod_merge:
	for line in rmod_merge:
		line = line.strip().split()
		# print( int(line[2])-int(line[1]) )
		if line[0] not in putative_chimeras:
			putative_chimeras[ line[0] ] = [ int(line[2])-int(line[1]), line[3] ]
			continue
		if line[0] in putative_chimeras:
			putative_chimeras[ line[0] ][0] += int(line[2])-int(line[1])
			if line[3] not in putative_chimeras[ line[0] ][1]:
				putative_chimeras[ line[0] ][1] += ","
				putative_chimeras[ line[0] ][1] += line[3]

# pprint( putative_chimeras )

# bash( "rm rmod_merged.bed" )

# add bp of alignments that overlap (same RMOD TE, same motif)
rmask_collapse_overlaps = {}

for alignment in rmask_data:
	if alignment[0] not in rmask_data:
		rmask_collapse_overlaps[ alignment[0] ] = alignment[1]
	elif alignment[0] in rmask_data:
		rmask_collapse_overlaps[ alignment[0] ] += alignment[1]

print( "len(rmask_collapse_overlaps):", len(rmask_collapse_overlaps) )

# removes lesser/keeps better alignments
rmask_collapse_better = {}

for key, val in rmask_collapse_overlaps.items():
	temp = key.split( "_|_|_" )
	if temp[0] not in rmask_collapse_better:
		rmask_collapse_better[ temp[0] ] = [ temp[1], val ]
	elif temp[0] in rmask_collapse_better:
		temp2 = [ temp[1], val ]
		if rmask_collapse_better[ temp[0] ][1] > temp2[1]:
			rmask_collapse_better[ temp[0] ] = [ temp[1], val ]

print( "len(rmask_collapse_better):", len(rmask_collapse_better) )


# print( putative_chimeras[ "rnd-5_family-574#DNA/hAT-Ac" ] )

# sys.exit( "subtest" )

# list of unknown elements from uclust_no_genes
summary_file = [ [ "## RMOD", "aligned_to", "percent_align" ] ]

# list of shown chimeric sequences (sequences that matched significantly with more than one type of element)
list_of_chimeras = []

# check if the sequence in uclust_no_genes and if the alignment is high, then rename; otherwise, skip
with open( uclust_postmatch, "w" ) as postmatch:
	for record in SeqIO.parse(uclust_no_genes, "fasta"):
		put_chi = True
		truly_unknown = True
		legacy_id = record.id
		summary = [ record.id ]
		record.description = ""
		if record.id in rmask_collapse_better:
			alignment = rmask_collapse_better[ record.id ][1] / len( record.seq )
			if alignment > threshold:
				# if "Unknown" in record.id:
				# 	print( record.id )
				val = rmask_collapse_better[ record.id ]
				record.id = record.id.split("#")[0] + "#" + val[0]
				summary.extend([ val[0], str(alignment) ])
				put_chi = False
			truly_unknown = False
			# check for chimeras:
			chimeric_align = putative_chimeras[ legacy_id ][0] / len( record.seq )
			if chimeric_align > threshold and put_chi:
				# print( legacy_id, putative_chimeras[ legacy_id ], chimeric_align )
				record.id = legacy_id.split("#")[0] + "#Chimera/" + legacy_id.split("#")[1]
				list_of_chimeras.append( legacy_id )
			# checks if the element did not match RMASK, and thus unknown.
			elif chimeric_align < threshold:
				record.id = legacy_id.split( "#" )[0] + "#Unknown"
		SeqIO.write( record, postmatch, "fasta")
		print( file=postmatch )
		if len( summary ) == 1:
			summary.extend([ "unidentified", "-" ])
		summary_file.append( summary )


# write all the chimeras to a file:
with open( "chimeras.lst", "w" ) as chimeras:
	print( "##RMOD_seq\tRMASK_align", file=chimeras )
	for chimer in list_of_chimeras:
		alignments = putative_chimeras[ chimer ][1].split( "," )
		for i in range( len( alignments ) ):
			print( "{0}\t{1}".format( chimer, alignments[i] ), file=chimeras )

# pprint( summary_file )

bash( "echo 'This works'" )

sys.exit()


# ==============================================================================================
# CODE THE BELOW NEATLY; ALSO comment
# ==============================================================================================




pinko = []

count_in_post_match = 0
tmp_counter = 0
aligned_counter = 0

unknowns = []

with open( uclust_postmatch, "w" ) as f:
	for record in SeqIO.parse(uclust_no_genes, "fasta"):
		temp = []
		temp.append( record.id )
		tmp_counter+=1

		if "Unknown" in record.id:
			unknowns.append( record.id )

		for item in data:
			if record.id == item[0]:
				alignment = item[2]/len(record.seq)
				if alignment > threshold:
					# print( record.id, end="\t" )
					print( item[1][6:] )
					record.id = record.id + "::" + item[1][6:]
					temp.append( item[1] )
					temp.append( str(alignment) )
					aligned_counter+=1
		SeqIO.write( record, f, "fasta")
		print( file=f ) # prints a newline between sequences
		count_in_post_match+=1
		# print("temp: ", temp)
		if len(temp) == 3:
			pinko.append( temp )
		elif len(temp) == 1:
			temp.append("unidentified")
			temp.append("-")
			pinko.append(temp)


print( "len(unknowns):", len(unknowns) )
print( "len(list(set(unknowns)))", len(list(set(unknowns))) )

print("aligned_counter: ", aligned_counter)

print( "count_in_post_match = ", count_in_post_match )
print( "tmp_counter = ", tmp_counter )

with open( summary_file, "w" ) as f:
	for item in pinko:
		print( "\t".join(item), file=f )


	# if record.id in genes:
	# 	# IF IT IS NOT, THEN WRITES THE SEQUENCE TO THE FILE
	# 	print(record.format("fasta"), file=f)
	# 	adding_to_genes.append(record.id)
	# 	genes_count+=1

 #
 #                                               ,,
 #   mm                                          db
 #   MM
 # mmMMmm   .gP"Ya  `7Mb,od8 `7MMpMMMb.pMMMb.  `7MM  `7MMpMMMb.  `7MM  `7MM  ,pP"Ybd
 #   MM    ,M'   Yb   MM' "'   MM    MM    MM    MM    MM    MM    MM    MM  8I   `"
 #   MM    8M//////   MM       MM    MM    MM    MM    MM    MM    MM    MM  `YMMMa.
 #   MM    YM.    ,   MM       MM    MM    MM    MM    MM    MM    MM    MM  L.   I8
 #   `Mbmo  `Mbmmd' .JMML.   .JMML  JMML  JMML..JMML..JMML  JMML.  `Mbod"YML.M9mmmP'
 #
 #
