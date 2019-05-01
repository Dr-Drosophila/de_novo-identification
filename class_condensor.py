"""
========= FILE INFORMATION =========================================================
File name:			class_condensor.py
Run using:			class_condensory.py -a <assembly size> -i <input GFF> -o <output file> -s <species>
Input(s):			-i <input GFF (file) [uclust_no_genes.fasta.out (within repeatMasker directory)]>
Output(s):			-o <output file>
Description:
	ouput file format

	Class 	Repeat Count  	Total len of class in genome	Percent composition of assembly.

Author: 	Chinmay Rele
Date:		2019/02/15

========= IMPORTS =================================================================
"""

import sys # to import a file from the commmand line
import getopt # TO PASS IN INPUT AND OUTPUT PARAMETERS
from Bio import SeqIO # Need Seq.IO for parsing and making FASTA file
import collections # TEMP TO SEE IF THERE ARE REPEATS IN ANY LISTS
import re # USED TO SPLIT STRINGS TO SEPARATE THE NAME AND ALIGNMENT OF THE REPEATMASKER LIBRARY ANNOTATION.
import math # for absolute methods
import pickle # to export the dictionary to a file
import pprint # to print dictionaries

"""
========= FUNCTIONS ===============================================================
"""

if __name__ == "__main__":
    # print( "from main" )
    try:
        opts, args = getopt.getopt( sys.argv[1:],"hi:o:s:",["help", "input_file=", "output_file=", "species_short=" ] )
    except getopt.GetoptError:
        print( 'python3 class_condensory.py -s <species_short> -i <input GFF> -o <output file> -s <species_short>' )
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print( 'python3 class_condensory.py -s <species_short> -i <input GFF> -o <output file>' )
            sys.exit()
        elif opt in ("-i", "--input_file"):
            input_file = arg
        elif opt in ("-o", "--output_file"):
            output_file = arg
        elif opt in ("-s", "--species_short"):
            species_short = arg
    print( 'input_file is:', input_file )
    print( 'output_file is:', output_file )
    print( 'species_short is:', species_short )
    print(  )

"""
========= CODE ====================================================================
"""

# pickle_in = open("dict.pickle","rb")
# dict_annotations = pickle.load(pickle_in)

# print(dict_annotations.keys())

# create an empty list to hold most of the data
gff = []

for line in open( input_file, "r" ):
	if "##" not in line:
		line = line.strip()
		temp = []
		line = line.split("\t")
		for item in line:
			if len(item) != 0:
				temp.append( item )
		gff.append(temp)

gff = gff[3:]

# repeat_info [nested list]: condensed information about each identified repeat
# each sublist:
# 		[ repeat class, length ]

repeat_info = []

gff = gff[:15]

for item in gff:
	temp = []
	temp.append( item[10] )
	item[11] = item[11].replace( "(", "" ).replace( ")", "" )
	# finding the length of the repeat (using absolute value because might be on reverse strand)
	temp.append( math.fabs(int(item[6])-int(item[5])) )
	repeat_info.append( temp )


# repeat_count [dictionary]: of the repeat and its count
repeat_count = {}

# repeat_bp [dictionary]: of the repeat and the number of base pairs it accomodates
repeat_bp = {}

for item in repeat_info:
	if item[0] not in repeat_count:
		repeat_count[ item[0] ] = 1
		repeat_bp[ item[0] ] = item[1]
	elif item[0] in repeat_count:
		repeat_count[ item[0] ] += 1
		repeat_bp[ item[0] ] += item[1]

# pprint.pprint( repeat_count )
# pprint.pprint( repeat_bp )

repeat_classes = []
for item in repeat_count:
	repeat_classes.append( item )

temp = []
for key,val in repeat_bp.items():
	temp.append(val)
sum_bp = sum( temp )

test = []


# this repeat_info[] has nothing to do with the older list.
# we do not need to store the older list anymore; overwriting
repeat_info = []
for item in repeat_classes:
	temp = []
	temp.append( item )
	temp.append( repeat_count.get( item ) )
	temp.append( repeat_bp.get( item ) )
	temp.append( temp[2]/sum_bp )
	test.append( temp[2]/sum_bp )
	repeat_info.append( temp )

print(repeat_info)

# need to write the list to an output
with open( output_file, "w" ) as output:
	# header separated by "##"
	print( "## Frequency of classes within assembly.\n## Encapsulates class, frequency, base pairs and the percent the class covered within all identified TEs (out of 1).\n##\n## Species\tClass\tcount\tbps\tpercent_of_TEs", file=output )
	for item in repeat_info:
		temp = []
		temp.append( species_short )
		temp.append( item[0] )
		temp.append( str(item[1]) )
		temp.append( str(item[2]) )
		temp.append( str(item[3]) )
		print( "\t".join(temp), file=output )

"""
---------------------------
	END OF FILE
---------------------------
"""
