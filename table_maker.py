"""
================================= FILE INFORMATION ========================================
File name:			table_maker.py
Run using:			python3 ../table_maker.py consensi.fa.classified $uclust_postmatch $family_table
Input(s):			$consensi_classified( File that has BLASTX data. Minor editing to remove commented lines. ); $uclust_post_match( Raw UCLUSTed Dmel data. )
Output(s):			$family_table.txt( Dataset \t Family \t Frequency )
Description:
	Takes in consensi.fa.cassified and $uclust_post_match and finds the number of counts of the repeats in each family.
	Will be done for RepeatModeler data as well as RepeatMasker data.

Author: 	Chinmay Rele
Date:		date

===================================== IMPORTS =============================================
"""

import sys # to import a file from the commmand line
import getopt # TO PASS IN INPUT AND OUTPUT PARAMETERS
from Bio import SeqIO # Need Seq.IO for parsing and making FASTA file
import collections # TEMP TO SEE IF THERE ARE REPEATS IN ANY LISTS
from collections import Counter # Counts frequencies


"""
========= FUNCTIONS ===============================================================
"""

if __name__ == "__main__":
	# print( "from main" )
	try:
		opts, args = getopt.getopt( sys.argv[1:],"hc:p:f:",["help", "consensi=", "post_match=", "family_table="] )
	except getopt.GetoptError:
		print( 'python3 blastx_parser.py -c <consensi> -p <post_match> -f <family_table>' )
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			print( 'python3 blastx_parser.py -c <consensi> -p <post_match> -f <family_table>' )
			sys.exit()
		elif opt in ("-c", "--consensi"):
			consensi = arg
		elif opt in ("-p", "--post_match"):
			post_match = arg
		elif opt in ("-f", "--family_table"):
			family_table = arg
	print( 'no main input_file is:', consensi )
	print( 'no main post_match is:', post_match )
	print( 'no main family_table is:', family_table )


"""
======================================== CODE =============================================
"""

final_data = []

r_modeler = []
for record in SeqIO.parse(consensi, "fasta"):
	r_modeler.append(record.id)

temp = []

for item in r_modeler:
	item = item.split("#")[1]
	temp.append(item)

modeler_dict = Counter(temp)
modeler_list = []
for key,value in modeler_dict.items():
	item = ["RepeatModeler", key, value]
	modeler_list.append(item)

for item in modeler_list:
	item[2] = str( item[2] )
	temp = "\t".join( item )
	final_data.append( temp )

# made a list of repeat_modeler_data

r_masker = []

for record in SeqIO.parse(post_match, "fasta"):
	r_masker.append(record.id)

temp = []

for item in r_masker:
	item = item.split("#")[1]
	item = item.split("::")[-1]
	temp.append(item)

masker_dict = Counter(temp)
masker_list = []
for key,value in masker_dict.items():
	item = ["RepeatMasker", key, value]
	masker_list.append(item)

print(len(masker_list))

for item in masker_list:
	item[2] = str( item[2] )
	temp = "\t".join( item )
	final_data.append( temp )

print( masker_list[0:3] )

with open( family_table, "w" ) as f:
	for item in final_data:
		f.write( item + "\n")

"""
---------------------------
	END OF FILE
---------------------------
"""
