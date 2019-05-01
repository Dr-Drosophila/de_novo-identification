"""
========= FILE INFORMATION =========================================================
File name:			class_name_condensory.py
Run using:			class_name_condensory.py -i <input_file> -o <output_file>
Input(s):			-i <input_file> tab-delimited file that has all species and all unaltered class names
Output(s):			-o <output_file> tab-delimited file that has all species and all altered/condensed class names
Description:
	Dictionarythat creates and shortens names, but retains values.

Author: 	Chinmay Rele
Date:		date

========= IMPORTS =================================================================
"""

import sys # to import a file from the commmand line
import getopt # TO PASS IN INPUT AND OUTPUT PARAMETERS
from Bio import SeqIO # Need Seq.IO for parsing and making FASTA file
import collections # TEMP TO SEE IF THERE ARE REPEATS IN ANY LISTS
import re # USED TO SPLIT STRINGS TO SEPARATE THE NAME AND ALIGNMENT OF THE REPEATMASKER LIBRARY ANNOTATION.
from math import fabs # for absolute methods


"""
========= FUNCTIONS ===============================================================
"""

## get information from bash
if __name__ == "__main__":
	try:
		opts, args = getopt.getopt( sys.argv[1:],"hi:o:",["help", "input_file=", "output_file="] )
	except getopt.GetoptError:
		print( 'python3 class_condensory.py  -i <input_file> -o <output file> ' )
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			print( 'python3 class_condensory.py -i <input_file> -o <output file>' )
			sys.exit()
		elif opt in ("-i", "--intput_file"):
			input_file = arg
		elif opt in ("-o", "--output_file"):
			output_file = arg

	print( 'input_file is:', input_file )
	print( 'output_file is:', output_file )
	print(  )



"""
========= CODE ====================================================================
"""

# read in files and create a list of classes.

uncondensed_classes = []
full_data = []
with open( input_file, "r" ) as input:
	for line in input:
		if "##" not in line : #and "Unknown" not in line
			full_data.append( line.strip().split("\t") )
			uncondensed_classes.append( line.strip().split("\t")[1] )

# get rid of all classes that are uncertain (contain a ? in name)
for i in range(5):
	for item in full_data:
		if "?" in item[1]:
			full_data.remove( item )
for i in range(5):
	for item in uncondensed_classes:
		if "?" in item:
			uncondensed_classes.remove( item )

## manual creation of dictionary
## for custom disctionaries
# class_dictionary = {
# 	"DNA": [],
# 	"DNA/hAT": [],
# 	"DNA/CMC": [],
# 	"DNA/Kolobok": [],
# 	"DNA/Maverick": [],
# 	"DNA/Merlin": [],
# 	"DNA/Mule": [],
# 	"DNA/P": [],
# 	"DNA/PIF": [],
# 	"DNA/PiggyBac": [],
# 	"DNA/Sola": [],
# 	"DNA/TcMar": [],
# 	"DNA/Zator": [],
# 	"LINE/CR1": [],
# 	"LINE/I": [],
# 	"LINE/Penelope": [],
# 	"LINE/L2": [],
# 	"LINE/R1": [],
# 	"LINE/R2": [],
# 	"LINE/RTE": [],
# 	"Low_complexity": [],
# 	"LTR/Copia": [],
# 	"LTR/Gypsy": [],
# 	"LTR/Pao": [],
# 	"Other": [],
# 	"RC/Helitron": [],
# 	"Satellite": [],
# 	"Simple_repeat": [],
# 	# "": [], # to add more if need be
# 	"Unknown": []
# }
# for item in uncondensed_classes:
# 	tmp = item.split( "-" )
# 	if tmp[0] in class_dictionary:
# 		if item not in class_dictionary[tmp[0]]:
# 			class_dictionary[tmp[0]].append( item )


## automatic creation of dictionary
class_dictionary = {}
for item in uncondensed_classes:
	tmp = item.split( "/" )
	if tmp[0] not in class_dictionary:
		class_dictionary[tmp[0]] = [ item ]
	elif tmp[0] in class_dictionary:
		if item not in class_dictionary[tmp[0]]:
			class_dictionary[tmp[0]].append( item )

classes = [ key for key in class_dictionary.keys() ]

if "Satellite" in classes and "Simple_Repeat" in classes:
	reps = class_dictionary["Satellite"] + class_dictionary["Simple_Repeat"]
	del class_dictionary["Satellite"]
	del class_dictionary["Simple_Repeat"]
	class_dictionary["Satellite/Simple_Repeat"] = reps

# [ print (k, v) for k, v in class_dictionary.items() ]


# create a list of the final data so if need to edit downstream, can
final_data = []
for row in full_data:
	sat_SR_count = 0
	temp = []
	temp.append( row[0] )
	unknown = False
	for compressed, uncompressed in class_dictionary.items():
		for item in uncompressed:
			if row[1] == item:
				temp.append(compressed)
	temp.append( str(fabs(float((row[2])))) )
	temp.append( str(fabs(float((row[3])))) )
	temp.append( str(fabs(float((row[4])))) )
	final_data.append(temp)

for row in final_data:
	if row[0] == "dere" and "Unknown" in row[1]:
		print(row)
print()
# split Rmask annotated repeats
for item in final_data:
	if "::" in item[1]:
		item[1] = item[1].split("::")[0]

for row in final_data:
	if row[0] == "dbia" and "Unknown" in row[1]:
		print(row)

# finds the total bp of Sat/SR sequences in assembly
sat_SR_count = {}
for row in final_data:
	if row[0] not in list( sat_SR_count.keys() ):
		sat_SR_count[ row[0] ] = 0
	if "Satellite/Simple_Repeat" in row:
		sat_SR_count[ row[0] ] += int(float( row[3] ))

# print( sat_SR_count )

# accounts for unknown sequences within the rows and assigns to satellite/SR
for row in final_data:
	if "Unknown" in row:
		temp = int(float( row[3] ))
		if temp > sat_SR_count[ row[0] ]:
			row[3] = str(fabs(temp - sat_SR_count[ row[0] ]))

temp_dic = {}
for item in final_data:
	name_rep = item[0] + "_" + item[1]
	if name_rep not in list( temp_dic.keys() ):
		temp_dic[ name_rep ] = 0
	if name_rep in list( temp_dic.keys() ):
		temp_dic[ name_rep ] += int(float( item[3] ))

# for key, value in temp_dic.items():
# 	print( key, "\t", value )

# find total percent of of each repeat within data
species_reps = {}
for item in final_data:
	if item[0] not in list( species_reps.keys() ):
		species_reps[ item[0] ] = 0
	if item[0] in list( species_reps.keys() ):
		species_reps[ item[0] ] += int(float( item[3] ))
for item in final_data:
	if item[0] in list( species_reps.keys() ):
		item[4] = str(float(item[3])/species_reps[ item[0] ])

[ print( key, "\t", value ) for key, value in species_reps.items() ]

# [ print(row) for row in final_data ]

# write output data in TIDY format
# encasulates no header for easier parsing in R
with open( output_file, "w" ) as output:
	print( "species\tclass\tcopy_number\tbp\tpercent", file=output )
	for item in final_data:
		print( "\t".join(item), file=output )

"""
---------------------------
	END OF FILE
---------------------------
"""
