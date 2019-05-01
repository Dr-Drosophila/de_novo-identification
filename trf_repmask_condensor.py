"""
========= FILE INFORMATION =========================================================
File name:			trf_repmask_condensor.py
Run using:			python3 ../trf_repmask_condensor.py -t <trf_output> -r <repeat_masker_output> -o <edited ouput>
Input(s):			-t <trf_output>
					-r <repeat_masker_output>
Output(s):			-o <edited ouput>
					trf_dat.bed = trf in bed format
					rmask_out.bed = repeatmasker in bed format

Description:
	Takes in the TRF and the RepeatMasker output that identified TEs and combines them.

Author: 	Chinmay Rele
Date:		2019/02/26

========= IMPORTS =================================================================
"""

import sys # to import a file from the commmand line
import getopt # TO PASS IN INPUT AND OUTPUT PARAMETERS
from Bio import SeqIO # Need Seq.IO for parsing and making FASTA file
import collections # TEMP TO SEE IF THERE ARE REPEATS IN ANY LISTS
import re # USED TO SPLIT STRINGS TO SEPARATE THE NAME AND ALIGNMENT OF THE REPEATMASKER LIBRARY ANNOTATION.
import math # for absolute methods

"""
========= FUNCTIONS ===============================================================
"""

## get information from bash
if __name__ == "__main__":
	try:
		opts, args = getopt.getopt( sys.argv[1:],"ht:r:o:",["help", "trf_output=", "repeat_masker_output=", "edited_ouput="] )
	except getopt.GetoptError:
		print( 'python3 ../trf_repmask_condensor.py \n\t -t <trf_output> \n\t -r <repeat_masker_output> \n\t -o <edited_ouput> ' )
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			print( 'python3 ../trf_repmask_condensor.py \n\t -t <trf_output> \n\t -r <repeat_masker_output> \n\t -o <edited_ouput>' )
			sys.exit()
		elif opt in ("-t", "--trf_output"):
			trf_output = arg
		elif opt in ("-r", "--repeat_masker_output"):
			repeat_masker_output = arg
		elif opt in ("-o", "--edited_ouput"):
			edited_ouput = arg

	print( 'trf_output is:', trf_output )
	print( 'repeat_masker_output is:', repeat_masker_output )
	print( 'edited_ouput is:', edited_ouput )
	print(  )
"""
========= CODE ====================================================================
"""

sum = 0

# read in the TRF file and rename items
with open( trf_output, "r" ) as trf, open( edited_ouput, "w" ) as out, open( "trf_dat.bed", "w" ) as trf_bed:
	assembly = ""
	for line in trf:
		if "Sequence" in line:
			assembly = line.strip().split(" ")[1]
		elif len(line.strip()) != 0:
			temp = line.strip().split(" ")
			if len(temp) == 15:
				temp.insert( 0, assembly )
				rep = temp[-2]
				temp = temp[0:3]
				sum += int( temp[2] ) - int(temp[1])
				temp.append( "Simple_repeat" )
				print( "\t".join(temp) , file=out )
				print( "\t".join(temp) , file=trf_bed )


print( "sum = %d" %sum )

# empty list to add all lines of repeat_masker_output to, and then to put it back in the file.
rmask_file = []

# read in and parse the rmasker output
with open( repeat_masker_output, "r" ) as rmask, open( edited_ouput, "a" ) as out, open( "rmask_out.bed", "w" ) as rmask_bed:
	for line in rmask:
		if "SW" not in line and "score" not in line and len(line.strip()) > 0:
			if "Satellite" in line or "Simple_repeat" in line:
				temp = line.strip().split()
				temp2 = [ temp[4], temp[5], temp[6], temp[10] ]
				print( "\t".join(temp2) , file=out)
				print( "\t".join(temp2), file=rmask_bed)
			else:
				rmask_file.append( line.strip() )


with open( "shortened_rmask_file.out", "w" ) as out:
	[ print( item, file=out ) for item in rmask_file ]
