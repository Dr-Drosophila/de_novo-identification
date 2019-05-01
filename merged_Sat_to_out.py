"""
========= FILE INFORMATION =========================================================
File name:			merged_Sat_to_out.py
Run using:			python3 ../merged_Sat_to_out.py -p <merged_satellites> -s <shortened_out_file> -o <final_output>
Input(s):			-p <merged_satellites>
					-s <shortened_out_file>
Output(s):			-o <final_output>
Description:
	Appends merged satellite data to the shortened_out_file, with correct formatting.

Author: 	Chinmay Rele
Date:		2019/02/27

========= IMPORTS =================================================================
"""

import sys # to import a file from the commmand line
import getopt # TO PASS IN INPUT AND OUTPUT PARAMETERS
from Bio import SeqIO # Need Seq.IO for parsing and making FASTA file
import collections # TEMP TO SEE IF THERE ARE REPEATS IN ANY LISTS
import re # USED TO SPLIT STRINGS TO SEPARATE THE NAME AND ALIGNMENT OF THE REPEATMASKER LIBRARY ANNOTATION.
import math # for absolute method
import random # for random choice

"""
========= FUNCTIONS ===============================================================
"""

if __name__ == "__main__":
    try:
        opts, args = getopt.getopt( sys.argv[1:],"hp:s:o:",["help", "merged_satellites=", "shortened_out_file=", "final_output="] )
    except getopt.GetoptError:
        print( 'python3 ../merged_Sat_to_out.py -p <merged_satellites> -s <shortened_out_file> -o <final_output>' )
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print( 'python3 ../merged_Sat_to_out.py -p <merged_satellites> -s <shortened_out_file> -o <final_output>' )
            sys.exit()
        elif opt in ("-p", "--merged_satellites"):
            merged_satellites = arg
        elif opt in ("-s", "--shortened_out_file"):
            shortened_out_file = arg
        elif opt in ("-o", "--final_output"):
            final_output = arg

    print( 'merged_satellites is:', merged_satellites )
    print( 'shortened_out_file is:', shortened_out_file )
    print( 'final_output is:', final_output )
    print(  )

"""
========= CODE ====================================================================
"""

# open and read merged satellite data into a list so can edit and append to file later.
merged_Sat_list = []
with open( merged_satellites, "r" ) as file:
	for line in file:
		merged_Sat_list.append( line.strip().split() )

# a list of formatted elements
formatted = []
for item in merged_Sat_list:
	if item[3] not in ["Satellite", "Simple_Repeat"]:
		item[3] = random.choice(["Satellite", "Simple_Repeat"])
	temp = [ "000", "0.0", "0.0", "0.0", item[0], item[1], item[2], "(00000)", "+", "rnd-n_family-000", item[3], item[1], item[2], "(000)", "0000", "*" ]
	formatted.append( temp )

# appending the two files (to include Sattelites in a single file)
with open( shortened_out_file, "r" ) as in_file, open( final_output, "w" ) as out_file:
	for item in in_file:
		print( "\t".join( item.strip().split() ), file=out_file )
	for item in formatted:
		print( "\t".join(item), file=out_file )
