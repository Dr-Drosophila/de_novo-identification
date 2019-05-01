"""
========= FILE INFORMATION =========================================================
File name:			nanopore_id_condensor.py
Run using:			python3 nanopore_id_condensor.py -i $fasta_file -o $output_fasta
Input(s):			-
Output(s):			dsim_nanopore (FASTA)
Description:
	Condenses seq ids of nanopore data
 	Consensus_Consensus_Consensus_utg000004l_pilon_pilon_pilon
		to
	utg000004l

Author: 	Chinmay Rele
Date:		2019/01/30

========= IMPORTS =================================================================
"""

import sys # to import a file from the commmand line
import getopt # TO PASS IN INPUT AND OUTPUT PARAMETERS
from Bio import SeqIO # Need Seq.IO for parsing and making FASTA file
import collections # TEMP TO SEE IF THERE ARE REPEATS IN ANY LISTS
import re # USED TO SPLIT STRINGS TO SEPARATE THE NAME AND ALIGNMENT OF THE REPEATMASKER LIBRARY ANNOTATION.

"""
========= FUNCTIONS ===============================================================
"""

if __name__ == "__main__":
    # print( "from main" )
    try:
        opts, args = getopt.getopt( sys.argv[1:],"hi:o:",["help", "input_file=", "output_file="] )
    except getopt.GetoptError:
        print( 'python3 nanopore_id_condensor.py -i <input_file> -o <output_file>' )
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print( 'python3 nanopore_id_condensor.py -i <input_file> -o <output_file>' )
            sys.exit()
        elif opt in ("-i", "--output_file"):
            input_file = arg
        elif opt in ("-o", "--output_file"):
            output_file = arg
    print( 'no main input_file is:', input_file )
    print( 'no main output_file is:', output_file )

"""
========= CODE ====================================================================
"""

# open input and output files:
#     parse
#     for each record:
#         shorten name
#         write to output_file

with open( input_file , "r" ) as original, open( output_file , "w" ) as corrected:
	records = SeqIO.parse(original, 'fasta')
	for record in records:
		record.id = record.id.rsplit('::', 1)[0]
		record.description = ""
		if len( record.id ) >= 50:
			print( "{} still has length >= 50.".format( record.id ) )
        # print( record.id )
		SeqIO.write( record, corrected, 'fasta' )


"""
---------------------------
	END OF FILE
---------------------------
"""
