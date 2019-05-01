"""
========= FILE INFORMATION =========================================================
File name:			sequence_getter.py
Run using:			python3 sequence_getter.py
Input(s):			input file name, seq_name(s)
Output(s):			FASTA with specified sequence(s)
Description:
	from a file that has many FASTA sequences, gets a single sequence and outputs to a file.

Author: 	Chinmay Rele
Date:		2019/01/30

run using:
	python3 ../sequence_getter.py input_fasta seq_name(s) output_fasta

========= IMPORTS =================================================================
"""

import sys # to import a file from the commmand line
from Bio import SeqIO # Need Seq.IO for parsing and making FASTA file
import collections # TEMP TO SEE IF THERE ARE REPEATS IN ANY LISTS
import re # USED TO SPLIT STRINGS TO SEPARATE THE NAME AND ALIGNMENT OF THE REPEATMASKER LIBRARY ANNOTATION.

"""
========= FUNCTIONS ===============================================================
"""

# fucntions



"""
========= CODE ====================================================================
"""

input_file_name = sys.argv[-3]
sequence_name = sys.argv[-2]
output_file_name = sys.argv[-1]

print( "input_file_name:\t", input_file_name )
print( "sequence_name:\t\t", sequence_name )
print( "output_file_name:\t", output_file_name )

sequence_list = sequence_name.split(',')
print( "sequence_list: ", sequence_list )


with open(input_file_name) as input, open(output_file_name, 'a') as ouput:
	records = SeqIO.parse(input, 'fasta')
	for record in records:
		for seq in sequence_list:
			if seq in record.description:
				SeqIO.write(record, ouput, 'fasta')
				print(file=output)


"""
---------------------------
	END OF FILE
---------------------------
"""
