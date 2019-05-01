"""
========= FILE INFORMATION =========================================================
File name:			postmatch_name_adder.py
Run using:			python3 postmatch_name_adder.py
Input(s):			-
Output(s):			-
Description:
	Takes the file and adds the name of the species to it

Author: 	Chinmay Rele
Date:		2019/04/24
Version:	0.0.3

========= IMPORTS =================================================================
"""

import sys						# system usage
import os 						# to create folders
import shutil, sys  			# to remove directories
from tqdm import tqdm as pbar	# to print progress bar on terminal
import time
from pprint import pprint 		# to print items to screen enatly
from Bio import SeqIO			# Parsing and making fasta_files

"""
========= FUNCTIONS ===============================================================
"""

# functions

"""
========= CODE ====================================================================
"""

## get path
path = os.getcwd() + "/og_files"
# print( path )

## create and populate list of fasta files
base_postmatch_files = []

# grab all file names from the path
for root, dirs, files in os.walk( path ):
	for filename in files:
		if "postmatch" in filename and ".py" not in filename and filename != "all_postmatch.fasta":
			base_postmatch_files.append( "./og_files/" + filename )

# pprint( base_postmatch_files )

with open( "all_postmatch.fasta", "w" ) as all_fasta:
	for file in pbar( base_postmatch_files, desc="Species: " ):
		species_name = file.split( "_" )[1].split("/")[1] + "_"
		records = SeqIO.parse( file, "fasta" )
		for record in pbar( records, desc="Each Record" ):
			record.id = species_name + record.id
			record.description = ""
			SeqIO.write( record, all_fasta, "fasta" )
