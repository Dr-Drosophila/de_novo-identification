"""
========= FILE INFORMATION =========================================================
File name:			postmatch_name_adder.py
Run using:			python3 postmatch_name_adder.py
Input(s):			-
Output(s):			-
Description:
	Takes the file and adds the name of the species to it

Author: 	Chinmay Rele
Date:		2019/04/18
Version:	0.0.1

========= IMPORTS =================================================================
"""

import requests				# get images
import urllib
import os 					# to create folders
from fpdf import FPDF 		# to create PDF
import shutil, sys  		# to remove directories
import tqdm					# to print progress bar on terminal
from pprint import pprint 	# to print items to screen enatly

"""
========= FUNCTIONS ===============================================================
"""

# functions

"""
========= CODE ====================================================================
"""

path = os.getcwd()
# print( path )

for root, dirs, files in os.walk( path ):
	for filename in files:
		print( filename )
