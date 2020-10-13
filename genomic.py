#!/usr/bin/env python
# -*- coding: iso-8859-1 -*-

#### author: Nicolas Radomski ####
# usefull functions from a personal module
# the present module genomic.py (version 20201006, Septembre 2020) was prepared and tested with Python and Conda packages below (Name/Version/Build/Channel)
#- python/3.7.8/h6f2ec95_1_cpython/conda-forge
#- biopython/1.78/py38h1e0a361_0/conda-forge
# for instance the function congratulation of the module genomic.py must be called in the main script like below:
# import genomic
# genomic.congratulation()

# import Python modules
import os
import sys
import re
import string
from Bio import SeqIO

# function printing a congratulation message
def congratulation():
	print('#### Congratulation .... The workflow finished successfully ... Thank you for your trust')

# function creating a directory if it does not exist or returning an error message and exiting from the systhem if it exists
def nodir_makedir_error(directory):
	if not os.path.exists(directory):
		os.mkdir(directory)
		print("#### The directory %s is successfully created" %directory)
	else:    
		print("#### ERROR: The directory %s already exists and should be removed of renamed" %directory)
		sys.exit("#### ERROR: Please, remove or rename the directory %s before to recall the script" %directory)

# function creating a directory if it does not exist or returning a warning message if it exists
def nodir_makedir_warning(directory):
	if not os.path.exists(directory):
		os.mkdir(directory)
		print("#### The directory %s is successfully created" %directory)
	else:    
		print("#### WARNING: The directory %s already exists" %directory)

# function returning an error message and exiting from the systhem if a directory is emptyness (e.g. usually a failed previous step of the workflow), or a successfull message if this last is full
def emptydir_error_success(directory):
	if len(os.listdir(directory)) == 0:
		print("#### ERROR: Results are not produced in %s" %directory)
		sys.exit("#### ERROR: Please, check standard error (stderr.log) and output (stdout.log) before to recall the script" %directory)
	else:
		wd = os.getcwd()
		output = wd + '/' + directory
		print("#### Results are successfully produced in {} and can be found in {}" .format(directory, output))

# function returning an warning message from the systhem if a file does not exist (e.g. usually a failed previous step of the workflow), or a successfull message if this last exist
def absentefile_warning_success(expectedfile):
	if os.path.exists(expectedfile) is False:
		print("#### WARNING: The expected file does not exist in %s" %expectedfile)
	else:
		wd = os.getcwd()
		output = wd + '/' + expectedfile
		print("#### The expected file is successfully produced in {} and can be found in {}" .format(expectedfile, output))

# function returning an error message from the systhem and exiting if a file does not exist (e.g. usually a failed previous step of the workflow), or a successfull message if this last exist
def absentefile_error_success(expectedfile):
	if os.path.exists(expectedfile) is False:
		print("#### ERROR: The expected file does not exist in %s" %expectedfile)
		sys.exit("#### ERROR: Please, check standard error (stderr.log) and output (stdout.log) before to recall the script" %directory)
	else:
		wd = os.getcwd()
		output = wd + '/' + expectedfile
		print("#### The expected file is successfully produced in {} and can be found in {}" .format(expectedfile, output))

# function adding a prefix (e.g. usually ID sample) to mutiple files of a specific directory
def prefix_files_dir(directory, prefix):
	wd = os.getcwd()
	path = wd + '/' + directory
	files = [ f for f in os.listdir(path) if os.path.isfile(os.path.join(path,f)) ]
	for f in files:
		os.rename(path + '/' + f, path + '/' + prefix + '_' + f)
	print("#### ID prefix is added in files from %s" %directory)
	print("#### Renamed files:", files)

# replace headers of fasta files with a specific string (i.e. ID sample)
def string_in_fasta(infilepath, outfilepath, ID):
	sedcmd = 'sed "s@>.*@>' + ID + '@" ' + infilepath + ' > ' + outfilepath
	print("#### Execute %s" %sedcmd)
	os.system(sedcmd)
	print("#### The new fasta file with ID sample is in %s" % outfilepath)

# add incrementing numbers in the headers of a fasta file
def number_in_fasta(infilepath, outfilepath):
	inf = open(infilepath)
	data = inf.readlines()
	inf.close()
	outf = open(outfilepath, "w")
	c = 1
	l = 1
	for i in data:
		i = re.sub("\n|\r", "", i)
		if c%2 != 0:
			outf.write(i+"_" +str(l) +"\n")
			l+=1
		else:
			outf.write(i +"\n")
		c += 1
	outf.close()
	print("#### The new fasta file with contig number is in %s" % outfilepath)

# convert a fasta file from multi to single lines
def multi_single_fasta(infilepath, outfilepath):
	with open(infilepath, "rU") as input_handle, open(outfilepath, "w") as output_handle:
		sequences = SeqIO.parse(input_handle, "fasta")
		count = SeqIO.write(sequences, output_handle, "fasta-2line")
		print("#### %s contigs are converted in a single line fasta file" % count)
		print("#### The new single line fasta file is in %s" % outfilepath)

# convert a fasta file from single to multi lines
def single_multi_fasta(infilepath, outfilepath):
	with open(infilepath, "rU") as input_handle, open(outfilepath, "w") as output_handle:
		sequences = SeqIO.parse(input_handle, "fasta-2line")
		count = SeqIO.write(sequences, output_handle, "fasta")
		print("#### %s contigs are converted in a multi lines fasta file" % count)
		print("#### The new multi line fasta file is in %s" % outfilepath)

# remove small contigs from a fasta file
def parse_contigs_fasta(infilepath, outfilepath, length):
	input_seq_iterator = SeqIO.parse(infilepath, "fasta-2line")
	long_contigs = [record for record in input_seq_iterator if len(record.seq) > length]
	print("#### {} contigs are higher than {} bases" .format(len(long_contigs), length))
	SeqIO.write(long_contigs, outfilepath, "fasta-2line")
	print("#### The new single line fasta file is in %s" % outfilepath)
