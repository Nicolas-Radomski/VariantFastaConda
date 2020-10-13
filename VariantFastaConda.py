#!/usr/bin/env python
# -*- coding: iso-8859-1 -*-

#### author: Nicolas Radomski ####
# version for conda in a cluster: conda enviroment PairedEndVariantCalling
# transform a vcf file from the main scripts PairedEndVariant.py and VariantAlignmentConda.py into a multi-fasta file for downstream phylogenomic analyses, including unvariable sites from the single contig reference genome, as well as single nucleotide polymorphisms (SNPs) and/or small Insertions/Deletions (InDels) from each samples
# the present main script VariantFastaConda.py corresponds to an adaptation with Python3 of the main script VCFtoPseudoGenome.py that Arnaud Felten developped with Python2
# the module genomic.py has to be with the present main script PairedEndAssemblyConda.py to lunch it properly
# the present main script VariantFasta.py and module genomic.py (version 20201006, Octobre 2020) were prepared and tested with Python and Conda packages below (Name/Version/Build/Channel)
#- python/3.8.5/h1103e12_9_cpython/conda-forge
#- biopython/1.78/py38h1e0a361_0/conda-forge

'''
#### exemple of Bash command (bash_VariantFastaConda.sh) ####
#!/bin/bash
#SBATCH -p Research
#SBATCH -o %x.%N.%j.out
#SBATCH -e %x.%N.%j.err
#SBATCH --cpus-per-task=1
#SBATCH --job-name=test-20201012
source /global/conda/bin/activate;conda activate PairedEndVariantCalling; \
python /global/bio/projets/GAMeR/Nicolas-Radomski/PairedEndVariant/VariantFastaConda.py \
	-i /global/bio/projets/GAMeR/Nicolas-Radomski/PairedEndVariant/VariantCalling/7_alignment/filtered.snps.indels.vcf \
	-ref /global/bio/projets/GAMeR/Nicolas-Radomski/PairedEndVariant/VariantCalling/1_reference/Enteritidis_P125109.fasta \
	-o /global/bio/projets/GAMeR/Nicolas-Radomski/PairedEndVariant/VariantCalling/8_matrix/results \
	--NoINDELs

#### exemple of Bash command execution ####
sbatch bash_VariantFastaConda.sh
'''

import os, sys, csv, argparse, genomic
from Bio import SeqIO
from itertools import chain

# parse arguments
def get_parser():

	# function asking arguments
	parser = argparse.ArgumentParser(description='transform a vcf file (e.g. from the main scripts PairedEndVariant.py and VariantAlignmentConda.py) into a multi-fasta file for downtream phylogenomic analyses, including unvariable sites from the single contig reference genome, as well as single nucleotide polymorphisms (SNPs) and/or small Insertions/Deletions (InDels) from each samples')

	# setting of arguments

	parser.add_argument('-i', action="store", dest='VCF', 
						type=str, required=True, help='vcf file (REQUIRED)')

	parser.add_argument('-ref', dest='ref', 
						type=str, required=True, help='reference genome in FASTA (REQUIRED)')

	parser.add_argument('-o', action="store", dest='output', 
						type=str, default='output.fasta', help='output (default:output.fasta)')

	parser.add_argument('--NoINDELs', dest='NoINDELs', action='store_true', help='remove INDELs position (default:False)', default=False)

	return parser

# write a multi Fasta file from a dictionnary (outputFile is the name of the Fasta output / dicoResult is a dictionnary with Fasta headers as keys and sequences as values)
def write_fasta(outputFile, dicoResult):
	outTab = open(outputFile, "w")
	for element in dicoResult :
		header = element.split('/')[-1]
		outTab.write('>' + header + '\n')
		sequence = ''.join(dicoResult[element])
		i = 0
		for nt in sequence :
			outTab.write(nt)
			i+=1
			if i%70==0 :
				outTab.write('\n')
		outTab.write('\n\n')
	outTab.close()	

# stock the sequences from a Fasta file into a dictionnary (FASTAfile is the name of the Fasta file / record_dict is a dictionnary with Fasta headers as keys and sequences as values)
def fasta_to_dico(FASTAfile):
	handle = open(FASTAfile, "r")
	record_dict = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
	handle.close()
	return record_dict

# check if the reference Fasta file contain only one contig (dicoRef is a dictionnary with Fasta headers as keys and sequences as values)
def check_ref(dicoRef):
	if len(dicoRef.keys())>1 :
		print ("#### ERROR: There are more than one sequence in the reference Fasta file .... Please, change the reference Fasta file")
		sys.exit(1)

#main function	
def main():
	# get parser object
	parser=get_parser()
		
	# print parser.help if there are no arguments in the command
	if len(sys.argv)==1:
		parser.print_help()
		sys.exit(1)

	# extract arguments from parser
	Arguments=parser.parse_args()

	# stock reference in a dictionary
	ref = fasta_to_dico(Arguments.ref)
	check_ref(ref)
	seqRef = str(ref[list(ref.keys())[0]].seq)

	firstLine = True
	listGenome = []
	dicoResult_ALL = {}
	pos = 0

	if Arguments.NoINDELs :
		VCFfile = open(Arguments.VCF, 'r')
		lines = VCFfile.readlines()
		VCFfile.close()
		tmpFile = open("tmp_pseudo.vcf",'w')

		for line in lines :
			row = line.split('\t')
			if line[0]=='#':
				tmpFile.write(line)
			elif(len(row[3])==1):
				alternatif = line.split('\t')[4]
				if not(',' in alternatif) and len(alternatif)==1:
					tmpFile.write(line)	
				else :
					alternatifs = alternatif.split(',')	
					flag = True
					for element in alternatifs :
						if len(element)>1:
							flag = False
							break
					if flag : 	
						tmpFile.write(line)		
		tmpFile.close()
		Arguments.VCF = "tmp_pseudo.vcf"

	with open(Arguments.VCF, 'r') as csvfile:
		TSVreader = csv.reader(csvfile, delimiter='\t')
		for row in TSVreader :
			if("##" not in row[0]):
				if firstLine :
					nbGenomes = len(row) - 9					
					for i in row[9:] :
						dicoResult_ALL[i] = []
						listGenome.append(i)
					firstLine = False					
				else :				
					variant = []
					variant.append(row[3])
					variantModif = []
					if(row[4] != '.' and '*' not in row[4]):					
						variantModif.append(row[4].split(','))
						variantModif = list(chain.from_iterable(variantModif))
						variant = variant + variantModif
						newPos = int(row[1]) 
						if(newPos > pos + 1):
							for genome in dicoResult_ALL :
								dicoResult_ALL[genome].append(seqRef[pos:newPos-1])
						pos = newPos
						y = 0
						j = 0	
						for i in row[9:] :
							var = i.split(":")[0]
							maxLen = 0
							for baseVar in variant :
								if len(baseVar) > maxLen :
									maxLen = len(baseVar)
							if maxLen > 1 :
								if j == 0 :
									variantGap = []
									for baseVar in variant :
										variantGap.append(baseVar + '-'*(maxLen-len(baseVar)))
									variant = variantGap							
							if(var=='.'):
								dicoResult_ALL[listGenome[y]].append('N')
							else:	
								dicoResult_ALL[listGenome[y]].append(variant[int(var)])
							y += 1
							j += 1
	if(pos < len(seqRef)):
		for genome in dicoResult_ALL :
			dicoResult_ALL[genome].append(seqRef[pos:len(seqRef)])	

	# print dicoResult_ALL							
	write_fasta(Arguments.output, dicoResult_ALL)
	if Arguments.NoINDELs :
		os.remove("tmp_pseudo.vcf")

	# congtratulate users with the function congratulation of the module genomic.py
	genomic.congratulation()

# driver code: if the code above is a script, call  main() function, rather than to considere it as a module
if __name__ == "__main__":
	# calling main() function
	main() 
