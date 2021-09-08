# !/usr/bin/env python3
# Eleanor Robin Heuchan 
# Version 1.0.0
# Last edit date: 16.08.21

instructions = ['''
Usage: $ python3 fasta-search.py <mode> <input_file> ...
	''', '''
is - INTERRUPTION SEARCH
     Searches an input fasta file for possible one, two and three base 
     interruptions and gives as outputs: search sequence (optional); total 
     appearances of each interruption; number of reads in which each 
     interruption appears; allele frequency.
     Bases immediately preceding and following possible interruptions must be
     given as optional arguments, however these can be left blank by inputting
     empty quotes (''). 
     The preceding and following bases can be modified to search for 
     interruptions of greater than three bases, in this instance -v should 
     always be used for clarity. 

     Usage: $ python3 fasta-search.py is <input_file> <preceding_bases> 
	        <following_bases> 
            
            optional: <minimum counts to display (default=1)> 
                      <[-rev] also search for reverse complements>
                      <[-v] display full search sequences>
                      <[-ord] sort the output into descending order by reads>
                      <[-o outfilename] outputs to file, if no file name is
                      given the default file name is [input_file].csv>
	''', '''
ss - SEQUENCE SEARCH
     Searches an input fasta file for a given sequence and outputs: total reads
     in file; search sequence; total appearances of sequence; number of reads
     in which the sequence appears; allele frequency. 

     Usage: $ python3 fasta-search.py ss <input_file> <search_sequence>
         
            optional: <[-rev] also search for reverse complement>
	''', '''
sf - SEQUENCE FILTER
     Searches an input fasta file for a given sequence and outputs any reads 
     containing the sequence. 
     Use '> output.fasta' on the command line to write to a new file. 

     Usage: $ python3 fasta-search.py ss <input_file> <search_sequence>
            
              optional: <[-rev] also search for reverse complements>
	''']

from sys import argv
import pandas as pd
from Bio import SeqIO
import csv

# Assign mode and input file

mode = str(argv[1])

# Check if in help mode before proceeding

if mode == '-h' or mode == '-help':
	for i in instructions:
		print(i)
	exit()

# Assign input file to variable 

f1 = argv[2]

# Set optional input defaults

minimum = 1
revcomp = False
verbose = False
ordered = False
out_to_file = False

# Assign inputs to variables according to mode 

if mode.lower() == 'is':

	if str(f1) == '-h' or str(f1) == '-help':
		print(instructions[1])
		exit()

	preceding = str(argv[3])
	following = str(argv[4])

	if len(argv) == 6:

		if argv[5] == '-rev':
			revcomp = True
		elif argv[5] == '-v':
			verbose = True
		elif argv[5] == '-ord':
			ordered = True
		elif argv[5] == '-o':
			out_to_file = True 
			outFile = (str(f1)+ '.csv') 
		else:
			minimum = int(argv[5])

	elif len(argv) == 7:

		if argv[5] == '-rev' or argv[6] == '-rev':
			revcomp = True
		if argv[5] == '-v' or argv[6] == '-v':
			verbose = True
		if argv[5] == '-ord' or argv[6] == '-ord':
			ordered = True
		if argv[6] == '-o':
			out_to_file = True
			outFile = (str(f1)+ '.csv')
		elif argv[5] == '-o': 
			out_to_file = True 
			outFile = (str(argv[6]) + '.csv')
		if argv[5] != '-rev' and argv[5] != '-v' and argv[5] != '-ord' and argv[5] != '-o':
			minimum = int(argv[5])

	elif len(argv) == 8:

		if argv[5] == '-rev' or argv[6] == '-rev':
			revcomp = True 
		if argv[5] == '-v' or argv[6] == '-v' or argv[7] == '-v':
			verbose = True
		if argv[5] == '-ord' or argv[6] == '-ord' or argv[7] == '-ord':
			ordered = True
		if argv[7] == '-o':
			out_to_file = True
			outFile = (str(f1)+ '.csv') 
		elif argv[6] == '-o':
			out_to_file = True 
			outFile = (str(argv[7]) + '.csv')
		if argv[5] != '-rev' and argv[5] != '-v' and argv[5] != '-ord':
			minimum = int(argv[5])

	elif len(argv) == 9:

		if argv[5] == '-rev' or argv[6] == '-rev':
			revcomp = True 
		if argv[5] == '-v' or argv[6] == '-v' or argv[7] == '-v':
			verbose = True
		if argv[6] == '-ord' or argv[7] == '-ord' or argv[8] == '-ord':
			ordered = True	
		if argv[8] == '-o':
			out_to_file = True
			outFile = (str(f1)+ '.csv')
		elif argv[7] == '-o':
			out_to_file = True 
			outFile = (str(argv[8]) + '.csv')
		if argv[5] != '-rev' and argv[5] != '-v':
			minimum = int(argv[5])

	elif len(argv) == 10:

		if argv[5] == '-rev' or argv[6] == '-rev':
			revcomp = True 		
		if argv[6] == '-v' or argv[7] == '-v':
			verbose = True
		if argv[7] == '-ord' or argv[8] == '-ord':
			ordered = True	
		if argv[9] == '-o':
			out_to_file = True
			outFile = (str(f1)+ '.csv') 
		elif argv[8] == '-o':
			out_to_file = True 
			outFile = (str(argv[9]) + '.csv')
		if argv[5] != '-rev':
			minimum = int(argv[5])

	elif len(argv) == 11:
		
		if argv[6] == '-rev':
			revcomp = True
		if argv[7] == '-v':
			verbose = True
		if argv[8] == '-ord':
			ordered = True
		if argv[9] == '-o':
			out_to_file = True
			outFile = (str(argv[10]) + '.csv')
		minimum = int(argv[5])


elif mode.lower() == 'ss':

	if str(f1) == '-h' or str(f1) == '-help':
		print(instructions[2]) #CHANGE 
		exit()

	query = str(argv[3]).upper()

	if len(argv) == 5:
		if argv[4] == '-rev':
			revcomp = True
				
elif mode.lower() == 'sf':

	if str(f1) == '-h' or str(f1) == '-help':
		print(instructions[3]) #CHANGE
		exit()

	query = str(argv[3]).upper()

	if len(argv) == 5:
		if argv[4] == '-rev':
			revcomp = True

# List of interruptions to search for 

interruptions = [
	'C', 'G', 'A', 'T',
	'CC', 'CG', 'CA', 'CT', 
	'GC', 'GG', 'GA', 'GT',
	'AC', 'AG', 'AA', 'AT',
	'TC', 'TG', 'TA', 'TT',
	'TTT', 'TTC', 'TTA', 'TTG', 
	'TCT', 'TCC', 'TCA', 'TCG', 
	'TAT', 'TAC', 'TAA', 'TAG', 
	'TGT', 'TGC', 'TGA', 'TGG', 
	'CTT', 'CTC', 'CTA', 'CTG', 
	'CCT', 'CCC', 'CCA', 'CCG', 
	'CAT', 'CAC', 'CAA', 'CAG', 
	'CGT', 'CGC', 'CGA', 'CGG', 
	'ATT', 'ATC', 'ATA', 'ATG', 
	'ACT', 'ACC', 'ACA', 'ACG',
	'AAT', 'AAC', 'AAA', 'AAG', 
	'AGT', 'AGC', 'AGA', 'AGG', 
	'GTT', 'GTC', 'GTA', 'GTG', 
	'GCT', 'GCC', 'GCA', 'GCG', 
	'GAT', 'GAC', 'GAA', 'GAG', 
	'GGT', 'GGC', 'GGA', 'GGG'
	]

# Dictionary for finding reverse complements 

complements = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}

# Reads input fasta file to pandas dataframe 

seqList = []
sequences = SeqIO.parse(open(f1), 'fasta')
for record in sequences:
	des = record.description
	des = str('>' + des)
	seq = str(record.seq)
	seqList.append([des, seq])
seq_df = pd.DataFrame(seqList, columns=['ID','Seq'])

# IS

def interruption_search():

	outList = []

	for i in interruptions:

		total_app = 0
		total_reads = 0 

		query = preceding + i + following 

		if revcomp == True:
			query_rev = ''.join([complements[base] for base in query[::-1]])
			for index in seq_df['Seq']:
				count = index.count(query) + index.count(query_rev)
				total_app += count
				if count > 0:
					total_reads += 1 
		else:
			for index in seq_df['Seq']:
				count = index.count(query)
				total_app += count
				if count > 0: 
					total_reads += 1

		af = total_reads / (len(seq_df.index))
		af_short = round(af, 5)

		if total_reads >= minimum:
			if verbose == True and revcomp == True:
				outList.append([i, query + ' (' + query_rev + ')' , total_app, total_reads, af_short])
			elif verbose == True: 
				outList.append([i, query, total_app, total_reads, af_short])
			else:
				outList.append([i, total_app, total_reads, af_short])

	if ordered == True:
		outList = sorted(outList, key=lambda x:x[2], reverse=True)

	if out_to_file == True:
		if verbose == True:
			headers = ['Interruption', 'Query', 'Count', 'Reads', 'AF']
		else:
			headers = ['Interruption', 'Count', 'Reads', 'AF']
		
		with open(outFile, 'w') as outcsv:
			writer = csv.writer(outcsv)
			writer.writerow(headers)
			writer.writerows(outList)
	else:
		out_df = pd.DataFrame(outList, columns=['Interruption', 'Count', 'Reads', 'AF'])
		print(out_df.to_string(index=False))

# SS 

def sequence_search():

	total_app = 0
	total_reads = 0

	if revcomp == True:
		query_rev = ''.join([complements[base] for base in query[::-1]])
		for index in seq_df['Seq']:
			count = index.count(query) + index.count(query_rev)
			total_app += count
			if count > 0:
				total_reads += 1
	else:
		for index in seq_df['Seq']:
			count = index.count(query)
			total_app += count	
			if count > 0:
				total_reads += 1

	af = total_reads / (len(seq_df.index))
	af_short = round(af, 5)
	
	print('\n' + f1)
	print('Total Reads in File: ' + str(len(seq_df.index)))
	print('Search Sequence: ' + str(query))
	print('Total Appearances of Seq.: ' + str(total_app))
	print('Reads Containing Seq.: ' + str(total_reads))
	print('Allele Frequency: ' + str(af_short))

# Search for a given sequence and output reads in which it appears 

def sequence_filter():
	
	outList = []
	if revcomp == True:
		query_rev = ''.join([complements[base] for base in query[::-1]])
		for index in seq_df.index:
			seq = seq_df['Seq'][index]
			count = seq.count(query) + seq.count(query_rev)
			if count > 0:
				outList.append(seq_df['ID'][index])
				outList.append(seq)
	else:
		for index in seq_df.index:
			seq = seq_df['Seq'][index]
			count = seq.count(query)
			if count > 0:
				outList.append(seq_df['ID'][index])
				outList.append(seq)
	for i in outList:
		print(i)

# Run appropriate function according to mode 

if mode.lower() == 'is':
	interruption_search()
elif mode.lower() == 'ss':
	sequence_search()
elif mode.lower() == 'sf':
	sequence_filter()