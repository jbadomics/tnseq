#!/usr/bin/env python

import sys
import os
import csv
import Bio
from Bio import SeqIO, SeqFeature
from Bio.SeqRecord import SeqRecord
import argparse
from argparse import RawTextHelpFormatter

# silence Biopython warnings
import warnings
from Bio import BiopythonParserWarning
from Bio import BiopythonWarning
warnings.simplefilter('ignore', BiopythonParserWarning)
warnings.simplefilter('ignore', BiopythonWarning)

# script help and usage
# parser=argparse.ArgumentParser(
    # description='This script parses a Genbank file and writes FASTA amino acid sequences (.faa) of putative multiheme\nc-type cytochromes (3 or more CXXCH motifs; if this qualification met, counts also CXXXCH motifs)\n\nNOTE: Organisms with multiple Genbank records (e.g. those with multiple chromosomes or plasmids)\nshould be concatenated into a single .gbk file before executing this script. For example:\n% cat NC_000001.gbk NC_000002.gbk [...NC_00000n.gbk] > concatenated.gbk\n\nRequires BioPython v. 1.65 or later (http://biopython.org/wiki/Download)', 
    # epilog='Author: Jon Badalamenti, Bond Lab, University of Minnesota (http://www.thebondlab.org)\nJune 2015\n \n', formatter_class=RawTextHelpFormatter)
# parser.add_argument('[GENBANK FILE]', help='Genbank file containing translated amino acid sequences. Requires SOURCE annotation\nto be present in the header. For example:\n \nLOCUS       NC_002939            3814128 bp    DNA     circular CON 16-MAY-2014\nDEFINITION  Geobacter sulfurreducens PCA chromosome, complete genome.\nACCESSION   NC_002939\nVERSION     NC_002939.5  GI:400756305\nDBLINK      BioProject: PRJNA57743\nKEYWORDS    RefSeq.\nSOURCE      Geobacter sulfurreducens PCA <--- *MUST BE PRESENT* \n  ORGANISM  Geobacter sulfurreducens PCA <--- *MUST BE PRESENT*')
# args=parser.parse_args()

# define input file handle
genbankFile = open(sys.argv[1], 'r')

insertionEvents=[]

print "reading insertion data..."

with open(sys.argv[2], 'r') as f:
	insertionPointsFile = csv.reader(f,delimiter='\t')	
	for chromosome,coordinate,insertion in insertionPointsFile:
		insertionEvents.append([chromosome, coordinate, insertion])

print "done"
print "tabulating insertions by locus tag..."

outputFileName = "%s.tabulated_insertions.coding.txt" % sys.argv[3]
outputFile = open(outputFileName, 'wb')

noHitsFileName = "%s.nohits.txt" % sys.argv[3]
noHitsFile = open(noHitsFileName, 'wb')

intergenicHitsFileName = "%s.tabulated_insertions.intergenic.txt" % sys.argv[3]
intergenicHitsFile = open(intergenicHitsFileName, 'wb')

CDS_list = []

NtermTrim = int(sys.argv[4])
CtermTrim = float(sys.argv[5])

# parse input genbank file with SeqIO
for sequenceRecord in SeqIO.parse(genbankFile, "genbank"):
	
	# loop over features in genbank file
	for feature in sequenceRecord.features:
	
		if feature.type == 'CDS' or feature.type == 'tRNA' or feature.type == 'rRNA':

			locusTag = ''.join(feature.qualifiers["locus_tag"])
			product = ''.join(feature.qualifiers["product"])
			startCoord = int(feature.location.start.position)
			endCoord = int(feature.location.end.position)
			realStartCoord = int(startCoord + ((endCoord - startCoord) * NtermTrim))
			realEndCoord = int(endCoord - ((endCoord - startCoord) * CtermTrim))
			print '%i\t%i\t%i' % (startCoord, endCoord, realEndCoord)
			
			CDS_list.append((startCoord, endCoord))

			chromosome = []
			insertionSites = []
			hits = []
			
			for insertionEvent in insertionEvents:
				if realStartCoord <= int(insertionEvent[1]) <= realEndCoord:
					chromosome = insertionEvent[0]
					hits.append(int(insertionEvent[2]))
					insertionSites.append(insertionEvent[1])
					coordinate = insertionEvent[1]
				
			#if sum(hits) > 100:
			if hits:
				codingHits = "%s\t%s\t%i\t%s" % (chromosome, locusTag, sum(hits), len(insertionSites))
				outputFile.write(codingHits+"\n")
			
			if not chromosome:
				print "no transposon insertions in %s\t%s" % (locusTag, product)
				noHits = "%s\t%s" % (locusTag, product)
				noHitsFile.write(noHits+"\n")
	
	print 'tabulating intergenic insertions...'
			
	for i,pospair in enumerate(CDS_list[1:]):
		last_end = CDS_list[i][1]
		this_start = pospair[0]
		if this_start - last_end >= 1:
			chromosome = []
			insertionSites = []
			hits = []
			intergenicHits = []
			for insertionEvent in insertionEvents:
				if last_end < int(insertionEvent[1]) < this_start:
					chromosome = insertionEvent[0]
					hits.append(int(insertionEvent[2]))
					insertionSites.append(insertionEvent[1])
					coordinate = insertionEvent[1]
					intergenicHits = "%s\t%s\t%i" % (chromosome, coordinate, sum(hits))
					intergenicHitsFile.write(intergenicHits+"\n")
			
print "...done"
print "insertions in coding regions written to %s" % outputFileName
print "intergenic insertions written to %s" % intergenicHitsFileName

intergenicHitsFile.close()
noHitsFile.close()
outputFile.close()
genbankFile.close()
