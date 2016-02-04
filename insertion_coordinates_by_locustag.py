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
parser=argparse.ArgumentParser(
    description='Given a tab-delimited file providing genomic coordinates of transposon insertions, this script\nparses a Genbank file and tabulates the total number of transposon insertions and insertion sites\nper locus tag. It also tabulates intergenic insertions and locus tags without any insertions. \n\nNOTE: Organisms with multiple Genbank records (e.g. those with multiple chromosomes or plasmids)\nshould be concatenated into a single .gbk file before executing this script. For example:\n% cat NC_000001.gbk NC_000002.gbk [...NC_00000n.gbk] > concatenated.gbk\n\nRequires BioPython v. 1.65 or later (http://biopython.org/wiki/Download)', 
    epilog='Author: Jon Badalamenti, Bond Lab, University of Minnesota (http://www.thebondlab.org)\nhttp://github.com/jbadomics/tnseq\nFebruary 2016\n \n', formatter_class=RawTextHelpFormatter)
parser.add_argument('[GENBANK FILE]', help='Genbank file containing, at a minimum, locus tags and corresponding genomic\ncoordinates')
parser.add_argument('[DATA FILE]', help='insertion data as three-column tab-delimited file containing chromosome,\ngenomic coordinate, and total number of insertions at that coordinate, e.g.\nchr1\t327491\t1639')
parser.add_argument('[LOCUS TAG]', help='the locus tag for which to report insertion sites; output file names will be prepended with the locus tag as a prefix')
parser.add_argument('[N-terminal trim]', type=float, help='percent (expressed as a decimal) of the total gene length to trim from the\nN-terminus of the translated protein')
parser.add_argument('[C-terminal trim]', type=float, help='percent (expressed as a decimal) of the total gene length to trim from the\nC-terminus of the translated protein')
args=parser.parse_args()

# define input file handle
genbankFile = open(sys.argv[1], 'r')

insertionEvents=[]

print "reading insertion data..."

with open(sys.argv[2], 'r') as f:
	insertionPointsFile = csv.reader(f,delimiter='\t')	
	for chromosome,coordinate,insertion in insertionPointsFile:
		insertionEvents.append([chromosome, coordinate, insertion])

outputFileName = "%s.insertions.txt" % sys.argv[3]
outputFile = open(outputFileName, 'wb')

# CDS_list = []

locusTag = sys.argv[3]

NtermTrim = float(sys.argv[4])
CtermTrim = float(sys.argv[5])

# parse input genbank file with SeqIO
for sequenceRecord in SeqIO.parse(genbankFile, "genbank"):
	
	# loop over features in genbank file
	for feature in sequenceRecord.features:
	
		if feature.type == 'CDS' or feature.type == 'tRNA' or feature.type == 'rRNA' or feature.type == 'ncRNA':
		
			if ''.join(feature.qualifiers["locus_tag"]) == locusTag:

				matchedLocusTag = ''.join(feature.qualifiers["locus_tag"])
				product = ''.join(feature.qualifiers["product"])
				strand = int(feature.location.strand)
				if feature.type == 'CDS':
					aaLength = len(''.join(feature.qualifiers["translation"]))
					ntLength = aaLength * 3
				else:
					print '%s is not a protein coding sequence' % locusTag

				if strand == 1:
					startCoord = int(feature.location.start.position)
					endCoord = int(feature.location.end.position)
					realStartCoord = int(round(startCoord + ((endCoord - startCoord) * NtermTrim)))
					realEndCoord = int(round(endCoord - ((endCoord - startCoord) * CtermTrim)))
					print 'tabulating insertions in %s %s (+ strand, %i nt, %i aa)...' % (locusTag, product, ntLength, aaLength)
				
				if strand == -1:
					startCoord = int(feature.location.end.position) #this is not a typo
					endCoord = int(feature.location.start.position) #this is not a typo
					realStartCoord = int(round(startCoord - ((startCoord - endCoord) * NtermTrim)))
					realEndCoord = int(round(endCoord + ((startCoord - endCoord) * CtermTrim)))
					print 'tabulating insertions in %s %s (- strand, %i nt, %i aa)...' % (locusTag, product, ntLength, aaLength)

try: matchedLocusTag
except NameError: matchedLocusTag = None

if matchedLocusTag == None:
	print 'ERROR: %s does not exist in the provided genome' % locusTag
	print 'check the locus tag and try again'
	sys.exit()

totalSites = []

for insertionEvent in insertionEvents:
	if strand and strand == 1:
		if realStartCoord <= int(insertionEvent[1]) <= realEndCoord:
			totalSites.append(insertionEvent[1])
			chromosome = insertionEvent[0]
			numHits = insertionEvent[2]
			insertionCoord = int(insertionEvent[1])
			outputString = "%s\t%s\t%i\t%s" % (chromosome, locusTag, insertionCoord, numHits)
			if numHits:
				print outputString
				outputFile.write(outputString+"\n")
		
	if strand and strand == -1:
		if realEndCoord <= int(insertionEvent[1]) <= realStartCoord:
			totalSites.append(insertionEvent[1])
			chromosome = insertionEvent[0]
			numHits = insertionEvent[2]
			insertionCoord = int(insertionEvent[1])
			outputString = "%s\t%s\t%i\t%s" % (chromosome, locusTag, insertionCoord, numHits)
			if numHits:
				print outputString
				outputFile.write(outputString+"\n")

if totalSites:
	print '%i total unique insertion sites in %s' % (len(totalSites), locusTag)
	print "insertions written to %s" % outputFileName
if not totalSites:
	print "no transposon insertions in %s" % locusTag

outputFile.close()
genbankFile.close()
