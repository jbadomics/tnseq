#!/usr/bin/env python

import sys
import os
from subprocess import call
import csv
import Bio
import re
from Bio import SeqIO, SeqFeature
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import GC
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
    epilog='Author: Jon Badalamenti, Bond Lab, University of Minnesota (http://www.thebondlab.org)\nhttp://github.com/jbadomics/tnseq\nJanuary 2016\n \n', formatter_class=RawTextHelpFormatter)
parser.add_argument('[GENBANK FILE]', help='Genbank file containing, at a minimum, locus tags and corresponding genomic\ncoordinates')
parser.add_argument('[DATA FILE]', help='insertion data as three-column tab-delimited file containing chromosome,\ngenomic coordinate, and total number of insertions at that coordinate, e.g.\nchr1\t327491\t1639')
parser.add_argument('[prefix]', help='output file names will be prepended with this prefix')
parser.add_argument('[N-terminal trim]', type=float, help='percent (expressed as a decimal) of the total gene length to trim from the\nN-terminus of the translated protein')
parser.add_argument('[C-terminal trim]', type=float, help='percent (expressed as a decimal) of the total gene length to trim from the\nC-terminus of the translated protein')
args=parser.parse_args()

insertionEvents=[]

print "reading insertion data..."

# read insertion data into list
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

# define list of CDS coordinates for calculating intergenic regions
CDS_list = [[0,0]]

NtermTrim = float(sys.argv[4])
CtermTrim = float(sys.argv[5])

# parse input genbank file with SeqIO and loop over all sequence records
genbankFile = open(sys.argv[1], 'r')
for sequenceRecord in SeqIO.parse(genbankFile, "genbank"):
	genome = sequenceRecord.seq
	genomeSequence = str(sequenceRecord.seq)
	TAcoordinates = []
	totalTAsites = genome.count('TA')
	
	# collect coordinates of all TA sites in the genome
	for m in re.finditer('TA', genomeSequence):
         TAcoordinates.append(m.end())
	print "working on %s %s (%i bp, %2.2f%% GC, %i TA sites)...\n" % (''.join(sequenceRecord.id), ''.join(sequenceRecord.description).rstrip('.').replace(', complete genome', ''), len(genome), GC(genome), totalTAsites)
	
	for feature in sequenceRecord.features:
	
		if feature.type == 'CDS' or feature.type == 'tRNA' or feature.type == 'rRNA' or feature.type == 'ncRNA':

			locusTag = ''.join(feature.qualifiers["locus_tag"])
			product = ''.join(feature.qualifiers["product"])
			strand = int(feature.location.strand)
			
			geneTAcoordinates = []
			TAinsertionSites = []
			nonTAsites = []
			hits = []
			
			# deal with strand-specific coordinate definitions
			if strand == 1:
				strandSign = '+'
				startCoord = int(feature.location.start.position)
				endCoord = int(feature.location.end.position)
				# correct coordinates based on positional arguments
				realStartCoord = int(round(startCoord + ((endCoord - startCoord) * NtermTrim)))
				realEndCoord = int(round(endCoord - ((endCoord - startCoord) * CtermTrim)))
				geneSequence = str(genome[realStartCoord:realEndCoord])
				gene_TA_count = geneSequence.count('TA')
				# collect all TA sites in within corrected coordinates
				for n in re.finditer('TA', geneSequence):
					geneTAcoordinates.append(n.end() + realStartCoord)
				# collect uncorrected coordinates in CDS_list for determining intergenic regions later
				CDS_list.append([startCoord, endCoord])
				for insertionEvent in insertionEvents:
					chromosome = insertionEvent[0]
					# essential test to match coordinates of seqRecord with the same chromosome in the input hits file
					if chromosome == ''.join(sequenceRecord.id):
						if realStartCoord <= int(insertionEvent[1]) <= realEndCoord:
							hits.append(int(insertionEvent[2]))
							coordinate = int(insertionEvent[1])
							if coordinate in geneTAcoordinates:
								TAinsertionSites.append(coordinate)
							else:
								nonTAsites.append(coordinate)
				
			if strand == -1:
				strandSign = '-'
				startCoord = int(feature.location.end.position) #this is not a typo - biopython defines start and end coordinates from left to right regardless of strandedness
				endCoord = int(feature.location.start.position) #this is not a typo
				# correct coordinates based on positional arguments
				realStartCoord = int(round(startCoord - ((startCoord - endCoord) * NtermTrim)))
				realEndCoord = int(round(endCoord + ((startCoord - endCoord) * CtermTrim)))
				geneSequence = str(genome[realEndCoord:realStartCoord])
				gene_TA_count = geneSequence.count('TA')
				# collect all TA sites in within corrected coordinates
				for n in re.finditer('TA', geneSequence):
					geneTAcoordinates.append(n.end() + realEndCoord)
				# collect uncorrected coordinates in CDS_list for determining intergenic regions later
				CDS_list.append([endCoord, startCoord])
				for insertionEvent in insertionEvents:
					chromosome = insertionEvent[0]
					# essential test to match coordinates of seqRecord with the same chromosome in the input hits file
					if chromosome == ''.join(sequenceRecord.id):
						if realEndCoord <= int(insertionEvent[1]) <= realStartCoord:
							hits.append(int(insertionEvent[2]))
							coordinate = int(insertionEvent[1])
							if coordinate in geneTAcoordinates:
								TAinsertionSites.append(coordinate)
							else:
								nonTAsites.append(coordinate)
			
			if hits:
				if len(geneTAcoordinates) > 0:
					hitPercent = (float(len(TAinsertionSites)) / len(geneTAcoordinates)) * 100
					codingHits = "%s\t%s\t%s\t%i\t%s\t%i\t%i\t%2.1f\t%s" % (''.join(sequenceRecord.id), locusTag, strandSign, sum(hits), len(TAinsertionSites), len(nonTAsites), len(geneTAcoordinates), hitPercent, product)
				else:
					hitPercent = '*'
					codingHits = "%s\t%s\t%s\t%i\t%s\t%i\t%i\t%s\t%s" % (''.join(sequenceRecord.id), locusTag, strandSign, sum(hits), len(TAinsertionSites), len(nonTAsites), len(geneTAcoordinates), hitPercent, product)
				outputFile.write(codingHits+"\n")
			
			if not hits:
				print "no transposon insertions in %s\t%s" % (locusTag, product)
				noHits = "%s\t%s" % (locusTag, product)
				noHitsFile.write(noHits+"\n")
	
	print 'tabulating intergenic insertions in %s %s...\n' % (''.join(sequenceRecord.id), ''.join(sequenceRecord.description).rstrip('.').replace(', complete genome', ''))
	
	# define intergenic coordinates from CDS_list
	intergenicCoords = []			
	for i in range(1, len(CDS_list) - 1):
		last_end = CDS_list[i-1][1]
		this_start = CDS_list[i][0]
		intergenicCoords.append([last_end, this_start])
	
	# loop over intergenic regions and report insertion statistics if length > 0
	for coordPair in intergenicCoords:
		start = int(coordPair[0])
		end = int(coordPair[1])
		if end - start >= 1:
			for insertionEvent in insertionEvents:
				chromosome = insertionEvent[0]
				# essential test to match coordinates of seqRecord with the same chromosome in the input hits file
				if chromosome == ''.join(sequenceRecord.id):
					if start < int(insertionEvent[1]) < end:
						hits = int(insertionEvent[2])
						coordinate = insertionEvent[1]
						# not sum(hits) here since intergenic ranges are not defined by locus tags
						# rather, tabulate only by coordinate
						# future improvement: collect and report information about flanking genes and their orientation
						intergenicHits = "%s\t%s\t%i" % (chromosome, coordinate, hits)
						intergenicHitsFile.write(intergenicHits+"\n")
			
print "...done"
print "insertions in coding regions written to %s" % outputFileName
print "intergenic insertions written to %s" % intergenicHitsFileName

intergenicHitsFile.close()
noHitsFile.close()
outputFile.close()
genbankFile.close()
