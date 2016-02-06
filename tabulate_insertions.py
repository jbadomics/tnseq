#!/usr/bin/env python

import sys
import os
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

NtermTrim = float(sys.argv[4])
CtermTrim = float(sys.argv[5])

# parse input genbank file with SeqIO

for sequenceRecord in SeqIO.parse(genbankFile, "genbank"):
	genome = sequenceRecord.seq
	genomeSequence = str(sequenceRecord.seq)
	TAcoordinates = []
	totalTAsites = genome.count('TA')
	for m in re.finditer('TA', genomeSequence):
         TAcoordinates.append(m.end())
	print "working on %s %s (%i bp, %2.2f%% GC, %i TA sites)" % (''.join(sequenceRecord.id), ''.join(sequenceRecord.description).rstrip('.').replace(', complete genome', ''), len(genome), GC(genome), totalTAsites)
	
	for feature in sequenceRecord.features:
	
		if feature.type == 'CDS' or feature.type == 'tRNA' or feature.type == 'rRNA' or feature.type == 'ncRNA':

			locusTag = ''.join(feature.qualifiers["locus_tag"])
			product = ''.join(feature.qualifiers["product"])
			strand = int(feature.location.strand)
			
			geneTAcoordinates = []

			chromosome = []
			insertionSites = []
			hits = []
			
			if strand == 1:
				startCoord = int(feature.location.start.position)
				endCoord = int(feature.location.end.position)
				realStartCoord = int(round(startCoord + ((endCoord - startCoord) * NtermTrim)))
				realEndCoord = int(round(endCoord - ((endCoord - startCoord) * CtermTrim)))
				geneSequence = str(genome[realStartCoord:realEndCoord])
				gene_TA_count = geneSequence.count('TA')
				for n in re.finditer('TA', geneSequence):
					geneTAcoordinates.append(n.end() + realStartCoord)
				CDS_list.append((startCoord, endCoord))
				for insertionEvent in insertionEvents:
					if realStartCoord <= int(insertionEvent[1]) <= realEndCoord:
						chromosome = insertionEvent[0]
						hits.append(int(insertionEvent[2]))
						insertionSites.append(insertionEvent[1])
						coordinate = insertionEvent[1]
				
			if strand == -1:
				startCoord = int(feature.location.end.position) #this is not a typo
				endCoord = int(feature.location.start.position) #this is not a typo
				realStartCoord = int(round(startCoord - ((startCoord - endCoord) * NtermTrim)))
				realEndCoord = int(round(endCoord + ((startCoord - endCoord) * CtermTrim)))
				geneSequence = str(genome[realEndCoord:realStartCoord])
				gene_TA_count = geneSequence.count('TA')
				for n in re.finditer('TA', geneSequence):
					geneTAcoordinates.append(n.end() + realEndCoord)
				CDS_list.append((endCoord, startCoord))
				for insertionEvent in insertionEvents:
					if realEndCoord <= int(insertionEvent[1]) <= realStartCoord:
						chromosome = insertionEvent[0]
						hits.append(int(insertionEvent[2]))
						insertionSites.append(insertionEvent[1])
						coordinate = insertionEvent[1]
				
			#if sum(hits) > 100:
			if hits:
				codingHits = "%s\t%s\t%i\t%s\t%s\t%i" % (chromosome, locusTag, sum(hits), len(insertionSites), product, len(geneTAcoordinates))
				#print codingHits
				outputFile.write(codingHits+"\n")
			
			if not chromosome:
				print "no transposon insertions in %s\t%s" % (locusTag, product)
				noHits = "%s\t%s" % (locusTag, product)
				noHitsFile.write(noHits+"\n")
	
	print 'tabulating intergenic insertions in %s %s' % (''.join(sequenceRecord.id), ''.join(sequenceRecord.description))...
			
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
