#!/usr/bin/env python

import sys
import os
import csv
import re
import Bio
from Bio import SeqIO, SeqFeature
from Bio.SeqUtils import GC
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
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
    description='Given 1) a tab-delimited file providing genomic coordinates of transposon insertions and\n2) a locus tag of interest, this script parses a Genbank file and tabulates the total number of\ntransposon insertions for each insertion site within the genomic coordinates of the given locus tag.\n\nNOTE: Organisms with multiple Genbank records (e.g. those with multiple chromosomes or plasmids)\nshould be concatenated into a single .gbk file before executing this script. For example:\n% cat NC_000001.gbk NC_000002.gbk [...NC_00000n.gbk] > concatenated.gbk\n\nRequires BioPython v. 1.65 or later (http://biopython.org/wiki/Download)', 
    epilog='Author: Jon Badalamenti, Bond Lab, University of Minnesota (http://www.thebondlab.org)\nhttp://github.com/jbadomics/tnseq\nFebruary 2016\n \n', formatter_class=RawTextHelpFormatter)
parser.add_argument('[GENBANK FILE]', help='Genbank file containing, at a minimum, locus tags and corresponding genomic\ncoordinates')
parser.add_argument('[DATA FILE]', help='insertion data as three-column tab-delimited file containing chromosome,\ngenomic coordinate, and total number of insertions at that coordinate, e.g.\nchr1\t327491\t1639')
parser.add_argument('[LOCUS TAG]', help='the locus tag for which to report insertion sites; output files will be named\n<locus tag>.insertions.txt')
parser.add_argument('[N-terminal trim]', type=float, help='percent (expressed as a decimal) of the total gene length to trim from the\nN-terminus of the translated protein')
parser.add_argument('[C-terminal trim]', type=float, help='percent (expressed as a decimal) of the total gene length to trim from the\nC-terminus of the translated protein')
args=parser.parse_args()

genbankFile = open(sys.argv[1], 'r')

insertionEvents=[]

print "reading insertion data..."
# read insertion data into list
with open(sys.argv[2], 'r') as f:
	insertionPointsFile = csv.reader(f,delimiter='\t')	
	for chromosome,coordinate,insertion in insertionPointsFile:
		insertionEvents.append([chromosome, coordinate, insertion])

locusTag = sys.argv[3]
NtermTrim = float(sys.argv[4])
CtermTrim = float(sys.argv[5])

outputFileName = "%s.insertions.txt" % sys.argv[3]
outputFile = open(outputFileName, 'wb')

totalGenomeLength = []
totalTAsites = []

# parse input genbank file with SeqIO and loop over all sequence records
for sequenceRecord in SeqIO.parse(genbankFile, "genbank"):
	# collect genome information
	genome = sequenceRecord.seq
	genomeSequence = str(sequenceRecord.seq)
	totalGenomeLength.append(len(genomeSequence))
	seqID = ''.join(sequenceRecord.id)
	OrganismName = ''.join(sequenceRecord.annotations["source"])
	totalTAsites.append(genomeSequence.count('TA'))
	
	# collect coordinates of all TA sites in the genome
	TAcoordinates = []	
	for m in re.finditer('TA', genomeSequence):
         TAcoordinates.append(m.end())
	
	# loop over features in genbank file
	for feature in sequenceRecord.features:
	
		if feature.type == 'CDS' or feature.type == 'tRNA' or feature.type == 'rRNA' or feature.type == 'ncRNA':
		
			# find locus tag that matches sys.argv[3]
			if ''.join(feature.qualifiers["locus_tag"]) == locusTag:

				matchedLocusTag = ''.join(feature.qualifiers["locus_tag"])
				product = ''.join(feature.qualifiers["product"])
				strand = int(feature.location.strand)
				
				geneTAcoordinates = []
				
				# collect locus tag statistics for CDS and RNAs
				if feature.type == 'CDS':
					aaLength = len(''.join(feature.qualifiers["translation"]))
					ntLength = aaLength * 3
				else:
					print '%s is not a protein coding sequence' % locusTag
					ntLength = abs(int(feature.location.start.position) - int(feature.location.end.position))
					aaLength = 0

                # deal with strand-specific coordinate definitions
				if strand == 1:
					startCoord = int(feature.location.start.position)
					endCoord = int(feature.location.end.position)
					# correct coordinates based on positional arguments
					realStartCoord = int(round(startCoord + ((endCoord - startCoord) * NtermTrim)))
					realEndCoord = int(round(endCoord - ((endCoord - startCoord) * CtermTrim)))
					geneSequence = str(genome[realStartCoord:realEndCoord])
					gene_TA_count = geneSequence.count('TA')
					# collect all TA sites in within corrected coordinates and report TA density
					for n in re.finditer('TA', geneSequence):
						geneTAcoordinates.append(n.end() + realStartCoord)
					TAdensity = gene_TA_count / (float(realEndCoord - realStartCoord) / 1000)
					print '%s %s\t%i..%i(+), %i nt, %i aa, %i TA sites, %2.1f TA sites/kbp\n' % (locusTag, product, startCoord, endCoord, ntLength, aaLength, gene_TA_count, TAdensity)
				
				if strand == -1:
					startCoord = int(feature.location.end.position) #this is not a typo - biopython defines start and end coordinates from left to right regardless of strandedness
					endCoord = int(feature.location.start.position) #this is not a typo
					# correct coordinates based on positional arguments
					realStartCoord = int(round(startCoord - ((startCoord - endCoord) * NtermTrim)))
					realEndCoord = int(round(endCoord + ((startCoord - endCoord) * CtermTrim)))
					geneSequence = str(genome[realEndCoord:realStartCoord])
					# rc_geneSequence = str(genome[realEndCoord:realStartCoord].reverse_complement())
					gene_TA_count = geneSequence.count('TA')
					# collect all TA sites in within corrected coordinates and report TA density
					for n in re.finditer('TA', geneSequence):
						geneTAcoordinates.append(n.end() + realEndCoord)
					TAdensity = gene_TA_count / (float(realStartCoord - realEndCoord) / 1000)
					print '%s %s\t%i..%i(-), %i nt, %i aa, %i TA sites, %2.1f TA sites/kbp\n' % (locusTag, product, endCoord, startCoord, ntLength, aaLength, gene_TA_count, TAdensity)

	totalGoodSites = []
	totalBadSites = []

	# set strand if locus tag not found in sequenceRecord and advance for loop to next sequenceRecord if strand = None
	try: strand
	except NameError: strand = None
	if strand != None:
		for insertionEvent in insertionEvents:
			chromosome = insertionEvent[0]
			# essential test to match coordinates of seqRecord with the same chromosome in the input hits file
			if chromosome == seqID:
				# deal with strand-specific coordinate definitions
				if strand == 1:
					if realStartCoord <= int(insertionEvent[1]) <= realEndCoord:
						numHits = insertionEvent[2]
						insertionCoord = int(insertionEvent[1])
						if insertionCoord in geneTAcoordinates:
							# report statistics for TA sites
							totalGoodSites.append(insertionCoord)
							outputString = "%s\t%s\t%i\t%s\tTA site" % (chromosome, locusTag, insertionCoord, numHits)
						else:
							# report statistics for non-TA sites
							totalBadSites.append(insertionCoord)
							outputString = "%s\t%s\t%i\t%s\t*" % (chromosome, locusTag, insertionCoord, numHits)
						if numHits:
							print outputString
							outputFile.write(outputString+"\n")
		
			if strand == -1:
				if chromosome == seqID:
				# essential test to match coordinates of seqRecord with the same chromosome in the input hits file
					if realEndCoord <= int(insertionEvent[1]) <= realStartCoord:
						numHits = insertionEvent[2]
						insertionCoord = int(insertionEvent[1])
						if insertionCoord in geneTAcoordinates:
							# report statistics for TA sites
							totalGoodSites.append('TA')
							outputString = "%s\t%s\t%i\t%s\tTA site" % (chromosome, locusTag, insertionCoord, numHits)
						else:
							# report statistics for non-TA sites
							outputString = "%s\t%s\t%i\t%s\t*" % (chromosome, locusTag, insertionCoord, numHits)
							totalBadSites.append(insertionCoord)
						if numHits:
							print outputString
							outputFile.write(outputString+"\n")

		# print results for each possible outcome: 1) TA sites hit; 2) only non-TA sites hit; 3) no hits at all
		if totalGoodSites:
			print '%s (%i bp, %2.1f%% GC) contains %i total theoretical TA insertion sites' % (OrganismName, sum(totalGenomeLength), GC(genome), sum(totalTAsites))
			print '%i total theoretical TA insertion sites are present in %s' % (gene_TA_count, locusTag)
			print '%i of %i sites (%3.1f %%) in %s were hit' % (len(totalGoodSites), gene_TA_count, (len(totalGoodSites) / float(gene_TA_count)) * 100,locusTag)
			print "insertion positions written to %s" % outputFileName
			sys.exit() # exit before looping over all sequence records
		if totalBadSites and not totalGoodSites:
			print '%i total theoretical TA insertion sites are present in %s' % (gene_TA_count, locusTag)
			print 'insertions in %s occur only in non-TA sites' % locusTag
			sys.exit() # exit before looping over all sequence records
		if not totalGoodSites and not totalBadSites:
			print "no transposon insertions found in %s" % locusTag
			sys.exit() # exit before looping over all sequence records
		
# raise error and exit script if provided locus tag does not exist
try: matchedLocusTag
except NameError: matchedLocusTag = None
if matchedLocusTag == None:
	print 'ERROR: %s does not exist in the provided genome' % locusTag
	print 'check the locus tag and try again'
	sys.exit()

genbankFile.close()
