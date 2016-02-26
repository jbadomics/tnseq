#!/usr/bin/env python

from __future__ import print_function
import sys
import Bio
from Bio import SeqIO, SeqFeature
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

# Make file write mode compatible for Python 2 and Python 3
writemode = 'wb' if sys.version_info < (3,) else 'w'


# script help and usage
parser=argparse.ArgumentParser(
    description='extracts full nucleotide sequence record in FASTA format from a (multi)-Genbank file.',
    epilog='Author: Jon Badalamenti, Bond Lab, University of Minnesota (http://www.thebondlab.org)\nhttp://github.com/jbadomics/tnseq\nFebruary 2016\n \n', formatter_class=RawTextHelpFormatter)
parser.add_argument('[GENBANK FILE]', help='Genbank file from which to extract genome sequence')
args=parser.parse_args()

outputFileName = sys.argv[1].replace(".gbk", ".fasta")

with open(outputFileName, writemode) as outputFile:
	with open(sys.argv[1], 'r') as genbankFile:
		for sequenceRecord in SeqIO.parse(genbankFile, "genbank"):
			seqHeader = ''.join(sequenceRecord.id)
			genomeSequence = str(sequenceRecord.seq)
			print(seqHeader, len(genomeSequence))
			SeqIO.write(sequenceRecord, outputFile, "fasta")
