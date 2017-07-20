#! /cm/shared/apps/python/2.7.6/bin/python
##! /usr/bin/python

# Summary: for a file containing genomic regions, this script determines the DNA sequence and various attributes relating to base frequencies
#
# example call: 
# python lookupDnaSequence.py --infile=mm9_regions.txt --chrom=chrom_mm9 --chromstart=chromstart_mm9 --chromend=chromend_mm9 --genome=mm9 --genomeDir=. --includeSequence --includeFrequencies --includeCgiParameters --annotatePatternPositions=CG --outfile=mm9_regions_seq.txt
# python bisulfiteReadFiltering.py --infile=test1_bismark_mm10.deduplicated.sam --outfile=test1_bismark_mm10.deduplicated.filtered.sam --genome=mm10 --genomeDir=.
# python bisulfiteReadFiltering.py --infile=test1_bismark_mm10.deduplicated.sam --outfile=test1_bismark_mm10.deduplicated.filtered.sam --genome=mm10 --genomeDir=~/genome
#
# interactive analysis:
# import os
# os.chdir("D:/Epigenomics/Software/CeMM_scripts/")
# import bisulfiteReadFiltering
# bisulfiteReadFiltering.lookupSequence("mm10", "chr2", 5972848, 5972948)
# import imp; imp.reload(bisulfiteReadFiltering) # reload library after saving changes

import os
import math
import re
from Bio import SeqIO
import sys
#from downloadGenome import downloadGenome

# definition of global variables
chromFiles = {} # will be initialized on the first call to lookupSequence()

global genome
global debug
global rnaMode
debug = True


# removes newline characters at the end of a string (similar to rstrip() but leaves tabs intact)
def removeLinefeed(s):
	while s[-1] in ["\n", "\r"]: s = s[0:-1]
	return s

# obtain the DNA sequence for a given genomic region

#JK: New function lookupSequence()
#Major changes:
#1. load chromosome sequences from an indexed genome fasta (genome)
#2. check CIGAR for insertions or deletions: if present use this information to retrieve the correct genomic sequence
#if not retrieve genomic sequence as was before.
def lookupSequence(genome, chrom, chromstart, chromend, cigar,lowercaseRepeats=False):
	global chromFiles

	# load the corresponding chromosome file
	try:
		chromFile = chromFiles[chrom]
	except KeyError:
	   # chromFilename = genomePath + os.sep + chrom + ".fa"
		if not options.rnaMode:
			print "Reading chromosome file: "+ chrom
		chromFiles[chrom] = genome[chrom].seq
		if debug: print "Reading chromosome file completed"

	if "I" in cigar or "D" in cigar:
		deletion = 0
		insertion = 0
		pos = 0
		cigar_str = re.findall('([MID])', cigar) #JK: array containing the CIGAR IDs
		cigar_int = re.findall('(\d+)', cigar)   #JK: corresponding array containing the counts

		#JK: first calculate the new chromend by including insertions and deletions
		for i in range(len(cigar_str)):
			if cigar_str[i] == "D":
				deletion = deletion + int(cigar_int[i])
			if cigar_str[i] == "I":
				insertion = insertion + int(cigar_int[i])
		chromend = chromend + deletion - insertion
		#JK: get sequence using the updated chromend
		seq = list(str(chromFiles[chrom][chromstart:chromend]))
		#JK: perform insertions (as "N") and deletions in the retrieved genomic sequence according to the CIGAR info
		for i in range(len(cigar_str)):
			if cigar_str[i] == "D":
				for j in range(int(cigar_int[i])):
					seq.pop(pos)
			elif cigar_str[i] == "I":
				for j in range(int(cigar_int[i])):
					seq.insert(pos, "N")
				pos = pos + int(cigar_int[i])
			elif cigar_str[i] == "M":
				pos = pos + int(cigar_int[i])
		seq = "".join(seq)
	else:
		#retrieve genomic sequence with no modifications
		seq = str(chromFiles[chrom][chromstart:chromend])

	if not lowercaseRepeats: seq = seq.upper()
	if debug: print "Sequence of %s:%i-%i (%s) is: %s" % (chrom, chromstart, chromend, genome, seq)
	return seq

# def lookupSequence(genomePath, chrom, chromstart, chromend, lowercaseRepeats=False):
#     global chromFiles
#     # load the corresponding chromosome file
#     try:
#         chromFile = chromFiles[chrom]
#     except KeyError:
#         chromFilename = genomePath + os.sep + chrom + ".fa"
#         print "Reading chromosome file: "+ chromFilename
#         chromFiles[chrom] = SeqIO.read(open(chromFilename), "fasta")
#         if debug: print "Reading chromosome file completed"
#     # retrieve genomic sequence
#     seq = chromFiles[chrom].seq[chromstart:chromend]._data
#
#
#     if not lowercaseRepeats: seq = seq.upper()
#     if debug: print "Sequence of %s:%i-%i (%s) is: %s" % (chrom, chromstart, chromend, genomePath, seq)
#     return seq


# obtain the DNA sequence for a given genomic region (used to work, but now buggy? Needs testing)
def lookupSequence_python2(genomePath, chrom, chromstart, chromend, lowercaseRepeats=False):
	global chromFiles
	global offset
	try:
		chromFile = chromFiles[chrom]
	except KeyError:
		chromFilename = genomePath + os.sep + chrom + ".fa"
		chromFile = open(chromFilename,'r')
		# determine offset
		firstByte = chromFile.read(1)
		offset = 0
		if firstByte == ">": offset = len(chromFile.readline()) + 1 # +1 for first byte
		# make sure that the row length is 50 bp and the line separator is "\n"
		chromFile.seek(offset+49,0)
		seq = chromFile.read(3)
		validLetters = ["A","C","G","T","N"]
		if not seq[0].upper() in validLetters or not seq[1] == "\n" or not seq[2].upper() in validLetters:
			print "Incorrect file format detected for "+chrom+". Sequence around position 50: "+seq
			raise SystemExit
		chromFiles[chrom] = chromFile
	# retrieve genomic sequence
	pos = chromstart / 50 * 51 + chromstart % 50
	length = (chromend - chromstart + 50) / 50 * 51
	chromFile.seek(offset+pos,0)
	seq = chromFile.read(length)
	seq = seq.replace("\n","")
	seq = seq[0:(chromend-chromstart)]
	if not lowercaseRepeats: seq = seq.upper()
	if debug: print "Sequence of %s:%i-%i (%s) is: %s" % (chrom, chromstart, chromend, genomePath, seq)
	return seq


# obtain the DNA sequence for a given genomic region (marginally buggy)
def lookupSequence_python3(genomePath, chrom, chromstart, chromend, lowercaseRepeats=False):
	global chromFiles
	global offset
	try:
		chromFile = chromFiles[chrom]
	except KeyError:
		chromFilename = genomePath + os.sep + chrom + ".fa"
		chromFile = open(chromFilename, 'r')
		# determine offset
		firstByte = chromFile.read(1)
		offset = 0
		if firstByte == ">": offset = len(chromFile.readline()) + 1 # +1 for first byte
		# make sure that the row length is 50 bp and the line separator is "\n"
		chromFile.seek(offset + 49, 0)
		seq = chromFile.read(3)
		validLetters = ["A", "C", "G", "T", "N"]
		if not seq[0].upper() in validLetters or not seq[1] == "\n" or not seq[2].upper() in validLetters:
			print
			"Incorrect file format detected for " + chrom + ". Sequence around position 50: " + seq
			raise SystemExit
		chromFiles[chrom] = chromFile
	# retrieve genomic sequence by
	pos = int(math.floor(float(chromstart) / 50) * 51 + chromstart % 50) # this code is compatible with both Python 2 and Python 3
	length = int(math.floor(float((chromend - chromstart + 50)) / 50) * 51)
	chromFile.seek(offset + pos, 0)
	seq = chromFile.read(length)
	seq = seq.replace("\n", "")
	seq = seq[0:(chromend - chromstart)]
	if not lowercaseRepeats: seq = seq.upper()
	if debug: print("Sequence of %s:%i-%i (%s) is: %s" % (chrom, chromstart, chromend, genomePath, seq))
	return seq


# check whether a read shows sufficient evidence of correct bisulfite conversion
def assessBisulfiteConversionStatus(bisSeq,refSeq,g2a=False):
	length = len(bisSeq)
	nonContextTotal = 0
	nonContextConverted = 0
	nonContextUnconverted = 0
	if not g2a:
		for i in range(0,len(bisSeq)-1):
			if refSeq[i:i+2] == "CA" or refSeq[i:i+2] == "CC" or refSeq[i:i+2] == "CT":
				nonContextTotal += 1
				if bisSeq[i] == "T": nonContextConverted += 1
				if bisSeq[i] == "C": nonContextUnconverted += 1
	else:
		for i in range(0,len(bisSeq)-1):
			if refSeq[i:i+2] == "AG" or refSeq[i:i+2] == "GG" or refSeq[i:i+2] == "TG":
				nonContextTotal += 1
				if bisSeq[i+1] == "A": nonContextConverted += 1
				if bisSeq[i+1] == "G": nonContextUnconverted += 1
	conversionRatio = 0
	if (nonContextConverted+nonContextUnconverted) > 0: conversionRatio = float(nonContextConverted) / (nonContextConverted+nonContextUnconverted)
	retainRead = True
	# discard reads that align to fewer than N cytosines in the reference genome that are not in a CpG context
	if nonContextTotal < options.minNonCpgSites: retainRead = False
	# discard all reads with fewer than X% of cytosines outside of a CpG context showing up as thymines in the reads.
	if rnaMode:
		if conversionRatio > options.maxConversionRate: retainRead = False   #JK: keep methylated (unconverted) reads
	else:
		if conversionRatio < options.minConversionRate: retainRead = False



	if debug:
		print("ref: "+refSeq)
		print("bis: "+bisSeq)
		print("g2a: "+str(g2a))
		print("nonContextTotal: "+str(nonContextTotal))
		print("nonContextConverted: "+str(nonContextConverted))
		print("nonContextUnconverted: "+str(nonContextUnconverted))
		print("conversionRatio: "+str(conversionRatio))
		print("length: "+str(length))
		print("retainRead: "+str(retainRead))
		print("\n")
	return retainRead


def passFilter(line, genome):
	# parsing the current line of the SAM file
	cells = line.split("\t")
	if len(cells) < 11:
		raise Exception("WARNING: Skipping a line that is not consistent with BAM format: "+line)
	chrom = cells[2]
	chromstart = int(cells[3]) - 1 # correct for the 0-based vs. 1-based difference between BAM and lookupSequence()
	cigar = cells[5]
	bisSeq = cells[9]
	chromend = chromstart + len(bisSeq)

	#JK:set passfilter to false for unmapped reads and reads that have strange CIGAR strings
	if int(cells[1]) & 4:   #JK: Exclude discard unmapped reads
		return(False)

	if len(re.findall('[SHPX=]',cells[5]))!= 0 :  # JK skip if any other possible CIGAR characters except for MID (seldom and therefore not handled)
		return(False)

	if "N" in cells[5] and rnaMode:  # JK: Accept spliced reads (CIGAR N) without further checking
		return(True)
	elif "N" in cells[5] and not(rnaMode):
		return(False)



	# obtaining strand conversion information
	# if len(cells) > 15 and cells[14][0:5] == "XR:Z:" and cells[15][0:5] == "XG:Z:":
	#      # valid Bismark read, set conversion information accordingly
	#      if cells[15][0:7] == "XG:Z:CT": g2a = False
	#      elif cells[15][0:7] == "XG:Z:GA": g2a = True
	#      else: raise Exception("Unexpected SAM file format")
	# elif len(cells) > 12 and cells[12][0:5] == "ZS:Z:":
	#      if cells[12][0:6] == "ZS:Z:+": g2a = False
	#      elif cells[12][0:6] == "ZS:Z:-": g2a = True
	#      else: raise Exception("Unexpected SAM file format")
	# else:
	#   raise Exception("WARNING: Skipping line that is not consistent with any of the supported Bismark or BSMAP alignment formats: "+line)
	#JK: Check if read is mapped to the reverse strand if first in pair (not second in pair) or fwd strand if second in pair and thereby set g2a (ZS/XG fields are not available in bt)
	if (int(cells[1]) & 16 and not(int(cells[1]) & 128)) or (not(int(cells[1]) & 16) and int(cells[1]) & 128) : g2a = True
	else: g2a = False




	# obtaining the reference genome sequence
	try:
		#refSeq = lookupSequence(options.genomeDir + os.sep + options.genome, chrom, chromstart, chromend)
		refSeq = lookupSequence(genome, chrom, chromstart, chromend, cigar)
	except (Exception) as ex:
		raise Exception('Could not retrieve sequence for current line ("' + line[0:100] + '") due to the following exception: ' + str(type(ex)) + ": " + str(ex))

	# testing whether the current read shows sufficient evidence of correct bisulfite conversion
	return assessBisulfiteConversionStatus(bisSeq,refSeq,g2a)


# main analysis procedure
def performAnalysis(options):

	# prepare output file
	outfile = open(options.outfile, 'w')
	if options.skipped != "": skipped = open(options.skipped, 'w')

	# download genome sequence (if necessary)
	if not os.path.exists(options.genomeDir + os.sep + options.genome+ os.sep + options.genome + ".fa"):
		print("Genome assembly not found for '" + options.genome )
	#JK: gave error and not really required if one provides the genome
	# if not os.path.exists(options.genomeDir + os.sep + options.genome):
	#     print("Genome assembly not found for '" + options.genome + "', trying to obtain genome from UCSC Genome Browser")
	#     downloadGenome(options.genome, options.genomeDir)
	# if not os.path.exists(options.genomeDir + os.sep + options.genome):
	#     print("Could not obtain required genome assembly for '" + options.genome + "'")
		raise SystemExit

	# JK Initialize genome database (using SeqIO.index_db makes it possible to use a single fasta file instead of aving all chroms in separate fasta files)
	genome = SeqIO.index_db(os.path.join(options.genomeDir, options.genome, options.genome + ".idx"),
							os.path.join(options.genomeDir, options.genome, options.genome + ".fa"),
							"fasta")


	print("Genome " + options.genome + " indexed")

	# iterate through all lines of the SAM file
	infile = open(options.infile, 'r')
	count = 0
	warnings = 0
	acceptedRecords = 0
	unmappedRecords = 0
	skippedRecords = 0
	splicedRecords = 0
	cigarskippedRecords = 0
	unmatchedReads_accepted = 0
	unmatchedReads_skipped = 0
	isMate1 = True
	prevLine = ""
	print("Starting the read filtering")
	for line in infile:
		# process current line from the SAM file

		if count % 1000000 == 0: sys.stdout.write(".")
		if count % 100000000 == 0: sys.stdout.write("\n.")

		if line.split("\t")[0] == "HWI-D00693_0018:1:1208:19326:20646#K562_1_7/1":
			print("stop")


		if warnings > 25:
			raise Exception("Maximum number of warning messages reached. Aborting.")

		count = count + 1
		line = removeLinefeed(line)
		if debug: print(line)

		# optionally process header lines
		if count <= options.headerLines:
			if debug: print("Header line: "+line)
			outfile.write(line + "\n")
			if options.skipped != "": skipped.write(line + "\n")  #JK: write header lines in out AND skipped file

			continue

		if int(line.split("\t")[1]) & 4:   #JK: Exclude discard unmapped reads
			#skipped.write(line + "\n")
			unmappedRecords += 1
			#skippedRecords += 1
			#continue

		if len(re.findall('[SHPX=]',line.split("\t")[5]))!= 0 :  # JK skip if any other possible CIGAR characters except for MID (seldom and therefore not handled)
			#skipped.write(line + "\n")
			#skippedRecords += 1
			cigarskippedRecords += 1
			#continue

		if "N" in line.split("\t")[5] and rnaMode:  # JK: Accept spliced reads (CIGAR N) without further checking
			#outfile.write(line + "\n")
			#acceptedRecords += 1
			splicedRecords += 1
			#continue
		elif "N" in line.split("\t")[5] and not(rnaMode):
			#skipped.write(line + "\n")
			#skippedRecords += 1
			cigarskippedRecords += 1
			#continue


		# testing whether the current read or read pair shows sufficient evidence of correct bisulfite conversion
		try:
			if options.pairedEnd:
				if isMate1:
					isMate1 = False
					prevLine = line
					continue
				else:
					if debug: print("Processing paired-end read")
					if prevLine.split("\t")[0] != line.split("\t")[0]:
						# process unmatched read as a singleton and proceed
						if passFilter(prevLine, genome):
							outfile.write(prevLine + "\n")
							acceptedRecords += 1
							unmatchedReads_accepted += 1
						else:
							if options.skipped != "": skipped.write(prevLine + "\n")
							skippedRecords += 1
							unmatchedReads_skipped += 1
						isMate1 = False
						prevLine = line
					else:
						# process valid paired-end read pair
						if passFilter(prevLine, genome) and passFilter(line, genome):
							outfile.write(prevLine + "\n" + line + "\n")
							acceptedRecords += 2
						else:
							if options.skipped != "": skipped.write(prevLine + "\n" + line + "\n")
							skippedRecords += 2
						isMate1 = True
						prevLine = ""
			else:
				if passFilter(line, genome):
					outfile.write(line + "\n")
					acceptedRecords += 1
				else:
					if options.skipped != "": skipped.write(line + "\n")
					skippedRecords += 1
		except (Exception) as ex:
			warnings += 1
			print(str(ex))

	# print summary statistics
	print("Total number of processed lines:   "+str(count))
	print("Total number of accepted records: "+str(acceptedRecords))
	print("Total number of spliced records: "+str(splicedRecords))
	print("Total number of CIGAR failed records: "+str(cigarskippedRecords))
	print("Total number of unmapped records: "+str(unmappedRecords))
	print("Total number of skipped records:   "+str(skippedRecords))
	if options.pairedEnd:
		print("Total number of accepted unmatched reads: "+str(unmatchedReads_accepted))

	# cleaning up
	infile.close()
	outfile.close()
	if options.skipped != "": skipped.close()


if __name__ == '__main__':
	print
	"Starting program..."
	# constructing command line parser
	import optparse

	parser = optparse.OptionParser()
	parser.add_option('--infile', action='store', type='string', dest='infile',help='Specify the name of the input file',default="")
	parser.add_option('--outfile', action='store', type='string', dest='outfile',help='Specify the name of the output file', default="")
	parser.add_option('--skipped', action='store', type='string', dest='skipped',help='Specify the name of the file to store skipped reads (optional)', default="")
	parser.add_option('--skipHeaderLines', action='store', type='int', dest='headerLines',help='Specify the number of header lines that should be skipped (default: 0)',default=0)
	parser.add_option('--pairedEnd', action='store_true', dest='pairedEnd',help='Specify paired-end mode (default: False)', default=False)
	parser.add_option('--minNonCpgSites', action='store', type='int', dest='minNonCpgSites',help='Specify the minimum number of genomic cytosines outside of a CpG context required for estimating the conversion rate', default=5)
	parser.add_option('--minConversionRate', action='store', type='float', dest='minConversionRate',help='Specify the minimum conversion rate for cytosines outside of a CpG context', default=0.9)
	#Add option maxConversionRate for CORE-Seq RNA
	parser.add_option('--maxConversionRate', action='store', type='float', dest='maxConversionRate',help='Specify the maximum conversion rate for cytosines outside of a CpG context (CORE-seq)', default=0.1)
	parser.add_option('--genome', action='store', type='string', dest='genome',help='Specify the genome for which the DNA sequence should be retrieved')
	parser.add_option('--genomeDir', action='store', type='string', dest='genomeDir',help='Specify the name of the directory containing the genome sequence data (will be downloaded from USCS Genome Browser if not available locally)',default=".")
	parser.add_option('-v', '--verbose', action='store_true', dest='verbose',help='Print debugging information (default: False)', default=False)
	#Add option to activate RNA-mode
	parser.add_option('-r', '--rna', action='store_true', dest='rnaMode',help='Activate RNA mode (default: False)', default=False)

	(options, args) = parser.parse_args()
	debug = options.verbose
	rnaMode = options.rnaMode
	if options.outfile == "": options.outfile = options.infile + ".output.txt"
	if not options.infile or not options.genome:
		print("Mandatory parameters missing. Program will terminate now.")
		print("\nYour parameter settings:")
		print(options)
		raise SystemExit
	try:
		performAnalysis(options)
	except Exception as e:
		print e
		exit(1)

	print
	"Program successfully terminating...."
