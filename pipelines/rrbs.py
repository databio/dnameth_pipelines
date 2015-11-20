#!/usr/bin/env python

"""
RRBS pipeline

"""

from argparse import ArgumentParser
import os, re
import os.path
import sys
import subprocess
import yaml

# Argument Parsing
################################################################################
parser = ArgumentParser(description='Pypiper arguments.')

parser.add_argument('-P', '--pypiper', dest='pypiper_dir',
					default=os.path.dirname(os.path.dirname(os.path.dirname(os.path.realpath(__file__)))) + "/pypiper",
					type=str)
# Set up a pointer to the pypiper code you want to use:
# Just grab the single pypiper arg, and add it to the path here; leaving all other args for later parsing.
args = parser.parse_known_args()[0]
os.sys.path.insert(0, args.pypiper_dir)
from pypiper import Pypiper
from pypiper import ngstk

# Add Pypiper arguments to the arguments list and then reparse.
parser = Pypiper.add_pypiper_args(parser)

# Pipeline options
## Generic options for all pipelines
parser.add_argument('-c','--config', dest='config', type=str, required=True, help='Path to YAML configuration file')
parser.add_argument('-i', '--unmapped-bam', nargs="+", dest='unmapped_bam', required=True,
					help="Input unmapped bam file(s))")
parser.add_argument('-s', '--sample-name', default="default", dest='sample_name', type=str,
					help='Sample Name') # default means deduction from filename, except .bam extension
parser.add_argument('--flowcell', default="", dest='flowcell', type=str, help='Flowcell number')
parser.add_argument('--lane', default="", dest='lane', type=str, help='Lane number')
parser.add_argument('--instr', default="", dest='instr', type=str, help='Instrument name')
parser.add_argument('-r', '--project-root', default="", dest='project_root', type=str,
					help='Directory in which the project will reside.')
parser.add_argument('-g', '--genome', dest='genome_assembly', type=str, required=True, help='Genome Assembly')
parser.add_argument('-C', '--no-checks', dest='no_check', action='store_true',
						default=False, help='Skip sanity checks')
parser.add_argument('-p', '--paired-end', dest='paired_end', action='store_true', help='Paired End Mode')
parser.add_argument('-q', '--single-end', dest='paired_end', action='store_false', help='Single End Mode')
# AS: RRBS as used by us is mainly single-end... maybe it should be the default? Will break with other pipeline's defaults?
# AS: I got the code from Angelo to run it in paired end mode

## Pipeline-specific options
parser.add_argument('-t', '--trimgalore', dest='trimmomatic', action="store_false", default=True, help='Use trimgalore instead of trimmomatic?')
args = parser.parse_args()

if not args.unmapped_bam:
	parser.print_help() #or, print_usage() for less verbosity
	raise SystemExit

# Merging
################################################################################
# If 2 unmapped bam files are given, then these are to be merged.
# Must be done here to initialize the sample name correctly

merge = False
if len(args.unmapped_bam) > 1:
	merge = True
	if args.sample_name == "default":
		args.sample_name = "merged"
else:
	if args.sample_name == "default":
		args.sample_name = os.path.splitext(os.path.basename(args.unmapped_bam[0]))[0]

# Set up environment path variables
################################################################################
# Set up an container class to hold paths
class Container:
	pass

paths = Container()
paths.scripts_dir = os.path.dirname(os.path.realpath(__file__))

# import the yaml config
config = yaml.load(open(args.config, 'r'))

# Resources
paths.resources_dir = config["resources"]["resources"]
paths.adapter_file = os.path.join(paths.resources_dir, "adapters", "epignome_adapters_2_add.fa")
paths.rrbs_adapter_file = os.path.join(paths.resources_dir, "adapters", "RRBS_adapters.fa")
paths.ref_genome = os.path.join(paths.resources_dir, "genomes")
paths.ref_genome_fasta = os.path.join(paths.resources_dir, "genomes", args.genome_assembly, args.genome_assembly + ".fa")
paths.chrom_sizes = os.path.join(paths.resources_dir, "genomes", args.genome_assembly, args.genome_assembly + ".chromSizes")

paths.bismark_indexed_genome = os.path.join(paths.resources_dir, "genomes", args.genome_assembly, "indexed_bismark_bt2")
paths.bismark_spikein_genome = os.path.join(paths.resources_dir, "genomes", "meth_spikein_k1_k3", "indexed_bismark_bt1")

# Tools
paths.python = config["tools"]["python"]
paths.java = config["tools"]["java"]
paths.picard_jar = config["tools"]["picard"]
paths.trimmomatic_jar = config["tools"]["trimmomatic"]
paths.trimmomatic_epignome_jar = config["tools"]["trimmomatic_epignome"]
paths.trimgalore = config["tools"]["trimgalore"]
paths.bowtie1 = config["tools"]["bowtie1"]
paths.bowtie2 = config["tools"]["bowtie2"]
paths.bed2bigBed = config["tools"]["bed2bigBed"]
paths.bed2bigWig = config["tools"]["bed2bigWig"]
paths.bismark = config["tools"]["bismark"]
paths.deduplicate_bismark = config["tools"]["deduplicate_bismark"]
paths.bsmap = config["tools"]["bsmap"]
paths.samtools = config["tools"]["samtools"]

#Output
paths.pipeline_outfolder_abs = os.path.abspath(os.path.join(args.project_root, args.sample_name))

#Biseq paths:
paths.biseq_tools_dir = config["tools"]["biseqMethCalling_tools"]
paths.genomes_split = os.path.join(paths.resources_dir, "genomes_split")

# Create a Pypiper object, and start the pipeline (runs initial setting and logging code to begin)
mypiper = Pypiper(name="RRBS", outfolder=paths.pipeline_outfolder_abs, args=args)

# Create a ngstk object
myngstk = ngstk.NGSTk(paths)

myngstk.make_sure_path_exists(os.path.join(paths.pipeline_outfolder_abs, "unmapped_bam"))

if merge:
	#raise NotImplementedError("Sample merging currently not implemented for RRBS")  # TODO AS: merge currently deactivated for RRBS

	# inactive code
	# There is no reason to deactivate the merge...
	# But currently we can only merge bams, not fastq.gz inputs.
	input_bams = args.unmapped_bam
	merge_folder = os.path.join(paths.pipeline_outfolder_abs, "unmapped_bam")
	sample_merged_bam = args.sample_name + ".merged.bam"
	output_merge = os.path.join(merge_folder, sample_merged_bam)
	cmd = myngstk.merge_bams(input_bams, output_merge)

	mypiper.call_lock(cmd, output_merge)
	args.unmapped_bam = os.path.join(paths.pipeline_outfolder_abs, "unmapped_bam", sample_merged_bam)  #update unmapped bam reference
	local_unmapped_bam_abs = os.path.join(paths.pipeline_outfolder_abs, "unmapped_bam", sample_merged_bam)
	input_ext = ".bam"
else:
	# Link the file into the unmapped_bam directory
	args.unmapped_bam = args.unmapped_bam[0]

	if not os.path.isabs(args.unmapped_bam):
		args.unmapped_bam = os.path.abspath(args.unmapped_bam)
		print args.unmapped_bam

	if args.unmapped_bam.endswith(".bam"):
		input_ext = ".bam"
	elif args.unmapped_bam.endswith(".fastq.gz"):
		input_ext = ".fastq.gz"
	else:
		raise NotImplementedError("This pipeline can only deal with .bam or .fastq.gz files")

	local_unmapped_bam_abs = os.path.join(paths.pipeline_outfolder_abs, "unmapped_bam", args.sample_name + input_ext)
	mypiper.callprint("ln -sf " + args.unmapped_bam + " " + local_unmapped_bam_abs, shell=True)

#check for file exists:
if not os.path.exists(local_unmapped_bam_abs):
	raise Exception(local_unmapped_bam_abs + " is not a file")

# Record file size of input file

cmd = "stat -Lc '%s' " + local_unmapped_bam_abs
input_size = mypiper.checkprint(cmd)
input_size = float(input_size.replace("'",""))

mypiper.report_result("File_mb", round((input_size/1024)/1024,2))
mypiper.report_result("Instrument",args.instr)
mypiper.report_result("Flowcell",args.flowcell)
mypiper.report_result("Lane",args.lane)
if args.paired_end:
	mypiper.report_result("Type","paired")
else:
	mypiper.report_result("Type","single")
mypiper.report_result("Genome",args.genome_assembly)


# Fastq conversion
################################################################################
mypiper.timestamp("### Fastq conversion: ")
fastq_folder = os.path.join(paths.pipeline_outfolder_abs, "fastq")
out_fastq_pre = os.path.join(fastq_folder, args.sample_name)
unaligned_fastq = out_fastq_pre + "_R1.fastq"

if input_ext ==".bam":
	print("Found bam file")
	cmd = myngstk.bam_to_fastq(local_unmapped_bam_abs, out_fastq_pre, args.paired_end)
	mypiper.call_lock(cmd, unaligned_fastq)
	if not args.no_check:
		raw_reads = myngstk.count_reads(local_unmapped_bam_abs, args.paired_end)
		mypiper.report_result("Raw_reads", str(raw_reads))
		fastq_reads = myngstk.count_reads(unaligned_fastq, args.paired_end)
		mypiper.report_result("Fastq_reads", fastq_reads)
		fail_filter_reads = myngstk.count_fail_reads(local_unmapped_bam_abs, args.paired_end)
		pf_reads = int(raw_reads) - int(fail_filter_reads)
		mypiper.report_result("PF_reads", str(pf_reads))
		if fastq_reads != int(raw_reads):
			raise Exception("Fastq conversion error? Number of reads doesn't match unaligned bam")

elif input_ext == ".fastq.gz":
	print("Found gz fastq file")
	cmd = "gunzip -c " + local_unmapped_bam_abs + " > " + unaligned_fastq
	myngstk.make_sure_path_exists(fastq_folder)
	mypiper.call_lock(cmd, unaligned_fastq, shell=True)
	if not args.no_check:
		# Can't make a check here like we can for .bam files. oh well.
		fastq_reads = myngstk.count_reads(unaligned_fastq, args.paired_end)
		mypiper.report_result("Fastq_reads", fastq_reads)


# Adapter trimming (Trimmomatic)
################################################################################
mypiper.timestamp("### Adapter trimming: ")

# We need to detect the quality encoding type of the fastq.
cmd = paths.python + " -u " + os.path.join(paths.scripts_dir, "detect_quality_code.py") + " -f " + unaligned_fastq

encoding_string = mypiper.checkprint(cmd)

if encoding_string.find("phred33") != -1:
	encoding = "phred33"
elif encoding_string.find("phred64") != -1:
	encoding = "phred64"
else:
	raise Exception("Unknown quality encoding type: "+encoding_string)

if args.trimmomatic:

	trimmed_fastq = out_fastq_pre + "_R1_trimmed.fq"
	trimmed_fastq_R2 = out_fastq_pre + "_R2_trimmed.fq"

	memory = str(config["parameters"]["trimmomatic"]["memory"])
	threads = str(config["parameters"]["trimmomatic"]["threads"])
	illuminaclip = str(config["parameters"]["trimmomatic"]["illuminaclip"])

	# REMARK AS: instead of trim_galore we try to use Trimmomatic for now
	# - we are more compatible with the other pipelines
	# - better code base, not a python wrapper of a perl script (as trim_galore)
	# - rrbs-mode not needed because biseq has the same functionality

	# REMARK NS:
	# The -Xmx4000m restricts heap memory allowed to java, and is necessary
	#  to prevent java from allocating lots of memory willy-nilly
	# if it's on a machine with lots of memory, which can lead
	# to jobs getting killed by a resource manager. By default, java will
	# use more memory on systems that have more memory, leading to node-dependent
	# killing effects that are hard to trace.

	cmd = paths.java + " -Xmx" + memory + "g -jar " + paths.trimmomatic_epignome_jar
	if args.paired_end:
		cmd += " PE"
	else:
		cmd += " SE"
	cmd += " -" + encoding
	cmd += " -threads " + threads + " "
	#cmd += " -trimlog " + os.path.join(fastq_folder, "trimlog.log") + " "
	if args.paired_end:
		cmd += out_fastq_pre + "_R1.fastq "
		cmd += out_fastq_pre + "_R2.fastq "
		cmd += out_fastq_pre + "_R1_trimmed.fq "
		cmd += out_fastq_pre + "_R1_unpaired.fq "
		cmd += out_fastq_pre + "_R2_trimmed.fq "
		cmd += out_fastq_pre + "_R2_unpaired.fq "
	else:
		cmd += out_fastq_pre + "_R1.fastq "
		cmd += out_fastq_pre + "_R1_trimmed.fq "
	cmd += "ILLUMINACLIP:" + paths.rrbs_adapter_file + illuminaclip

else: # use trim_galore
	# Trim galore requires biopython, cutadapt modules. RSeQC as well (maybe?)
	#   --- $trim_galore -q $q --phred33 -a $a --stringency $s -e $e --length $l --output_dir $output_dir $input_fastq

	raise NotImplementedError("TrimGalore no longer supported")

	if args.paired_end:
		raise NotImplementedError("TrimGalore for PE RRBS not implemented")
	input_fastq = out_fastq_pre + "_R1.fastq "

	# With trimgalore, the output file is predetermined.
	trimmed_fastq = out_fastq_pre + "_R1_trimmed.fq"

	output_dir=fastq_folder

	#Adapter
	a="AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"

	cmd = paths.trimgalore
	cmd += " -q 20" 	#quality trimming
	cmd += " --" + encoding
	cmd += " -a " + a
	cmd += " --stringency 1"		#stringency: Overlap with adapter sequence required to trim a sequence
	cmd += " -e 0.1" 	#Maximum allowed error rate
	cmd += " --length 16"	#Minimum Read length
	# by unchangeable default Trimmomatic discards reads of lenth 0 (produced by ILLUMINACLIP):
	cmd += " --output_dir " + output_dir + " " + input_fastq

# Trimming command has been constructed, using either trimming options.
# The code to run it is the same either way:

mypiper.call_lock(cmd, trimmed_fastq)   # TODO AS: maybe another lock_name?
trimmed_reads_count = myngstk.count_reads(trimmed_fastq, args.paired_end)

mypiper.clean_add(os.path.join(fastq_folder, "*.fastq"), conditional=True)
mypiper.clean_add(os.path.join(fastq_folder, "*.fq"), conditional=True)
mypiper.clean_add(os.path.join(fastq_folder, "*.log"), conditional=True)
mypiper.clean_add(fastq_folder, conditional=True)

if not args.no_check:
	mypiper.report_result("Trimmed_reads", trimmed_reads_count)

# RRBS alignment with BSMAP.
################################################################################
mypiper.timestamp("### BSMAP alignment: ")
bsmap_folder = os.path.join(paths.pipeline_outfolder_abs, "bsmap_" + args.genome_assembly)  # e.g. bsmap_hg19
myngstk.make_sure_path_exists(bsmap_folder)
# no tmp folder needed for BSMAP alignment

out_bsmap = os.path.join(bsmap_folder, args.sample_name + ".bam")
#original filename in bash version of pipeline:
#bsmap_aligned_bam=$align_dir/$sample.all.aligned.bsmap.$mismatches.mism.$allow_multimapper.r.bam


# REMARK NS: In previous versions of the pipeline, TRIMGALORE used a .fq
# extension, while trimmomatic used .fastq.
# I updated it so that the trimmmomatic path also outputs a .fq file, so this
# doesn't have to vary based on trimmer.

cmd = paths.bsmap
cmd += " -a " + out_fastq_pre + "_R1_trimmed.fq"
if args.paired_end:
	cmd += " -b " + out_fastq_pre + "_R2_trimmed.fq"
cmd += " -d " + paths.ref_genome_fasta
cmd += " -o " + out_bsmap
cmd += " -D " + str(config["parameters"]["bsmap"]["rrbs_mapping_mode"])
cmd += " -w " + str(config["parameters"]["bsmap"]["equal_best_hits"])
cmd += " -v " + str(config["parameters"]["bsmap"]["mismatch_rate"])
cmd += " -r " + str(config["parameters"]["bsmap"]["report_repeat"])
cmd += " -p " + str(config["parameters"]["bsmap"]["processors"])
cmd += " -n " + str(config["parameters"]["bsmap"]["map_to_strands"])
cmd += " -s " + str(config["parameters"]["bsmap"]["seed_size"])
cmd += " -S " + str(config["parameters"]["bsmap"]["random_number_seed"])
cmd += " -f " + str(config["parameters"]["bsmap"]["filter"])
cmd += " -q " + str(config["parameters"]["bsmap"]["quality_threshold"])
cmd += " -u"       # report unmapped reads (into same bam file)
if args.paired_end:
	cmd += " -m " + str(config["parameters"]["bsmap"]["minimal_insert_size"])
	cmd += " -x " + str(config["parameters"]["bsmap"]["maximal_insert_size"])

mypiper.call_lock(cmd, out_bsmap)

if not args.no_check:
	# BSMap apparently stores all the reads (mapped and unmapped) in
	# its output bam; to count aligned reads, then, we have to use
	# a -F4 flag (with count_mapped_reads instead of count_reads).
	x = myngstk.count_mapped_reads(out_bsmap, args.paired_end)
	mypiper.report_result("Aligned_reads", x)
	mypiper.report_result("Aligned_rate", str(float(x)/float(trimmed_reads_count)))
	# In addition, BSMap can (if instructed by parameters) randomly assign
	# multimapping reads. It's useful to know how many in the final bam were such.
	x = myngstk.count_multimapping_reads(out_bsmap, args.paired_end)
	mypiper.report_result("Multimap_reads", x)
	mypiper.report_result("Multimap_rate", str(float(x)/float(trimmed_reads_count)))

# Clean up big intermediate files:
mypiper.clean_add(os.path.join(bsmap_folder, "*.fastq"))
mypiper.clean_add(os.path.join(bsmap_folder, "*.fq"))

# Run biseq-methcalling:
################################################################################
mypiper.timestamp("### biseqMethCalling: ")

# Python Software Requirements for biseqMethCalling
# REMARK AS: all packages are available via "easy_install --user <lib>"
# pip is also a possibility if available (currently not on CeMM infrastructure)
#
# Direct links just in case:
# - biopython: wget https://pypi.python.org/pypi/biopython or wget http://biopython.org/DIST/biopython-1.63.zip
# - bitarray: wget https://pypi.python.org/packages/source/b/bitarray/bitarray-0.8.1.tar.gz
# - guppy: wget https://pypi.python.org/packages/source/g/guppy/guppy-0.1.10.tar.gz
# - pysam: wget https://code.google.com/p/pysam/downloads/detail?name=pysam-0.7.5.tar.gz

biseq_output_path = os.path.join(paths.pipeline_outfolder_abs, "biseqMethcalling_" + args.genome_assembly)
biseq_output_path_web = os.path.join(biseq_output_path, "web")
biseq_output_path_temp = os.path.join(biseq_output_path, "temp")

myngstk.make_sure_path_exists (biseq_output_path)

cmd = paths.python + " -u " + os.path.join(paths.scripts_dir, "biseqMethCalling.py")
cmd += " --sampleName=" + args.sample_name
cmd += " --alignmentFile=" + out_bsmap      # this is the absolute path to the bsmap aligned bam file
cmd += " --methodPrefix=RRBS"
cmd += " --rrbsMode"
cmd += " --checkRestriction"
cmd += " --minFragmentLength=" + str(config["parameters"]["biseqMethCalling"]["minFragmentLength"])
cmd += " --maxFragmentLength=" + str(config["parameters"]["biseqMethCalling"]["maxFragmentLength"])
cmd += " --pfStatus=" + str(config["parameters"]["biseqMethCalling"]["pfStatus"])
cmd += " --maxMismatches=" + str(config["parameters"]["biseqMethCalling"]["maxMismatches"])
cmd += " --baseQualityScoreC=" + str(config["parameters"]["biseqMethCalling"]["baseQualityScoreC"])
cmd += " --baseQualityScoreNextToC=" + str(config["parameters"]["biseqMethCalling"]["baseQualityScoreNextToC"])
cmd += " --laneSpecificStatistics"
cmd += " --bigBedFormat"
cmd += " --deleteTemp"
cmd += " --toolsDir=" + paths.biseq_tools_dir
cmd += " --outputDir=" + biseq_output_path
cmd += " --webOutputDir=" + biseq_output_path_web
cmd += " --tempDir=" + biseq_output_path_temp
cmd += " --timeDelay=" + str(config["parameters"]["biseqMethCalling"]["timeDelay"])
cmd += " --genomeFraction=" + str(config["parameters"]["biseqMethCalling"]["genomeFraction"])
cmd += " --smartWindows=" + str(config["parameters"]["biseqMethCalling"]["smartWindows"])
cmd += " --maxProcesses=" + str(config["parameters"]["biseqMethCalling"]["maxProcesses"])
cmd += " --genomeDir=" + paths.genomes_split
cmd += " --inGenome=" + args.genome_assembly
cmd += " --outGenome=" + args.genome_assembly
# TODO AS: Investigate what happens with biseq in the case of paired-end data

# The dog genome has 38 chromosomes (plus one X chromosome). It's probably best to check here for these rarely used
# reference genomes:
# The default value for includedChromosomes is chr1-30, X, Y, Z (sufficient for human and mouse genomes)
# REMARK NS: This is a hack to account for the way biseqMethCalling restricts to
# default chroms. THis should be fixed in biseq in the future, but for now, this
# lets us run dog samples using the default pipeline. hack!
if args.genome_assembly == "canFam3":
	cmd += ' --includedChromosomes="chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,' \
	       'chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chr23,chr24,chr25,chr26,chr27,chr28,chr29,chr30,chrX,' \
	       'chrY,chrZ,chr31,chr32,chr33,chr34,chr35,chr36,chr37,chr38"'

# Deactivated options:
#cmd += " --appendStatisticsOutput=" + stat_output  # TODO AS: I disable this option for now. This is an analysis-global file where every biseq run writes to
#stat_output = os.path.join(biseq_output_path, "RRBS_biseq_statistics.txt")  # general stats file independent of sample

biseq_finished_helper = os.path.join(biseq_output_path, "biseq.completed")
cmd2 = "touch " + biseq_finished_helper

mypiper.call_lock([cmd, cmd2], target=biseq_finished_helper)

# Now parse some results for pypiper result reporting.
read_variables = ['uniqueSeqMotifCount','totalSeqMotifCount','bisulfiteConversionRate','globalMethylationMean']
totalSeqMotifCount = 0.0
uniqueSeqMotifCount = 0.0
for var in read_variables:
	cmd = paths.python + " -u " + os.path.join(paths.scripts_dir, "tsv_parser.py")
	cmd += " -i " + biseq_output_path + "/RRBS_statistics_" + args.sample_name + ".txt"
	cmd += " -c " + var
	x = mypiper.checkprint(cmd, shell=True)

	if var == 'totalSeqMotifCount':
		totalSeqMotifCount = float(x)
	if var == 'uniqueSeqMotifCount':
		uniqueSeqMotifCount = float(x)

	if var == 'uniqueSeqMotifCount':
		mypiper.report_result('Unique_CpGs', x)
	elif var == 'totalSeqMotifCount':
		mypiper.report_result('Total_CpGs', x)
		mypiper.report_result('meanCoverage', str( totalSeqMotifCount/uniqueSeqMotifCount ))
	else:
		mypiper.report_result(var, x)

################################################################################
mypiper.timestamp("### Make bigbed: ")
# REMARK AS: Make bigwig uses a bismark output file. For RRBS we don't have the bismark cov file
# (essentially a bedgraph file) which the tool bed2bigWig would need
# REMARK AS: UCSC tracks are generated by biseq-methcalling

# First, convert the bed format into the bigBed input style.
# This is how biseq did it, but it's actually unnecessary; instead we can just go straight off the output file.
# Left command here for posterity.
# awk '{ printf "%s\t%s\t%s\t\047%s%[\04720\047]\047\t%s\t%s\n", $1, $2, $3, $5/10, $5, $6 }' RRBS_cpgMethylation_01_2276TU.bed > f


# bigbed conversion input file is the biseq methylation calls output file
biseq_methcall_file = biseq_output_path + "/RRBS_cpgMethylation_" + args.sample_name + ".bed"

bigbed_output_path = os.path.join(paths.pipeline_outfolder_abs, "bigbed_" + args.genome_assembly)
bigwig_output_path = os.path.join(paths.pipeline_outfolder_abs, "bigwig_" + args.genome_assembly)


# bedToBigBed RRBS_cpgMethylation_01_2276TU.bed ~/linkto/resources/genomes/hg19/hg19.chromSizes RRBS_cpgMethylation_test2.bb

myngstk.make_sure_path_exists (bigbed_output_path)
myngstk.make_sure_path_exists (bigwig_output_path)
bigbed_output_file = os.path.join(bigbed_output_path,"RRBS_" + args.sample_name + ".bb")
out_bedGraph = os.path.join(bigwig_output_path,"RRBS_" + args.sample_name + ".bedGraph")
out_bigwig = os.path.join(bigwig_output_path, "RRBS_" + args.sample_name + ".bw")

cmd = paths.bed2bigBed
cmd += " " + biseq_methcall_file
cmd += " " + paths.chrom_sizes
cmd += " " + bigbed_output_file

# REMARK NS: As of June 2015, IGV will load bigBed files for methylation
# in a unique format if the *filename contains  "RRBS_cpgMethylation" -- see
# https://github.com/igvteam/igv/blob/master/src/org/broad/igv/methyl/MethylTrack.java
# This is obviously not ideal, but I will create a link with this filename
# to the original file (even for WGBS tracks) so that you could load these into
# IGV if you want:

filename_hack_link_file = os.path.join(bigbed_output_path,"RRBS_cpgMethylation_" + args.sample_name + ".bb")
cmd2 = "ln -sf " + os.path.relpath(bigbed_output_file, bigbed_output_path) + " " + filename_hack_link_file

mypiper.call_lock([cmd, cmd2], bigbed_output_file)

# Let's also make bigwigs:

# First convert to bedGraph
# hard coding tabs doesn't seem to work:
#cmd = "awk '{ printf \\\"%s\t%s\t%s\t%s\n\\\", $1, $2, $3, $5/10 }'"
cmd = "awk -v OFS='\t' '{ print $1, $2, $3, $5/10 }'"
cmd += " " + biseq_methcall_file
cmd += " > " + out_bedGraph

cmd2 = paths.bed2bigWig
cmd2 += " " + out_bedGraph
cmd2 += " " + paths.chrom_sizes
cmd2 += " " + out_bigwig

mypiper.call_lock([cmd, cmd2], out_bigwig, shell=True)

################################################################################
# Calculate neighbor methylation matching
mypiper.timestamp("### Neighbor Methylation Matching: ")
nmm_output_dir = os.path.join(paths.pipeline_outfolder_abs, "nmm_" + args.genome_assembly)
myngstk.make_sure_path_exists (nmm_output_dir)
nmm_outfile=os.path.join(nmm_output_dir, args.sample_name + ".nmm.bed")

cmd = paths.python + " -u " + os.path.join(paths.scripts_dir, "methylMatch.py")
cmd += " --inFile=" + out_bsmap      # this is the absolute path to the bsmap aligned bam file
cmd += " --methFile=" + biseq_methcall_file
cmd += " --outFile=" + nmm_outfile
cmd += " --cores=4"
cmd += " -q"

mypiper.call_lock(cmd, nmm_outfile)

################################################################################
mypiper.timestamp("### Bismark spike-in alignment: ")
# currently using bowtie1 instead of bowtie2

# get unaligned reads out of BSMAP bam
bsmap_unalignable_bam = os.path.join(bsmap_folder, args.sample_name + "_unalignable.bam")
mypiper.call_lock(paths.samtools + " view -bh -f 4 -F 128 "+out_bsmap+" > " + bsmap_unalignable_bam, bsmap_unalignable_bam, shell=True)

# Re-flag the unaligned paired-end reads to make them look like unpaired for Bismark
if args.paired_end:
	bsmap_unalignable_bam_output = os.path.join(bsmap_folder, args.sample_name + "_unalignable_reflagged.bam")
	cmd = paths.python + " -u " + os.path.join(paths.scripts_dir, "pe_flag_changer.py")
	cmd += " -i " + bsmap_unalignable_bam
	cmd += " -o " + bsmap_unalignable_bam_output
	mypiper.call_lock(cmd, bsmap_unalignable_bam_output)
	mypiper.clean_add(bsmap_unalignable_bam, conditional=True)
	bsmap_unalignable_bam = bsmap_unalignable_bam_output

# convert BAM to fastq
bsmap_fastq_unalignable_pre = os.path.join(bsmap_folder, args.sample_name + "_unalignable")
bsmap_fastq_unalignable = bsmap_fastq_unalignable_pre  + "_R1.fastq"
cmd = myngstk.bam_to_fastq(bsmap_unalignable_bam, bsmap_fastq_unalignable_pre, args.paired_end)
mypiper.call_lock(cmd, bsmap_fastq_unalignable)

# actual spike-in analysis
spikein_folder = os.path.join(paths.pipeline_outfolder_abs, "bismark_spikein")
myngstk.make_sure_path_exists(spikein_folder)
spikein_temp = os.path.join(spikein_folder, "bismark_temp")
myngstk.make_sure_path_exists(spikein_temp)
out_spikein_base = args.sample_name + ".spikein.aln"

out_spikein = os.path.join(spikein_folder, out_spikein_base + ".bam")
cmd = paths.bismark + " " + paths.bismark_spikein_genome + " "
cmd += bsmap_fastq_unalignable_pre + "_R1.fastq"
cmd += " --bam --unmapped"
cmd += " --path_to_bowtie " + paths.bowtie1
#	cmd += " --bowtie2"
cmd += " --temp_dir " + spikein_temp
cmd += " --output_dir " + spikein_folder
cmd += " --basename=" + out_spikein_base
#cmd += " -p 4"
cmd += " -n 0" #allow no mismatches

mypiper.call_lock(cmd, out_spikein, nofail=True)

# Clean up the unmapped file which is copied from the parent
# bismark folder to here:
mypiper.clean_add(os.path.join(spikein_folder,"*.fastq"), conditional=True)
mypiper.clean_add(os.path.join(spikein_folder,"*.fq"), conditional=True)
mypiper.clean_add(spikein_temp) # For some reason, the temp folder is not deleted.



################################################################################
mypiper.timestamp("### PCR duplicate removal (Spike-in): ")
# Bismark's deduplication forces output naming, how annoying.
#out_spikein_dedup = spikein_folder + args.sample_name + ".spikein.aln.deduplicated.bam"
out_spikein_dedup = re.sub(r'.bam$', '.deduplicated.bam', out_spikein)
if not args.paired_end:
	cmd = paths.deduplicate_bismark + " --single "    # TODO: needs module load bismark or absolute path to this tool
	cmd += out_spikein
	cmd += " --bam"
else:
	cmd = paths.deduplicate_bismark + " --paired "
	cmd += out_spikein
	cmd += " --bam"

out_spikein_sorted = out_spikein_dedup.replace('.deduplicated.bam', '.deduplicated.sorted')
cmd2 = paths.samtools + " sort " + out_spikein_dedup + " " + out_spikein_sorted
cmd3 = paths.samtools + " index " + out_spikein_sorted + ".bam"
cmd4 = "rm " + out_spikein_dedup
mypiper.call_lock([cmd, cmd2, cmd3, cmd4], out_spikein_sorted + ".bam.bai", nofail=True)

# Spike-in methylation calling
################################################################################
mypiper.timestamp("### Methylation calling (testxmz) Spike-in: ")

cmd1 = paths.python + " -u " + os.path.join(paths.scripts_dir, "testxmz.py")
cmd1 += " " + out_spikein_sorted + ".bam" + " " + "K1_unmethylated"
cmd1 += " >> " + mypiper.pipeline_stats_file
cmd2 = cmd1.replace("K1_unmethylated", "K3_methylated")
mypiper.callprint(cmd1, shell=True, nofail=True)
mypiper.callprint(cmd2, shell=True, nofail=True)

# TODO: check if results file exists already...


# PDR calculation:
################################################################################

# PDR not applied to PE case because bisulfiteReadConcordanceAnalysis.py crashes
if not args.paired_end:

	mypiper.timestamp("### PDR (Partial Disordered Methylation) analysis")

	pdr_output_dir = os.path.join(paths.pipeline_outfolder_abs, "pdr_" + args.genome_assembly)
	myngstk.make_sure_path_exists (pdr_output_dir)

	# convert aligned bam to sam

	pdr_in_samfile = os.path.join(pdr_output_dir, args.sample_name + ".aligned.sam") # gets deleted after, see some lines below
	mypiper.call_lock(paths.samtools + " view " + out_bsmap + " > " + pdr_in_samfile, pdr_in_samfile, shell=True)

	# PDR calculation:
	#
	# output files:
	pdr_bedfile=os.path.join(pdr_output_dir, args.sample_name + ".pdr.bed")

	produce_sam = True  # TODO AS: make this an option somewhere
	concordsam=os.path.join(pdr_output_dir, args.sample_name + ".concordant.sam")
	discordsam=os.path.join(pdr_output_dir, args.sample_name + ".discordant.sam")

	# command::
	cmd1 = paths.python + " -u " + os.path.join(paths.scripts_dir, "bisulfiteReadConcordanceAnalysis.py")
	cmd1 += " --infile=" + pdr_in_samfile
	cmd1 += " --outfile=" + pdr_bedfile
	cmd1 += " --skipHeaderLines=0"
	cmd1 += " --genome=" + args.genome_assembly
	cmd1 += " --genomeDir=" + paths.ref_genome
	cmd1 += " --minNonCpgSites=3"   # These two parameters are not relevant for PDR analysis
	cmd1 += " --minConversionRate=0.9"

	if produce_sam == True:
		cmd1 += " --concordantOutfile=" + concordsam
		cmd1 += " --discordantOutfile=" + discordsam
		#TODO: perhaps convert them to bam *cough*

	#call:
	mypiper.call_lock(cmd1, pdr_bedfile)

	# delete huge input SAM file
	mypiper.clean_add(pdr_in_samfile)

# Final sorting and indexing
################################################################################
# create sorted and indexed BAM files for visualization and analysis
# bsmap already outputs a sorted and indexed bam file

# Cleanup
################################################################################
mypiper.stop_pipeline()
