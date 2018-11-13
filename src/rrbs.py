#!/usr/bin/env python

"""
RRBS pipeline
"""

__author__ = "Nathan Sheffield"
__email__ = "nathan@code.databio.org"
__credits__ = ["Charles Dietz", "Johanna Klughammer", "Christoph Bock", "Andreas Schoenegger"]
__license__ = "GPL3"
__version__ = "0.3.0-dev"


import os
import re
import pypiper
from helpers import get_epi_cmd


def _parse_args(cmdl):
	from argparse import ArgumentParser

	parser = ArgumentParser(description='Pipeline')

	# First, add arguments from Pypiper, including
	# 1. pypiper options, 2. looper connections, 3. common options,
	# using the all_args=True flag (you can customize this).
	# Adds options including; for details, see:
	# http://github.com/epigen/pypiper/command_line_args.md
	parser = pypiper.add_pypiper_args(parser, all_args=True)

	# Add any pipeline-specific arguments
	parser.add_argument("-t", "--trimgalore", dest="trimgalore", action="store_true",
		help='Use trimgalore instead of trimmomatic?')
	parser.add_argument("-e", "--epilog", dest='epilog', action="store_true",
		help='Use epilog for meth calling?')
	parser.add_argument("--pdr", dest="pdr", action="store_true",
		help='Calculate Proportion of Discordant Reads (PDR)?')
	parser.add_argument("--rrbs-fill", dest="rrbs_fill", type=int, default=4,
		help="Number of bases from read end to regard as unreliable and ignore due to RRBS chemistry")
	parser.add_argument("--dark-bases", type=int, default=None,
		help="Number of bases from to prepend to R1 from R2 for dark sequencing")

	args = parser.parse_args(cmdl)

	# Translate pypiper method of read type specification into flag-like option.
	if args.single_or_paired == "paired":
		args.paired_end = True
	else:
		args.paired_end = False

	# Input is required.
	if not args.input:
		parser.print_help()
		raise SystemExit

	return args


def main(cmdl):
	args = _parse_args(cmdl)
	# Create a PipelineManager object and start the pipeline
	outfolder = os.path.abspath(os.path.join(args.output_parent, args.sample_name))
	pm = pypiper.PipelineManager(name="RRBS", outfolder=outfolder, args=args, version=__version__)

	# Set up a few additional paths not in the config file
	pm.config.tools.scripts_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)), "tools")
	pm.config.resources.ref_genome_fasta = os.path.join(pm.config.resources.genomes, args.genome_assembly, args.genome_assembly + ".fa")
	pm.config.resources.chrom_sizes = os.path.join(pm.config.resources.genomes, args.genome_assembly, args.genome_assembly + ".chromSizes")
	pm.config.resources.genomes_split = os.path.join(pm.config.resources.resources, "genomes_split")
	pm.config.resources.bismark_spikein_genome = os.path.join(pm.config.resources.genomes, pm.config.resources.spikein_genome, "indexed_bismark_bt1")

	# Epilog indexes
	pm.config.resources.methpositions = os.path.join(pm.config.resources.genomes, args.genome_assembly, "indexed_epilog", args.genome_assembly + "_cg.tsv.gz")
	pm.config.resources.spikein_methpositions = os.path.join(pm.config.resources.genomes, pm.config.resources.spikein_genome, "indexed_epilog", pm.config.resources.spikein_genome + "_index.tsv.gz")

	pm.config.parameters.pipeline_outfolder = outfolder

	print(pm.config)
	tools = pm.config.tools  # Convenience alias
	param = pm.config.parameters
	resources = pm.config.resources

	# Create a ngstk object
	ngstk = pypiper.NGSTk(pm=pm)

	raw_folder = os.path.join(param.pipeline_outfolder, "raw/")
	fastq_folder = os.path.join(param.pipeline_outfolder, "fastq/")

	# Merge/Link sample input and Fastq conversion
	# These commands merge (if multiple) or link (if single) input files,
	# then convert (if necessary, for bam, fastq, or gz format) files to fastq.
	################################################################################
	pm.timestamp("### Merge/link and fastq conversion: ")

	local_input_files = ngstk.merge_or_link([args.input, args.input2], raw_folder, args.sample_name)
	cmd, out_fastq_pre, unaligned_fastq = ngstk.input_to_fastq(local_input_files, args.sample_name, args.paired_end, fastq_folder)
	pm.run(cmd, unaligned_fastq,
		follow=ngstk.check_fastq(local_input_files, unaligned_fastq, args.paired_end))
	pm.clean_add(out_fastq_pre + "*.fastq", conditional=True)

	pm.report_result("File_mb", ngstk.get_file_size(local_input_files))
	pm.report_result("Read_type", args.single_or_paired)
	pm.report_result("Genome", args.genome_assembly)



	if args.dark_bases:
		pm.timestamp("### Dark sequencing mode: ")
		cmd = tools.scripts_dir + "/darkSeqCombineReads.pl " + \
			out_fastq_pre + "_R1.fastq " +\
			out_fastq_pre + "_R2.fastq " +\
			out_fastq_pre + "_undark_R1.fastq " +\
			str(args.dark_bases)
		out_fastq_pre = out_fastq_pre + "_undark"
		unaligned_fastq = out_fastq_pre + "_R1.fastq"
		pm.run(cmd, unaligned_fastq)
		args.paired_end = False



	################################################################################
	pm.timestamp("### Adapter trimming: ")

	# We need to detect the quality encoding type of the fastq.

	if args.paired_end: 
		# Just look at the first read
		cmd = tools.python + " -u " + os.path.join(tools.scripts_dir,
		"detect_quality_code.py") + " -f " + unaligned_fastq[0]
	else:
		cmd = tools.python + " -u " + os.path.join(tools.scripts_dir,
		"detect_quality_code.py") + " -f " + unaligned_fastq
		
	encoding_string = pm.checkprint(cmd)
	if encoding_string.find("phred33") != -1:
		encoding = "phred33"
	elif encoding_string.find("phred64") != -1:
		encoding = "phred64"
	else:
		raise Exception("Unknown quality encoding type: "+encoding_string)

	if args.trimgalore:
		# Trim galore requires biopython, cutadapt modules. RSeQC as well (maybe?)
		#   --- $trim_galore -q $q --phred33 -a $a --stringency $s -e $e --length $l --output_dir $output_dir $input_fastq

		raise NotImplementedError("TrimGalore no longer supported")

		if args.paired_end:
			raise NotImplementedError("TrimGalore for PE RRBS not implemented")
		input_fastq = out_fastq_pre + "_R1.fastq "

		# With trimgalore, the output file is predetermined.
		trimmed_fastq = out_fastq_pre + "_R1_trimmed.fq"

		output_dir = fastq_folder

		# Adapter
		a = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"

		cmd = tools.trimgalore
		cmd += " -q 20"  # quality trimming
		cmd += " --" + encoding
		cmd += " -a " + a
		cmd += " --stringency 1"  # stringency: Overlap with adapter sequence required to trim a sequence
		cmd += " -e 0.1"  # Maximum allowed error rate
		cmd += " --length 16"  # Minimum Read length
		# by unchangeable default Trimmomatic discards reads of lenth 0 (produced by ILLUMINACLIP):
		cmd += " --output_dir " + output_dir + " " + input_fastq

	else:
		# Trimmomatic

		trimmed_fastq = out_fastq_pre + "_R1_trimmed.fq"
		trimmed_fastq_R2 = out_fastq_pre + "_R2_trimmed.fq"

		# REMARK AS: instead of trim_galore we try to use Trimmomatic for now
		# - we are more compatible with the other pipelines
		# - better code base, not a python wrapper of a perl script (as trim_galore)
		# - rrbs-mode not needed because biseq has the same functionality

		# REMARK NS:
		# The -Xmx4000m restricts heap memory allowed to java, and is necessary
		# to prevent java from allocating lots of memory willy-nilly
		# if it's on a machine with lots of memory, which can lead
		# to jobs getting killed by a resource manager. By default, java will
		# use more memory on systems that have more memory, leading to node-dependent
		# killing effects that are hard to trace.

		cmd = tools.java + " -Xmx" + str(pm.mem) + " -jar " + tools.trimmomatic_epignome
		if args.paired_end:
			cmd += " PE"
		else:
			cmd += " SE"
		cmd += " -" + encoding
		cmd += " -threads " + str(pm.cores) + " "
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
		cmd += "ILLUMINACLIP:" + resources.adapter_file + param.trimmomatic.illuminaclip


	# Trimming command has been constructed, using either trimming options.
	# The code to run it is the same either way:

	pm.run(cmd, trimmed_fastq,
		follow = ngstk.check_trim(trimmed_fastq, args.paired_end, trimmed_fastq_R2,
			fastqc_folder = os.path.join(param.pipeline_outfolder, "fastqc/")))

	pm.clean_add(os.path.join(fastq_folder, "*.fastq"), conditional = True)
	pm.clean_add(os.path.join(fastq_folder, "*.fq"), conditional = True)
	pm.clean_add(os.path.join(fastq_folder, "*.log"), conditional = True)
	pm.clean_add(fastq_folder, conditional = True)


	# RRBS alignment with BSMAP.
	################################################################################
	pm.timestamp("### BSMAP alignment: ")
	bsmap_folder = os.path.join(param.pipeline_outfolder, "bsmap_" + args.genome_assembly)  # e.g. bsmap_hg19
	ngstk.make_sure_path_exists(bsmap_folder)
	# no tmp folder needed for BSMAP alignment

	out_bsmap = os.path.join(bsmap_folder, args.sample_name + ".bam")

	cmd = tools.bsmap
	cmd += " -a " + out_fastq_pre + "_R1_trimmed.fq"
	if args.paired_end:
		cmd += " -b " + out_fastq_pre + "_R2_trimmed.fq"
	cmd += " -d " + resources.ref_genome_fasta
	cmd += " -o " + out_bsmap
	cmd += " " + str(param.bsmap.rrbs_mapping_mode)
	cmd += " -w " + str(param.bsmap.equal_best_hits)
	cmd += " -v " + str(param.bsmap.mismatch_rate)
	cmd += " -r " + str(param.bsmap.report_repeat)
	cmd += " -p " + str(param.bsmap.processors)
	cmd += " -n " + str(param.bsmap.map_to_strands)
	cmd += " -s " + str(param.bsmap.seed_size)
	cmd += " -S " + str(param.bsmap.random_number_seed)
	cmd += " -f " + str(param.bsmap.filter)
	cmd += " -q " + str(param.bsmap.quality_threshold)
	cmd += " -u"       # report unmapped reads (into same bam file)
	cmd += " -V 1"     # set verbosity level
	if args.paired_end:
		cmd += " -m " + str(param.bsmap.minimal_insert_size)
		cmd += " -x " + str(param.bsmap.maximal_insert_size)

	def check_bsmap():
		# BSMap apparently stores all the reads (mapped and unmapped) in
		# its output bam; to count aligned reads, then, we have to use
		# a -F4 flag (with count_mapped_reads instead of count_reads).
		ar = ngstk.count_mapped_reads(out_bsmap, args.paired_end)
		pm.report_result("Aligned_reads", ar)
		rr = float(pm.get_stat("Raw_reads"))
		tr = float(pm.get_stat("Trimmed_reads"))
		pm.report_result("Alignment_rate", round(float(ar) *
	 100 / float(tr), 2))
		pm.report_result("Total_efficiency", round(float(ar) * 100 / float(rr), 2))

		# In addition, BSMap can (if instructed by parameters) randomly assign
		# multimapping reads. It's useful to know how many in the final bam were such.
		mr = ngstk.count_multimapping_reads(out_bsmap, args.paired_end)
		pm.report_result("Multimap_reads", mr)
		pm.report_result("Multimap_rate", round(float(mr) *
	 100 / float(tr), 2))

	pm.run(cmd, out_bsmap, follow = check_bsmap)

	# bsmap2.90 requires that
	cmd2 = tools.samtools + " sort -o " + out_bsmap + " " + out_bsmap
	cmd3 = tools.samtools + " index " + out_bsmap
	pm.run([cmd2, cmd3], out_bsmap + ".bai", nofail=True)

	# Clean up big intermediate files:
	pm.clean_add(os.path.join(bsmap_folder, "*.fastq"))
	pm.clean_add(os.path.join(bsmap_folder, "*.fq"))

	# Run biseq-methcalling:
	################################################################################
	pm.timestamp("### Biseq methylation calling: ")

	# Python Software Requirements for biseq
	# REMARK AS: all packages are available via "easy_install --user <lib>"
	# pip is also a possibility if available (currently not on CeMM infrastructure)
	#
	# Direct links just in case:
	# - biopython: wget https://pypi.python.org/pypi/biopython or wget http://biopython.org/DIST/biopython-1.63.zip
	# - bitarray: wget https://pypi.python.org/packages/source/b/bitarray/bitarray-0.8.1.tar.gz
	# - guppy: wget https://pypi.python.org/packages/source/g/guppy/guppy-0.1.10.tar.gz
	# - pysam: wget https://code.google.com/p/pysam/downloads/detail?name=pysam-0.7.5.tar.gz

	biseq_output_path = os.path.join(param.pipeline_outfolder, "biseq_" + args.genome_assembly)
	biseq_output_path_web = os.path.join(biseq_output_path, "web")
	biseq_output_path_temp = os.path.join(biseq_output_path, "temp")

	ngstk.make_sure_path_exists (biseq_output_path)

	cmd = tools.python + " -u " + os.path.join(tools.scripts_dir, "biseqMethCalling.py")
	cmd += " --sampleName=" + args.sample_name
	cmd += " --alignmentFile=" + out_bsmap      # this is the absolute path to the bsmap aligned bam file
	cmd += " --methodPrefix=RRBS"
	cmd += " --rrbsMode"
	cmd += " --restrictionSite=" + str(param.biseq.restrictionSite) # specify the pattern of restriction sites
	cmd += " --checkRestriction"
	cmd += " --minFragmentLength=" + str(param.biseq.minFragmentLength)
	cmd += " --maxFragmentLength=" + str(param.biseq.maxFragmentLength)
	cmd += " --pfStatus=" + str(param.biseq.pfStatus)
	cmd += " --maxMismatches=" + str(param.biseq.maxMismatches)
	cmd += " --baseQualityScoreC=" + str(param.biseq.baseQualityScoreC)
	cmd += " --baseQualityScoreNextToC=" + str(param.biseq.baseQualityScoreNextToC)
	cmd += " --laneSpecificStatistics"
	cmd += " --bigBedFormat"
	cmd += " --deleteTemp"
	cmd += " --toolsDir=" + tools.biseq_tools
	cmd += " --outputDir=" + biseq_output_path
	cmd += " --webOutputDir=" + biseq_output_path_web
	cmd += " --tempDir=" + biseq_output_path_temp
	cmd += " --timeDelay=" + str(param.biseq.timeDelay)
	cmd += " --genomeFraction=" + str(param.biseq.genomeFraction)
	cmd += " --smartWindows=" + str(param.biseq.smartWindows)
	cmd += " --maxProcesses=" + str(param.biseq.maxProcesses)
	cmd += " --genomeDir=" + resources.genomes_split
	cmd += " --inGenome=" + args.genome_assembly
	cmd += " --outGenome=" + args.genome_assembly
	# TODO AS: Investigate what happens with biseq in the case of paired-end data

	# The dog genome has 38 chromosomes (plus one X chromosome). It's probably best to check here for these rarely used
	# reference genomes:
	# The default value for includedChromosomes is chr1-30, X, Y, Z (sufficient for human and mouse genomes)
	# REMARK NS: This is a hack to account for the way biseq restricts to
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

	pm.run([cmd, cmd2], target = biseq_finished_helper)

	# Now parse some results for pypiper result reporting.
	read_variables = ['uniqueSeqMotifCount','totalSeqMotifCount','bisulfiteConversionRate','globalMethylationMean']
	totalSeqMotifCount = 0.0
	uniqueSeqMotifCount = 0.0
	for var in read_variables:
		cmd = tools.python + " -u " + os.path.join(tools.scripts_dir, "tsv_parser.py")
		cmd += " -i " + os.path.join(biseq_output_path, "RRBS_statistics_" + args.sample_name + ".txt")
		cmd += " -c " + var
		x = pm.checkprint(cmd, shell = True)

		if var == 'totalSeqMotifCount':
			totalSeqMotifCount = float(x)
		if var == 'uniqueSeqMotifCount':
			uniqueSeqMotifCount = float(x)

		if var == 'uniqueSeqMotifCount':
			pm.report_result('Unique_CpGs', x)
		elif var == 'totalSeqMotifCount':
			pm.report_result('Total_CpGs', x)
			pm.report_result('meanCoverage', str( totalSeqMotifCount/uniqueSeqMotifCount ))
		else:
			pm.report_result(var, x)

	################################################################################
	pm.timestamp("### Make bigbed: ")
	# REMARK AS: Make bigwig uses a bismark output file. For RRBS we don't have the bismark cov file
	# (essentially a bedgraph file) which the tool bedGraphToBigWig would need
	# REMARK AS: UCSC tracks are generated by biseq-methcalling

	# First, convert the bed format into the bigBed input style.
	# This is how biseq did it, but it's actually unnecessary; instead we can just go straight off the output file.
	# Left command here for posterity.
	# awk '{ printf "%s\t%s\t%s\t\047%s%[\04720\047]\047\t%s\t%s\n", $1, $2, $3, $5/10, $5, $6 }' RRBS_cpgMethylation_01_2276TU.bed > f


	# bigbed conversion input file is the biseq methylation calls output file
	biseq_methcall_file = os.path.join(biseq_output_path, "RRBS_cpgMethylation_" + args.sample_name + ".bed")

	bigbed_output_path = os.path.join(param.pipeline_outfolder, "bigbed_" + args.genome_assembly)
	bigwig_output_path = os.path.join(param.pipeline_outfolder, "bigwig_" + args.genome_assembly)

	ngstk.make_sure_path_exists (bigbed_output_path)
	ngstk.make_sure_path_exists (bigwig_output_path)
	bigbed_output_file = os.path.join(bigbed_output_path,"RRBS_" + args.sample_name + ".bb")
	out_bedGraph = os.path.join(bigwig_output_path,"RRBS_" + args.sample_name + ".bedGraph")
	out_bigwig = os.path.join(bigwig_output_path, "RRBS_" + args.sample_name + ".bw")

	cmd = tools.bedToBigBed
	cmd += " " + biseq_methcall_file
	cmd += " " + resources.chrom_sizes
	cmd += " " + bigbed_output_file

	# REMARK NS: As of June 2015, IGV will load bigBed files for methylation
	# in a unique format if the *filename contains  "RRBS_cpgMethylation" -- see
	# https://github.com/igvteam/igv/blob/master/src/org/broad/igv/methyl/MethylTrack.java
	# This is obviously not ideal, but I will create a link with this filename
	# to the original file (even for WGBS tracks) so that you could load these into
	# IGV if you want:

	filename_hack_link_file = os.path.join(bigbed_output_path,"RRBS_cpgMethylation_" + args.sample_name + ".bb")
	cmd2 = "ln -sf " + os.path.relpath(bigbed_output_file, bigbed_output_path) + " " + filename_hack_link_file

	pm.run([cmd, cmd2], bigbed_output_file)

	# Let's also make bigwigs:

	# First convert to bedGraph
	cmd = "awk -v OFS='\t' '{ print $1, $2, $3, $5/10 }'"
	cmd += " " + biseq_methcall_file
	cmd += " > " + out_bedGraph

	pm.clean_add(out_bedGraph, conditional=True)

	cmd2 = tools.bedGraphToBigWig
	cmd2 += " " + out_bedGraph
	cmd2 += " " + resources.chrom_sizes
	cmd2 += " " + out_bigwig

	pm.run([cmd, cmd2], out_bigwig, shell=True)

	################################################################################

	def build_epilog_command(readsfile, sitesfile, context, outdir, skip_epis=False):
		ngstk.make_sure_path_exists(outdir)
		outfile = os.path.join(outdir, "all_calls.txt")
		epis_file = os.path.join(outdir, "all_epialleles.txt") \
			if param.epilog.epialleles and not skip_epis else None
		process_logfile = os.path.join(outdir, "processing_statistics.txt") \
			if param.epilog.track_process_stats else None
		return outfile, get_epi_cmd(tools.epilog, readsfile, sitesfile, outfile,
			min_rlen=param.epilog.read_length_threshold, min_qual=param.epilog.qual_threshold,
			strand_method=param.epilog.strand_method, rrbs_fill=args.rrbs_fill,
			mem_gig=param.epilog.mem_gig, context=context, cores=pm.cores,
			keep_chrom_files=param.epilog.keep_chrom_files,
			epis_file=epis_file, process_logfile=process_logfile)

	if args.epilog:
		pm.timestamp("### Epilog methylation calling: ")
		epilog_output_dir = os.path.join(
				param.pipeline_outfolder, "epilog_" + args.genome_assembly)
		epi_tgt, epi_cmd = build_epilog_command(out_bsmap, resources.methpositions,
			context=param.epilog.context, outdir=epilog_output_dir)
		pm.run(epi_cmd, target=epi_tgt, lock_name="epilog", nofail=True)

		"""
		epilog_outfile = os.path.join(
				epilog_output_dir, args.sample_name + "_epilog.bed")
		epilog_summary_file = os.path.join(
				epilog_output_dir, args.sample_name + "_epilog_summary.bed")
	
		cmd = tools.epilog
		cmd += " call"
		cmd += " --infile=" + out_bsmap  # absolute path to the bsmap aligned bam
		cmd += " --positions=" + resources.methpositions
		cmd += " --outfile=" + epilog_outfile
		cmd += " --summary-filename=" + epilog_summary_file
		cmd += " --cores=" + str(pm.cores)
		cmd += " --qual-threshold=" + str(param.epilog.qual_threshold)
		cmd += " --read-length-threshold=" + str(param.epilog.read_length_threshold)
		cmd += " --rrbs-fill=" + str(args.rrbs_fill)
		cmd += " --use-strand"    # Strand mode required because this isn't a bismark alignment.
	
		pm.run(cmd, epilog_outfile, nofail=True)
		"""


	################################################################################
	pm.timestamp("### Bismark alignment (spike-in): ")
	# currently using bowtie1 instead of bowtie2

	# get unaligned reads out of BSMAP bam
	bsmap_unalignable_bam = os.path.join(bsmap_folder, args.sample_name + "_unalignable.bam")
	pm.run(tools.samtools + " view -bh -f 4 -F 128 "+out_bsmap+" > " + bsmap_unalignable_bam, bsmap_unalignable_bam, shell=True)

	# Re-flag the unaligned paired-end reads to make them look like unpaired for Bismark
	if args.paired_end:
		bsmap_unalignable_bam_output = os.path.join(bsmap_folder, args.sample_name + "_unalignable_reflagged.bam")
		cmd = tools.python + " -u " + os.path.join(tools.scripts_dir, "pe_flag_changer.py")
		cmd += " -i " + bsmap_unalignable_bam
		cmd += " -o " + bsmap_unalignable_bam_output
		pm.run(cmd, bsmap_unalignable_bam_output)
		pm.clean_add(bsmap_unalignable_bam, conditional=True)
		bsmap_unalignable_bam = bsmap_unalignable_bam_output

	# convert BAM to fastq
	bsmap_fastq_unalignable_pre = os.path.join(bsmap_folder, args.sample_name + "_unalignable")
	bsmap_fastq_unalignable = bsmap_fastq_unalignable_pre  + "_R1.fastq"
	cmd = ngstk.bam_to_fastq(bsmap_unalignable_bam, bsmap_fastq_unalignable_pre, args.paired_end)
	pm.run(cmd, bsmap_fastq_unalignable)

	# actual spike-in analysis
	spikein_folder = os.path.join(param.pipeline_outfolder, "bismark_spikein")
	ngstk.make_sure_path_exists(spikein_folder)
	spikein_temp = os.path.join(spikein_folder, "bismark_temp")
	ngstk.make_sure_path_exists(spikein_temp)
	out_spikein_base = args.sample_name + ".spikein.aln"

	out_spikein = os.path.join(spikein_folder, out_spikein_base + ".bam")
	cmd = tools.bismark + " " + resources.bismark_spikein_genome + " "
	cmd += bsmap_fastq_unalignable_pre + "_R1.fastq"
	cmd += " --bam --unmapped"
	if (os.path.isdir(tools.bowtie1)):
		# If tools.bowtie1 is not a directory, assume owtie is in the path,
		# in which case bismark doesn't need it.
		cmd += " --path_to_bowtie " + tools.bowtie1
	#	cmd += " --bowtie2"
	cmd += " --temp_dir " + spikein_temp
	cmd += " --output_dir " + spikein_folder
	cmd += " --basename=" + out_spikein_base
	#cmd += " -p 4"
	cmd += " -n 0" #allow no mismatches

	pm.run(cmd, out_spikein, nofail=True)

	# Clean up the unmapped file which is copied from the parent
	# bismark folder to here:
	pm.clean_add(os.path.join(spikein_folder,"*.fastq"), conditional=True)
	pm.clean_add(os.path.join(spikein_folder,"*.fq"), conditional=True)
	pm.clean_add(out_spikein, conditional=True)
	pm.clean_add(spikein_temp)



	################################################################################
	pm.timestamp("### PCR duplicate removal (spike-in): ")
	# Bismark's deduplication forces output naming, how annoying.
	#out_spikein_dedup = spikein_folder + args.sample_name + ".spikein.aln.deduplicated.bam"
	out_spikein_dedup = re.sub(r'.bam$', '.deduplicated.bam', out_spikein)
	if not args.paired_end:
		cmd = tools.deduplicate_bismark + " --single "    # TODO: needs module load bismark or absolute path to this tool
		cmd += out_spikein
		cmd += " --bam"
	else:
		cmd = tools.deduplicate_bismark + " --paired "
		cmd += out_spikein
		cmd += " --bam"

	out_spikein_sorted = re.sub(r'.deduplicated.bam$', '.deduplicated.sorted.bam', out_spikein_dedup)
	cmd2 = tools.samtools + " sort " + out_spikein_dedup + " -o " + out_spikein_sorted
	cmd3 = tools.samtools + " index " + out_spikein_sorted
	pm.run([cmd, cmd2, cmd3], out_spikein_sorted + ".bai", nofail=True)
	pm.clean_add(out_spikein_dedup, conditional=False)

	# Spike-in methylation calling
	################################################################################
	pm.timestamp("### Testxmz methylation calling (spike-in): ")
	spike_chroms = ngstk.get_chrs_from_bam(out_spikein_sorted)

	for chrom in spike_chroms:
		cmd1 = tools.python + " -u " + os.path.join(tools.scripts_dir, "testxmz.py")
		cmd1 += " " + out_spikein_sorted + " " + chrom
		cmd1 += " >> " + pm.pipeline_stats_file
		pm.run(cmd1, lock_name="spikein", nofail=True)

	if args.epilog:
		# spike in conversion efficiency calculation with epilog
		pm.timestamp("### Epilog methylation calling (spike-in): ")
		epi_tgt, epi_cmd = build_epilog_command(out_bsmap, resources.spikein_methpositions,
			context="C", outdir=spikein_folder, skip_epis=True)
		pm.run(epi_cmd, target=epi_tgt, lock_name="epilog", nofail=True)

	"""
	epilog_spike_outfile=os.path.join(
			spikein_folder, args.sample_name + "_epilog.bed")
	epilog_spike_summary_file=os.path.join(
			spikein_folder, args.sample_name + "_epilog_summary.bed")
	
	cmd = tools.epilog
	cmd += " call"
	cmd += " --infile=" + out_spikein_sorted # absolute path to the bsmap aligned bam
	cmd += " --positions=" + resources.spikein_methpositions
	cmd += " --outfile=" + epilog_spike_outfile
	cmd += " --summary-filename=" + epilog_spike_summary_file
	cmd += " --cores=" + str(pm.cores)
	cmd += " --qual-threshold=30"    # quality_threshold
	cmd += " --read-length-threshold=30"    # read length cutoff
	cmd += " --rrbs-fill=0"
	
	pm.run(cmd, epilog_spike_outfile, nofail=True)
	
	# Now parse some results for pypiper result reporting.
	for chrom in spike_chroms:
		cmd = tools.python + " -u " + os.path.join(tools.scripts_dir, "tsv_parser.py")
		cmd += " -i " + os.path.join(spikein_folder, epilog_spike_summary_file)
		cmd += " -r context=C chr=" + chrom
	
		cmd_total = cmd + " -c " + "total"
		x = pm.checkprint(cmd_total, shell=True)
		pm.report_result(chrom+'_count_EL', x)
		cmd_rate = cmd + " -c " + "rate"
		x = pm.checkprint(cmd_rate, shell=True)
		pm.report_result(chrom+'_meth_EL', x)
	"""

	# PDR calculation:
	################################################################################

	# PDR not applied to PE case because bisulfiteReadConcordanceAnalysis.py is single-end only
	if not args.paired_end and args.pdr:

		pm.timestamp("### PDR (Partial Disordered Methylation) analysis")

		pdr_output_dir = os.path.join(param.pipeline_outfolder, "pdr_" + args.genome_assembly)
		ngstk.make_sure_path_exists (pdr_output_dir)

		# convert aligned bam to sam

		pdr_in_samfile = os.path.join(pdr_output_dir, args.sample_name + ".aligned.sam") # gets deleted after, see some lines below
		pm.run(tools.samtools + " view " + out_bsmap + " > " + pdr_in_samfile, pdr_in_samfile, shell=True)

		# PDR calculation:
		#
		# output files:
		pdr_bedfile=os.path.join(pdr_output_dir, args.sample_name + ".pdr.bed")

		produce_sam = False  # TODO AS: make this an option somewhere
		concordsam=os.path.join(pdr_output_dir, args.sample_name + ".concordant.sam")
		discordsam=os.path.join(pdr_output_dir, args.sample_name + ".discordant.sam")

		# command::
		cmd1 = tools.python + " -u " + os.path.join(tools.scripts_dir, "bisulfiteReadConcordanceAnalysis.py")
		cmd1 += " --infile=" + pdr_in_samfile
		cmd1 += " --outfile=" + pdr_bedfile
		cmd1 += " --skipHeaderLines=0"
		cmd1 += " --genome=" + args.genome_assembly
		cmd1 += " --genomeDir=" + resources.genomes
		cmd1 += " --minNonCpgSites=3"   # These two parameters are not relevant for PDR analysis
		cmd1 += " --minConversionRate=0.9"

		if produce_sam:
			cmd1 += " --produce_sam"
			cmd1 += " --concordantOutfile=" + concordsam
			cmd1 += " --discordantOutfile=" + discordsam
			#TODO: perhaps convert them to bam *cough*

		#call:
		pm.run(cmd1, pdr_bedfile, nofail = True)

		# delete huge input SAM file
		pm.clean_add(os.path.join(pdr_output_dir,"*.sam"), conditional = True)
		pm.clean_add(pdr_output_dir, conditional = True)

		if os.path.isfile(os.path.join(tools.scripts_dir, "extractPDR.pl")):

			pm.timestamp("### PDR (Perl version by Kendell)")
			pdr_out = os.path.join(pdr_output_dir, args.sample_name + ".pdr")

			cmd = "perl " + os.path.join(tools.scripts_dir, "extractPDR.pl")
			cmd += " " + os.path.join(pdr_output_dir, args.sample_name) + " " + args.genome_assembly + ""
			cmd += " " + out_bsmap

			pm.run(cmd, target = pdr_out, nofail = True)



	# Final sorting and indexing
	################################################################################
	# create sorted and indexed BAM files for visualization and analysis
	# bsmap already outputs a sorted and indexed bam file

	# Cleanup
	################################################################################
	pm.stop_pipeline()


if __name__ == "__main__":
	import sys
	main(sys.argv[1:])
