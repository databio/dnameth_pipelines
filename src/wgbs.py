#!/usr/bin/env python

"""
WGBS pipeline
"""

__author__ = "Nathan Sheffield"
__email__ = "nathan@code.databio.org"
__credits__ = ["Charles Dietz", "Johanna Klughammer", "Christoph Bock", "Andreas Schoenegger"]
__license__ = "GPL3"
__version__ = "0.2.0-dev"

from argparse import ArgumentParser
import os
import re
import subprocess
import pypiper


parser = ArgumentParser(description='Pipeline')

# First, add arguments from Pypiper, including
# 1. pypiper options, 2. looper connections, 3. common options,
# using the all_args=True flag (you can customize this).
# Adds options including; for details, see:
# http://github.com/epigen/pypiper/command_line_args.md
parser = pypiper.add_pypiper_args(parser, all_args=True)


parser.add_argument('-e', '--epilog', dest='epilog', action="store_true", default=False,
	help='Use epilog for meth calling?')

parser.add_argument('--single2', dest='single2', action="store_true", default=False,
	help='Single secondary mode: any reads not mapping in paired-end mode will \
			be aligned using single-end mode, and then analyzed. Only valid for \
			paired-end mode. ')

args = parser.parse_args()

if args.single_or_paired == "paired":
	args.paired_end = True
else:
	args.paired_end = False


# Merging
################################################################################
# If 2 input files are given, then these are to be merged.
# Must be done here to initialize the sample name correctly

merge = False
if len(args.input) > 1:
	merge = True
	if args.sample_name == "default":
		args.sample_name = "merged"
else:
	if args.sample_name == "default":
		# Default sample name is derived from the input file
		args.sample_name = os.path.splitext(os.path.basename(args.input[0]))[0]

if not args.input:
	parser.print_help()
	raise SystemExit

# Create a PipelineManager object and start the pipeline
outfolder = os.path.abspath(os.path.join(args.output_parent, args.sample_name))
pm = pypiper.PipelineManager(name="WGBS", outfolder=outfolder, args=args, version=__version__)

# Set up a few additional paths not in the config file
pm.config.tools.scripts_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)), "tools")
pm.config.resources.ref_genome_fasta = os.path.join(pm.config.resources.genomes, args.genome_assembly, args.genome_assembly + ".fa")
pm.config.resources.chrom_sizes = os.path.join(pm.config.resources.genomes, args.genome_assembly, args.genome_assembly + ".chromSizes")
pm.config.resources.genomes_split = os.path.join(pm.config.resources.resources, "genomes_split")
try:
	pm.config.resources.bismark_spikein_genome = os.path.join(pm.config.resources.genomes, pm.config.resources.spikein_genome, "indexed_bismark_bt1")
except:
	pm.config.resources.bismark_spikein_genome = None

pm.config.resources.bismark_indexed_genome = os.path.join(pm.config.resources.genomes, args.genome_assembly, "indexed_bismark_bt2")

# Epilog indexes
pm.config.resources.methpositions = os.path.join(pm.config.resources.genomes, args.genome_assembly, "indexed_epilog", args.genome_assembly + "_cg.tsv.gz")

if pm.config.resources.bismark_spikein_genome:
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

# Adapter trimming
################################################################################
pm.timestamp("### Adapter trimming: ")

# We need to detect the quality encoding type of the fastq.
cmd = tools.python + " -u " + os.path.join(tools.scripts_dir, "detect_quality_code.py") + " -f " + unaligned_fastq
encoding_string = pm.checkprint(cmd)
if encoding_string.find("phred33") != -1:
	encoding = "phred33"
elif encoding_string.find("phred64") != -1:
	encoding = "phred64"
else:
	raise Exception("Unknown quality encoding type: "+encoding_string)


trimmed_fastq = out_fastq_pre + "_R1_trimmed.fq"
trimmed_fastq_R2 = out_fastq_pre + "_R2_trimmed.fq"

cmd = tools.java + " -Xmx" + str(pm.mem) + " -jar " + tools.trimmomatic
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
cmd += " " + param.trimmomatic.trimsteps
cmd += " ILLUMINACLIP:" + resources.adapter_file + param.trimmomatic.illuminaclip

pm.run(cmd, trimmed_fastq, 
	follow = ngstk.check_trim(trimmed_fastq, trimmed_fastq_R2, args.paired_end,
		fastqc_folder = os.path.join(param.pipeline_outfolder, "fastqc/")))

pm.clean_add(os.path.join(fastq_folder, "*.fastq"), conditional=True)
pm.clean_add(os.path.join(fastq_folder, "*.fq"), conditional=True)
pm.clean_add(os.path.join(fastq_folder, "*.log"), conditional=True)
pm.clean_add(fastq_folder, conditional=True)


# WGBS alignment with bismark.
################################################################################
pm.timestamp("### Bismark alignment: ")
# Bismark will start multiple instances of bowtie, so we have to split
# the alotted cores among the instances. Otherwise we will use 2x or 4x the number
# of cores that we aresupposed to. It will start 2 threads in
# normal mode, and 4 in --non-directional mode.

if param.bismark.nondirectional:
	bismark_bowtie_threads = 4
else:
	bismark_bowtie_threads = 2

bismark_cores = int(pm.cores) // bismark_bowtie_threads

if int(pm.cores) % bismark_bowtie_threads != 0:
	print("inefficient core request; make divisible by " + 	str(bismark_bowtie_threads))

bismark_folder = os.path.join(param.pipeline_outfolder, "bismark_" + args.genome_assembly )
ngstk.make_sure_path_exists(bismark_folder)
bismark_temp = os.path.join(bismark_folder, "bismark_temp" )
ngstk.make_sure_path_exists(bismark_temp)

if args.paired_end:
	out_bismark = os.path.join(bismark_folder, args.sample_name + "_pe.bam")
else:
	out_bismark = os.path.join(bismark_folder, args.sample_name + ".bam")

cmd = tools.bismark + " " + resources.bismark_indexed_genome + " "
if args.paired_end:
	cmd += " --1 " + out_fastq_pre + "_R1_trimmed.fq"
	cmd += " --2 " + out_fastq_pre + "_R2_trimmed.fq"
else:
	cmd += out_fastq_pre + "_R1_trimmed.fq"
cmd += " --bam --unmapped"
# Bowtie may be specified in raw form to indicate presence on path.
if tools.bowtie2 != "bowtie2":
	cmd += " --path_to_bowtie " + tools.bowtie2
cmd += " --bowtie2"
cmd += " --temp_dir " + bismark_temp
cmd += " --output_dir " + bismark_folder
if args.paired_end:
	cmd += " --minins 0"
	cmd += " --maxins " + str(param.bismark.maxins)
cmd += " -p " + str(bismark_cores) # Number of processors
cmd += " --basename=" + args.sample_name

# By default, BS-seq libraries are directional, but this can be turned off
# in bismark for non-directional protocols
if param.bismark.nondirectional:
	cmd += " --non_directional"

def check_bismark():
	ar = ngstk.count_mapped_reads(out_bismark, args.paired_end)
	pm.report_result("Aligned_reads", ar)
	rr = float(pm.get_stat("Raw_reads"))
	tr = float(pm.get_stat("Trimmed_reads"))
	pm.report_result("Alignment_rate", round(float(ar) *
 100 / float(tr), 2))
	pm.report_result("Total_efficiency", round(float(ar) * 100 / float(rr), 2))

	mr = ngstk.count_multimapping_reads(out_bismark, args.paired_end)
	pm.report_result("Multimap_reads", mr)
	pm.report_result("Multimap_rate", round(float(mr) *
 100 / float(tr), 2))

pm.run(cmd, out_bismark, follow = check_bismark)


# Secondary single mode:
# align unmapped in single end mode?
if args.paired_end and args.single2:
	pm.timestamp("### Bismark secondary single-end alignment: ")
	out_bismark_se =[]
	for read_n in ["1", "2"]:  # Align each read in single end mode
		read_string = "R" + str(read_n)
		bismark2_folder = os.path.join(bismark_folder, "se" + str(read_string))
		ngstk.make_sure_path_exists(bismark2_folder)
		bismark2_temp = os.path.join(bismark2_folder, "bismark2_temp" )
		ngstk.make_sure_path_exists(bismark2_temp)
		out_bismark2 = os.path.join(bismark2_folder, args.sample_name + read_string +  ".bam")

		unmapped_reads_pre = os.path.join(bismark_folder, args.sample_name)

		cmd = tools.bismark + " " + resources.bismark_indexed_genome + " "
		cmd += unmapped_reads_pre + "_unmapped_reads_" + str(read_n) + ".fq"
		cmd += " --bam --unmapped"
		# Bowtie may be specified in raw form to indicate presence on path.
		if tools.bowtie2 != "bowtie2":
			cmd += " --path_to_bowtie " + tools.bowtie2
		cmd += " --bowtie2"
		cmd += " --temp_dir " + bismark2_temp
		cmd += " --output_dir " + bismark2_folder
		cmd += " --basename="  + args.sample_name + read_string
		cmd += " -p " + str(bismark_cores)
		if param.bismark.nondirectional:
			cmd += " --non_directional"

		pm.run(cmd, out_bismark2)
		out_bismark_se.append(out_bismark2)

	# Now merge, sort, and analyze the single-end data
	merged_bismark = args.sample_name + "_SEmerged.bam"
	output_merge = os.path.join(bismark_folder, merged_bismark)
	cmd = ngstk.merge_bams(out_bismark_se, output_merge, in_sorted = "FALSE", 
		tmp_dir = resources.tmp_dir)

	pm.run(cmd, output_merge)
	# Sort by read name
	sorted_bismark = args.sample_name + "_SEsorted.bam"
	output_sort = os.path.join(bismark_folder, sorted_bismark)

	cmd = tools.samtools + " sort -n -o " + output_merge + " " + output_sort
	pm.run(cmd, output_sort)

	cmd = tools.python + " -u " + os.path.join(tools.scripts_dir, "rematch_pairs.py")
	cmd += " -i " + output_sort

	pm.run(cmd, lock_name="rematch")




pm.timestamp("### PCR duplicate removal: ")
# Bismark's deduplication forces output naming, how annoying.
out_dedup = os.path.join(bismark_folder, args.sample_name + "_pe.deduplicated.bam")
out_dedup = re.sub(r'.bam$', '.deduplicated.bam', out_bismark)

cmd = tools.deduplicate_bismark
if args.paired_end:
	cmd += " --paired "
else:
	cmd += " --single "
cmd += out_bismark
cmd += " --bam"

pm.run(cmd, out_dedup, follow = lambda:
	pm.report_result("Deduplicated_reads", ngstk.count_reads(out_dedup, args.paired_end)))


pm.timestamp("### Aligned read filtering: ")

# convert bam file into sam file and sort again to
# compensate for a sorting issue of "deduplicate_bismark"
sam_temp = os.path.join(bismark_folder, "sam_temp")
ngstk.make_sure_path_exists(sam_temp)
out_sam = os.path.join(bismark_folder, args.sample_name + ".aln.deduplicated.sam")
#Is this an old version of samtools?
#cmd = tools.samtools + " sort -n -o " + out_dedup + " " + out_dedup.replace(".bam", "_sorted") + " | " + tools.samtools + " view -h - >" + out_sam
#cmd = tools.samtools + " sort -n " + out_dedup + " " + " | " + tools.samtools + " view -h - >" + out_sam
cmd = tools.samtools + " sort -n " + out_dedup + " -o " + out_sam


pm.run(cmd, out_sam, shell=True)

#sorted file same size as presorted?
#pm.report_result("Filtered_reads", ngstk.count_reads(out_sam_filter, args.paired_end)) = ngstk.count_reads(out_sam, args.paired_end)
#if sorted_reads != deduplicated_reads:
#	raise Exception("Sorted size doesn't match deduplicated size.")

out_sam_filter = os.path.join(bismark_folder, args.sample_name + ".aln.dedup.filt.sam")

headerLines = subprocess.check_output(tools.samtools + " view -SH " + out_sam + "| wc -l", shell=True).strip()
cmd = tools.python + " " + os.path.join(tools.scripts_dir, "bisulfiteReadFiltering_forRNA.py")
cmd += " --infile=" + out_sam
cmd += " --outfile=" + out_sam_filter
cmd += " --skipHeaderLines=" + headerLines
cmd += " --genome=" + args.genome_assembly
cmd += " --genomeDir=" + resources.genomes
cmd += " --minNonCpgSites=3"
cmd += " --minConversionRate=0.9"
if args.paired_end:
	cmd = cmd + " --pairedEnd"

pm.run(cmd, out_sam_filter, follow=lambda:
	pm.report_result("Filtered_reads", ngstk.count_reads(out_sam_filter, args.paired_end)))

# Clean up all intermediates
pm.clean_add(out_bismark) # initial mapped bam file
pm.clean_add(os.path.join(bismark_folder, "*.fastq"))
pm.clean_add(os.path.join(bismark_folder, "*.fq"))
pm.clean_add(out_dedup) # deduplicated bam file
pm.clean_add(out_sam) # dedup conversion to sam
pm.clean_add(out_sam_filter) # after filtering


# Methylation extractor
################################################################################
# REMARK NS:
# Bismark methylation extractor produces various outpus, but unfortunately none
# are great. The default "coverage" (.bismark.cov) file is thus:
# chr	start	stop	meth	methylated	unmethylated
# chr17	4890653	4890653	100	1	0
# chr17	5334751	5334751	100	1	0
# This output lacks strand information, so you don't know if the coordinate is
# pointing to a C or G on the + strand unless you look it up in the reference genome.
# The "cytosine_report" file has all the info, but includes an entry for every
# CpG, covered or not:
# chr17	3000204	+	0	0	CG	CGT
# chr17	3000205	-	0	0	CG	CGA
# chr17	4890653	-	1	0	CG	CGA
# Solution: Use the cytosine_report file, and filter out any uncovered reads.

pm.timestamp("### Methylation calling (bismark extractor): ")

extract_dir = os.path.join(bismark_folder, "extractor")
ngstk.make_sure_path_exists(extract_dir)
out_extractor = os.path.join(extract_dir, re.sub(r'.sam$', '.bismark.cov', os.path.basename(out_sam_filter)))
out_cpg_report = re.sub(r'.bismark.cov$', '.CpG_report.txt.gz', out_extractor)

cmd = tools.bismark_methylation_extractor
if args.paired_end:
	cmd += " --paired-end --no_overlap"
else:
	cmd += " --single-end"
cmd += " --report"
cmd += " --bedGraph"
cmd += " --merge_non_CpG"
cmd += " --cytosine_report"
cmd += " --genome_folder " + resources.bismark_indexed_genome
cmd += " --gzip"
cmd += " --output " + extract_dir
cmd += " " + out_sam_filter

pm.run(cmd,  out_cpg_report)

# TODO: make these boolean flags options to the pipeline
keep_bismark_report = True
keep_non_standard_chromosomes = False
adjust_minus_strand = True

# prepare outputs:
out_cpg_report_filt = re.sub(r'.CpG_report.txt.gz$', '.CpG_report_filt.txt', out_cpg_report)
out_cpg_report_filt_cov = re.sub(r'.CpG_report.txt.gz$', '.CpG_report_filt.cov', out_cpg_report)

# remove uncovered regions:
# Update to Bismark version 17 now gzips this output.
cmd = "gzip -c -d "
cmd += " " + out_cpg_report
cmd = " | awk '{ if ($4+$5 > 0) print; }'"
cmd += " > " + out_cpg_report_filt
pm.run(cmd,  out_cpg_report_filt, shell=True)

# convert the bismark report to the simpler coverage format and adjust the coordinates
# of CpG's on the reverse strand while doing so (by substracting 1 from the start):
cmd = tools.Rscript + " " + os.path.join(os.path.dirname(os.path.dirname(os.path.realpath(__file__))), "scripts", "convertBismarkReport.R") # disable coverage filter, because we have already used `awk` to achieve this result
cmd += " --formats=cov,min,gibberish"
cmd += " --noCovFilter"
if keep_non_standard_chromosomes:
	cmd += " --noChromFilter"
if not adjust_minus_strand:
	cmd += " --noAdjustMinusStrand"
cmd += " -i " + out_cpg_report_filt
pm.run(cmd,  out_cpg_report_filt_cov)


# tidy up:
if not keep_bismark_report:
	pm.clean_add(out_cpg_report_filt)


# Make bigwig
################################################################################
pm.timestamp("### Make bigwig: ")

bedGraph = re.sub(".bismark.cov$", ".bedGraph", out_extractor)
sort_bedGraph = re.sub(".bedGraph$", ".sort.bedGraph", bedGraph)
out_bigwig = re.sub(".bedGraph$", ".bw", bedGraph)
cmd1 = "sed '1d' " + bedGraph + " | LC_COLLATE=C sort -k1,1 -k2,2n - " + " > " + sort_bedGraph
cmd2 = tools.bedGraphToBigWig + " " + sort_bedGraph + " " + resources.chrom_sizes
cmd2 += " " + out_bigwig

pm.run([cmd1, cmd2], out_bigwig)

################################################################################

if args.epilog:
	# out_bismark must be indexed in order for epilog to use it
	# should we do this on out_ded
	out_dedup_sorted = re.sub(r'.bam$',"_sort.bam", out_dedup)
	cmd2 = tools.samtools + " sort -@ " + str(pm.cores) + " -o " + out_dedup_sorted + " " + out_dedup
	cmd3 = tools.samtools + " index " + out_dedup_sorted
	pm.run([cmd2, cmd3], out_dedup_sorted + ".bai")

	pm.timestamp("### Epilog Methcalling: ")
	epilog_output_dir = os.path.join(
			param.pipeline_outfolder, "epilog_" + args.genome_assembly)
	ngstk.make_sure_path_exists(epilog_output_dir)
	epilog_outfile = os.path.join(
			epilog_output_dir, args.sample_name + "_epilog.bed")
	epilog_summary_file = os.path.join(
			epilog_output_dir, args.sample_name + "_epilog_summary.bed")

	cmd = tools.epilog
	cmd += " --infile=" + out_dedup_sorted  # absolute path to the aligned bam
	cmd += " --positions=" + resources.methpositions
	cmd += " --outfile=" + epilog_outfile
	cmd += " --summary-filename=" + epilog_summary_file
	cmd += " --cores=" + str(pm.cores)
	cmd += " --rrbs-fill-count=0"    # Turn off RRBS mode

	pm.run(cmd, epilog_outfile, nofail=True)

# Spike-in alignment
################################################################################
# currently using bowtie1 instead of bowtie2
if resources.bismark_spikein_genome:
	pm.timestamp("### Bismark spike-in alignment: ")
	spikein_folder = os.path.join(param.pipeline_outfolder, "bismark_spikein" )
	ngstk.make_sure_path_exists(spikein_folder)
	spikein_temp = os.path.join(spikein_folder, "bismark_temp" )
	ngstk.make_sure_path_exists(spikein_temp)
	out_spikein_base = args.sample_name + ".spikein.aln"

	#out_spikein = spikein_folder + args.sample_name + "_R1_trimmed.fastq_unmapped_reads_1.fq_bismark_pe.bam"

	unmapped_reads_pre = os.path.join(bismark_folder, args.sample_name)
	if args.paired_end:
		out_spikein = os.path.join(spikein_folder, out_spikein_base + "_pe.bam")
	else:
		out_spikein = os.path.join(spikein_folder, out_spikein_base + ".bam")
	cmd = tools.bismark + " " + resources.bismark_spikein_genome + " "
	if args.paired_end:
		cmd += " --1 " + unmapped_reads_pre + "_unmapped_reads_1.fq"
		cmd += " --2 " + unmapped_reads_pre + "_unmapped_reads_2.fq"
	else:
		cmd += unmapped_reads_pre + "_unmapped_reads.fq"
	cmd += " --bam --unmapped"
	# Bowtie may be specified in raw form to indicate presence on path.
	if tools.bowtie1 != "bowtie":
		cmd += " --path_to_bowtie " + tools.bowtie1
	#cmd += " --bowtie2"
	cmd += " --temp_dir " + spikein_temp
	cmd += " --output_dir " + spikein_folder
	if args.paired_end:
		cmd += " --minins 0"
		cmd += " --maxins " + str(param.bismark.maxins)
	cmd += " --basename="  + out_spikein_base
	if param.bismark.nondirectional:
		cmd += " --non_directional"


	pm.run(cmd, out_spikein, nofail=True)
	# Clean up the unmapped file which is copied from the parent
	# bismark folder to here:
	pm.clean_add(os.path.join(spikein_folder, "*.fq"), conditional=False)
	pm.clean_add(spikein_temp)


	pm.timestamp("### PCR duplicate removal (Spike-in): ")
	# Bismark's deduplication forces output naming, how annoying.
	#out_spikein_dedup = spikein_folder + args.sample_name + ".spikein.aln.deduplicated.bam"
	out_spikein_dedup = re.sub(r'.bam$', '.deduplicated.bam', out_spikein)


	cmd = tools.deduplicate_bismark
	if args.paired_end:
		cmd += " --paired "
	else:
		cmd += " --single "
	cmd += out_spikein
	cmd += " --bam"

	out_spikein_sorted = out_spikein_dedup.replace('.deduplicated.bam', '.deduplicated.sorted')
	cmd2 = tools.samtools + " sort " + out_spikein_dedup + " -o " + out_spikein_sorted
	cmd3 = tools.samtools + " index " + out_spikein_sorted + ".bam"
	cmd4 = "rm " + out_spikein_dedup
	pm.run([cmd, cmd2, cmd3, cmd4], out_spikein_sorted +".bam.bai", nofail=True)

	# Spike-in methylation calling
	################################################################################
	pm.timestamp("### Methylation calling (testxmz) Spike-in: ")
	spike_chroms = ngstk.get_chrs_from_bam(out_spikein_sorted + ".bam")

	for chrom in spike_chroms:
		cmd1 = tools.python + " -u " + os.path.join(tools.scripts_dir, "testxmz.py")
		cmd1 += " " + out_spikein_sorted + ".bam" + " " + chrom
		cmd1 += " >> " + pm.pipeline_stats_file
		pm.callprint(cmd1, shell=True, nofail=True)


	# spike in conversion efficiency calculation with epilog
	epilog_output_dir = os.path.join(
			param.pipeline_outfolder, "epilog_" + args.genome_assembly)
	ngstk.make_sure_path_exists (epilog_output_dir)
	epilog_spike_outfile=os.path.join(
			spikein_folder, args.sample_name + "_epilog.bed")
	epilog_spike_summary_file=os.path.join(
			spikein_folder, args.sample_name + "_epilog_summary.bed")


	cmd = tools.epilog
	cmd += " --infile=" + out_spikein_sorted + ".bam"  # absolute path to the bsmap aligned bam
	cmd += " --positions=" + resources.spikein_methpositions
	cmd += " --outfile=" + epilog_spike_outfile
	cmd += " --summary=" + epilog_spike_summary_file
	cmd += " --cores=" + str(pm.cores)
	cmd += " --qual-threshold=30"
	cmd += " --read-length-threshold=30"
	cmd += " --rrbs-fill-count=0"    # no rrbs mode for WGBS pipeline

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


# Final sorting and indexing
################################################################################
# create sorted and indexed BAM files for visualization and analysis
pm.timestamp("### Final sorting and indexing: ")

#out_header = bismark_folder + args.sample_name + ".reheader.bam"
out_final = os.path.join(bismark_folder, args.sample_name + ".final.bam")
# temp_folder = os.path.join(bismark_folder, "tmp")

# # Sort
# cmd = tools.java + " -Xmx" + str(pm.mem)
# # This sort can run out of temp space on big jobs; this puts the temp to a
# # local spot.
# cmd += " -Djava.io.tmpdir=" + str(temp_folder)
# cmd += " -jar " + tools.picard + " SortSam"
# cmd += " I=" + out_sam_filter
# cmd += " O=" + out_final
# cmd += " SORT_ORDER=coordinate"
# cmd += " VALIDATION_STRINGENCY=SILENT"
# cmd += " CREATE_INDEX=true"
# pm.run(cmd, out_final, lock_name="final_sorting")

cmd = tools.samtools + " sort -@ " + str(pm.cores) + " " + out_sam_filter + " -o " + out_final
cmd2 = tools.samtools + " index " + out_final
pm.run([cmd, cmd2], out_final + ".bai")

# Cleanup
################################################################################
# remove temporary folders
pm.clean_add(bismark_temp)
pm.clean_add(sam_temp)
pm.stop_pipeline()
