#!/usr/bin/env python

"""
WGBS pipeline
"""

__author__ = "Nathan Sheffield"
__email__ = "nathan@code.databio.org"
__credits__ = ["Charles Dietz", "Johanna Klughammer", "Christoph Bock", "Andreas Schoenegger"]
__license__ = "GPL3"
__version__ = "0.1"
__status__ = "Development"

from argparse import ArgumentParser
import os, re
import sys
import subprocess
import yaml
import pypiper

parser = ArgumentParser(description='Pipeline')

# First, add arguments from Pypiper, including
# 1. pypiper options, 2. looper connections, 3. common options,
# using the all_args=True flag (you can customize this).
# Adds options including; for details, see:
# http://github.com/epigen/pypiper/command_line_args.md
parser = pypiper.add_pypiper_args(parser, all_args=True)

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


# Create a PipelineManager object and start the pipeline
pm = pypiper.PipelineManager(name = "RRBS", outfolder = os.path.abspath(os.path.join(args.output_parent, args.sample_name)), args = args)

# Set up a few additional paths not in the config file
pm.config.tools.scripts_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)), "tools")
pm.config.resources.adapter_file = os.path.join(pm.config.resources.resources, "adapters", "epignome_adapters_2_add.fa")
pm.config.resources.rrbs_adapter_file = os.path.join(pm.config.resources.resources, "adapters", "RRBS_adapters.fa")
pm.config.resources.ref_genome = os.path.join(pm.config.resources.resources, "genomes")
pm.config.resources.ref_genome_fasta = os.path.join(pm.config.resources.resources, "genomes", args.genome_assembly, args.genome_assembly + ".fa")
pm.config.resources.chrom_sizes = os.path.join(pm.config.resources.resources, "genomes", args.genome_assembly, args.genome_assembly + ".chromSizes")
pm.config.resources.genomes_split = os.path.join(pm.config.resources.resources, "genomes_split")

pm.config.resources.methpositions = os.path.join(pm.config.resources.resources, "regions", "cgs", args.genome_assembly + ".cgs.txt")

pm.config.resources.bismark_spikein_genome = os.path.join(pm.config.resources.resources, "genomes", "meth_spikein_k1_k3", "indexed_bismark_bt1")
pm.config.resources.bismark_indexed_genome = os.path.join(pm.config.resources.resources, "genomes", args.genome_assembly, "indexed_bismark_bt2")

pm.config.parameters.pipeline_outfolder = os.path.abspath(os.path.join(args.output_parent, args.sample_name))

print(pm.config)
tools = pm.config.tools  # Convenience alias
param = pm.config.parameters
resources = pm.config.resources

# Create a ngstk object
myngstk = pypiper.NGSTk(args.config_file)

myngstk.make_sure_path_exists(os.path.join(param.pipeline_outfolder, "unmapped_bam"))

if merge:
	if not args.input.endswith(".bam"):
		raise NotImplementedError("Currently we can only merge bam inputs")
	merge_folder = os.path.join(param.pipeline_outfolder, "unmapped_bam")
	sample_merged_bam = args.sample_name + ".merged.bam"
	output_merge = os.path.join(merge_folder, sample_merged_bam)
	cmd = myngstk.merge_bams(args.input, output_merge)

	pm.run(cmd, output_merge)
	args.input = os.path.join(param.pipeline_outfolder, "unmapped_bam", sample_merged_bam)  #update unmapped bam reference
	local_unmapped_bam_abs = os.path.join(param.pipeline_outfolder, "unmapped_bam", sample_merged_bam)
	input_ext = ".bam"
else:
	# Link the file into the unmapped_bam directory
	args.input = args.input[0]

	if not os.path.isabs(args.input):
		args.input = os.path.abspath(args.input)
		print args.input

	if args.input.endswith(".bam"):
		input_ext = ".bam"
	elif args.input.endswith(".fastq.gz"):
		input_ext = ".fastq.gz"
	else:
		raise NotImplementedError("This pipeline can only deal with .bam or .fastq.gz files")

	local_unmapped_bam_abs = os.path.join(param.pipeline_outfolder, "unmapped_bam", args.sample_name + input_ext)
	pm.callprint("ln -sf " + args.input + " " + local_unmapped_bam_abs, shell=True)

#check for file exists:
if not os.path.exists(local_unmapped_bam_abs):
	raise Exception(local_unmapped_bam_abs + " is not a file")

# Record file size of input file

cmd = "stat -Lc '%s' " + local_unmapped_bam_abs
input_size = pm.checkprint(cmd)
input_size = float(input_size.replace("'",""))

pm.report_result("File_mb", round((input_size/1024)/1024,2))
pm.report_result("Read_type",args.single_or_paired)
pm.report_result("Genome",args.genome_assembly)

# Fastq conversion
################################################################################
pm.timestamp("### Fastq conversion: ")
fastq_folder = os.path.join(param.pipeline_outfolder, "fastq")
out_fastq_pre = os.path.join(fastq_folder, args.sample_name)
unaligned_fastq = out_fastq_pre + "_R1.fastq"

def check_fastq():
	raw_reads = myngstk.count_reads(local_unmapped_bam_abs, args.paired_end)
	pm.report_result("Raw_reads", str(raw_reads))
	fastq_reads = myngstk.count_reads(unaligned_fastq, args.paired_end)
	pm.report_result("Fastq_reads", fastq_reads)
	fail_filter_reads = myngstk.count_fail_reads(local_unmapped_bam_abs, args.paired_end)
	pf_reads = int(raw_reads) - int(fail_filter_reads)
	pm.report_result("PF_reads", str(pf_reads))
	if fastq_reads != int(raw_reads):
		raise Exception("Fastq conversion error? Number of reads doesn't match unaligned bam")

if input_ext ==".bam":
	print("Found bam file")
	cmd = myngstk.bam_to_fastq(local_unmapped_bam_abs, out_fastq_pre, args.paired_end)
	pm.run(cmd, unaligned_fastq, follow=check_fastq)
elif input_ext == ".fastq.gz":
	print("Found gz fastq file")
	cmd = "gunzip -c " + local_unmapped_bam_abs + " > " + unaligned_fastq
	myngstk.make_sure_path_exists(fastq_folder)
	pm.run(cmd, unaligned_fastq, shell=True, follow=lambda:
		pm.report_result("Fastq_reads",  myngstk.count_reads(unaligned_fastq, args.paired_end)))

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

cmd = tools.java + " -Xmx" + str(param.trimmomatic.memory) + "g -jar " + tools.trimmomatic_epignome
if args.paired_end:
	cmd += " PE"
else:
	cmd += " SE"
cmd += " -" + encoding
cmd += " -threads " + str(param.trimmomatic.threads) + " "
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
cmd += " HEADCROP:6 ILLUMINACLIP:" + resources.adapter_file + param.trimmomatic.illuminaclip

pm.run(cmd, trimmed_fastq, follow=lambda:
	pm.report_result("Trimmed_reads",  myngstk.count_reads(trimmed_fastq, args.paired_end)))

pm.clean_add(os.path.join(fastq_folder, "*.fastq"), conditional=True)
pm.clean_add(os.path.join(fastq_folder, "*.fq"), conditional=True)
pm.clean_add(os.path.join(fastq_folder, "*.log"), conditional=True)
pm.clean_add(fastq_folder, conditional=True)


# WGBS alignment with bismark.
################################################################################
pm.timestamp("### Bismark alignment: ")
bismark_folder = os.path.join(param.pipeline_outfolder, "bismark_" + args.genome_assembly )
myngstk.make_sure_path_exists(bismark_folder)
bismark_temp = os.path.join(bismark_folder, "bismark_temp" )
myngstk.make_sure_path_exists(bismark_temp)

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
cmd += " --path_to_bowtie " + tools.bowtie2
cmd += " --bowtie2"
cmd += " --temp_dir " + bismark_temp
cmd += " --output_dir " + bismark_folder
if args.paired_end:
	cmd += " --minins 0"
	cmd += " --maxins 5000"
cmd += " -p 8 " # Number of processors
cmd += " --basename=" +args.sample_name

pm.run(cmd, out_bismark)

def check_bismark():
	x = myngstk.count_mapped_reads(out_bsmap, args.paired_end)
	pm.report_result("Aligned_reads", x)
	x = myngstk.count_multimapping_reads(out_bsmap, args.paired_end)
	pm.report_result("Multimap_reads", x)


pm.run(cmd, out_bismark, follow=check_bismark)

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

pm.run(cmd, out_dedup, follow=lambda:
	pm.report_result("Deduplicated_reads", myngstk.count_reads(out_dedup, args.paired_end)))


pm.timestamp("### Aligned read filtering: ")

# convert bam file into sam file and sort again to
# compensate for a sorting issue of "deduplicate_bismark"
sam_temp = os.path.join(bismark_folder, "sam_temp")
myngstk.make_sure_path_exists(sam_temp)
out_sam = os.path.join(bismark_folder, args.sample_name + ".aln.deduplicated.sam")
cmd = tools.samtools + " sort -n -o " + out_dedup + " " + out_dedup.replace(".bam", "_sorted") + " | " + tools.samtools + " view -h - >" + out_sam

pm.run(cmd, out_sam, shell=True)

#sorted file same size as presorted?
#pm.report_result("Filtered_reads", myngstk.count_reads(out_sam_filter, args.paired_end)) = myngstk.count_reads(out_sam, args.paired_end)
#if sorted_reads != deduplicated_reads:
#	raise Exception("Sorted size doesn't match deduplicated size.")

out_sam_filter = os.path.join(bismark_folder, args.sample_name + ".aln.dedup.filt.sam")

headerLines = subprocess.check_output(tools.samtools + " view -SH " + out_sam + "| wc -l", shell=True).strip()
cmd = tools.python + " " + os.path.join(tools.scripts_dir, "bisulfiteReadFiltering_forRNA.py")
cmd += " --infile=" + out_sam
cmd += " --outfile=" + out_sam_filter
cmd += " --skipHeaderLines=" + headerLines
cmd += " --genome=" + args.genome_assembly
cmd += " --genomeDir=" + resources.ref_genome
cmd += " --minNonCpgSites=3"
cmd += " --minConversionRate=0.9"
if args.paired_end:
	cmd = cmd + " --pairedEnd"

pm.run(cmd, out_sam_filter, follow=lambda:
	pm.report_result("Filtered_reads", myngstk.count_reads(out_sam_filter, args.paired_end)))

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
# pointing to a C or G on the + strand unles you look it up in the reference genome.
# The "cytosine_report" file has all the info, but includes an entry for every
# CpG, covered or not:
# chr17	3000204	+	0	0	CG	CGT
# chr17	3000205	-	0	0	CG	CGA
# chr17	4890653	-	1	0	CG	CGA
# Solution: Use the cytosine_report file, and filter out any uncovered reads.

pm.timestamp("### Methylation calling (bismark extractor): ")

extract_dir = os.path.join(bismark_folder, "extractor")
myngstk.make_sure_path_exists(extract_dir)
out_extractor = os.path.join(extract_dir, re.sub(r'.sam$', '.bismark.cov', os.path.basename(out_sam_filter)))
out_cpg_report = re.sub(r'.bismark.cov$', '.CpG_report.txt', out_extractor)

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
out_cpg_report_filt = re.sub(r'.CpG_report.txt$', '.CpG_report_filt.txt', out_cpg_report)
out_cpg_report_filt_cov = re.sub(r'.CpG_report.txt$', '.CpG_report_filt.cov', out_cpg_report)

# remove uncovered regions:
cmd = "awk '{ if ($4+$5 > 0) print; }'"
cmd += " " + out_cpg_report
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

bedGraph = out_extractor.replace(".bismark.cov",".bedGraph")
out_bigwig = bedGraph.replace(".bedGraph", ".bw")
cmd = tools.bed2bigWig + " " + bedGraph + " " + resources.chrom_sizes
cmd += " " + out_bigwig

pm.run(cmd, out_bigwig, shell=False)


# Spike-in alignment
################################################################################
# currently using bowtie1 instead of bowtie2
pm.timestamp("### Bismark spike-in alignment: ")
spikein_folder = os.path.join(param.pipeline_outfolder, "bismark_spikein" )
myngstk.make_sure_path_exists(spikein_folder)
spikein_temp = os.path.join(spikein_folder, "bismark_temp" )
myngstk.make_sure_path_exists(spikein_temp)
out_spikein_base = args.sample_name + ".spikein.aln"

#out_spikein = spikein_folder + args.sample_name + "_R1_trimmed.fastq_unmapped_reads_1.fq_bismark_pe.bam"

unmapped_reads_pre = os.path.join(bismark_folder, args.sample_name)

out_spikein = os.path.join(spikein_folder, out_spikein_base + ".bam")
cmd = tools.bismark + " " + resources.bismark_spikein_genome + " "
if args.paired_end:
	cmd += " --1 " + unmapped_reads_pre + "_unmapped_reads_1.fq"
	cmd += " --2 " + unmapped_reads_pre + "_unmapped_reads_2.fq"
else:
	cmd += unmapped_reads_pre + "_unmapped_reads.fq"
cmd += " --bam --unmapped"
cmd += " --path_to_bowtie " + tools.bowtie1
#cmd += " --bowtie2"
cmd += " --temp_dir " + spikein_temp
cmd += " --output_dir " + spikein_folder
if args.paired_end:
	cmd += " --minins 0"
	cmd += " --maxins 5000"
cmd += " --basename="  + out_spikein_base
#cmd += " -p 4"

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
cmd2 = tools.samtools + " sort " + out_spikein_dedup + " " + out_spikein_sorted
cmd3 = tools.samtools + " index " + out_spikein_sorted + ".bam"
cmd4 = "rm " + out_spikein_dedup
pm.run([cmd, cmd2, cmd3, cmd4], out_spikein_sorted +".bam.bai", nofail=True)

# Spike-in methylation calling
################################################################################
pm.timestamp("### Methylation calling (testxmz) Spike-in: ")

cmd1 = tools.python + " -u " + os.path.join(tools.scripts_dir, "testxmz.py")
cmd1 += " " + out_spikein_sorted + ".bam" + " " + "K1_unmethylated"
cmd1 += " >> " + pm.pipeline_stats_file
cmd2 = cmd1.replace("K1_unmethylated", "K3_methylated")
pm.callprint(cmd1, shell=True, nofail=True)
pm.callprint(cmd2, shell=True, nofail=True)


# Final sorting and indexing
################################################################################
# create sorted and indexed BAM files for visualization and analysis
pm.timestamp("### Final sorting and indexing: ")

#out_header = bismark_folder + args.sample_name + ".reheader.bam"
out_final = os.path.join(bismark_folder, args.sample_name + ".final.bam")


# Sort
cmd = tools.java + " -Xmx4g -jar"
# This sort can run out of temp space on big jobs; this puts the temp to a
# local spot.
#cmd += " -Djava.io.tmpdir=`pwd`/tmp"
cmd += " " + tools.picard + " SortSam"
cmd += " I=" + out_sam_filter
cmd += " O=" + out_final
cmd += " SORT_ORDER=coordinate"
cmd += " VALIDATION_STRINGENCY=SILENT"
cmd += " CREATE_INDEX=true"
pm.run(cmd, out_final, lock_name="final_sorting")

# Cleanup
################################################################################
# remove temporary folders
pm.clean_add(bismark_temp)
pm.clean_add(sam_temp)
pm.stop_pipeline()
