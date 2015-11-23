#!/usr/bin/env python
"""
WGBS pipeline
documentation.
"""
from argparse import ArgumentParser
import os, re
import os.path
import sys
from subprocess import call
import subprocess
from datetime import datetime
import yaml
# Argument Parsing
################################################################################
parser = ArgumentParser(description='Pypiper arguments.')
# Set up a pointer to the pypiper code you want to use:
# Just grab the single pypiper arg, and add it to the path here; leaving all other args for later parsing.
parser.add_argument('-P', '--pypiper', dest='pypiper_dir',
					default=os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.realpath(__file__)))), "pypiper"),
					type=str)

args = parser.parse_known_args()[0]
os.sys.path.insert(0, args.pypiper_dir)
from pypiper import Pypiper
from pypiper import ngstk

# Add Pypiper arguments to the arguments list
parser = Pypiper.add_pypiper_args(parser)
parser.add_argument('-c','--config', dest='config', type=str, required=True, help='Path to YAML configuration file')
parser.add_argument('-i', '--unmapped-bam', nargs="+", dest='unmapped_bam', required=True, help="Input unmapped bam file(s))")
parser.add_argument('-s', '--sample-name', default="default", dest='sample_name', type=str, help='Sample Name') # default means deduction from filename, except .bam extension
parser.add_argument('-r', '--project-root', default="", dest='project_root', type=str, help='Directory in which the project will reside.')
parser.add_argument('-g', '--genome', dest='genome_assembly', type=str, required=True, help='Genome Assembly')
parser.add_argument('-C', '--no-checks', dest='no_check', action='store_true', default=False, help='Skip sanity checks')
parser.add_argument('-p', '--paired-end', dest='paired_end', action='store_true', help='Paired End Mode')
parser.add_argument('-q', '--single-end', dest='paired_end', action='store_false', help='Single End Mode')
args = parser.parse_args()

# Merging
################################################################################
# If 2 unmapped bam files are given, then these are to be merged.
# Must be done here to initialize the sample name correctly
merge = False
if (len(args.unmapped_bam) > 1):
	merge = True
	if (args.sample_name == "default"):
		args.sample_name = "merged";
else:
	if (args.sample_name == "default"):
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
paths.ref_genome = os.path.join(paths.resources_dir, "genomes")
paths.bismark_indexed_genome = os.path.join(paths.resources_dir, "genomes", args.genome_assembly, "indexed_bismark_bt2")
paths.bismark_spikein_genome = os.path.join(paths.resources_dir, "genomes", "meth_spikein_k1_k3", "indexed_bismark_bt1")
paths.chrom_sizes = os.path.join(paths.resources_dir, "genomes", args.genome_assembly, args.genome_assembly + ".chromSizes")

# Tools
paths.python = config["tools"]["python"]
paths.java = config["tools"]["java"]
paths.Rscript = config["tools"]["Rscript"]
paths.picard_jar = config["tools"]["picard"]
paths.trimmomatic_jar = config["tools"]["trimmomatic"]
paths.trimmomatic_epignome_jar = config["tools"]["trimmomatic_epignome"]
paths.bowtie1 = config["tools"]["bowtie1"]
paths.bowtie2 = config["tools"]["bowtie2"]
paths.bed2bigBed = config["tools"]["bed2bigBed"]
paths.bed2bigWig = config["tools"]["bed2bigWig"]
paths.bismark = config["tools"]["bismark"]
paths.deduplicate_bismark = config["tools"]["deduplicate_bismark"]
paths.bismark_methylation_extractor = config["tools"]["bismark_methylation_extractor"]
paths.samtools = config["tools"]["samtools"]

#Output
paths.pipeline_outfolder_abs = os.path.abspath(os.path.join(args.project_root, args.sample_name))

# Create a Pypiper object, and start the pipeline
mypiper = Pypiper(name="WGBS", outfolder=paths.pipeline_outfolder_abs, args=args)

# Create a ngstk object
myngstk = ngstk.NGSTk(paths)

unmapped_bam_folder = os.path.join(paths.pipeline_outfolder_abs, "unmapped_bam")
myngstk.make_sure_path_exists(unmapped_bam_folder)

if merge:
	input_bams = args.unmapped_bam
	merge_folder = unmapped_bam_folder
	sample_merged_bam = args.sample_name + ".merged.bam"
	output_merge = os.path.join(merge_folder, sample_merged_bam)
	cmd =  myngstk.merge_bams(input_bams, output_merge)

	mypiper.call_lock(cmd, output_merge)
	args.unmapped_bam = os.path.join(unmapped_bam_folder, sample_merged_bam)  #update unmapped bam reference
	local_unmapped_bam_abs = os.path.abspath(args.unmapped_bam)
else:
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

	local_unmapped_bam_abs = os.path.abspath(os.path.join(unmapped_bam_folder, args.sample_name + input_ext))
	mypiper.callprint("ln -sf " + args.unmapped_bam + " " + local_unmapped_bam_abs, shell=True)


print("Input Unmapped Bam: " + args.unmapped_bam)
#check for file exists:
if not os.path.isfile(local_unmapped_bam_abs):
	raise Exception(local_unmapped_bam_abs + " is not a file")


# Fastq conversion
################################################################################
mypiper.timestamp("### Fastq conversion: ")

fastq_folder = os.path.join(paths.pipeline_outfolder_abs, "fastq")
out_fastq_pre = os.path.join(fastq_folder, args.sample_name)

cmd = myngstk.bam_to_fastq(local_unmapped_bam_abs, out_fastq_pre, args.paired_end)
mypiper.call_lock(cmd, out_fastq_pre + "_R1.fastq")
mypiper.clean_add(os.path.join(out_fastq_pre, "*.fastq"), conditional=True)

# Sanity checks:
if not args.no_check:
	raw_reads = myngstk.count_reads(local_unmapped_bam_abs, args.paired_end)
	mypiper.report_result("Raw_reads", str(raw_reads))
	fastq_reads = myngstk.count_reads(out_fastq_pre + "_R1.fastq", args.paired_end)
	mypiper.report_result("Fastq_reads", fastq_reads)
	if (fastq_reads != int(raw_reads)):
		raise Exception("Fastq conversion error? Size doesn't match unaligned bam")

# Adapter trimming
################################################################################
mypiper.timestamp("### Adapter trimming: ")

memory = str(config["parameters"]["trimmomatic"]["memory"])
threads = str(config["parameters"]["trimmomatic"]["threads"])
illuminaclip = str(config["parameters"]["trimmomatic"]["illuminaclip"])
encoding = "phred33"

trimmed_fastq = out_fastq_pre + "_R1_trimmed.fq"
trimmed_fastq_R2 = out_fastq_pre + "_R2_trimmed.fq"

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
cmd += " HEADCROP:6 ILLUMINACLIP:" + paths.adapter_file + illuminaclip

mypiper.clean_add(os.path.join(fastq_folder, "*.fastq"), conditional=True)
mypiper.clean_add(os.path.join(fastq_folder, "*.fq"), conditional=True)
mypiper.clean_add(os.path.join(fastq_folder, "*.log"), conditional=True)
mypiper.clean_add(fastq_folder, conditional=True)

mypiper.call_lock(cmd, trimmed_fastq)   # TODO AS: maybe another lock_name?
trimmed_reads_count = myngstk.count_reads(trimmed_fastq, args.paired_end)
if not args.no_check:
	mypiper.report_result("Trimmed_reads", trimmed_reads_count)

# WGBS alignment with bismark.
################################################################################
mypiper.timestamp("### Bismark alignment: ")
bismark_folder = os.path.join(paths.pipeline_outfolder_abs, "bismark_" + args.genome_assembly )
myngstk.make_sure_path_exists(bismark_folder)
bismark_temp = os.path.join(bismark_folder, "bismark_temp" )
myngstk.make_sure_path_exists(bismark_temp)

if args.paired_end:
	out_bismark = os.path.join(bismark_folder, args.sample_name + "_pe.bam")
else:
	out_bismark = os.path.join(bismark_folder, args.sample_name + ".bam")

cmd = paths.bismark + " " + paths.bismark_indexed_genome + " "
if args.paired_end:
	cmd += " --1 " + out_fastq_pre + "_R1_trimmed.fq"
	cmd += " --2 " + out_fastq_pre + "_R2_trimmed.fq"
else:
	cmd += out_fastq_pre + "_R1_trimmed.fq"
cmd += " --bam --unmapped"
cmd += " --path_to_bowtie " + paths.bowtie2
cmd += " --bowtie2"
cmd += " --temp_dir " + bismark_temp
cmd += " --output_dir " + bismark_folder
if args.paired_end:
	cmd += " --minins 0"
	cmd += " --maxins 5000"
cmd += " -p 8 " # Number of processors
cmd += " --basename=" +args.sample_name

mypiper.call_lock(cmd, out_bismark)

if not args.no_check:
	x = myngstk.count_reads(out_bismark, args.paired_end)
	mypiper.report_result("Aligned_reads", x)


mypiper.timestamp("### PCR duplicate removal: ")
# Bismark's deduplication forces output naming, how annoying.
out_dedup = os.path.join(bismark_folder, args.sample_name + "_pe.deduplicated.bam")
out_dedup = re.sub(r'.bam$', '.deduplicated.bam', out_bismark)

cmd = paths.deduplicate_bismark
if args.paired_end:
	cmd += " --paired "
else:
	cmd += " --single "
cmd += out_bismark
cmd += " --bam"

mypiper.call_lock(cmd, out_dedup)

if not args.no_check:
	deduplicated_reads = myngstk.count_reads(out_dedup, args.paired_end)
	mypiper.report_result("Deduplicated_reads", deduplicated_reads)


mypiper.timestamp("### Aligned read filtering: ")

# convert bam file into sam file and sort again to
# compensate for a sorting issue of "deduplicate_bismark"
sam_temp = os.path.join(bismark_folder, "sam_temp")
myngstk.make_sure_path_exists(sam_temp)
out_sam = os.path.join(bismark_folder, args.sample_name + ".aln.deduplicated.sam")
cmd = paths.samtools + " sort -n -o " + out_dedup + " " + out_dedup.replace(".bam", "_sorted") + " | " + paths.samtools + " view -h - >" + out_sam

mypiper.call_lock(cmd, out_sam, shell=True)

if not args.no_check:
	#sorted file same size as presorted?
	sorted_reads = myngstk.count_reads(out_sam, args.paired_end)
	if sorted_reads != deduplicated_reads:
		raise Exception("Sorted size doesn't match deduplicated size.")

out_sam_filter = os.path.join(bismark_folder, args.sample_name + ".aln.dedup.filt.sam")

headerLines = subprocess.check_output(paths.samtools + " view -SH " + out_sam + "| wc -l", shell=True).strip()
cmd = paths.python + " " + os.path.join(paths.scripts_dir, "bisulfiteReadFiltering_forRNA.py")
cmd += " --infile=" + out_sam
cmd += " --outfile=" + out_sam_filter
cmd += " --skipHeaderLines=" + headerLines
cmd += " --genome=" + args.genome_assembly
cmd += " --genomeDir=" + paths.ref_genome
cmd += " --minNonCpgSites=3"
cmd += " --minConversionRate=0.9"
if args.paired_end:
	cmd = cmd + " --pairedEnd"

mypiper.call_lock(cmd, out_sam_filter)

if not args.no_check:
	x = myngstk.count_reads(out_sam_filter, args.paired_end)
	mypiper.report_result("Filtered_reads", x)


# Clean up all intermediates
mypiper.clean_add(out_bismark) # initial mapped bam file
mypiper.clean_add(os.path.join(bismark_folder, "*.fastq"))
mypiper.clean_add(os.path.join(bismark_folder, "*.fq"))
mypiper.clean_add(out_dedup) # deduplicated bam file
mypiper.clean_add(out_sam) # dedup conversion to sam
mypiper.clean_add(out_sam_filter) # after filtering


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

mypiper.timestamp("### Methylation calling (bismark extractor): ")

extract_dir = os.path.join(bismark_folder, "extractor")
myngstk.make_sure_path_exists(extract_dir)
out_extractor = os.path.join(extract_dir, re.sub(r'.sam$', '.bismark.cov', os.path.basename(out_sam_filter)))
out_cpg_report = re.sub(r'.bismark.cov$', '.CpG_report.txt', out_extractor)

cmd = paths.bismark_methylation_extractor
if args.paired_end:
	cmd += " --paired-end --no_overlap"
else:
	cmd += " --single-end"
cmd += " --report"
cmd += " --bedGraph"
cmd += " --merge_non_CpG"
cmd += " --cytosine_report"
cmd += " --genome_folder " + paths.bismark_indexed_genome
cmd += " --gzip"
cmd += " --output " + extract_dir
cmd += " " + out_sam_filter

mypiper.call_lock(cmd,  out_cpg_report)

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
mypiper.call_lock(cmd,  out_cpg_report_filt, shell=True)

# convert the bismark report to the simpler coverage format and adjust the coordinates
# of CpG's on the reverse strand while doing so (by substracting 1 from the start):
cmd = paths.Rscript + " " + os.path.join(os.path.dirname(os.path.dirname(os.path.realpath(__file__))), "bin", "convertBismarkReport.R") # disable coverage filter, because we have already used `awk` to achieve this result
cmd += " --formats=cov,min,gibberish"
cmd += " --noCovFilter"
if keep_non_standard_chromosomes:
	cmd += " --noChromFilter"
if not adjust_minus_strand:
	cmd += " --noAdjustMinusStrand"
cmd += " -i " + out_cpg_report_filt
mypiper.call_lock(cmd,  out_cpg_report_filt_cov)


# tidy up:
if not keep_bismark_report:
	mypiper.clean_add(out_cpg_report_filt)


# Make bigwig
################################################################################
mypiper.timestamp("### Make bigwig: ")

bedGraph = out_extractor.replace(".bismark.cov",".bedGraph")
out_bigwig = bedGraph.replace(".bedGraph", ".bw")
cmd = paths.bed2bigWig + " " + bedGraph + " " + paths.chrom_sizes
cmd += " " + out_bigwig

mypiper.call_lock(cmd, out_bigwig, shell=False)


# Spike-in alignment
################################################################################
# currently using bowtie1 instead of bowtie2
mypiper.timestamp("### Bismark spike-in alignment: ")
spikein_folder = os.path.join(paths.pipeline_outfolder_abs, "bismark_spikein" )
myngstk.make_sure_path_exists(spikein_folder)
spikein_temp = os.path.join(spikein_folder, "bismark_temp" )
myngstk.make_sure_path_exists(spikein_temp)
out_spikein_base = args.sample_name + ".spikein.aln"

#out_spikein = spikein_folder + args.sample_name + "_R1_trimmed.fastq_unmapped_reads_1.fq_bismark_pe.bam"

unmapped_reads_pre = os.path.join(bismark_folder, args.sample_name)

out_spikein = os.path.join(spikein_folder, out_spikein_base + ".bam")
cmd = paths.bismark + " " + paths.bismark_spikein_genome + " "
if args.paired_end:
	cmd += " --1 " + unmapped_reads_pre + "_unmapped_reads_1.fq"
	cmd += " --2 " + unmapped_reads_pre + "_unmapped_reads_2.fq"
else:
	cmd += unmapped_reads_pre + "_unmapped_reads.fq"
cmd += " --bam --unmapped"
cmd += " --path_to_bowtie " + paths.bowtie1
#cmd += " --bowtie2"
cmd += " --temp_dir " + spikein_temp
cmd += " --output_dir " + spikein_folder
if args.paired_end:
	cmd += " --minins 0"
	cmd += " --maxins 5000"
cmd += " --basename="  + out_spikein_base
#cmd += " -p 4"

mypiper.call_lock(cmd, out_spikein, nofail=True)
# Clean up the unmapped file which is copied from the parent
# bismark folder to here:
mypiper.clean_add(os.path.join(spikein_folder, "*.fq"), conditional=False)
mypiper.clean_add(spikein_temp)


mypiper.timestamp("### PCR duplicate removal (Spike-in): ")
# Bismark's deduplication forces output naming, how annoying.
#out_spikein_dedup = spikein_folder + args.sample_name + ".spikein.aln.deduplicated.bam"
out_spikein_dedup = re.sub(r'.bam$', '.deduplicated.bam', out_spikein)


cmd = paths.deduplicate_bismark
if args.paired_end:
	cmd += " --paired "
else:
	cmd += " --single "
cmd += out_spikein
cmd += " --bam"

out_spikein_sorted = out_spikein_dedup.replace('.deduplicated.bam', '.deduplicated.sorted')
cmd2 = paths.samtools + " sort " + out_spikein_dedup + " " + out_spikein_sorted
cmd3 = paths.samtools + " index " + out_spikein_sorted + ".bam"
cmd4 = "rm " + out_spikein_dedup
mypiper.call_lock([cmd, cmd2, cmd3, cmd4], out_spikein_sorted +".bam.bai", nofail=True)

# Spike-in methylation calling
################################################################################
mypiper.timestamp("### Methylation calling (testxmz) Spike-in: ")

cmd1 = paths.python + " -u " + os.path.join(paths.scripts_dir, "testxmz.py")
cmd1 += " " + out_spikein_sorted + ".bam" + " " + "K1_unmethylated"
cmd1 += " >> " + mypiper.pipeline_stats_file
cmd2 = cmd1.replace("K1_unmethylated", "K3_methylated")
mypiper.callprint(cmd1, shell=True, nofail=True)
mypiper.callprint(cmd2, shell=True, nofail=True)


# Final sorting and indexing
################################################################################
# create sorted and indexed BAM files for visualization and analysis
mypiper.timestamp("### Final sorting and indexing: ")

#out_header = bismark_folder + args.sample_name + ".reheader.bam"
out_final = os.path.join(bismark_folder, args.sample_name + ".final.bam")


# Sort
cmd = paths.java + " -Xmx4g -jar"
# This sort can run out of temp space on big jobs; this puts the temp to a
# local spot.
#cmd += " -Djava.io.tmpdir=`pwd`/tmp"
cmd += " " + paths.picard_jar + " SortSam"
cmd += " I=" + out_sam_filter
cmd += " O=" + out_final
cmd += " SORT_ORDER=coordinate"
cmd += " VALIDATION_STRINGENCY=SILENT"
cmd += " CREATE_INDEX=true"
mypiper.call_lock(cmd, out_final, lock_name="final_sorting")

# Cleanup
################################################################################
# remove temporary folders
mypiper.clean_add(bismark_temp)
mypiper.clean_add(sam_temp)
mypiper.stop_pipeline()
