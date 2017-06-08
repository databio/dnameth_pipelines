#!/usr/bin/python

import pysam
from argparse import ArgumentParser

argument_parser = ArgumentParser(
    description='Script to change paired read SAM flag into unpaired for Bismark processing')

argument_parser.add_argument(
    '-i', '--filein',
    help='input file name',
    required=True,
    type=str)

argument_parser.add_argument(
    '-o','--fileout',
    help='output file name',
    required=False,
    type=str)

name_space = argument_parser.parse_args()

samfile_in_name = ''

if name_space.filein:
    samfile_in_name = name_space.filein

if name_space.fileout:
    samfile_out_name = name_space.fileout
else:
    samfile_out_name = samfile_in_name[:-4] + '_REFLAGGED.bam'

if not (samfile_in_name.endswith('.bam') and samfile_out_name.endswith('.bam')):
    raise Exception("Input and output files must be in BAM format")

samfile_in = pysam.AlignmentFile(samfile_in_name, "rb", check_sq=False)
samfile_out = pysam.AlignmentFile(samfile_out_name, "wb", template=samfile_in)

for read in samfile_in.fetch(until_eof=True):

    new_segment = read
    new_segment.flag = read.flag & ~235
    new_segment.mrnm = -1
    new_segment.next_reference_start = -1
    samfile_out.write(read)

samfile_in.close()
samfile_out.close()
