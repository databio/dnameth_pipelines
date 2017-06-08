#!/usr/bin/env python
"""
This is just a simple parser to grab, from a tsv file, a certain column.
You can also specify to grab this column only from a certain row, which you can
key with any number of column key specificities.
I wrote it to pull output from biseq methcalling results, but it could
be used for anything... just give a tsv file and a key column to pull, and any
rows to restrict to if desired.
"""

__author__ = "Nathan C. Sheffield"

import sys, os, csv
import subprocess
from argparse import ArgumentParser


parser = ArgumentParser(description='TSV Parser')

parser.add_argument('-i', '--input-file', dest='input_file',
	help="Supply tsv file.",
	required=True)

parser.add_argument('-c', '--col', dest='column',
	default="uniqueSeqMotifCount", nargs="+",
	help="Choose column(s) to extract",
	required=True)

parser.add_argument('-r', '--rowmatch', dest='rowmatch',
	default="", nargs="+", type=str,
	help="Choose row(s) to extract; format: key=value. Multiple, space-separated \
	key=value pairs are accepted.",
	required=False)

parser.add_argument('-k', '--keys', dest='include_keys',
 	default=False, action="store_true",
	help="Print out keys? By default, No.",
	required=False)

parser.add_argument('-H', '--no-header', dest='header',
 	default=True, action="store_false",
	help="Turns off header row parsing (column names).",
	required=False)

args = parser.parse_args()

input_open = open(args.input_file, 'rb')
if args.header:
	f = csv.DictReader(input_open, delimiter="\t")
else:
	f = csv.reader(input_open, delimiter="\t")

fail = False

for row in f:  # iterates the rows of the file in order
	#print(row)
	if args.rowmatch:
		fail = False
		for rowmatch_string in args.rowmatch:
			rowmatch_list = rowmatch_string.split("=")
			# print(rs)
			# rowmatch_list[0]  # key
			# rowmatch_list[1]  # value
			if not args.header:
				rowmatch_list[0] = int(rowmatch_list[0])
			if not row[rowmatch_list[0]] == rowmatch_list[1]:
				fail = True
	if fail:
		# The first column value did not match any of the keyed rows;
		# So we don't print out this row.
		continue

	for key in args.column:
		if args.include_keys:
			print("\t".join([key, row[key]]))
		else:
			if args.header:
				print(row[key])
			else:
				# maybe the key was an integer?
				print(row[int(key)])
