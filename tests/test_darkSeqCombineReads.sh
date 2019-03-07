#!/usr/bin/env bash
# script to test darkSeqCombineReads.pl
perl ../src/tools/darkSeqCombineReads.pl ./testR1.fastq ./testR2.fastq ./actual3.fastq 3
perl ../src/tools/darkSeqCombineReads.pl ./testR1.fastq ./testR2.fastq ./actual4.fastq 4
#there should be no difference
echo $(diff ./expected3.fastq ./actual3.fastq)
echo $(diff ./expected4.fastq ./actual4.fastq)
rm ./actual3.fastq
rm ./actual4.fastq
