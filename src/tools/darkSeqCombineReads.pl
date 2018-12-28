#! /usr/bin/env perl 
# count mismatches between read1 and read2 for dark sequencing
# inputs: read1file, read2file, outputfile, number of bases from the beginning
# of read2 to be added to the beginning of read1
# usage: perl darkSeqCombineReads.pl read1.fq read2.fq combined.fq 3

# input arguments
$file1 = shift;  # read1 from @ARGV
$file2 = shift;
$fileout = shift;
# $startind = shift; # which base in read 2 to start on when comparing to read1
$addlength = shift; # add the first x bases from read 2 to read1

open(my $fh, "<", $file1);
open(my $fh2, "<", $file2);
open(my $fhout, ">", $fileout);

# Loop through reads
while($seq1header = <$fh>) {
    
    $seq1 = <$fh>;
    $seq2header = <$fh2>;
    $seq2 = <$fh2>;
    $seq2add = substr($seq2, 0, $addlength); #bases that will be added to seq1
    <$fh2>; # skip third line
    $linethree = <$fh>;
    
    # add this read to updated fastq file only if it passed filter 
    # (filtered=N means it was not filtered so we should keep it)
    if (($seq1header =~ /:N:/) and ($seq2header =~ /:N:/)) {
        print $fhout $seq1header;
        print $fhout ($seq2add . $seq1);
        print $fhout $linethree; # line three of seq1
        print $fhout (substr(<$fh2>, 0, $addlength) . <$fh>); # concatenate quality scores
    } else {
        # skip quality scores to move to next read
        <$fh>;
        <$fh2>;
    }
}







# output file

# assign file for read1 and read2


# count number of total reads
# count number of reads from R2 that do not 100% completely match R1 (bases 4-9)
# write this to a file





######## Other script, run once for each sample
# read file that has counts for each substr and compile that into totals for a single sample

#opendir my $dir, "/some/path" or die "Cannot open directory: $!";
#my @files = readdir $dir;
#closedir $dir;