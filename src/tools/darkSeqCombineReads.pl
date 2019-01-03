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
    chomp($seq2 = <$fh2>);
    $seq2length = length $seq2;
    $seq2add = substr($seq2, 0, $addlength); #bases that will be added to seq1
    <$fh2>; # skip third line
    $linethree = <$fh>;
    
    # add this read to updated fastq file only if it passed filter 
    # (filtered=N means it was not filtered so we should keep it)
    if (($seq1header =~ /:N:/) and ($seq2header =~ /:N:/)) {
                
        $newseq = ($seq2add . $seq1);
        $qual1 = <$fh>;
        $qual2 = <$fh2>;
        $newqual = (substr($qual2, 0, $addlength) . $qual1);
        
        # if you are not concerned with the overlapping sequence between R1 and R2
        # you could delete this for loop to save computational time
        ## for the bases that overlap, keep whichever has the best quality score
        ## loop through overlapping bases, add whichever has the best score to $seq2add
        for (my $i=$addlength; i < $seq2length; $i++) {
        
            # only compare quality scores if bases don't match
            # perl is zero indexed so $addength + 1 (1 index) base is the first to be compared
            if (substr($seq2, i, 1) ne substr($seq1, i-$addlength, 1)) {
                $thisbase1 = substr($seq1, i-$addlength, 1);
                $thisbase2 = substr($seq2, i, 1);
                
                # compare quality scores, if seq2 is better then change
                # otherwise keep current base and quality score
                # higher number/ASCII character is better phred score
                if (ord(substr($qual2, i, 1)) > ord(substr($qual1, i-$addlength, 1))) {
                    substr($newseq, i, 1) = thisbase2;
                    substr($newqual, i, 1) = substr($qual2, i, 1);
                }
            }
            # else, make no changes to existing base
        
        }
        
        print $fhout $seq1header;
        print $fhout $newseq;
        print $fhout $linethree; # line three of seq1
        print $fhout $newqual; # concatenate quality scores
        
        
        
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