# count mismatches between read1 and read2 for dark sequencing
# inputs: read1file, read2file, number to start at in read2 file (1 indexed)
# you always start with the first base on read1
# the comparison is by exact text matching so comparing N and A, C, G, or T would give a mismatch

$file1 = shift;  # read1 from @ARGV
$file2 = shift;
$startind = shift; # which base in read 2 to start on when comparing to read1 (1 indexed)
open(my $fh, "<", $file1);
open(my $fh2, "<", $file2);


# Loop through reads
my $mismatch = 0;
my $seqcount = 0;
while($header1 = <$fh>) {
    $seqcount += 1;
    <$fh2>; # get rid of header2
    my $seq1 = <$fh>;
    my $seq2 = <$fh2>;
    chomp $seq1;
    chomp $seq2;
    # $seq2 =~ s/^\s+//;
    # chomp($seq2);
    $seq2length = length($seq2);
    $comparelength = $seq2length - ($startind - 1);
    # substr EXPR,OFFSET,LENGTH
    $firstfewseq1 = substr($seq1, 0, $comparelength);
    # I made startind 1 indexed so subtract 1
    $matchingseq2 = substr($seq2, $startind - 1, $comparelength);
    
    #if ($seqcount%500000 == 0) {
    #    print("firstfewseq1 = [$firstfewseq1]");
    #    print("matchingseq2 = [$matchingseq2]");
    #}
    
    
	<$fh>;<$fh>;<$fh2>;<$fh2>;
    
    if ($firstfewseq1 ne $matchingseq2) {
        $mismatch += 1;

    } 
    #if ($seqcount%500000 == 0) {
    #    print("seq1 = $firstfewseq1, seq2 = $matchingseq2");
    #    print("seq2length = $seq2length\n");
    #    print("length firstfewseq1 = " . length(firstfewseq1). "\n");
    #}   
}
# print $seq2length . "\n";
# print("comparelength = $comparelength\n");
# print("seq2length = $seq2length\n");

# the output. Proportion of reads where the compared portions do not match
print STDOUT ($mismatch / $seqcount);
