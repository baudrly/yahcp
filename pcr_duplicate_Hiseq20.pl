#!/usr/bin/env perl

# 19/03/2015
# Take two fastq files together and revove the PCR duplicates i.e whose ends have exactly the same sequences.
# and remove tags of 10 bp 

if ($#ARGV != 3) {
    print "usage: ./pcr_duplicate.pl leftsource rightsource leftout rightout\n";
    exit;
}

my ($leftsource, $rightsource, $leftout, $rightout) = @ARGV;

# open input files
if ($leftsource =~ /\.gz$/) {
    open(LEFT, "pigz -d -c $leftsource |");
} else {
    open(LEFT, "< $leftsource");
}
if ($rightsource =~ /\.gz$/) {
    open(RIGHT, "pigz -d -c $rightsource |");
} else {
    open(RIGHT, "< $rightsource");
}

# open output files
if  ($leftout =~ /\.gz$/) {
    open(LOUT, "| pigz > $leftout");
} else {
    open(LOUT, "> $leftout");
}
if  ($rightout =~ /\.gz$/) {
    open(ROUT, "| pigz > $rightout");
} else {
    open(ROUT, "> $rightout");
}

my $f = 0;
my %group;
while (1) {
    my $line_before1 = <LEFT> or last; 
    my $line_before2 = <RIGHT> or die "Unexpected line ending";

    my $line1 = <LEFT> or die "Unexpected line ending";
    my $line2 = <RIGHT> or die "Unexpected line ending";
    my $word1 = substr $line1, 0, 20;
    my $word2 = substr $line2, 0, 20;
	
    if (exists($group{$word1.$word2}) || exists($group{$word2.$word1})) {
        <LEFT> or die "Unexpected line ending";
        <RIGHT> or die "Unexpected line ending";
        <LEFT> or die "Unexpected line ending";
        <RIGHT> or die "Unexpected line ending";
    } else {
	$group{$word1.$word2} = 1;
      
	print LOUT $line_before1, substr($line1, 10);
	print ROUT $line_before2, substr($line2, 10);
	$line1 = <LEFT> or die "Unexpected line ending";
	$line2 = <RIGHT> or die "Unexpected line ending";
	print LOUT $line1;
	print ROUT $line2;
        $line1 = <LEFT> or die "Unexpected line ending";
        $line2 = <RIGHT> or die "Unexpected line ending";
	print LOUT substr($line1, 10);
	print ROUT substr($line2, 10);
    }
    $f++;
}
print STDERR "There are $f reads parsed.\n";
close(LEFT); close(RIGHT); close(LOUT); close(ROUT);

0;

