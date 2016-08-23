#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use FindBin qw($RealBin);
use lib "$RealBin/lib/dapPerlGenomicLib";
use EnsemblRestQuery;

my $build = "GRCh38";
my %opts = (s => 'human', l => 60);
GetOptions
(
    \%opts,
    'h|help',
    'g|grch37',
    'l|line_length=i',
    'f|flanks=i',
    's|species=s',
) or usage("Syntax error!\n");

usage() if $opts{h};
usage("No coordinates provided!") if not @ARGV;

my $restQuery = new EnsemblRestQuery();
if ($opts{g}){
    $restQuery->useGRCh37Server();
}

while (my $region = shift){
    dnaFromRegion($region);
}

sub dnaFromRegion{
    my $region = shift;
    my $strand = "1";
    if ($region =~ /\:(-?1)$/){
        $strand =  $1;
    }else{
        $region .= ':1';
    }
    if ($opts{f}){
        my @spl = split(/[:\-]/, $region); 
        $spl[1] -= $opts{f};
        $spl[2] += $opts{f};
        $region = "$spl[0]:$spl[1]-$spl[2]:$strand";
    }
    my $endpoint =  "/sequence/region/$opts{s}/$region";
    my $seq = $restQuery->queryEndpoint($endpoint);
    if ($seq){
        print ">$seq->{id}\n";
        printDna($seq->{seq});
    }
}
 
sub printDna{
    my $dna = shift;
    if ($opts{f}){
        my $s = lc(substr($dna, 0, $opts{f})); 
        $s .= uc(substr($dna, $opts{f}, length($dna) - 2 * $opts{f})); 
        $s .= lc(substr($dna, length($dna) - $opts{f})); 
        $dna = $s;
    }
    if ($opts{l}){
        for (my $i = 0; $i < length($dna); $i += $opts{l}){
            if (length($dna) > $i + $opts{l}){
                print substr($dna, $i, $opts{l});
            }else{
                print substr($dna, $i,);
            }
            print "\n";
        }
    }else{
        print "$dna\n";
    }
}

sub usage{
    my $msg = shift;
    print "\n$msg\n" if $msg;
    my $exe = fileparse($0);
    print <<EOT

Usage: $0 chr1:100000-100200 [chr2:123456-234567 ...] [options]

Options: 
    -g, --grch37
        Use the GRCh37 REST server instead of the default

    -s, --species STRING
        Retrieve DNA for this species (default = human)
    
    -l,  --line_length INT
        Length of sequence lines in output. Default = 60. 
        Set to 0 to output DNA as a single line.

    -f,  --flanks INT
        Number of bases to add either side of the target region.
        Target bases will be in upper case while flanking bases 
        will be in lower case.

    -h, --help
        Show this message and exit

Description:

Retrives DNA from one or more regions of the genome using Ensembl's REST API. 
DNA will be printed in FASTA format. You may optionally specify the strand by
adding ":1" (+ strand) or ":-1" (- strand) at the end of your region (e.g. 
chr1:100000-100200:1). 

Examples:

    $exe 1:100000-100200 
    (Outputs DNA for coordinates 100,000 to 100,200 of human chromosome one, build GRCh38)

    $exe 1:100000-100200:-1 
    (As above but outputs the - strand)

    $exe 1:100000-100200 -s mouse
    (Retrieve mouse DNA)

    $exe 1:100000-100200 -g 
    (Retrieve DNA from human GRCh37)

    $exe 1:100000-100200 -f 50 
    (Get DNA for the flanking 50 bp as well)

    $exe 1:100000-100200 -l 100
    (Output lines of 100 letters in length)

    $exe 1:100000-100200 -l 0
    (Output a single line for DNA sequence regardless of length)

    $exe chr1:100000-100200 chrX:1234567-2345678
    (more than one region can be specified. Leading 'chr' is optional)

EOT
;
    exit 1 if $msg;
    exit;
}


