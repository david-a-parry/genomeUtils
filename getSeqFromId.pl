#!/usr/bin/env perl
use warnings;
use strict;
use Getopt::Long;
use Data::Dumper;
use Pod::Usage;
use Getopt::Long;
use FindBin qw($RealBin);
use lib "$RealBin/lib/dapPerlGenomicLib";
use EnsemblRestQuery;
use IdParser;

my $restQuery = new EnsemblRestQuery();
my @ids = ();
my %opts = 
(
    i => \@ids,
    s => "human",
    n => 60,
);
GetOptions
(
    \%opts,
    'l|list=s',
    'n|line_length=i',
    'i|ids=s{,}',
    's|species=s',
    'r|grch37',
    't|type=s',
    'h|?|help',
    'm|manual',
) or pod2usage( -message => "Syntax error.", -exitval => 2 );
pod2usage( -verbose => 2 ) if ($opts{m});
pod2usage( -verbose => 1 ) if ($opts{h});

pod2usage
(
     -message => "ERROR: Either -l/--list or -i/--ids argument is required.\n", 
     -exitval => 2 
) if ( not $opts{l} and not @ids);

checkType();

if ($opts{r}){
    $restQuery->useGRCh37Server();
}

if ($opts{l}){
    open (my $GENES, $opts{l}) or die "Could not open -list '$opts{l}' for reading: $!\n";
    while (my $line = <$GENES>){
        my @s = split(/\s+/, $line); 
        push @ids, $s[0];
    }
}
if (not @ids){
    die "No gene IDs identified!\n";
}

my @cross_ref = ();
my %dbs       = ();
foreach my $g (@ids){
    getSeq($g);
}


#########################################################
sub checkType{
    return if not $opts{t};
    my @valid = qw /
        genomic
        cdna
        cds
        protein
    /;
    if (not grep { $_ eq lc($opts{t}) } @valid){
        die "Unrecognised --type '$opts{t}' - valid types are:\n" .
            join("\n", @valid) . "\n";
    }
    $opts{t} = lc($opts{t});# REST API is case-sensitive
}


#########################################################
sub getSeq{
    my $id = shift;
    my $id_parser = new IdParser();
    $id_parser->parseId($id);
    my $ensid;
    if ($id_parser->get_isEnsemblId()){
         $ensid = $id;
    }else{
        if ($id_parser->get_isTranscript() or $id_parser->get_isProtein() ) {
            my $o = $id_parser->get_isTranscript ? 'transcript' : 'translation';
            my $r = $restQuery->getViaXreg($id, $opts{s}, $o);
            if ($r and ref $r eq 'ARRAY'){
                if (@$r > 1){
                    warn "WARNING: Multiple objects identified by ".
                      "cross-reference search for $id - picking the first.\n";
                }
                foreach my $tr (@$r){
                    if (exists $tr->{id}){
                        $ensid = $tr->{id};
                    }
                }
            }
        }else{
            my $gene = $restQuery->getGeneDetails($id, $opts{s});
            if ($gene){
                $ensid = $gene->{id};
            }
        }
    }
    if (not $ensid){
        warn "WARNING: No object identified for ID \"$id\"\n";
        return;
    }
    my $endpoint =  "/sequence/id/$ensid?multiple_sequences=1";
    $endpoint .= ";type=$opts{t}" if $opts{t};
    my $seq = $restQuery->queryEndpoint($endpoint);
    if ($seq){
        foreach my $s (@$seq){
            print ">$s->{id}\n";
            printDna($s->{seq});
        }
    }
}

#########################################################
sub printDna{
    my $dna = shift;
    if ($opts{n}){
        for (my $i = 0; $i < length($dna); $i += $opts{n}){
            if (length($dna) > $i + $opts{n}){
                print substr($dna, $i, $opts{n});
            }else{
                print substr($dna, $i,);
            }
            print "\n";
        }
    }else{
        print "$dna\n";
    }
    print "\n";
}

#########################################################

=head1 NAME

getSeqFromId.pl - retrieve IDs for given genes

=head1 SYNOPSIS
 
        getSeqFromId.pl -g ABCD1 
        getSeqFromId.pl -l gene_list.txt 
        getSeqFromId.pl -h (display help message)
        getSeqFromId.pl -m (display manual page)

=cut

=head1 ARGUMENTS

=over 8

=item B<-i    --ids>

One or more gene symbols or IDs to search and retrieve IDs for.

=item B<-l    --list>

A list of gene symbols or IDs to use instead of or in conjunction with -g/--gene_ids. 

=item B<-s    --species>

Species genome to use. Defaults = human. 

When using gene symbols or other identifiers such as RefSeq or Uniprot, this program will attempt to identify the corresponding Ensembl gene but will only search this species. Therefore, if you are attempting to find the coordinates of a mouse RefSeq transcript you must use this option to specify 'mouse' as the species or the corresponding Ensembl gene will not be found. 

This option is ignored if Ensembl identifiers are used because they are already species specific.

=item B<-r    --grch37>

Use this option if you want to use the GRCh37 Ensembl REST server for gene queries rather than the current REST server.

=item B<-t    --type>

Sequence type to output. By default it will be decided by the type of identifier given. Valid options are: 

        genomic
        cdna
        cds
        protein

=item B<-n  --line_length>

Length of sequence lines in output. Default = 60. Set to 0 to output DNA as a single line.

=item B<-h    --help>

Show this script's help information.

=item B<-m    --manual>

Show this script's manual page.

=back

=cut

=head1 DESCRIPTION

This program uses Ensembl's REST API to output sequences for given gene/transcript/protein IDs. 

=cut

=head1 AUTHOR

David A. Parry


=head1 COPYRIGHT AND LICENSE

Copyright 2016 David A. Parry

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with this program. If not, see <http://www.gnu.org/licenses/>.

=cut




