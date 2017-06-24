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

my $restQuery = new EnsemblRestQuery();
my @var_ids = ();
my %opts = 
(
    v => \@var_ids,
    s => "human",
);
GetOptions
(
    \%opts,
    'l|list=s',
    'v|var_ids=s{,}',
    's|species=s',
    'r|grch37',
    'h|?|help',
    'm|manual',
) or pod2usage( -message => "Syntax error.", -exitval => 2 );
pod2usage( -verbose => 2 ) if ($opts{m});
pod2usage( -verbose => 1 ) if ($opts{h});

pod2usage
(
     -message => "ERROR: Either -l/--list or -v/--var_ids argument is required.\n", 
     -exitval => 2 
) if ( not $opts{l} and not @var_ids);

if ($opts{r}){
    $restQuery->useGRCh37Server();
}

if ($opts{l}){
    open (my $VARS, $opts{l}) or die "Could not open --list '$opts{l}' for reading: $!\n";
    while (my $line = <$VARS>){
        my @s = split(/\s+/, $line); 
        push @var_ids, $s[0] if $s[0];
    }
}
if (not @var_ids){
    die "No variant IDs identified!\n";
}

print join("\t", qw/
    ID
    Chr
    Start
    End
    Assembly
    Alleles
    MAF
    Class
    Minor
    Ancestral
/). "\n";
foreach my $id (@var_ids){
    my $var_hash = $restQuery->queryEndpoint("/variation/$opts{s}/$id");
    next if not $var_hash;
    foreach my $m (@{$var_hash->{mappings}}){
        print join( "\t", 
                    $id, 
                    $m->{seq_region_name},
                    $m->{start},
                    $m->{end},
                    $m->{assembly_name},
                    $m->{allele_string},
                    $var_hash->{MAF} || "NA",
                    $var_hash->{var_class} || "NA",
                    $var_hash->{minor_allele} || "NA",
                    $var_hash->{ancestral_allele} || "NA",
                    
        ) . "\n";
    }
}

=head1 NAME

varInfo.pl - retrieve genomic coordinates for given genes

=head1 SYNOPSIS
 
        varInfo.pl -v rs1111 
        varInfo.pl -l var_list.txt 
        varInfo.pl -h (display help message)
        varInfo.pl -m (display manual page)

=cut

=head1 ARGUMENTS

=over 8

=item B<-v    --var_ids>

One or more variant symbols or IDs to search for and retrieve.

=item B<-l    --list>

A list of variant IDs to use instead of or in conjunction with -v/--var_ids. 

=item B<-s    --species>

Species genome to use. Defaults = human. 

=item B<-r    --grch37>

Use this option if you want to use the GRCh37 Ensembl REST server for gene queries rather than the current REST server.

=item B<-h    --help>

Show this script's help information.

=item B<-m    --manual>

Show this script's manual page.

=back

=cut

=head1 EXAMPLES

 varInfo.pl -v rs1111
 #get info for human dbSNP variant rs1111

 varInfo.pl -v rs48520904  -s mouse
 #get info for mouse dbSNP variant rs48520904

 varInfo.pl -l var_list.txt
 #get information for a list of human variant IDs

 
=cut

=head1 DESCRIPTION

This program uses Ensembl's REST API to retrieve information including genomic coordinates, alleles and MAF (human only) from variant IDs.

Output is tab-delimited text.

=cut

=head1 AUTHOR

David A. Parry


=head1 COPYRIGHT AND LICENSE

Copyright 2017 David A. Parry

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with this program. If not, see <http://www.gnu.org/licenses/>.

=cut

