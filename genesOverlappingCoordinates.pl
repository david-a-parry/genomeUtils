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
my @coords = ();
my %opts = 
(
    c => \@coords,
    s => "human",
    f => 0,
);
GetOptions
(
    \%opts,
    'c|coords=s{,}',
    'l|list=s',
    's|species=s',
    'r|grch37',
    'f|flanks=i',
    'n|no_coordinates',
    'h|?|help',
    'm|manual',
) or pod2usage( -message => "Syntax error.", -exitval => 2 );
pod2usage( -verbose => 2 ) if ($opts{m});
pod2usage( -verbose => 1 ) if ($opts{h});

pod2usage
(
     -message => "ERROR: c-/--coords argument is required.\n", 
     -exitval => 2 
) if ( not $opts{l} and not @coords);

if ($opts{r}){
    $restQuery->useGRCh37Server();
}

if ($opts{l}){
    open (my $LIST, $opts{l}) or die "Could not open --gene_list '$opts{l}' for reading: $!\n";
    while (my $line = <$LIST>){
        my @s = split(/\s+/, $line); 
        push @coords, $s[0];
    }
}
if (not @coords){
    die "No coordinates identified!\n";
}
checkCoords();

foreach my $c (@coords){
    my $genes = getGenesFromCoords($c);
    if ($opts{n}){
        $genes =~ s/\|/\t/g;
        $genes =~ s/\,/\n/g;
        print "$genes\n";
    }else{
        $genes ||= 'Not Found';
        print "$c\t$genes\n";
    }
}

#########################################################
sub getGenesFromCoords{
    my $region = shift;
    my $endpoint = "overlap/region/$opts{s}/$region?feature=gene";
    my $genes = $restQuery->queryEndpoint($endpoint);
    my %ids = ();
    foreach my $g (@$genes){
        my $id = $g->{id};
        my $sy = $g->{external_name};
        $ids{"$id|$sy"} = undef;
    }
    return join(",", sort keys %ids);
}

#########################################################
sub checkCoords{
    foreach my $c (@coords){
        my $chr;
        my $start;
        my $end; 
        if ($c =~ /^(\w+):([\d\,]+)$/){
            $chr = $1;
            my $pos = $2;
            $pos =~ s/\,//g;
            $start = $pos;
            $end = $pos;
        }elsif ($c =~ /^(\w+):([\d\,]+)-([\d\,]+)$/){
            $chr = $1;
            $start = $2;
            $end = $3;
            $start =~ s/\,//g;
            $end =~ s/\,//g;
            
        }else{
            die "Invalid region: $c\n";
        }
        if ($opts{f}){
            $start -= $opts{f};
            $end += $opts{f};
            $start = $start > 0 ? $start : 1;
        }
        $c = "$chr:$start-$end";
    }
}

#########################################################
=head1 NAME

genesOverlappingCoordinates.pl - retrieve genomic coordinates for given genes

=head1 SYNOPSIS
 
        genesOverlappingCoordinates.pl -c chr1:100000-200000 
        genesOverlappingCoordinates.pl -l coord_list.intervals 
        genesOverlappingCoordinates.pl -h (display help message)
        genesOverlappingCoordinates.pl -m (display manual page)

=cut

=head1 ARGUMENTS

=over 8

=item B<-c    --coords>

One or more intervals to retrieve overlapping genes from in the format 'chr1:1000-2000'.

=item B<-l    --list>

A list of coordinates to use instead of or in conjunction with -c/--coords.

=item B<-n    --no_coordinates>

Use this flag to output a list of gene IDs and symbols without the input 
coordinates.

=item B<-s    --species>

Species genome to use. Defaults = human. 

=item B<-r    --grch37>

Use this option if you want to use the GRCh37 Ensembl REST server for gene queries rather than the current REST server.

=item B<-f    --flanks>

Add this many bp to each end of each interval and search genes overlapping this 
larger region.

=item B<-h    --help>

Show this script's help information.

=item B<-m    --manual>

Show this script's manual page.

=back

=cut

=head1 DESCRIPTION

TODO!

=cut

=head1 AUTHOR

David A. Parry


=head1 COPYRIGHT AND LICENSE

Copyright 2016 David A. Parry

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with this program. If not, see <http://www.gnu.org/licenses/>.

=cut


