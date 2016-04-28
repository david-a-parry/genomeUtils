#!/usr/bin/env perl
use warnings;
use strict;
use Getopt::Long;
use Data::Dumper;
use Pod::Usage;
use Getopt::Long;
use FindBin qw($RealBin);
use lib "$RealBin/lib/dapPerlGenomicLib";
use SortGenomicCoordinates;
use IdParser;
use EnsemblRestQuery;

my $id_parser = new IdParser();
my $restQuery = new EnsemblRestQuery();
my @gene_ids = ();
my %opts = 
(
    g => \@gene_ids,
    s => "human",
);
GetOptions
(
    \%opts,
    'l|list=s',
    'g|gene_ids=s{,}',
    's|species=s',
    'r|grch37',
    'j|join_regions',
    'q|quiet',
    'x|sort_regions',
    'h|?|help',
    'm|manual',
) or pod2usage( -message => "Syntax error.", -exitval => 2 );
pod2usage( -verbose => 2 ) if ($opts{m});
pod2usage( -verbose => 1 ) if ($opts{h});

pod2usage
(
     -message => "ERROR: Either -l/--list or -g/--gene_ids argument is required.\n", 
     -exitval => 2 
) if ( not $opts{l} and not @gene_ids);

if ($opts{r}){
    $restQuery->useGRCh37Server();
}

if ($opts{l}){
    open (my $GENES, $opts{l}) or die "Could not open --gene_list '$opts{l}' for reading: $!\n";
    while (my $line = <$GENES>){
        my @s = split(/\s+/, $line); 
        push @gene_ids, $s[0];
    }
}
if (not @gene_ids){
    die "No gene IDs identified!\n";
}

my @regions = ();
foreach my $g (@gene_ids){
    my $region = get_region_from_gene($g);
    push @regions, $region if defined $region;
}

if ($opts{j} or $opts{x}){
    my $reg_obj = SortGenomicCoordinates->new
    ( 
        array => \@regions, 
        type => 'bed', 
        col => 1, 
    );
    if ($opts{j}){
        my $region_ref = $reg_obj->prep();
        @regions = ();
        foreach my $reg (@$region_ref) {
            my %strand = ();
            my $str;
            my @names = ();
            foreach my $inf (@{$reg->{info}}){
                my @i = split("\t", $inf);
                push @names,  $i[4];
                $strand{$i[3]}++;
            }
            
            if (keys %strand > 1){
                my $n = 0;
                foreach my $k (keys %strand){
                    $strand{$k} > $n ? $str = $k : $str = $str;
                    $n = $strand{$k};
                }
            }else{
                $str = (keys %strand)[0];
            }
            my $out = join(
                "\t", 
                $reg->{chrom}, 
                $reg->{start},
                $reg->{end},
                $str,
                join("/", @names) ,
            );
            push @regions, $out;
        }
    }elsif($opts{x}){
        my $region_ref = $reg_obj->order();
        @regions = @$region_ref;
    }
}

foreach my $r (@regions){
    print "$r\n";
}

#########################################################
sub geneFromEnst{
    my $id = shift;
    if (not $opts{q}){
        print STDERR "Identifying parent gene from Ensembl transcript $id...\n";
    }
    return $restQuery->getParent($id, 1);
}

#########################################################
sub geneFromEnsp{
    my $id = shift;
    if (not $opts{q}){
        print STDERR "Identifying parent gene from Ensembl protein $id...\n";
    }
    my $par = $restQuery->getParent($id);
    if ($par){
        if (exists $par->{id}){
            return geneFromEnst($par->{id});
        }
    }
}


#########################################################
sub get_region_from_gene{
    my $id = shift;
    $id_parser->parseId($id);
    my $gene_hash; 
    my @lookups = ();
    if (not $opts{q}){
        print STDERR "Interpretting ID \"$id\" as of type \"" . 
          $id_parser->get_identifierType() . "\"...\n";
    }
    if ($id_parser->get_isEnsemblId()){
        if ( $id_parser->get_isTranscript() ){
            $gene_hash = geneFromEnst($id);
        }elsif( $id_parser->get_isProtein() ) {
            $gene_hash = geneFromEnsp($id);
        }else{
            $gene_hash = $restQuery->lookUpEnsId($id, 1);
        }
    }elsif($id_parser->get_isTranscript()  or $id_parser->get_isProtein() ) {
        if (not $opts{q}){
            print STDERR "Identifying Ensembl gene via transcript cross-reference...\n";
        }
        my $transcript = $restQuery->getTranscriptViaXreg($id, $opts{s});
        if ($transcript and ref $transcript eq 'ARRAY'){
            if (@$transcript > 1){
                print STDERR "WARNING: Multiple transcripts identified by ".
                  "cross-reference search for $id - picking the first.\n";
            }
            my $tr = $transcript->[0];
            if (exists $tr->{id}){
                $gene_hash = geneFromEnst($tr->{id});
            }
        }else{
            if (not $opts{s}){
                print STDERR "WARNING: No transcript identified for ID \"$id\"\n";
            }
        }
    }else{
        if (not $opts{q}){
            print STDERR "Identifying Ensembl gene via gene cross-reference...\n";
        }
        my $gene = $restQuery->getGeneViaXreg($id, $opts{s});
        if (ref $gene eq 'ARRAY'){
            foreach my $ge (@$gene){
                if ($ge->{id}){
                    my $ge_hash = $restQuery->lookUpEnsId($ge->{id}, 1);
                    if (uc($ge_hash->{display_name}) eq uc($id)){
                    #if gene symbol matches then we use this entry
                        $gene_hash = $ge_hash;
                        last;
                    }else{
                        push @lookups, $ge_hash;
                    }
                }
            }
            if (not $gene_hash){
                if (@lookups == 1){
                    $gene_hash = $lookups[0];
                }
            }
        }
    }
    if (not $gene_hash){
        print STDERR "WARNING: Could not identify gene for ID \"$id\"\n";
        if (@lookups){
            my $idstring = join("\n", map { $_->{display_name} } @lookups );
            print STDERR "Identified the following non-matching display names:\n".
                         "$idstring\n";
        }
        return;
    }
    my $strand = $gene_hash->{strand} > 0 ? "+" : "-";
    my $r = join
    (
        "\t", 
        $gene_hash->{seq_region_name},
        $gene_hash->{start},
        $gene_hash->{end},
        $strand,
        $gene_hash->{display_name} . "/" . $gene_hash->{id},
    );
    if (not $opts{q}){
        my $coord = "$gene_hash->{seq_region_name}:$gene_hash->{start}-$gene_hash->{end}";
        print STDERR "For gene ID '$id', found gene " . 
          $gene_hash->{display_name} . "/" . $gene_hash->{id} .
          " with coordintes $coord ($gene_hash->{assembly_name})\n";
    }
    return $r;
}

#!/usr/bin/env perl

=head1 NAME

coordinatesFromGenes.pl - retrieve genomic coordinates for given genes

=head1 SYNOPSIS
 
        coordinatesFromGenes.pl -g ABCD1 
        coordinatesFromGenes.pl -l gene_list.txt 
        coordinatesFromGenes.pl -h (display help message)
        coordinatesFromGenes.pl -m (display manual page)

=cut

=head1 ARGUMENTS

=over 8

=item B<-g    --gene_ids>

One or more gene symbols or IDs to search and retrieve regions.

=item B<-l    --list>

A list of gene symbols or IDs to use instead of or in conjunction with -g/--gene_ids. 

=item B<-s    --species>

Species genome to use. Defaults = human. 

When using gene symbols or other identifiers such as RefSeq or Uniprot, this program will attempt to identify the corresponding Ensembl gene but will only search this species. Therefore, if you are attempting to find the coordinates of a mouse RefSeq transcript you must use this option to specify 'mouse' as the species or the corresponding Ensembl gene will not be found. 

This option is ignored if Ensembl identifiers are used because they are already species specific.

=item B<-r    --grch37>

Use this option if you want to use the GRCh37 Ensembl REST server for gene queries rather than the current REST server.

=item B<-x    --sort_regions>

Use this flag to sort regions in coordinate order in the output.

=item B<-j    --join_regions>

Use this flag to merge overlapping regions in the output. Output will also be sorted.

=item B<-q    --quiet>

Use this flag to supress printing of information to STDERR.

=item B<-h    --help>

Show this script's help information.

=item B<-m    --manual>

Show this script's manual page.

=back

=cut

=head1 EXAMPLES

 coordinatesFromGenes.pl -g ABCD1
 #get human GRCh38 coordinates for ABCD1

 coordinatesFromGenes.pl -g ABCD1 -s mouse
 #get mouse genome coordinates for Abcd1

 coordinatesFromGenes.pl -l gene_list.txt
 #get human coordinates for a list of genes

 coordinatesFromGenes.pl -g ENST00000218104 
 #Ensembl identifiers can be used for genes, transcripts or proteins

 coordinatesFromGenes.pl -g NM_000033
 #cross reference searches are possible using RefSeq, CCDS, or uniprot identifiers

 coordinatesFromGenes.pl -g P48410 -s mouse
 #if searching for a cross-reference in a species other than human, you must specify the species

 
=cut

=head1 DESCRIPTION

This program uses Ensembl's REST API to identify gene coordinates for one or more genes provided by the user. It can identify genes from gene symbols, Ensembl gene, transcript or protein identifiers as well as RefSeq, CCDS, or Uniprot IDs.

Output is in BED format.

=cut

=head1 AUTHOR

David A. Parry


=head1 COPYRIGHT AND LICENSE

Copyright 2016 David A. Parry

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with this program. If not, see <http://www.gnu.org/licenses/>.

=cut


