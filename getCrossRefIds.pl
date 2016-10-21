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
my @gene_ids = ();
my @db = ();
my %opts = 
(
    g => \@gene_ids,
    s => "human",
    a => 0,
    d => \@db,
);
GetOptions
(
    \%opts,
    'l|list=s',
    'g|gene_ids=s{,}',
    's|species=s',
    'r|grch37',
    't|transcript_mode',
    'a|all_levels',
    'd|database_names=s{,}',
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

my @cross_ref = ();
my %dbs       = ();
foreach my $g (@gene_ids){
    getIds($g);
}
outputIds();

#########################################################
sub outputIds{
    if (not @cross_ref){
        die "No IDs found!\n";
    }
    if (not @db){
        outputAllDb();
    }else{
        outputUserDb();
    }
}

#########################################################
sub outputAllDb{
    print join
    (
        "\t", 
        "id",
        "display_name",
        sort keys %dbs,
    ) . "\n";
    foreach my $c (@cross_ref){
        print join("\t", map {$c->{$_} || "." } "id", "display_name", sort keys %dbs) . "\n";
    }   
}

#########################################################
sub outputUserDb{
    my $found_one = 0;
    foreach my $d (@db){
        if (not exists $dbs{$d}){
            warn "WARNING: No db named '$d' found from queries.\n";
        }else{
            $found_one++;
        }
    }
    exit if not $found_one;
    print join
    (
        "\t", 
        "id",
        "display_name",
        @db,
    ) . "\n";
    foreach my $c (@cross_ref){
        print join("\t", map {$c->{$_} || "." } "id", "display_name", @db) . "\n";
    }
}

#########################################################
sub getIds{
    my $id = shift;
    if ($opts{t}){
        foreach my $tr_hash ( $restQuery->getTranscriptDetails($id, $opts{s}) ){
            recordIds($tr_hash);
        }
    }else{
        my $gene_hash = $restQuery->getGeneDetails($id, $opts{s});
        recordIds($gene_hash);
    }
}

#########################################################
sub recordIds{
    my $gene_hash = shift;
    return if not $gene_hash;
    my $xrefs = $restQuery->getXrefs
    (
        id => $gene_hash->{id},
        all_levels => $opts{a},
    );
    my %ids = ();
    foreach my $x (@$xrefs){
        if (not exists $ids{$x->{dbname}}){
            $ids{$x->{dbname}} = $x->{primary_id};
        }else{
            $ids{$x->{dbname}} .= "/$x->{primary_id}";
        } 
        $dbs{$x->{dbname}} = undef;
    }
    $ids{id} = $gene_hash->{id};
    $ids{display_name} = $gene_hash->{display_name};
    push @cross_ref, \%ids;
}


#########################################################

=head1 NAME

getCrossRefIds.pl - retrieve IDs for given genes

=head1 SYNOPSIS
 
        getCrossRefIds.pl -g ABCD1 
        getCrossRefIds.pl -l gene_list.txt 
        getCrossRefIds.pl -h (display help message)
        getCrossRefIds.pl -m (display manual page)

=cut

=head1 ARGUMENTS

=over 8

=item B<-g    --gene_ids>

One or more gene symbols or IDs to search and retrieve IDs for.

=item B<-l    --list>

A list of gene symbols or IDs to use instead of or in conjunction with -g/--gene_ids. 

=item B<-s    --species>

Species genome to use. Defaults = human. 

When using gene symbols or other identifiers such as RefSeq or Uniprot, this program will attempt to identify the corresponding Ensembl gene but will only search this species. Therefore, if you are attempting to find the coordinates of a mouse RefSeq transcript you must use this option to specify 'mouse' as the species or the corresponding Ensembl gene will not be found. 

This option is ignored if Ensembl identifiers are used because they are already species specific.

=item B<-r    --grch37>

Use this option if you want to use the GRCh37 Ensembl REST server for gene queries rather than the current REST server.

=item B<-t    --transcript_mode>

Look up transcript IDs instead of gene IDs.

=item B<-a    --all_levels>

Lookup IDs for all levels of a gene/transcript - i.e. identify all transcript and protein ID cross-references associated with a gene.

=item B<-d    --database_names>

Only output cross-references from these databases.

=item B<-q    --quiet>

Use this flag to supress printing of information to STDERR.

=item B<-h    --help>

Show this script's help information.

=item B<-m    --manual>

Show this script's manual page.

=back

=cut

=head1 EXAMPLES

 getCrossRefIds.pl -g ABCD1
 #get gene cross-ref IDs for human gene ABCD1

 getCrossRefIds.pl -g ABCD1 -s mouse
 #as above for mouse gene Abcd1

 getCrossRefIds.pl -g ABCD1 -t
 #lookup transcript cross-ref IDs

 getCrossRefIds.pl -g ABCD1 -a
 #lookup cross-ref IDs for all transcripts and proteins associated with these gene

 getCrossRefIds.pl -l gene_list.txt
 #get gene cross-ref IDs for a list of human genes

 getCrossRefIds.pl -g ENST00000218104 
 #Ensembl identifiers can be used for genes, transcripts or proteins

 getCrossRefIds.pl -g NM_000033
 #cross reference searches are possible using RefSeq, CCDS, or uniprot identifiers

 getCrossRefIds.pl -g P48410 -s mouse
 #if searching for a cross-reference in a species other than human, you must specify the species

 getCrossRefIds.pl -g ABCD1 -d EntrezGene MIM_MORBID
 #you may specify database IDs to output rather than outputting all found

 
=cut

=head1 DESCRIPTION

This program uses Ensembl's REST API to output database gene/transcript IDs for one or more IDs provided by the user. It can identify genes from gene symbols, Ensembl gene, transcript or protein identifiers as well as RefSeq, CCDS, or Uniprot IDs. By default, all cross-ref IDs identified by Ensembl will be given.

=cut

=head1 AUTHOR

David A. Parry


=head1 COPYRIGHT AND LICENSE

Copyright 2016 David A. Parry

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with this program. If not, see <http://www.gnu.org/licenses/>.

=cut


