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
    's|species',
    'd|defaultRestServer',
    'j|join_regions',
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

$restQuery->useGRCh37Server() unless $opts{d};

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

if ($opts{j}){
    my $reg_obj = SortGenomicCoordinates->new
    ( 
        array => \@regions, 
        type => 'bed', 
        col => 1, 
    );
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
}
foreach my $r (@regions){
    print "$r\n";
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
