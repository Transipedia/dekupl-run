#! /usr/bin/perl
#
use strict;
use warnings;

use CracTools::Utils;

my $counts_file = shift;

my $counts_fh = CracTools::Utils::getReadingFileHandle($counts_file);
my $header = <$counts_fh>;

print STDERR "Loading kmers into memory\n";
# First loop over the file and store k-mers into a hash
my %kmers;
while(<$counts_fh>) {
  chomp;
  #my ($pvalue,$kmer,@counts) = split(/\s+/,$_);
  my ($kmer,$pvalue,@counts) = split(/\s+/,$_);
  $pvalue = sprintf("%.10g", $pvalue);
  $kmers{$kmer} = { pvalue => $pvalue, counts => \@counts };
}

my $pvalue_treshold = 0.001;

print STDERR "Merging k-mers\n";
print "nb_merged_kmers\tassembly\t$header";
# Then we loop over the kmers with k-1 overlapping sequence and close pvalues
my @DNA_alphabet = qw(A T G C);
foreach my $kmer (keys %kmers) {
  next if defined $kmers{$kmer}->{assembly_kmer};

  # K-mer used for the output of the aggregated contig
  my $output_kmer = $kmer;
  my $nb_merged   = 1;
  my $assembly    = $kmer;

  # Finding neighbors
  foreach my $kmer_ext ("left","right") {
    my $curr_kmer = $kmer;
    my $found_neighbor = 1;
    while($found_neighbor) {
      $found_neighbor = 0;
      my $neighbor_kmer;
      my $neighbor_orientation;
      foreach my $kmer_orientation ("normal","revcom") {
        my $curr_kmer_oriented = $kmer_orientation eq "normal"?
                                 $curr_kmer :
                                 CracTools::Utils::reverseComplement($curr_kmer);
        foreach my $nuc (@DNA_alphabet) {
          my $kmer_seq;
          if($kmer_orientation eq "normal") {
            if($kmer_ext eq "left") {
              $kmer_seq = $nuc.substr($curr_kmer_oriented, 0, -1);
            } else {
              $kmer_seq = substr($curr_kmer_oriented, 1).$nuc;
            }
          } else {
            if($kmer_ext eq "left") {
              $kmer_seq = substr($curr_kmer_oriented, 1).$nuc;
            } else {
              $kmer_seq = $nuc.substr($curr_kmer_oriented, 0, -1);
            }
          }
          if(defined $kmers{$kmer_seq}) {
            #abs($kmers{$kmer_dir}->{pvalue} - $kmers{$kmer}->{pvalue}) < $pvalue_treshold) {
            if($found_neighbor) {
              $found_neighbor++;
              last;
            } else {
              $found_neighbor = 1;
              $neighbor_kmer = $kmer_seq;
              $neighbor_orientation = $kmer_orientation;
            }
          }
          last if $found_neighbor > 1;
        }
        last if $found_neighbor > 1;
      }

      # Now we can assemble
      if($found_neighbor == 1 && !defined $kmers{$neighbor_kmer}->{assembly_kmer}) {
        my $neighbor_kmer_oriented  = $neighbor_orientation eq "normal"? 
                                      $neighbor_kmer : 
                                      CracTools::Utils::reverseComplement($neighbor_kmer);
        # Assemble this k-mer with the current assembly
        if($kmer_ext eq "left") {
          $assembly = substr($neighbor_kmer_oriented, 0, 1).$assembly;
        } else {
          $assembly = $assembly.substr($neighbor_kmer_oriented, -1);
        }
        # Set the output kmer has the one with the lowest pvalue
        if($kmers{$neighbor_kmer}->{pvalue} < $kmers{$output_kmer}->{pvalue}) {
          $output_kmer = $neighbor_kmer;
        }
        $nb_merged++;
        # Mark this kmer has already used
        $kmers{$neighbor_kmer}->{assembly_kmer} = $kmer;
        # Set the new k-mer to continue assembling
        $curr_kmer = $neighbor_kmer_oriented;
      } else {
        $found_neighbor = 0;
      }
    }
  }
  print join("\t",
    $nb_merged,
    $assembly,
    $output_kmer,
    $kmers{$output_kmer}->{pvalue},
    @{$kmers{$output_kmer}->{counts}}
  ),"\n";
}
