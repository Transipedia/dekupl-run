#! /usr/bin/perl
#
use strict;
use warnings;

use Pod::Usage;
use Getopt::Long;

use CracTools::Utils;

=head1 NAME

mergeCounts.pl - merge multiple counts files into a single counts table

=head1 SYNOPSIS

  mergeCounts.pl counts1.tsv counts2.tsv > merged-counts.tsv

  Options:

    --min-reccurence           INT  (default : 1)
    --min-reccurence-abundance INT  (default:  1)

=head1 AUTHOR

Jérôme Audoux C<jerome.audoux@inserm.fr>

=cut

my $min_recurrence = 1;
my $min_recurrence_abundance = 1;

GetOptions("min-recurrence=i"           => \$min_recurrence,
           "min-recurrence-abundance=i" => \$min_recurrence_abundance,
) or pod2usage(0);

my @count_files = @ARGV;

pod2usage(0) if !@count_files;

my @counts_files = @ARGV;
my @counts_it;


# Open all counts file
foreach my $f (@counts_files) {
  my $fh = CracTools::Utils::getReadingFileHandle($f);
  push(@counts_it,{fh => $fh});
}
my $min_kmer = undef;

# Init with the first kmer of each files
foreach my $it (@counts_it) {
  my $fh = $it->{fh};
  my $l = <$fh>;
  chomp $l;
  my ($kmer,$count) = split(/\s+/,$l,2);
  $it->{kmer} = $kmer;
  $it->{count} = $count;
  if(!defined $min_kmer || $kmer lt $min_kmer) {
    $min_kmer = $kmer;
  }
}

my $read = 1;
while($read) {
  $read = 0;

  # Get counts for the min_kmer
  my @counts;
  my $new_min_kmer;
  foreach my $it (@counts_it) {
    if(defined $it->{kmer} && $it->{kmer} eq $min_kmer) {
      push @counts, $it->{count};
      # Get next k-kmer in the file
      my $fh = $it->{fh};
      my $l = <$fh>;
      if($l) {
        chomp $l;
        my($kmer,$count) = split(/\s+/,$l,2);
        $it->{kmer} = $kmer;
        $it->{count} = $count;
        $read = 1;
      } else {
        $it->{kmer} = undef;
      }
    } else {
      push @counts, 0;
    }
    if(defined $it->{kmer} && (!defined $new_min_kmer || $it->{kmer} lt $new_min_kmer)) {
      $new_min_kmer = $it->{kmer};
    }
  }

  # Check the reccurence
  my $rec = 0;
  foreach my $c (@counts) {
    if($c >= $min_recurrence_abundance) {
      $rec++
    }
  }
  
  # Print the count
  if($rec >= $min_recurrence) {
    print join("\t",$min_kmer,@counts),"\n";
  }

  $min_kmer = $new_min_kmer;
}

