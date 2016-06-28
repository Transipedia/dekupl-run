#! /usr/bin/perl
#
use strict;
use warnings;

use CracTools::Utils;

my $diff_table = shift;
my $counts_table = shift;


# Index the diff table
my $diff_fh = CracTools::Utils::getReadingFileHandle($diff_table);
my %diff_hash;
while(<$diff_fh>) {
  my ($tag) = split /\s+/, $_;
  $diff_hash{$tag} = 1;
}

my $counts_fh = CracTools::Utils::getReadingFileHandle($counts_table);
my $header = <$counts_fh>;
print $header;
while(<$counts_fh>) {
  my ($tag) = split /\s+/, $_;
  my $rev_tag = CracTools::Utils::reverseComplement($tag);
  if(!defined $diff_hash{$tag} && !defined $diff_hash{$rev_tag}) {
    print $_;
  }
}
