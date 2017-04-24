#! /usr/bin/perl
#
use strict;
use warnings;

use CracTools::Utils;

my $diff_table = shift;
my $counts_table = shift;

my $diff_fh = CracTools::Utils::getReadingFileHandle($diff_table);
my $counts_fh = CracTools::Utils::getReadingFileHandle($counts_table);

# Init diff tag
my $diff_line = <$diff_fh>;
my ($diff_tag) = split /\s+/, $diff_line, 2;

# Print header
my $header = <$counts_fh>;
print $header;

while(my $counts_line = <$counts_fh>) {
  my ($tag) = split /\s+/, $counts_line, 2;
  while($diff_tag && $diff_tag lt $tag) {
    $diff_line = <$diff_fh>;
    if($diff_line) {
      ($diff_tag) = split /\s+/, $diff_line, 2;
    } else {
      $diff_tag = undef;
    }
  }
  if(!$diff_tag || $diff_tag ne $tag) {
    print $counts_line;
  }
}
