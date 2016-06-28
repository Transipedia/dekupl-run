#! /usr/bin/perl
#
use strict;
use warnings;

my $i = 0;
while (<>) {
  chomp;
  if($i % 4 == 1) {
    # complement the reversed DNA sequence
    $_ =~ tr/ACGTacgt/TGCAtgca/;
    $_ = reverse $_;
  } elsif($i % 3 == 3) {
    $_ = reverse $_; 
  }
  print $_,"\n";
  $i++;
}
