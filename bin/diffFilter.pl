#! /usr/bin/perl

#######################################################################
# The MIT License
#
# Copyright (c) 2017, Jérôme Audoux (jerome.audoux@inserm.fr)
#
# Permission is hereby granted, free of charge, to any person obtaining
# a copy of this software and associated documentation files 
# (the “Software”), to deal in the Software without restriction, 
# including without limitation the rights to use, copy, modify, merge,
# publish, distribute, sublicense, and/or sell copies of the Software,
# and to permit persons to whom the Software is furnished to do so, 
# subject to the following conditions:
#
# The above copyright notice and this permission notice shall be 
# included in all copies or substantial portions of the Software.
#
# The Software is provided “as is”, without warranty of any kind, 
# express or implied, including but not limited to the warranties of 
# merchantability, fitness for a particular purpose and 
# noninfringement. In no event shall the authors or copyright holders
# be liable for any claim, damages or other liability, whether in an 
# action of contract, tort or otherwise, arising from, out of or in 
# connection with the software or the use or other dealings in the 
# Software. 
#######################################################################

use strict;
use warnings;

sub getReadingFileHandle {
  my $file = shift;
  my $fh;
  if($file =~ /\.gz$/) {
    open($fh,"gunzip -c $file |") or die ("Cannot open $file");
  } else {
    open($fh,"< $file") or die ("Cannot open $file");
  }
  return $fh;
}

my $diff_table = shift;
my $counts_table = shift;

my $diff_fh   = getReadingFileHandle($diff_table);
my $counts_fh = getReadingFileHandle($counts_table);

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
