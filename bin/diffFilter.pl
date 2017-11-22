#! /usr/bin/perl

#######################################################################
# This file is part of Dekupl. 
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# Jérôme Audoux (jerome.audoux@inserm.fr), Copyright (C) 2017  
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
