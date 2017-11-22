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

use Pod::Usage;
use Getopt::Long;

use CracTools::Utils;
use CracTools::Output;

=head1 NAME

mergeCounts.pl - merge multiple counts files into a single counts table

=head1 SYNOPSIS

  mergeCounts.pl counts1.tsv counts2.tsv [--sep "\t"] [--field 1] [--no-header] > merged-counts.tsv

  Options:

    --sep               (default : "\t") Input separator
    --output-sep        (default:  "\t") Output separator
    --field             (default : 1)
    --no-header         Tell convertList.pl to keep the header (as-is)
    --min-counts=<int>  Filter-out rows that have no samples with a count higher to this threshold
                        (default : 0)

=head1 AUTHOR

Jérôme Audoux C<jerome.audoux@inserm.fr>

=cut

my $sep = "\t";
my $output_sep = "\t";
my $field = 1;
my $no_header = 0;
my $min_counts = 0;

GetOptions("sep=s"                => \$sep,
           "output-sep=s"         => \$output_sep,
           "field=i"              => \$field,
           "no-header"            => \$no_header,
           "min-counts=i"         => \$min_counts,
) or pod2usage(0);

my @count_files = @ARGV;

pod2usage(0) if !@count_files;

my %counts;
my @samples;

foreach my $file (@count_files) {
  my $fh = CracTools::Utils::getReadingFileHandle($file);
  my @headers;
  while(<$fh>) {
    next if $_ =~ /^#/;
    chomp;
    my($feat,@values) = split($sep,$_);

    # If this is the first line, the we can retrieve
    # the header line
    if(!@headers) {
      if($no_header) {
        @headers = map { "$file"."$_" } (1..scalar @values);
      } else {
        @headers = @values;
      }
      next;
    }
    
    for(my $i = 0; $i < @values; $i++) {
      $counts{$feat}{$headers[$i]} = $values[$i];
    }
  }
  push @samples,@headers;
}

my $output = CracTools::Output->new();
$output->printLine("feature",@samples) if !$no_header;
foreach my $feature (keys %counts) {
  my $max_row_count = 0;
  map { $max_row_count = $counts{$feature}{$_} if $counts{$feature}{$_} > $max_row_count } @samples;
  if($max_row_count >= $min_counts) {
    $output->printLine($feature,map { defined $counts{$feature}{$_}? $counts{$feature}{$_} : 0} @samples);
  }
}
