#! /usr/bin/perl
#
use strict;
use warnings;
use Pod::Usage;
use Getopt::Long;
use CracTools::Utils;
use CracTools::SAMReader;
use CracTools::Annotator;

=head1 NAME

annotateTags.pl - Annotate aligned tags (BAM) with a GFF file.

=head1 VERSION

version 0.001

=head1 SYNOPSIS

annotateTags.pl [options] file.bam annotations.gff

=head1 OPTIONS

  --max-dist-intergenic-tag-fusion INT  Max distance to consider two collinear
                                        intergenic tags as part of the same transcriptional 
                                        unit. [DEFAULT: 1000]
  --max-neighborhood-dist INT           Max distance to consider an intergenic tag as part 
                                        of an annotated gene in its neighborhood. [DEFAULT: 10000]
          
=head1 AUTHOR

Jérôme Audoux <jerome.audoux@gmail.com>

=head1 COPYRIGHT AND LICENSE

This software is Copyright (c) 2016 by Jérôme Audoux.

This is free software, licensed under:

  The GNU General Public License, Version 3, June 2007

=cut

my $max_dist_intergenic_tag_fusion = 1000;
my $max_neiborhood_dist            = 10000;

GetOptions( "max-dist-intergenic-tag-fusion=i"  => \$max_dist_intergenic_tag_fusion,
            "max-neighborhood-dist=i"           => \$max_neiborhood_dist,
          ) or pod2usage(-verbose => 1);

my $bam_file = shift @ARGV;
my $gff_file = shift @ARGV;

pod2usage(
  -message => "Mandatory argument 'file.bam' is missing",
  -verbose => 1,
) unless defined $bam_file;

pod2usage(
  -message => "Mandatory argument 'annotations.gff' is missing",
  -verbose => 1,
) unless defined $gff_file;

my $sam_reader = CracTools::SAMReader->new($bam_file);

print STDERR "Loading GFF annotations\n";
my $annotator  = CracTools::Annotator->new($gff_file, "fast");
#sleep 10;

print STDERR "Annotating tags\n";

print join("\t",qw(tag genomic_location annotation gene_name gene_id alignment)),"\n";

my $prev_tag_pos;
my $prev_tag_chr;
my $new_gene_id = 0;
my $sam_it     = $sam_reader->iterator();
while(my $line = $sam_it->()) {
  # Only consider primary alignements
  next if $line->isFlagged(256) || $line->isFlagged(2048);
  my $cigar_counts = $line->getCigarOperatorsCount();
  my $start = $line->pos;
  my $end   = $line->pos;
  map { $end += $cigar_counts->{$_} if defined $cigar_counts->{$_} } qw([M D N X =);
  
  # Retrive annotation from the alignment
  my $alignment_type = "NA";
  if(defined $cigar_counts->{D} && $cigar_counts->{D} > 0) {
    $alignment_type = "DELETION"; 
  } elsif(defined $cigar_counts->{I} && $cigar_counts->{I} > 0) {
    $alignment_type = "INSERTION";
  } elsif(defined $cigar_counts->{N} && $cigar_counts->{N} > 0) {
    $alignment_type = "SPLICE";
  } elsif(defined $line->getOptionalField('NM') && $line->getOptionalField('NM') > 0) {
    $alignment_type = "SUBSTITUTION";
  } elsif(scalar keys %{$cigar_counts} == 1 && defined $cigar_counts->{M}) {
    $alignment_type = "PERFECT_MATCH";
  } elsif(defined $cigar_counts->{S} && $cigar_counts->{S} > 0) {
    $alignment_type = "SOFTCLIPPED";
  }

  # Retrieve annotation the GFF
  my $annotation_type = "NA";
  my $gene_id   = "NA";
  my $gene_name = "NA";
  my $genomic_location = "NA";

  if($line->rname ne '*') {
    $genomic_location = $line->chr.":".$start."-".$end;
    
    # Look for annotations on the forward strand
    my ($annotation, $priority, $type) = $annotator->getBestAnnotationCandidate(
      $line->chr,
      $start,
      $end,
      1,
      \&getCandidatePriorityDefault
    );

    # Look for annotations on the reverse strand
    if(!defined $annotation) { 
      ($annotation, $priority, $type) = $annotator->getBestAnnotationCandidate(
        $line->chr,
        $start,
        $end,
        -1,
        \&getCandidatePriorityDefault
      );
      if(defined $annotation) {
      }
    }
    
    # Look for annotation of the 5prim neighborhood
    if(!defined $annotation) {
      my $annotation_forward = $annotator->getAnnotationNearestDownCandidates(
        $line->chr,
        $start,
        1,
      )->[0];
      my $annotation_reverse = $annotator->getAnnotationNearestDownCandidates(
        $line->chr,
        $start,
        -1,
      )->[0];

      # Decide wether forward or reverse annotation is closer to our tag
      if(defined $annotation_forward && defined $annotation_forward->{gene}) {
        if(defined $annotation_reverse && defined $annotation_reverse->{gene} &&
          $annotation_reverse->{gene}->{end} > $annotation_forward->{gene}->{end}) {
          $annotation = $annotation_reverse;
        } else {
          $annotation = $annotation_forward;
        }
      } else {
        $annotation = $annotation_reverse;
      }
      
      # if annotation is defined and close enough to our tag, we report it
      if(defined $annotation && defined $annotation->{gene} &&
        ($start - $annotation->{gene}->end) < $max_neiborhood_dist) {
        $type = "DOWNSTREAM_GENE";
      } else {
        $annotation = undef;
      }
    }
    
    # Look for annotation of the 3prim neighborhood
    if(!defined $annotation) {
      my $annotation_forward = $annotator->getAnnotationNearestUpCandidates(
        $line->chr,
        $end,
        1,
      )->[0];
      my $annotation_reverse = $annotator->getAnnotationNearestUpCandidates(
        $line->chr,
        $end,
        -1,
      )->[0];

      # Decide wether forward or reverse annotation is closer to our tag
      if(defined $annotation_forward && defined $annotation_forward->{gene}) {
        if(defined $annotation_reverse && defined $annotation_reverse->{gene} &&
          $annotation_reverse->{gene}->{start} < $annotation_forward->{gene}->{start}) {
          $annotation = $annotation_reverse;
        } else {
          $annotation = $annotation_forward;
        }
      } else {
        $annotation = $annotation_reverse;
      }
      if(defined $annotation && defined $annotation->{gene} &&
        ($annotation->{gene}->start - $end) < $max_neiborhood_dist) {
        $type = "UPSTREAM_GENE";
      } else {
        $annotation = undef;
      }
    }

    if(defined $annotation) {
      $annotation_type = $type;
      if(defined $annotation->{gene}) {
        $gene_name  = $annotation->{gene}->attribute('Name');
        ($gene_id)  = $annotation->{gene}->attribute('ID') =~ /gene:(.*)/;
      }
    } else {
      $annotation_type = "INTERGENIC";
      if(defined $prev_tag_pos &&
        ($prev_tag_chr ne $line->chr ||
        ($start - $prev_tag_pos) > $max_dist_intergenic_tag_fusion)) {
        $new_gene_id++;
      } 
      $gene_id = "TSPG$new_gene_id";
      $gene_name = $gene_id;
      $prev_tag_pos = $start;
      $prev_tag_chr = $line->chr;
    }
  }

  print join("\t",$line->qname,$genomic_location,$annotation_type,$gene_name,$gene_id,$alignment_type),"\n";
}

sub getCandidatePriorityDefault {
  my ($pos_start,$pos_end,$candidate) = @_;
  my ($priority,$type) = (-1,'');
  my ($gene,$exon) = ($candidate->{gene},$candidate->{exon});
  if(defined $gene) {
    if(defined $gene->attribute('biotype') && $gene->attribute('biotype') =~ /protein_coding/i) {
      if(defined $exon) {
        if(($exon->start <= $pos_start) && ($exon->end >= $pos_end)) {
          $priority = 1;
          if(defined $candidate->{three}) {
            $type = '3PRIM_UTR';
          } elsif(defined $candidate->{five}) {
            $type = '5PRIM_UTR';
          } elsif(defined $candidate->{CDS}) {
            $type = 'CDS';
          } else {
            $type = 'EXON';
          }
        } else {
          $priority = 2;
          $type = 'INXON';
        }
      } else {
        $priority = 4;
        $type = 'INTRON';
      }
    } else {
      if(defined $exon) {
        if(($exon->start <= $pos_start) && ($exon->end >= $pos_end)) {
          $priority = 3;
          $type = 'NON_CODING_EXON';
        }
      } else {
          $priority = 5;
          $type = 'NON_CODING_INTRON';
      }
    }
  }
  return ($priority,$type);
}
