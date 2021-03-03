package Genotypes;
use strict;
use warnings;
use Moose;
#use Mouse;
use namespace::autoclean;
use Carp;
use Scalar::Util qw (looks_like_number );
use List::Util qw ( min max sum );
use POSIX qw ( floor ceil );

# accession id and string of genotypes

has id => ( isa => 'Str',
	    is => 'ro',
	    required => 1,
	  );

# has genotypes => ( # e.g. '0100X021' i.e. only 0, 1, 2, X, and x with no separating characters.
#                   isa => 'Str',
#                   is => 'rw',
# 		  required => 1,
# 		 );

has genotypes01234 => ( # e.g. '01003021' i.e. only 0, 1, 2, 3, 4 with no separating characters.
		       isa => 'Str',
		       is => 'rw',
		       required => 1,
		      );

has quality_counts => ( # e.g. '1000 700 200 5 3' 4nd and 5th numbers are counts of gts with ambiguous gts.
		       isa => 'Str',
		       is => 'rw',
		       required => 1,
		      );

sub as_string{
  my $self = shift;
  my $s = $self->id() . "  " . $self->quality_counts() . "  " . $self->genotypes();
}

sub agmr_hgmr{
  my $self = shift;
  my $other = shift;
  my $gstr1 = $self->genotypes01234();
  my $gstr2 = $other->genotypes01234();
  my ($agmr, $hgmr) = (0, 0);
  my ($adenom, $hdenom) = (0, 0);
  for (my $i = 0; $i < length $gstr1; $i++) {
    my $g1 = substr($gstr1, $i, 1);
    my $g2 = substr($gstr2, $i, 1);
    next if($g1 >= 3  or $g2 >= 3); # count only if neither is bad
    $adenom++;
    $agmr++ if($g1 != $g2);
    next if($g1 == 1  or $g2 == 1);
    $hdenom++;
    $hgmr++ if($g1 != $g2);
  }
  $agmr = ($adenom > 0)? $agmr/$adenom : -1;
  $hgmr = ($hdenom > 0)? $hgmr/$hdenom : -1;
  return [$agmr, $hgmr];
}

sub remove_bad_markers{
  my $self = shift;
  my $marker_01s = shift; # array ref, with ids of 'bad' markers replaced with undef.
  my $new_genotypes_string = '';
  my @new_qual_counts = (0, 0, 0, 0, 0); # for updated accession gt quality counts
  my $L = length $self->genotypes01234();
  # print STDERR "L: $L  ", scalar @$marker_01s, "\n";
  die if(scalar @$marker_01s != $L);
  my $gtstr = $self->genotypes01234();
  for my $i (0..$L-1) {
    my $gt = substr($gtstr, $i, 1);
    if ($marker_01s->[$i] == 1) { # this marker is Ok
      $new_genotypes_string .= substr($gtstr, $i, 1);
      $new_qual_counts[$gt]++;
    } else {			# this marker is 'bad'
      # print STDERR "marker $i is bad.\n"; # die;
    }
  }
  $self->genotypes01234($new_genotypes_string);
  $self->quality_counts(join(" ", @new_qual_counts));
  # print STDERR "new quality counts:  ", join(" ", @new_qual_counts), "\n";
  die if(length $new_genotypes_string  !=  sum(@new_qual_counts));
}



1;
