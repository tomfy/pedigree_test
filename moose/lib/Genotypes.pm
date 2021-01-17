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

has genotypes => ( # e.g. '0100X021' i.e. only 0, 1, 2, X, and x with no separating characters.
                  isa => 'Str',
                  is => 'rw',
		  required => 1,
		 );

has quality_counts => ( # e.g. '1000 700 200 5 3' 4nd and 5th numbers are counts of gts mapping to X, x resp.
		       isa => 'Str',
		       is => 'rw',
		       required => 1,
		      );

sub as_string{
  my $self = shift;
  my $s = $self->id() . "  " . $self->quality_counts() . "  " . $self->genotypes();
}

# sub distances{
#   my $self = shift;
#   my $other = shift;
#   my $gts1 = $self->genotypes();
#   my $gts2 = $other->genotypes();
# die if(length $gts1 != length $gts2);
#   for my $i (0 .. length $gts1){
#     my $g1 = substr

sub agmr_hgmr{
  my $self = shift;
  my $other = shift;
  my $gstr1 = $self->genotypes();
  my $gstr2 = $other->genotypes();
  my ($agmr, $hgmr) = (0, 0);
  my ($adenom, $hdenom) = (0, 0);
  for (my $i = 0; $i < length $gstr1; $i++) {
    my $g1 = substr($gstr1, $i, 1);
    my $g2 = substr($gstr2, $i, 1);
    next if(uc $g1 eq 'X'  or uc $g2 eq 'X'); # count only if neither is X
    $adenom++;
    $agmr += ($g1 == $g2)? 0 : 1;
    if ($g1 == 0) {
      if ($g2 == 0) {
        $hdenom++;
      } elsif ($g2 == 2) {
        $hgmr++;
        $hdenom++;
      }
    } elsif ($g1 == 2) {
      if ($g2 == 2) {
        $hdenom++;
      } elsif ($g2 == 0) {
        $hgmr++;
        $hdenom++;
      }
    }
  }
  $agmr = ($adenom > 0)? $agmr/$adenom : -1;
  $hgmr = ($hdenom > 0)? $hgmr/$hdenom : -1;
  return [$agmr, $hgmr];
}

1
  ;
