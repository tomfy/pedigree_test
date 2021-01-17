package Pedigrees;
use strict;
use warnings;
use Moose;
#use Mouse;
use namespace::autoclean;
use Carp;
use Scalar::Util qw (looks_like_number );
use List::Util qw ( min max sum );
use POSIX qw ( floor ceil );
use Genotypes;

has pedigree_filename => (
			  isa => 'Str',
			  is => 'ro',
			  required => 1,
			 );

has accid_parents => (
		      isa => 'HashRef', # keys are accession ids, values: [matid, patid]
		      is => 'rw',
		      default => sub{ {} },
		     );

sub BUILD{
  my $self = shift;
  my $pedigree_filename = $self->pedigree_filename();

  my $accid_parentalids = {};
  open my $fh, "<", "$pedigree_filename" or die "couldn't open $pedigree_filename for reading.\n";
  my $first_line = <$fh>;
  die "pedigree file should have 'Accession' at start of first line.\n" if(! ($first_line =~ /^Accession/));
  while (my $line = <$fh>) {
    my @cols = split(" ", $line);
    my ($accid, $matid, $patid) = @cols[-3, -2, -1];
    next if($accid eq 'NA'  or  $matid eq 'NA'  or $patid eq 'NA'); # only store if both parents given.
    # print "$accid $matid $patid \n";
    my $parental_idpair = [$matid, $patid]; # PedigreeTest::order_idpair($matid, $patid);
    $accid_parentalids->{$accid} = $parental_idpair;
  }
  $self->accid_parents($accid_parentalids);
}

sub parents{
  my $self = shift;
  my $accid = shift;
  my $parents = $self->accid_parents()->{$accid} // undef;
  return (defined $parents)? @$parents : (undef, undef);
}

1;
