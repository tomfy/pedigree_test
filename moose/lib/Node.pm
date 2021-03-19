package Node;
use strict;
use warnings;
use Moose;
#use Mouse;
use namespace::autoclean;
use Carp;
use List::Util qw ( min max sum );

# A node representing part of a pedigree,
# it knows its Id, and optionally 

has id => (
	   isa => 'Str',
	   is => 'rw',
	   default => '',
	  );

has female_parent => (
		      isa => 'Maybe[Node]', # Node object
		      is => 'rw',
		      default => undef,
		     );

has male_parent => (
		    isa => 'Maybe[Node]',
		    is => 'rw',
		    default => undef,
		   );

has offspring => (
		  isa => 'Str',
		  is => 'rw',
		  default => '',
		  );

sub as_newick{
  my $self = shift;
  my $newick = '';
  my ($mindF, $maxdF, $mindM, $maxdM);
  if (defined $self->female_parent()) {
    (my $nwck, $mindF, $maxdF) = $self->female_parent()->as_newick();
    $newick .= $nwck;
    $mindF++; $maxdF++;
  } else {
    $newick .= 'NA';
    $mindF = 0; $maxdF = 0;
  }
  $newick .= ',';
  if (defined $self->male_parent()) {
    (my $nwck, $mindM, $maxdM) =  $self->male_parent()->as_newick();
    $newick .= $nwck;
       $mindM++; $maxdM++;
  } else {
    $newick .= 'NA';
       $mindM = 0; $maxdM = 0;
  }
  return ( '(' . $newick . ')' . $self->id, min($mindF, $mindM), max($maxdF, $maxdM) );
}

sub add_offspring{
  my $self = shift;
  my $new_offspring_id = shift;
  my $updated_offspring_id_string = $self->offspring() . "$new_offspring_id ";
  $self->offspring($updated_offspring_id_string);
  # print STDERR "offspringstring: ", $self->offspring(), "\n";
}

1;
