package PedigreesDAG;
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
use Node;
#use DAG;

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

has id_node => (
		isa => 'HashRef[Maybe[Str]]', # keys are accession ids, values are Node objects.
		is => 'rw',
		default => sub{ {} },
	       );

has childless_ids => (
		      isa => 'HashRef[Maybe[Str]]', #  key: id (string), value: Node Obj
		      is => 'rw',
		      default => sub{ {} },
		     );

has parentless_ids => (
		       isa => 'HashRef[Maybe[Str]]', #  key: id (string), value: Node Obj
		       is => 'rw',
		       default => sub{ {} },
		      );

# has the_dag => (
# 		isa => 'Maybe[DAG]',
# 		is => 'rw',
# 		default => undef,
# 		);

sub BUILD{
  my $self = shift;
  my $pedigree_filename = $self->pedigree_filename();

  my $accid_parentalids = {};
  open my $fh, "<", "$pedigree_filename" or die "couldn't open $pedigree_filename for reading.\n";
  my $first_line = <$fh>;
  die "pedigree file should have 'Accession' at start of first line.\n" if(! ($first_line =~ /^Accession/));
#  my $the_dag = DAG->new();
  while (my $line = <$fh>) {
    my @cols = split(" ", $line);
    my ($accid, $matid, $patid) = @cols[-3, -2, -1];
    next if($accid eq 'NA');

    
    my $matnode = ($matid eq 'NA')? undef : ( $self->id_node()->{$matid} // Node->new({id => $matid}) );
    my $patnode = ($patid eq 'NA')? undef : ( $self->id_node()->{$patid} // Node->new({id => $patid}) );
    $matnode->add_offspring( $accid ) if(defined $matnode);
    $patnode->add_offspring( $accid ) if(defined $patnode);

    my $anode = $self->id_node()->{$accid} // Node->new({id => $accid});
      $anode->female_parent( $matnode );
    $anode->male_parent( $patnode );
    $self->id_node()->{$accid} = $anode;

    $self->add_node($anode);

    next if($accid eq 'NA'  or  $matid eq 'NA'  or $patid eq 'NA'); # only store if both parents given.
    # print "$accid $matid $patid \n";
    my $parental_idpair = [$matid, $patid]; # PedigreeTest::order_idpair($matid, $patid);
    $accid_parentalids->{$accid} = $parental_idpair;
  }
  $self->accid_parents($accid_parentalids);
#  print STDERR "IS IT ACYCLIC?: ", $the_dag->is_it_acyclic(), "\n";
#  print STDERR "N ids: ", scalar keys %{$the_dag->id_node()}, "\n";

    print STDERR "IS IT ACYCLIC?: ", $self->is_it_acyclic(), "\n";
  print STDERR "N ids: ", scalar keys %{$self->id_node()}, "\n";
}

sub add_node{
  my $self = shift;
  my $node = shift;
  my $node_id = $node->id();
  $self->id_node()->{$node_id} = $node;
  if ($node->offspring() eq '') {
    $self->childless_ids()->{$node_id} = $node;
  }
  if (!defined $node->female_parent()  and  !defined $node->male_parent()) {
    $self->parentless_ids()->{$node_id} = $node;
  }
}

sub parents{
  my $self = shift;
  my $accid = shift;
  my $parents = $self->accid_parents()->{$accid} // undef;
  return (defined $parents)? @$parents : (undef, undef);
}

sub as_string{
  my $self = shift;
  my $string = '';
  my @ids = sort keys %{$self->id_node()};
  for my $id (@ids) {
    my $node = $self->id_node()->{$id};
    my ($nwck, $mind, $maxd) = $node->as_newick();
    my @offspring_array = split(" ", $node->offspring());
    $string .= "$id  $mind $maxd   $nwck  " . scalar @offspring_array . "\n";
  }
  return $string;
}

sub is_it_acyclic{
  my $self = shift;
  for my $childless_node ( map($self->childless_ids()->{$_}, sort keys %{$self->childless_ids()} ) ){
    my $res = $childless_node->ancestors_acyclic( {} );
    print STDERR "node ", $childless_node->id(), " and ancestors acyclic?  $res \n\n";
    return 0 if($res == 0);
  }
  return 1;
}

1;
