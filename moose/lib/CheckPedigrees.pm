package CheckPedigrees;
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
use GenotypesSet;
use Pedigrees;
use PedigreeCheck;

# accession id and string of genotypes

has genotypes_set => (
		     isa => 'Object',
		     is => 'ro',
		     required => 1,
		    );

has pedigrees => (
		  isa =>  'Object',
		  is => 'ro',
		  required => 1,
		 );

has pedigree_checks => (
			isa => 'HashRef',
			is => 'rw',
			default => sub { {} },
		       );

has matrand_checks => (
		       isa => 'HashRef',
		       is => 'rw',
		       default => sub { {} },
		      );

has patrand_checks => (
		       isa => 'HashRef',
		       is => 'rw',
		       default => sub { {} },
		      );

has randrand_checks => (
			isa => 'HashRef',
			is => 'rw',
			default => sub { {} },
		       );

sub BUILD{
  my $self = shift;

  my $accid_gts = $self->genotypes_set()->accid_genotypes(); # values are Genotypes objects
  my @accids = keys %$accid_gts;
  # while (my($accid, $gtobj) = each %$accid_gts) {
  #   print STDERR "$accid  ", $gtobj->id(), "\n";
  # }
  my $peds = $self->pedigrees();
  my $accid_parents = $peds->accid_parents();

  my $n_pedigree_checks = 0;
  while (my($accid, $accgtobj) = each %$accid_gts) {

    my ($matid, $patid) = $peds->parents($accid); # [matid, patid] if defined
    #   print STDERR "[$matid] [$patid] \n";
    next if(!defined $matid);

    my $matgtobj = $accid_gts->{$matid} // undef;
    next if (!defined $matgtobj);
    my $patgtobj = $accid_gts->{$patid} // undef;
    next if (!defined $patgtobj);

    $self->pedigree_checks()->{$accid} = PedigreeCheck->new({ mat_gtsobj => $matgtobj, pat_gtsobj => $patgtobj, acc_gtsobj => $accgtobj } );
    $n_pedigree_checks++;
    if($n_pedigree_checks % 100 == 0){
      print STDERR "# n pedigree checks created: $n_pedigree_checks \n";
    }
    #_get_a_pedcheck($accid_gts, $matid, $patid, $accid);

  #  $self->matrand_checks()->{$accid} = _get_a_pedcheck($accid_gts, $matid, $accids[int(rand(@accids))], $accid);

  #  $self->patrand_checks()->{$accid} = _get_a_pedcheck($accid_gts, $accids[int(rand(@accids))], $patid, $accid);

  #  $self->randrand_checks()->{$accid} = _get_a_pedcheck($accid_gts, $accids[int(rand(@accids))], $accids[int(rand(@accids))], $accid);

    # print $pedcheck->as_string(), "\n";
  }
}

sub as_string{
  my $self = shift;
   my $gtsetobj = $self->genotypes_set();
  my $the_string = '';
  $the_string .= "#n_accessions: " . $gtsetobj->n_accessions() . "\n";
  $the_string .= "#n_markers: " . $gtsetobj->n_markers() . "\n";
  $the_string .= "#delta: " . $gtsetobj->delta() . "\n";

  while (my ($accid, $pedchk) = each %{$self->pedigree_checks()}) {
 #   my $gtsobj = $gtsetobj->accid_genotypes()->{$accid};
#    $the_string .= $accid . "  " . $gtsobj->quality_counts() . "  ";
    $the_string .= $pedchk->as_string() . "\n";
  }
  return $the_string;
}



sub as_string_ns{
  my $self = shift;
  my $gtsetobj = $self->genotypes_set();
  my $the_string = '';
  $the_string .= "#n_accessions: " . $gtsetobj->n_accessions() . "\n";
  $the_string .= "#n_markers: " . $gtsetobj->n_markers() . "\n";
  $the_string .= "#delta: " . $gtsetobj->delta() . "\n";

  while (my ($accid, $pedchk) = each %{$self->pedigree_checks()}) {
    my $gtsobj = $gtsetobj->accid_genotypes()->{$accid};
    $the_string .= $accid . "  " . $gtsobj->quality_counts() . "  ";
    $the_string .= $pedchk->as_string_ns() . "\n";
  }
  return $the_string;
}

sub as_string_Ns{
  my $self = shift;
  my $gtsetobj = $self->genotypes_set();
  my $the_string = '';
  $the_string .= "# n_accessions: " . $gtsetobj->n_accessions() . "\n";
  $the_string .= "# n_markers: " . $gtsetobj->n_markers() . "\n";
  $the_string .= "# delta: " . $gtsetobj->delta() . "\n";

#  while (my ($accid, $pedchk) = each %{$self->pedigree_checks()}) {
    for my $accid (@{$gtsetobj->accession_ids()}){
      my $gtsobj = $gtsetobj->accid_genotypes()->{$accid};
      my $pedchk = $self->pedigree_checks()->{$accid} // undef;
    $the_string .= $accid . "  " . $gtsobj->quality_counts() . "  ";
    $the_string .= (defined $pedchk)? $pedchk->as_string_Ns() . "\n" : "  No pedigree \n";
  }
  return $the_string;
}

sub as_string_xs{
  my $self = shift;

  my $the_string = '';
  while (my ($accid, $pedchk) = each %{$self->pedigree_checks()}) {
    $the_string .= $pedchk->as_string_xs() . "   ";
    my $mrchk = $self->matrand_checks()->{$accid};
    $the_string .= $mrchk->as_string_xs() . "   ";
    my $prchk = $self->patrand_checks()->{$accid};
    $the_string .= $prchk->as_string_xs() . "   ";
    $the_string .= $self->randrand_checks()->{$accid}->as_string_xs() . "   ";
    $the_string .= "\n";
  }
  return $the_string;
}

sub as_string_zs{
  my $self = shift;

  my $the_string = '';
  while (my ($accid, $pedchk) = each %{$self->pedigree_checks()}) {
    $the_string .= $pedchk->as_string_zs() . "   ";
    my $mrchk = $self->matrand_checks()->{$accid};
    $the_string .= $mrchk->as_string_zs() . "   ";
    my $prchk = $self->patrand_checks()->{$accid};
    $the_string .= $prchk->as_string_zs() . "   ";
    $the_string .= $self->randrand_checks()->{$accid}->as_string_zs() . "   ";
    $the_string .= "\n";
  }
  return $the_string;
}

### non methods

# sub _nxyzs{	     # given mat, pat, child gts, get n000, n001, etc.
#   my $mat_gts = shift;
#   my $pat_gts = shift;
#   my $child_gts = shift;

#   my %n_c = ('000' => 0, '001' => 0, '002' => 0,   '010' => 0, '011' => 0, '012' => 0,
# 	     '020' => 0, '021' => 0, '022' => 0,   '110' => 0, '111' => 0, '112' => 0,
# 	     '120' => 0, '121' => 0, '122' => 0,   '220' => 0, '221' => 0, '222' => 0,
# 	     'X' => 0);
 
#   my @m_gts = split("", $mat_gts);
#   my @p_gts = split("", $pat_gts);
#   my @ch_gts = split("", $child_gts);
#   die "gt string length prob. \n" if (scalar @m_gts != scalar @p_gts  or  scalar @m_gts != scalar @ch_gts);
#   #  die "gt string length prob. " . length $mat_gts . " " . length $pat_gts . " " . length $child_gts . "\n"
#   #    if (length $mat_gts != length $pat_gts  or  length $mat_gts != length $child_gts);

#   #  for (my $i = 0; $i < length $mat_gts; $i++) {
#   while (my($i, $c) = each @ch_gts) {
#     my $m = $m_gts[$i];
#     my $p = $p_gts[$i];
#     my $tr; # e.g. '001'
#     #    print STDERR "$c  $m $p  $tr \n";
#     #   $tr .= $c;
#     if (uc $c eq 'X'  or  uc $m eq 'X'  or  uc $p eq 'X') { # missing data case
#       $tr = 'X';
#     } else {			# all 3 gts present. 
#       $tr = ($m < $p)? "$m$p$c" : "$p$m$c";
#     }
#     $n_c{$tr}++;
#   }
#   my @nxyzs = ($n_c{'000'}, $n_c{'001'}, $n_c{'002'},
# 	       $n_c{'010'}, $n_c{'011'}, $n_c{'012'},
# 	       $n_c{'020'}, $n_c{'021'}, $n_c{'022'},
# 	       $n_c{'110'}, $n_c{'111'}, $n_c{'112'},
# 	       $n_c{'120'}, $n_c{'121'}, $n_c{'122'},
# 	       $n_c{'220'}, $n_c{'221'}, $n_c{'222'}, $n_c{'X'}); #, scalar @ch_gts);
#   die if(scalar @ch_gts != sum(@nxyzs));
#   return \@nxyzs;
# }

sub _get_a_pedcheck{
  my $aid_gts = shift;
  my ($matid, $patid, $accid) = @_;
  my $accgtobj = $aid_gts->{$accid} // return undef;
  my $matgtobj = $aid_gts->{$matid} // return undef;
  my $patgtobj = $aid_gts->{$patid} // return undef;
  return  PedigreeCheck->new( { mat_gtsobj => $matgtobj, pat_gtsobj => $patgtobj, acc_gtsobj => $accgtobj } );
}

1;
