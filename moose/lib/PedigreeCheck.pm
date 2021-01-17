package PedigreeCheck;		# checking a single pedigree.
use strict;
use warnings;
use Moose;
#use Mouse;
use namespace::autoclean;
use Carp;
#use Scalar::Util qw (looks_like_number );
use List::Util qw ( min max sum );
#use POSIX qw ( floor ceil );
#use Genotypes;

# accession id and string of genotypes

# has accession_id => (
# 		     isa => 'Str',
# 		     is => 'ro',
# 		     required => 1,
# 		    );

# has maternal_id => (
# 		    isa => 'Str',
# 		    is => 'ro',
# 		    required => 1,
# 		   );

# has paternal_id => (
# 		    isa => 'Str',
# 		    is => 'ro',
# 		    required => 1,
# 		   );

has acc_gtsobj => (
		   isa => 'Object',
		   is => 'ro',
		   required => 1,
		  );
has mat_gtsobj => (
		   isa => 'Object',
		   is => 'ro',
		   required => 1,
		  );
has pat_gtsobj => (
		   isa => 'Object',
		   is => 'ro',
		   required => 1,
		  );

has am_distances => (
		     isa => 'ArrayRef', # [agmr, hgmr]
		     is => 'rw',
		     default => sub { [-1, -1] }, # agmr hgmr
		    );

has ap_distances => (
		     isa => 'ArrayRef',
		     is => 'rw',
		     default => sub { [-1, -1] },
		    );

has mp_distances => (
		     isa => 'ArrayRef',
		     is => 'rw',
		     default => sub { [-1, -1] },
		    );

has nxyzs => (
	      isa => 'ArrayRef[Int]',
	      is => 'rw',
	      required => 0,
	     );

has n00x => (
	     isa => 'ArrayRef[Int]',
	     is => 'rw',
	     default => sub{ [0, 0, 0, 0] },
	    );

has n01x => (
	     isa => 'ArrayRef[Int]',
	     is => 'rw',
	     default => sub{ [0, 0, 0, 0] },
	    );

has n02x => (
	     isa => 'ArrayRef[Int]',
	     is => 'rw',
	     default => sub{ [0, 0, 0, 0] },
	    );

has n11x => (
	     isa => 'ArrayRef[Int]',
	     is => 'rw',
	     default => sub{ [0, 0, 0, 0] },
	    );

has n12x => (
	     isa => 'ArrayRef[Int]',
	     is => 'rw',
	     default => sub{ [0, 0, 0, 0] },
	    );

has n22x => (
	     isa => 'ArrayRef[Int]',
	     is => 'rw',
	     default => sub{ [0, 0, 0, 0] },
	    );

has nX => (
	   isa => 'Int',
	   is => 'rw',
	   default => 0,
	  );

sub BUILD{
  my $self = shift;
  $self->triple_counts();
  my @nxyzs = @{$self->nxyzs()};

  my @n00xs = (@nxyzs[0..2], sum(@nxyzs[0..2]));
  $self->n00x(\@n00xs);

  my @n01xs = (@nxyzs[3..5], sum(@nxyzs[3..5]));
  $self->n01x(\@n01xs);

  my @n02xs = (@nxyzs[6..8], sum(@nxyzs[6..8]));
  $self->n02x(\@n02xs);

  my @n11xs = (@nxyzs[9..11], sum(@nxyzs[9..11]));
  $self->n11x(\@n11xs);

  my @n12xs = (@nxyzs[12..14], sum(@nxyzs[12..14]));
  $self->n12x(\@n12xs);

  my @n22xs = (@nxyzs[15..17], sum(@nxyzs[15..17]));
  $self->n22x(\@n22xs);

  $self->nX($nxyzs[18]);

  $self->distances();
}

sub as_string_ns{
  my $self = shift;
  my $the_string = sprintf("%s %s  ", $self->mat_gtsobj->id(), $self->pat_gtsobj->id());
  $the_string .= sprintf("%2i %2i %2i %2i   ", @{$self->n00x()});
  $the_string .= sprintf("%2i %2i %2i %2i   ", @{$self->n01x()});
  $the_string .= sprintf("%2i %2i %2i %2i   ", @{$self->n02x()});
  $the_string .= sprintf("%2i %2i %2i %2i   ", @{$self->n11x()});
  $the_string .= sprintf("%2i %2i %2i %2i   ", @{$self->n12x()});
  $the_string .= sprintf("%2i %2i %2i %2i   ", @{$self->n22x()});
  $the_string .= sprintf("%2i", $self->nX());

  $the_string .= sprintf("  %7.5f %7.5f", @{$self->am_distances()});
  $the_string .= sprintf("  %7.5f %7.5f", @{$self->ap_distances()});
  $the_string .= sprintf("  %7.5f %7.5f", @{$self->mp_distances()});
  return $the_string;
}

sub as_string_xs{
  my $self = shift;
  my @n_xyzs = @{$self->nxyzs()};
  my $the_string = sprintf("%s %s %s  ", $self->accession_id(), $self->maternal_id(), $self->paternal_id());
  for my $i (0..5) {
    my $j = 3*$i;
    my $ntot = sum(@n_xyzs[$j..$j+2]);

    $the_string .= ($ntot > 0)? sprintf("%5.4f %5.4f %5.4f %2i  ", map($_/$ntot, @n_xyzs[$j..$j+2]), $ntot ) : sprintf("-1 -1 -1 0  ");
  }
  return $the_string;
}

sub as_string_zs{		# lump together n001 n221, etc.
  my $self = shift;
  my $the_string = sprintf("%s %s %s  ", $self->acc_gtsobj->id(), $self->mat_gtsobj->id(), $self->pat_gtsobj->id());

  my $n00n22 = $self->n00x()->[3] + $self->n22x()->[3];
  my $z001221 = ($n00n22 > 0)? ($self->n00x->[1] + $self->n22x->[1])/$n00n22 : -1;
  my $z002220 = ($n00n22 > 0)? ($self->n00x->[2] + $self->n22x->[0])/$n00n22 : -1;
  
  my $n01n12 = $self->n01x()->[3] + $self->n12x()->[3];
  my $z012120 = ($n01n12 > 0)? ($self->n01x->[2] + $self->n12x->[0])/$n01n12 : -1;

  my $n02 = $self->n02x()->[3];
  my $z020022 = ($n02 > 0)? ($self->n02x->[0] + $self->n02x->[2])/$n02 : -1;
  $the_string .= sprintf("%5.4f %3i %5.4f %3i %5.4f %3i %5.4f %3i",
			 $z001221, $n00n22, $z002220, $n00n22, $z012120, $n01n12, $z020022, $n02);
  return $the_string;
}

sub triple_counts{   # given mat, pat, child gts, get n000, n001, etc.
  my $self = shift;
  my $mat_gts = $self->mat_gtsobj();
  my $pat_gts = $self->pat_gtsobj();
  my $child_gts = $self->acc_gtsobj();

  my %n_c = ('000' => 0, '001' => 0, '002' => 0,   '010' => 0, '011' => 0, '012' => 0,
	     '020' => 0, '021' => 0, '022' => 0,   '110' => 0, '111' => 0, '112' => 0,
	     '120' => 0, '121' => 0, '122' => 0,   '220' => 0, '221' => 0, '222' => 0,
	     'X' => 0);
 
  my @m_gts = split("", $mat_gts->genotypes());
  my @p_gts = split("", $pat_gts->genotypes());
  my @ch_gts = split("", $child_gts->genotypes());
  die "gt string length prob. \n" if (scalar @m_gts != scalar @p_gts  or  scalar @m_gts != scalar @ch_gts);
  #  die "gt string length prob. " . length $mat_gts . " " . length $pat_gts . " " . length $child_gts . "\n"
  #    if (length $mat_gts != length $pat_gts  or  length $mat_gts != length $child_gts);

  #  for (my $i = 0; $i < length $mat_gts; $i++) {
  while (my($i, $c) = each @ch_gts) {
    my $m = $m_gts[$i];
    my $p = $p_gts[$i];
    my $tr;			# e.g. '001'
    #    print STDERR "$c  $m $p  $tr \n";
    #   $tr .= $c;
    if (uc $c eq 'X'  or  uc $m eq 'X'  or  uc $p eq 'X') { # missing data case
      $tr = 'X';
    } else {			# all 3 gts present. 
      $tr = ($m < $p)? "$m$p$c" : "$p$m$c";
    }
    $n_c{$tr}++;
  }
  my @nxyzs = ($n_c{'000'}, $n_c{'001'}, $n_c{'002'},
	       $n_c{'010'}, $n_c{'011'}, $n_c{'012'},
	       $n_c{'020'}, $n_c{'021'}, $n_c{'022'},
	       $n_c{'110'}, $n_c{'111'}, $n_c{'112'},
	       $n_c{'120'}, $n_c{'121'}, $n_c{'122'},
	       $n_c{'220'}, $n_c{'221'}, $n_c{'222'}, $n_c{'X'}); #, scalar @ch_gts);
  die if(scalar @ch_gts != sum(@nxyzs));
  $self->nxyzs(\@nxyzs);
  # return \@nxyzs;
}

sub distances{
  my $self = shift;
  my $mat_gts = $self->mat_gtsobj();
  my $pat_gts = $self->pat_gtsobj();
  my $child_gts = $self->acc_gtsobj();
  $self->am_distances($child_gts->agmr_hgmr($mat_gts));
  $self->ap_distances($child_gts->agmr_hgmr($pat_gts));
  $self->mp_distances($mat_gts->agmr_hgmr($pat_gts));
}


1;
