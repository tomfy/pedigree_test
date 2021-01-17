package GenotypesSet;
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
use Marker;

# accession id and string of genotypes

has gt_matrix_filename => (
			   isa => 'Str',
			   is => 'ro',
			   required => 1,
			  );

has delta => (
	      isa => 'Num',
	      is => 'ro',
	      default => 0.5,
	     );

has n_accessions => (
	      isa => 'Num',
	      is => 'rw',
	      default => -1,
	     );

has n_markers => (
	      isa => 'Num',
	      is => 'rw',
	      default => -1,
	      );

has marker_ids => ( isa => 'ArrayRef[Str]', 
		    is => 'rw',
		    required => 0,
		    # default =>  sub { [] },
		  );

has markerid_marker => (
			isa => 'HashRef', # values are Marker objects
			is => 'rw',
			default => sub{ {} },
			);

has accid_genotypes => (	# 
			isa => 'HashRef', # values are Genotypes objects
			is => 'rw',
			default => sub { {} },
		       );

has gt_set_as_string => (
			 isa => 'Str',
			 is => 'rw',
			 default => '',
			);

sub BUILD {
  my $self = shift;
  my $filename = $self->gt_matrix_filename();

  open my $fh, "<", "$filename" or die "Couldn't open $filename for reading.\n";

  my $line1 = <$fh>;
  my @mrkr_ids = split(" ", $line1); # yes, input is whitespace separated.
  my $M = shift @mrkr_ids;
  die if($M ne 'MARKER');
  $self->marker_ids(\@mrkr_ids);
  while(my ($im, $marker_id) = each @mrkr_ids){
    $self->markerid_marker()->{$marker_id} = Marker->new({id => $marker_id});
  }

  my $out_string_012X = $line1; # if non-integer gts in input, set this to string with all 0,1,2, X gts
  my $out_string_counts012X = '';
  my $n_markers = scalar @mrkr_ids; print STDERR "n markers: $n_markers \n";
  my %id_gts = ();
  while (my $s = <$fh>) {
    next if($s =~ /^\s*$/); # skip whitespace-only lines
    my @gts = split(" ", $s);
    my $acc_id = shift @gts;
    #   print STDERR "acc id: [$acc_id] \n";
    my $q_counts_str = '';
    if (scalar @gts > 1) {	# 
      ($s, $q_counts_str) = $self->resolve_to_012X(\@gts);
      $out_string_012X .= "$acc_id $s\n";
    } else {
      $s = $gts[0];
      $out_string_012X .= "$acc_id $s\n";
    }
    die if(length $s != $n_markers);
    $id_gts{$acc_id} = Genotypes->new({id => $acc_id, genotypes => $s, quality_counts => $q_counts_str});
  }
  $self->n_accessions(scalar keys %id_gts);
  $self->n_markers($n_markers);
  $self->accid_genotypes(\%id_gts);
  $self->gt_set_as_string($out_string_012X);
}

sub resolve_to_012X{ # take gts which are non-integers in range [0,2], and round to 0, 1, 2, (or X if not withing delta of 0, 1, or 2)
  my $self = shift;
  my $gts = shift;		# array ref
  my $delta = $self->delta();
  my ($n0, $nX, $n1, $nx, $n2) = (0, 0, 0, 0, 0);

  my $gtstr012 = '';
  while (my($i, $agt) = each @$gts) {
    my $marker_id = $self->marker_ids()->[$i];
    my $marker_obj = $self->markerid_marker()->{$marker_id};
    my $integer_gt = undef;
    if ($agt <= $delta) {
      $integer_gt = '0'; $n0++;
    } elsif ($agt >= 1-$delta  and  $agt <= 1+$delta) {
      $integer_gt = '1'; $n1++;
    } elsif ($agt >= 2-$delta) {
      $integer_gt = '2'; $n2++;
    } else {
      if ($agt < 1) {
	$integer_gt = 'X'; $nX++;
      } else {
	$integer_gt = 'x'; $nx++;
      }
    }
    $gtstr012 .= $integer_gt;
    $marker_obj->increment_counts($integer_gt);
  }
  return ($gtstr012, "$n0 $n1 $n2 $nX $nx"); # e.g. 001012010201 etc.
}

sub as_string{
  my $self = shift;
  my @marker_ids =  @{$self->marker_ids()};
  my $s = "delta: " . $self->delta() . "\n";
  $s .= "MARKER " . join(" ", @marker_ids) . "\n";
  my @count_labels = ('n0', 'n1', 'n2', 'nX', 'nx');
  for my $j (0..4) {
    $s .= $count_labels[$j] . " ";
    for my $mid (@marker_ids) {
      $s .= sprintf(" %1i", $self->markerid_marker()->{$mid}->counts()->[$j] // 0);
    }
    $s .= "\n";
  }

  while (my ($accid, $gtsobj) = each %{$self->accid_genotypes()}) {
    $s .= $gtsobj->as_string . "\n";
  }
  return $s;
}

sub marker_qual_string{
  my $self = shift;
  my @marker_ids =  @{$self->marker_ids()};
  my $s = "# delta: " . $self->delta() . "\n";
  $s .= "# marker_id  n0 n1 n2  nX nx  hw\n";
  for my $mid (@marker_ids) {
    $s .= "$mid  ";
    my $n012 = $self->markerid_marker()->{$mid}->counts()->[5] // 0;
		for my $j (0..4) {
		  my $n = $self->markerid_marker()->{$mid}->counts()->[$j] // 0;
      $s .= sprintf(" %1i %7.5f", $n, ($n012>0)? $n/$n012 : -1);
    }
    $s .= "  " . $self->markerid_marker()->{$mid}->hardy_weinberg() . "\n";
  }

  # while (my ($accid, $gtsobj) = each %{$self->accid_genotypes()}) {
  #   $s .= $gtsobj->as_string . "\n";
  # }
  return $s;
}

1;
