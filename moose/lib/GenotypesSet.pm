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
# use Marker;

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

has accession_ids => (isa => 'ArrayRef[Str]',
		      is => 'rw',
		      required => 0,
		     );

has marker_ids => ( isa => 'ArrayRef[Str]', 
		    is => 'rw',
		    required => 0,
		    # default =>  sub { [] },
		  );

has marker_gt_counts  => ( isa => 'ArrayRef[ArrayRef[Int]]', 
			   is => 'rw',
			   required => 0,
			   # default =>  sub { [] },
			 );

# has markerid_marker => (
# 			isa => 'HashRef', # values are Marker objects
# 			is => 'rw',
# 			default => sub{ {} },
# 		       );

has accid_genotypes => (		  # 
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
  my @acc_ids = ();
  my @mrkr_ids = split(" ", $line1); # yes, input is whitespace separated.
  my $M = shift @mrkr_ids;
  die if($M ne 'MARKER');
  $self->marker_ids(\@mrkr_ids);
  my @mrkr_gtcounts = ();
  for (@mrkr_ids) {
    push @mrkr_gtcounts, [0, 0, 0, 0, 0];
  }
  $self->marker_gt_counts(\@mrkr_gtcounts);

  my $out_string_012Xx = $line1; # if non-integer gts in input, set this to string with all 0,1,2, X gts
  my $out_string_counts012X = '';
  my $n_markers = scalar @mrkr_ids; print STDERR "# n markers: $n_markers \n";
  my %id_gts = ();
  my $accessions_read = 0;
  while (my $s = <$fh>) {
    next if($s =~ /^\s*$/);	# skip whitespace-only lines
    my @gts = split(" ", $s);
    my $acc_id = shift @gts;
    push @acc_ids, $acc_id;
    #   print STDERR "acc id: [$acc_id] \n";
    next if (scalar @gts == 0);
    my ($s01234, $q_counts_str) = $self->resolve_to_01234(\@gts);
    my $s012Xx = $s01234;
    $s012Xx =~ s/3/X/g;
    $s012Xx =~ s/4/x/g;
    $out_string_012Xx .= "$acc_id $s012Xx\n"; 

    die if(length $s012Xx != $n_markers);
    $id_gts{$acc_id} = Genotypes->new({id => $acc_id, genotypes01234 => $s01234, genotypes => $s012Xx, quality_counts => $q_counts_str});
    $accessions_read++;
    print STDERR "# accessions read: $accessions_read \n" if(($accessions_read % 100) == 0);
  }
  $self->n_accessions(scalar keys %id_gts);
  $self->accession_ids(\@acc_ids);
  $self->n_markers($n_markers);
  $self->accid_genotypes(\%id_gts);
  $self->gt_set_as_string($out_string_012Xx);
}

sub resolve_to_01234{ # take gts which are non-integers in range [0,2], and round to 0, 1, 2, (or X if not withing delta of 0, 1, or 2)
  my $self = shift;
  my $gts = shift;		# array ref
  my $delta = $self->delta();

  my @gts01234 = ();
  for my $agt (@$gts) {
    if ($agt <= $delta) {	# round to 0
      push @gts01234, 0;
    } elsif (abs($agt - 1) <= $delta) { 
      push @gts01234, 1;
    } elsif ($agt >= 2-$delta) { # round to 2
      push @gts01234, 2;
    } else {
      if ($agt < 1) {
	push @gts01234, 3;
      } else {
	push @gts01234, 4;
      }
    }
  }

  $self->increment_marker_gt_counts(\@gts01234);

  my $gtstr01234 = join('', @gts01234);
  #   $self->increment_marker_gt_counts_alt($gtstr01234);
  my $n0 = ($gtstr01234 =~ tr/0//);
  my $n1 = ($gtstr01234 =~ tr/1//);
  my $n2 = ($gtstr01234 =~ tr/2//);
  my $nX = ($gtstr01234 =~ tr/3//);
  my $nx = ($gtstr01234 =~ tr/4//);
  return ($gtstr01234, "$n0 $n1 $n2 $nX $nx"); # e.g. 001012010201 etc.
}

sub increment_marker_gt_counts_alt{
  my $self = shift;
  my $igt_string = shift;
  my $mrkr_gt_counts = $self->marker_gt_counts();
  for my $i (0 .. (length $igt_string) -1) {
    my $igt = substr($igt_string, $i, 1);
    $mrkr_gt_counts->[$i]->[$igt]++;
  }
}

sub increment_marker_gt_counts{
  my $self = shift;
  my $igts = shift;		# array ref w 01234
  my $mrkr_gt_counts = $self->marker_gt_counts();
  while (my($i, $igt) = each @$igts) {
    $mrkr_gt_counts->[$i]->[$igt]++;
  }
}

sub as_string{
  my $self = shift;
  my @marker_ids =  @{$self->marker_ids()};
  my $s = "delta: " . $self->delta() . "\n";
  $s .= "MARKER " . join(" ", @marker_ids) . "\n";
  my @count_labels = ('n0', 'n1', 'n2', 'nX', 'nx');
  for my $j (0..4) {
    $s .= $count_labels[$j] . " ";
    #for my $mid (@marker_ids) {
    while (my($i, $mid) = each @marker_ids) {
      #    $s .= sprintf(" %1i", $self->markerid_marker()->{$mid}->counts()->[$j] // 0);
      $s .= sprintf(" %1i", $self->marker_gt_counts()->[$i]->[$j] // 0);
    }
    $s .= "\n";
  }

  #while (my ($accid, $gtsobj) = each %{$self->accid_genotypes()}) {
  for my $acc_id (@{$self->accession_ids()}) {
    my $gtsobj = $self->accid_genotypes()->{$acc_id};
    $s .= $gtsobj->as_string . "\n";
  }
  return $s;
}

sub marker_qual_string{
  my $self = shift;
  my @marker_ids =  @{$self->marker_ids()};
  my $s = "# delta: " . $self->delta() . "\n";
  $s .= "# marker_id  n0 n1 n2  nX nx  hw\n";
  while (my($i, $mid) = each @marker_ids) {
    $s .= "$mid  ";
    my @counts =  @{$self->marker_gt_counts()->[$i]};
    my $n012 = sum(@counts[0..2]);
    die if($self->n_accessions() != sum(@counts));
    for my $j (0..4) {
      my $n = $counts[$j] // 0;
      my $denom = ($j <= 2)? $n012 : sum(@counts);
      $s .= sprintf(" %1i %7.5f", $n, ($denom>0)? $n/$denom : -1);
    }
    $s .= "  " . hardy_weinberg(\@counts) . "\n";
  }
  return $s;
}

sub remove_bad_markers{
  my $self = shift;
  my $max_bad_gt_fraction = shift;
  my $min_hw_qual_param = shift;

  my @good_marker_ids = ();
  my @mrkr_ids = @{$self->marker_ids()};
  my @marker_thumbsupdown = (1) x scalar @mrkr_ids; # 1: good, 0: bad
  my $mrkr_gt_counts = $self->marker_gt_counts();
  while (my ($i, $markerid) = each @mrkr_ids) {
    my $marker_counts = $self->marker_gt_counts()->[$i];
    my $n_markers = $self->n_markers();
    my $n_accessions = $self->n_accessions();
    my $hw_qual_param = hardy_weinberg($marker_counts);
    my $bad_gt_fraction = ($marker_counts->[3] + $marker_counts->[4])/$n_accessions // 1;
    if ( ($bad_gt_fraction > $max_bad_gt_fraction) or ($hw_qual_param < $min_hw_qual_param) ) { # bad
      $marker_thumbsupdown[$i] = 0; # mark this marker as bad
    } else {			    # OK
      push @good_marker_ids, $markerid;
    }
  }
  $self->n_markers(scalar @good_marker_ids);
  $self->marker_ids(\@good_marker_ids);

  while ( my($accid, $gtsobj) = each %{$self->accid_genotypes()}) {
    $gtsobj->remove_bad_markers(\@marker_thumbsupdown);
  }
}

######### non-methods ##################

sub hardy_weinberg{ # if hardy-weinberg eq. frequencies, return value close to 1.
  my $counts = shift;		#
  my $n0 = $counts->[0] // 0;
  my $n1 = $counts->[1] // 0;
  my $n2 = $counts->[2] // 0;
  my $ntotal = $n0 + $n1 + $n2;
  my $hw_one = ($ntotal > 0)? ($n0**0.5 + $n2**0.5)/$ntotal**0.5 : -1;
  return $hw_one;
}

1;
