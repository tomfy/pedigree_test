package Cluster1d;
use strict;
use warnings;
use Math::Trig;			# so can use pi() for pi.
use List::Util qw(min max sum);

use constant PI => pi();


sub one_d_2cluster{
  my $xsar = shift;
  my $pow = shift // 'log';
  my $kernel_width = shift // undef;
  my @xs = sort {$a <=> $b} @$xsar; # sort low to high

  my $small = 1e-6;
  my $xsmall = undef;
  for my $x (@xs) {	    # find first (i.e. least) number >= $small
    if ($x >= $small) {
      $xsmall = $x;
      last;
    }
  }
  for my $x (@xs) {		#
    if ($x < $xsmall) {
      $x = $xsmall
    } else {
      last;
    }
  }

  if ($pow eq 'log') {
    my @logxs = map(log($_), @xs);
  #  my ($j_n_L, $j_h_opt, $j_lr_max) = jenks_2cluster(\@logxs);
    my ($km_n_L, $km_h_opt, $km_mom) = kmeans_2cluster(\@logxs);
    $km_h_opt = exp($km_h_opt);
    print STDERR $km_n_L-1, "  $kernel_width \n";
    my ($kde_n_L, $kde_h_opt, $min_kde_est) = kde_2cluster(\@logxs, $km_n_L-1, $kernel_width); #
    $kde_h_opt = exp($kde_h_opt);
    return ($km_n_L, $km_h_opt, $kde_n_L, $kde_h_opt);
  } else {
    my @txs = map($_**$pow, @xs);
  #  my ($j_n_L, $j_h_opt, $j_lr_max) = Cluster1d::jenks_2cluster(\@txs);
    my ($km_n_L, $km_h_opt, $km_mom) = kmeans_2cluster(\@txs);
    my ($kde_n_L, $kde_h_opt, $min_kde_est) = kde_2cluster(\@txs, $km_n_L-1, $kernel_width); #
    $km_h_opt = $km_h_opt**(1/$pow);
    $kde_h_opt = $kde_h_opt**(1/$pow);
    return ($km_n_L, $km_h_opt, $kde_n_L, $kde_h_opt);
  }
}


sub jenks_2cluster{ # divide into 2 clusters using jenks natural breaks
  my $xar = shift;  # array ref of real data values
  my @xs = sort {$a <=> $b} @$xar; # @xs is sorted small to large.

  my ($n_left_opt, $LR_max) = (1, -10000.0);
  my ($n, $sumx, $sumxsqr) = (scalar @xs, sum(@xs), 0); # sum(map($_*$_, @xs)));
  my ($n_left, $sumx_left, $sumxsqr_left) = (0, 0, 0);
  my ($n_right, $sumx_right, $sumxsqr_right) = ($n, $sumx, $sumxsqr);

  $LR_max *= scalar @xs;
  for my $x (@xs[0 .. $#xs-1]) {
    $n_left++; $n_right--;
    $sumx_left += $x; $sumx_right -= $x;
    #    $sumxsqr_left += $x*$x; $sumxsqr_right -= $x*$x;
    my $L = ($sumx_left)**2/$n_left;   # - $sumxsqr_left;
    my $R = ($sumx_right)**2/$n_right; # - $sumxsqr_right;
    my $LR = ($L + $R);
    if ($LR > $LR_max) {
      $LR_max = $LR;
      $n_left_opt = $n_left;
    }
  }
  my $h_opt = 0.5*($xs[$n_left_opt-1] + $xs[$n_left_opt]);
  return ($n_left_opt, $h_opt, $LR_max);
}


sub kmeans_2cluster{ # divide into 2 clusters by finding dividing value h s.t. h = 1/2 * (<x>_<h + <x>_>h)
  # however instead of guessing some initial cluster centers and iteratively refining, just 
  # consider break points with 1, 2, 3, etc. pts in the L-hand cluster,
  # until mean of means (i.e. mean of the mean of L cluster, and mean of R cluster)
  # lies between the two clusters.
  my $xar = shift;		   # array ref of real data values
  my @xs = sort {$a <=> $b} @$xar; # @xs is sorted small to large.

  my $h_opt;
  my ($n, $sumx, $sumxsqr) = (scalar @xs, sum(@xs), 0); # sum(map($_*$_, @xs)));
  my ($n_left, $sumx_left, $sumxsqr_left) = (0, 0, 0);
  my ($n_right, $sumx_right, $sumxsqr_right) = ($n, $sumx, $sumxsqr);

  my $mean_of_means;
  for my $x (@xs[0 .. $#xs-1]) {
    $n_left++; $n_right--;
    $sumx_left += $x; $sumx_right -= $x;
    $mean_of_means = 0.5*($sumx_left/$n_left + $sumx_right/$n_right);
    if ($mean_of_means < $xs[$n_left]  and  $mean_of_means >= $x) { # this is the place
      $h_opt = 0.5*($x + $xs[$n_left]);
      last;
    }
  }
#    print STDERR "kmeans  $n_left ", $xs[$n_left-1], "  ", $xs[$n_left], "  $mean_of_means \n"; # $p $LRp $LRp_max \n";
  return ($n_left, $h_opt, $mean_of_means);
}


sub kde_2cluster{
  # now refine using kernel density estimation.
  my $xar = shift;		# must be already sorted, low to high.
  my $i_opt = shift; # look for min of kde in neighborhood of $xar->[$i_opt]
  my $kernel_width = shift // undef;
  my @xs = @$xar;
  my $n = scalar @xs;
  my $n_left = $i_opt+1;
  my $n_right = $n - $n_left;
  if (!defined $kernel_width) { # consider 20% of pts near initial guess $i_opt
    my ($istart, $iend) = (max(0, $i_opt-int($n/10)), min($i_opt+int($n/10), $n-1));
    my ($max_delta, $n_pts) = (-1, 10); # find max difference between $xs[$i] and $xs[$i+n_pts-1]
    print STDERR "$istart $iend $n_pts \n";
    for (my $ii = $istart; ; $ii++) {
      my $jj = $ii + $n_pts - 1;
      last if($jj > $iend);
      my $delta = $xs[$jj] - $xs[$ii];
      $max_delta = $delta if ($delta > $max_delta);
    }
    $kernel_width = $max_delta/sqrt(2.0);
    print STDERR "### max_delta:  $max_delta \n";
  }
  print STDERR "# kernel width: $kernel_width \n";

  my $min_kde_est = 1e100;
  my $kde_x_est;
  my $kde_i_opt;
  for (my $j = max(0, $i_opt-int($n_left/6)); $j < min($i_opt+int($n_right/3), scalar @xs-1); $j++) {

    my $x =  0.5*($xs[$j] + $xs[$j+1]);
    my $kde_est = kde(\@xs, $x, $kernel_width, $j); # kde est of density at x

    if ($kde_est < $min_kde_est) {
      $min_kde_est = $kde_est;
      $kde_x_est = $x;
      $kde_i_opt = $j;
    }
    print STDERR "$j  $x   $kde_est \n";
  }
  my $kde_n_left = $kde_i_opt + 1;
  return ($kde_n_left, $kde_x_est, $min_kde_est);
}

sub kde{
  my $xs = shift;	     # set of sorted reals (low-to-hi)
  my $x = shift;	     # evaluate the kde at this x
  my $w = shift;	     # half-width of kernel
  my $i = shift // undef;    # largest index such that $xs->[$i] <= $x

  my $kde_sum = 0;

  for (my $j=$i; $j >= 0; $j--) { # go to smaller and smaller values of x until further than $w away
    my $arg = abs( $xs->[$j] - $x );
    last if ($arg >= $w);
    $kde_sum += kernel($w, $arg);
  }
  for (my $j=$i+1; $j < scalar @$xs; $j++) { # go to larger and larger values of x until further than $w away
    my $arg = abs( $xs->[$j] - $x );
    last if ($arg >= $w);
    $kde_sum += kernel($w, $arg);
  }
  return $kde_sum;
}

sub kernel{			# 2 at x=0, 0 at |x| >= w
  my $w = shift;
  my $x = shift;
  return (abs($x) >= $w)? 0 : (cos(PI*$x/$w) + 1)
}


1;
