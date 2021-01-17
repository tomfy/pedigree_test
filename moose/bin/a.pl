#!/usr/bin/perl -w
use strict;

my $line1 = <>;
my ($a, $n_acc) = split(" ", $line1);
die if($a ne '#n_accessions:');

my $line2 = <>;
my ($b, $n_markers) = split(" ", $line2);
die if($b ne '#n_markers:');

my $line3 = <>;
my ($c, $delta) = split(" ", $line3);
die if($c ne '#delta:');

while (<>) {
  my @cols = split(" ", $_);
  my ($accid, $matid, $patid) = @cols[0,6,7];
  my ($n0, $n1, $n2, $nX, $nx) = @cols[1..5];
  my $z001221 = z($cols[9] + $cols[29], $cols[11] + $cols[31]);
  my $z002220 = z($cols[10] + $cols[28], $cols[11] + $cols[31]);
  my $z012120 = z($cols[14] + $cols[24], $cols[15] + $cols[27]);
  my $z020022 = z($cols[16] + $cols[18], $cols[19]);
  my $z110 = z($cols[20], $cols[23]);
    my $z111 = z($cols[21], $cols[23]);
    my $z112 = z($cols[22], $cols[23]);
  # print "$accid $matid $patid  ", $nX+$nx, "  $z001221 $z002220 $z012120 $z020022 \n";
  chomp;
  print $_, "  ", $nX+$nx, "  $z001221 $z002220 $z012120 $z020022   $z110 $z111 $z112\n";
}

# 8   00 1,2
# 12  01 2
# 16  02 0,2
# 20  11
# 24  12 0
# 28  22 0,1

sub z{
  my $n = shift;
  my $d = shift;
  return ($d > 0)? $n/$d : -1;
}
