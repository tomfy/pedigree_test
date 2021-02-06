#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
# use Graphics::GnuplotIF qw(GnuplotIF);
# use Math::GSL::SF  qw( :all );

use File::Basename 'dirname';
use Cwd 'abs_path';
my ( $bindir, $libdir );
BEGIN {     # this has to go in Begin block so happens at compile time
  $bindir =
    dirname( abs_path(__FILE__) ) ; # the directory containing this script
  $libdir = $bindir . '/../lib';
  $libdir = abs_path($libdir);	# collapses the bin/../lib to just lib
}
use lib $libdir;
print STDERR "# libdir: $libdir \n";

use GenotypesSet;
use Pedigrees;
use CheckPedigrees;

my $gtfilename = undef;
my $pedigree_filename = undef;
my $delta = 0.25;
my $max_bad_gt_fraction = 1.0;
my $min_hw_qual_param = 0.0;
 GetOptions(
	    'gtsfile|gtfile|genotypesfile=s' => \$gtfilename,
	    'pedigreefile=s' => \$pedigree_filename,
	    'delta=f' => \$delta,
	    'max_bad_gt_fraction=f' => \$max_bad_gt_fraction,
	    'min_hw_qual=f' => \$min_hw_qual_param,
	   );

die "No genotypes matrix filename provided.\n" if(!defined $gtfilename);

print STDERR "# Read in the gts matrix file. Create GenotypesSet object.\n";
my $gtset = GenotypesSet->new({gt_matrix_filename => $gtfilename, delta => $delta});
print STDERR "# GenotypesSet object created.\n";

#print $gtset->as_string(), "\n";
print STDERR "# ", $gtset->marker_qual_string(), "\n";
#exit;
# print $gtset->gt_set_as_string(), "\n";
if( ($max_bad_gt_fraction < 1.0) or ($min_hs_qual_param > 0.0) ){
  $gtset->remove_bad_markers($max_bad_gt_fraction, $min_hw_qual_param);
}
print STDERR "# Create pedigrees object from file.\n";
my $pedigrees = Pedigrees->new({pedigree_filename => $pedigree_filename});
print STDERR "# Pedigrees object created.\n";

print STDERR "# Create CheckPedigrees object.\n";
my $check_peds = CheckPedigrees->new({genotypes_set => $gtset, pedigrees => $pedigrees});
print STDERR "# PedigreeChecks object created.\n";
print $check_peds->as_string_Ns();
#print $check_peds->as_string_zs();
