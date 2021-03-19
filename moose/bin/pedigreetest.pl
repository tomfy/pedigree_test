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
my $pedigree_table_filename = undef;
my $delta = 0.25; # real-number genotypes are rounded to nearest integer (0,1,2) if within +- $delta
my $max_bad_gt_fraction = 1.0;
my $min_hw_qual_param = 0.0;
my $n_random_parents = 0;
my $base_output_filename = 'out';
my $output_pedigrees = 0;
my $output_genotype_matrix = 0;

 GetOptions(
	    'gtsfile|gtfile|genotypesfile=s' => \$gtfilename,
	    'pedigreefile|pedtable=s' => \$pedigree_table_filename,
	    'delta=f' => \$delta,
	    'max_bad_gt_fraction=f' => \$max_bad_gt_fraction,
	    'min_hw_qual=f' => \$min_hw_qual_param,
	    'n_random_parents=i' => \$n_random_parents,
	    'out|basename=s' => \$base_output_filename,
	    'pedout!' => \$output_pedigrees,
	    'matrixout!' => \$output_genotype_matrix,
	   );

die "No genotypes matrix filename provided.\n" if(!defined $gtfilename);

# Read in the pedigree table:

print STDERR "# Creating pedigrees object from file: $pedigree_table_filename\n";
my $pedigrees = Pedigrees->new({pedigree_filename => $pedigree_table_filename});
print STDERR "# Pedigrees object created.\n";
if($output_pedigrees){
  my $pedigree_output_filename = $base_output_filename . '_pedigrees';
  open my $fhout, ">", "$pedigree_output_filename";
  print $fhout $pedigrees->as_string();
  close $fhout;
}

# Read in the genotype matrix file 

print STDERR "# Reading in the gts matrix file. Create GenotypesSet object.\n";
my $gtset = GenotypesSet->new({gt_matrix_filename => $gtfilename, delta => $delta});
print STDERR "# GenotypesSet object created.\n";

if( ($max_bad_gt_fraction < 1.0) or ($min_hw_qual_param > 0.0) ){
  print STDERR "# removing bad markers. max bad fraction: $max_bad_gt_fraction   min hw qual: $min_hw_qual_param \n";
  $gtset->clean_marker_set($max_bad_gt_fraction, $min_hw_qual_param);
}

if($output_genotype_matrix){
  my $gtmatrix_output_filename = $base_output_filename . '_genotypeset';
  open my $fhout, ">", "$gtmatrix_output_filename";
  print $fhout  $gtset->as_string();
  close $fhout;
}

print STDERR "# Create CheckPedigrees object.\n";
my $summary_output_filename = $base_output_filename . "_summary";
my $check_peds = CheckPedigrees->new({
				      genotypes_set => $gtset,
				      pedigrees => $pedigrees,
				      n_random_parents => $n_random_parents,
				      summary_filename => $summary_output_filename,
				     });
print STDERR "# PedigreeChecks object created.\n";


