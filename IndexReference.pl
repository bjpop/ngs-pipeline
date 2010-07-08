#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Std;
use IndexReference;
$| = 1;

my %OPTIONS = (i => "bwtsw");
getopts('i', \%OPTIONS);

die qq(
Usage:   IndexReference.pl [OPTIONS] <ref.fa|bfa>

OPTIONS: -i STR    Algorithm for constructing BWT index [$OPTIONS{i}]
) if(@ARGV < 1);

my $REFERENCE = shift @ARGV;
my @INDEX = ("-a $OPTIONS{i}"); 

IndexReference::mkRefDataBase($REFERENCE, @INDEX);
IndexReference::indexReference($REFERENCE);
