#!/usr/bin/env perl

use strict;
use warnings;
use VarFilter;

$| = 1;

my %OPTIONS = (D => 1000000);
getopts('D', \%OPTIONS);

die qq(
Usage:   VarFilter.pl [OPTIONS] <pileup file> 

OPTIONS: -D INT    maximum read depth for variant calling [$OPTIONS{D}]
) if(@ARGV < 1);


# XXX fixme, -p should be a separate option
my $VARFILTERFLAGS = "-p -D $OPTIONS{D}";
my $PILEUPFILE = shift @ARGV;

VarFilter::varFilter($VARFILTERFLAGS, $PILEUPFILE);
