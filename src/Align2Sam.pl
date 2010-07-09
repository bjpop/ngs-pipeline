#!/usr/bin/env perl

use strict;
use warnings;
use Align2Sam;

$| = 1;

die qq(
Usage:   Align2Sam.pl <ref.fa> <alignment.sai> 
) if(@ARGV < 2);

my $REFERENCE = shift @ARGV;
my $ALIGNMENT = shift @ARGV;
my @SAMFLAGS = ();

Align::align2Sam(@SAMFLAGS, $REFERENCE, $ALIGNMENT);
