#!/usr/bin/env perl

use strict;
use warnings;
use Pileup;

$| = 1;

die qq(
Usage:   Pileup.pl <dir> <reference> <alignment.bam>
) if(@ARGV < 3);

my ($DIR,$REFERENCE,$BAMALIGN) = @ARGV;

Pileup::pileup($DIR, "-c", $REFERENCE, $BAMALIGN);
