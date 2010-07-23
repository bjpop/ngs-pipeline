#!/usr/bin/env perl

use strict;
use warnings;
use MergeBamsAndIndex;

$| = 1;

die qq(
Usage:   Sam2BamSorted.pl <dir> <alignment_1.bam> .. <alignment_n.bam>
) if(@ARGV < 2);

my $DIR = shift @ARGV;

MergeBamsAndIndex::mergeBamsAndIndex($DIR, @ARGV);
