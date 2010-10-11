#!/usr/bin/env perl

use strict;
use warnings;
use Illumina2Sanger;
$| = 1;

die qq(
Usage:   Illumina2Sanger.pl <reads1.txt> <reads2.txt> ... <readsn.txt>
) if(@ARGV < 1);

Illumina2Sanger::convertManyFiles(@ARGV);
