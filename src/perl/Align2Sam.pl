#!/usr/bin/env perl

use strict;
use warnings;
use Align2Sam;

$| = 1;

die qq(
Usage:   Align2Sam.pl <ref.fa> <alignment.sai> 
) if(@ARGV < 2);

my ($REFERENCE, $SEQUENCE, $ALIGNMENT) = @ARGV;
my $SAMFLAG = "samse";

Align2Sam::align2Sam($SAMFLAG, $REFERENCE, $SEQUENCE, $ALIGNMENT);
