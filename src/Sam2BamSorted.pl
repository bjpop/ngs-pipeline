#!/usr/bin/env perl

use strict;
use warnings;
use Sam2BamSorted;;

$| = 1;

die qq(
Usage:   Sam2BamSorted.pl <ref.fa> <alignment.sam> 
) if(@ARGV < 2);

my ($REFERENCE, $SAMFILE) = @ARGV;
my $VIEWFLAG = "-b";

my $BAMFILE = Sam2BamSorted::sam2Bam($VIEWFLAG, $REFERENCE, $SAMFILE);
Sam2BamSorted::sortBam($BAMFILE);
