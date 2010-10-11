#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Std;
use Align;

$| = 1;

# Grab and set all options
my %OPTIONS = (a => "aln", d => "", t => 1);
getopts('a:d:t:', \%OPTIONS);

die qq(
Usage:   Align.pl [OPTIONS] <ref.fa> <reads.fastq> 

OPTIONS: -a STR	   algorithm to use for alignment [null]
         -d STR	   prefix for output directory [BWA]
         -t INT    max number of threads to use [$OPTIONS{t}]
) if(@ARGV < 2);

my $REFERENCE = shift @ARGV;
my $SEQUENCE = shift @ARGV;

Align::align($OPTIONS{a}, $OPTIONS{t}, $OPTIONS{d}, $REFERENCE, $SEQUENCE);
