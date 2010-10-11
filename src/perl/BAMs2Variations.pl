#!/usr/bin/env perl

use strict;
use warnings;
use Options;
use MergeBamsAndIndex;
use Pileup;
use VarFilter;
use Timer;

# Process command line arguments. Remove any arguments from command line vector.
my (%OPTIONS, $REFERENCE, @SORTED_BAMS) = Options::get_options();

# Steps 8-10. Merge BAM files, and align the result.
my $BAM_ALIGN = MergeBamsAndIndex::mergeBamsAndIndex($DIR, @SORTED_BAMS); 

Timer::timer('BAMs merged and indexed');

# Step 11. Construct pileup file. 
my $PILEUP = Pileup::pileup($DIR, "-c", $REFERENCE, $BAM_ALIGN);

Timer::timer('Pileup file created');

# Step 12-13. Run variation filter.
VarFilter::varFilter("-p -D $OPTIONS{D}", $PILEUP, ".all.snps");

Timer::timer('First VarFilter run');

VarFilter::varFilter("-D $OPTIONS{D}", $PILEUP, ".snps");

Timer::timer('Second VarFilter run');
