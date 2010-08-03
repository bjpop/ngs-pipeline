#!/usr/bin/env perl

use strict;
use warnings;
use Options;
use Logger;
use Align;
use IndexReference;
use Illumina2Sanger;
use OutputDir;
use Align2Sam;
use Sam2BamSorted;
use MergeBamsAndIndex;
use Pileup;
use VarFilter;
use Timer;

$| = 1;

Timer::init_timer('Timer.log');

Timer::timer('Start of job');

# Process command line arguments. Remove any arguments from command line vector.
my %OPTIONS = Options::get_options();

# Pick reference file and sequences from command line arguments vector.
my ($REFERENCE, @SEQUENCES) = @ARGV;

# Open the log file and redirect output to it
Logger::init_logger($OPTIONS{l});

# Step 1. Create reference database.
IndexReference::mkRefDataBase($REFERENCE, "-a $OPTIONS{i}");

Timer::timer('Reference database created');

# Step 2. Index the reference.
IndexReference::indexReference($REFERENCE);

Timer::timer('Reference database indexed');

# Create the output directory
my $DIR = OutputDir::mk_output_dir($SEQUENCES[0], $OPTIONS{d});

my $COUNT = 1;
my @SORTED_BAMS = ();
foreach (@SEQUENCES)
{
   # Step 3. Convert sequence to fastq format.
   my $SEQ_FASTQ = Illumina2Sanger::convertOneFile($_);

   Timer::timer('Sequence '.$COUNT.' converted to fastq format');

   # Step 4. Align sequence to the reference database.
   my $SEQ_ALIGN = Align::align($OPTIONS{a}, $OPTIONS{t}, $DIR, $REFERENCE, $SEQ_FASTQ);

   Timer::timer('Sequence '.$COUNT.' aligned to database');

   # Step 5. Convert alignment to SAM format.
   my $SEQ_SAM = Align2Sam::align2Sam("samse", $REFERENCE, $_, $SEQ_ALIGN);

   Timer::timer('Alignment for sequence '.$COUNT.' converted to SAM format');

   # Step 6. Convert SAM to BAM. 
   my $SEQ_BAM = Sam2BamSorted::sam2Bam("-b", $REFERENCE, $SEQ_SAM);

   Timer::timer('SAM file for sequence '.$COUNT.' converted to BAM format');

   # Step 7. Sort the BAM file.
   my $SEQ_BAM_SORTED = Sam2BamSorted::sortBam($SEQ_BAM);

   Timer::timer('BAM file for sequence '.$COUNT.' sorted');

   push @SORTED_BAMS, $SEQ_BAM_SORTED;
   $COUNT += 1;
}

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

Timer::stop_timer();
