#!/usr/bin/env perl

use strict;
use warnings;
use Options;
use Illumina2Sanger;
use Align2Sam;
use Sam2BamSorted;
use Timer;

# Process command line arguments. Remove any arguments from command line vector.
my (%OPTIONS, $REFERENCE, @SEQUENCES) = Options::get_options();

my $COUNT = 1;
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

   $COUNT += 1;
}
