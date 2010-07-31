#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Std;
use Align;
use File::Basename;
use IndexReference;
use Options;
use Illumina2Sanger;
use Logger;

$| = 1;

my %OPTIONS = Options::get_options();

# pick reference file and sequences from command line arguments vector
my ($REFERENCE, @SEQUENCES) = @ARGV;

# Open the log file and redirect stout to 
Logger::init_logger($OPTIONS{"l"});

# Step 1. Create reference database.
my @INDEX_FLAGS = ("-a $OPTIONS{i}");
IndexReference::mkRefDataBase($REFERENCE, @INDEX_FLAGS);
# Step 2. Index the reference.
IndexReference::indexReference($REFERENCE);


my $SEQ;
foreach $SEQ (@SEQUENCES)
{
   # Step 3. Convert sequence to fastq format.
   Illumina2Sanger::convertOneFile($SEQ);
}
