#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Std;
use Align;
use File::Basename;
use IndexReference;
use Options;
use Illumina2Sanger;

$| = 1;

my %OPTIONS = Options::get_options();

# pick reference file and sequences from command line arguments vector
my ($REFERENCE, @SEQUENCES) = @ARGV;

# variable initialisation
#my @ALN = ("-t $OPTIONS{t}");
#my @SAMSE = ();
#my @SAMPE = ();
#my @VIEW = ("-b");
#my @VREGION = ();
#my @SORT = ();
#my @PILEUP = ("-c");
#my @VARFILTER = ("-p", "-D $OPTIONS{D}");
#my $ALG = $OPTIONS{a};

# if log file specified
# XXX should push this to its own module.
if(defined $OPTIONS{l}) {
   open (FH, ">$OPTIONS{l}");
   close (FH);
   # Open the log file and redirect output to it
   open (STDERR, ">>$OPTIONS{l}");
   open (STDOUT, ">>$OPTIONS{l}");
   my $now = localtime time;
   print "Log File Created $now\n";
}

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
