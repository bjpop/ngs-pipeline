package Pileup;

use strict;
use warnings;


sub pileup {

   my ($DIR,$PILEUP_FLAGS,$REFERENCE,$BAMALIGN) = @_;
   my $PILEUPFILE = $DIR."reads_aligned.pileup";

   # construct a pileup file
   print "Constructing Pileup File\n";
   my $COMM = "samtools pileup $PILEUP_FLAGS -f $REFERENCE $BAMALIGN > $PILEUPFILE";
   print "$COMM\n";
   system($COMM);
   print "Pileup Finished\n\n";

   $PILEUPFILE;
}

1;
