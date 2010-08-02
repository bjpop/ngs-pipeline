package VarFilter;

use strict;
use warnings;


sub varFilter {

   my ($VARFILTERFLAGS,$PILEUPFILE,$OUT_EXTENSION) = @_;

   my $SNPSFILE = $PILEUPFILE.$OUT_EXTENSION;

   # run the variations filter
   print "Running Variations Filter\n";
   my $COMM = "samtools.pl varFilter $VARFILTERFLAGS $PILEUPFILE &> $SNPSFILE";
   print "$COMM\n";
   system($COMM);

   $SNPSFILE;
}

1;
