package VarFilter;

use strict;
use warnings;


sub varFilter {

   my ($VARFILTERFLAGS,$PILEUPFILE) = @_;

   my $SNPSFILE = $PILEUPFILE."all.snps"

   # run the variations filter
   print "Running Variations Filter\n";
   $COMM = "samtools.pl varFilter $VARFILTERFLAGS $PILEUPFILE &> $SNPSFILE;
   print "$COMM\n";
   system($COMM);

   $SNPSFILE;
}

1;
