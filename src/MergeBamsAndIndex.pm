package MergeBamsAndIndex;

use strict;
use warnings;

sub mergeBamsAndIndex {

   my $DIR = shift @_;
   my @BAMFILES = @_;
   my $COMM;
   my $BAMALIGN = $DIR."Binary/all_reads_aligned.bam";

   if($#BAMFILES > 0)
   {
      print "Merging Bam files\n";
      $COMM = "samtools merge $BAMALIGN @BAMFILES";
      print "$COMM\n";
      system($COMM);
      print "Bam files merged\n\n";
   }
   elsif ($#BAMFILES == 0)
   {
      $BAMALIGN = $BAMFILES[0];
   }
   else
   {
      die "Attempt to merge and index zero BAM files\n";
   }

   print "Indexing Alignment File\n";
   $COMM = "samtools index $BAMALIGN";
   print "$COMM\n";
   system($COMM);
   print "Indexing Finished\n\n";

   $BAMALIGN;
}

1;
