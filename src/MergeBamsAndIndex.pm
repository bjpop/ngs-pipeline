package MergeBamsAndIndex;

use strict;
use warnings;

sub mergeBamsAndIndex {

   my ($DIR, @BAMFILES) = @_;
   my $BAMALIGN = $DIR."Binary/all_reads_aligned.bam";
   my $num_bams = @BAMFILES;

   # If there is more than one BAM file then merge them.
   if($num_bams > 1)
   {
      print "Merging Bam files\n";
      my $COMM = "samtools merge $BAMALIGN @BAMFILES";
      print "$COMM\n";
      system($COMM);
      print "Bam files merged\n\n";
   }
   elsif ($num_bams == 1)
   {
      $BAMALIGN = $BAMFILES[0];
   }
   else
   {
      die "Attempt to merge and index zero BAM files\n";
   }

   print "Indexing Alignment File\n";
   my $COMM = "samtools index $BAMALIGN";
   print "$COMM\n";
   system($COMM);
   print "Indexing Finished\n\n";

   $BAMALIGN;
}

1;
