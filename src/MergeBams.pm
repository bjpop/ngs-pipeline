package MergeBams;

use strict;
use warnings;
use File::Basename;

sub mergeBams {
   my @BAMFILES = @_;

   if($#BAMFILES > 0)
   {
      print "Merging Bam files\n";
      $BAMALIGN = $DIR."Binary/all_reads_aligned.bam";
      $COMM = "samtools merge $BAMALIGN @BAMFILES";
      print "$COMM\n";
      system($COMM);
      print "Bam files merged\n\n";


#      print "Sorting Merged Alignments\n";
#      $SORTBAMALIGN = $DIR."Binary/all_reads_aligned.sorted";
#      my @TMP = @SORT;
#      push @TMP, $BAMALIGN;
#      push @TMP, $SORTBAMALIGN;
#      $COMM = "samtools sort @TMP";
#      print "$COMM\n";
#      system($COMM);
#      unlink($BAMALIGN);
#      print "Merged Alignment Sorting Finished\n\n";
   }
}
