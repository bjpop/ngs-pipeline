package MergeBams;

use strict;
use warnings;

sub mergeBams {

   my $DIR = shift @_;
   my @BAMFILES = @_;

   if($#BAMFILES > 0)
   {
      my $BAMALIGN = $DIR."Binary/all_reads_aligned.bam";
      print "Merging Bam files\n";
      my $COMM = "samtools merge $BAMALIGN @BAMFILES";
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

1;
