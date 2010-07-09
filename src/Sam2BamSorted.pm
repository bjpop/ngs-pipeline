package Sam2BamSorted;

use strict;
use warnings;
use File::Basename;

# Convert SAM to BAM
sub sam2Bam 
{
   my @VREGION = ();
   my ($VIEWFLAG, $REFERENCE, $SAMFILE) = @_;
   my ($name,$path,$suffix) = fileparse($SAMFILE, ".sam");
   my $BAMALIGN = $path."Binary/".$name.".bam";
   $COMM = "samtools view $VIEWFLAG -t $REFERENCE.fai -o $BAMALIGN $SAMFILE @VREGION";

   print "Converting to BAM\n";
   print "$COMM\n";
   system($COMM);
   print "BAM Conversion Finished\n\n";
}

# Sort a BAM file
# XXX fixme
sub sortBam
{
   # sort the alignments
   print "Sorting Alignments\n";
   $SORTBAMALIGN = $BAMALIGN;
   $SORTBAMALIGN =~ s/\.bam$/.sorted/;
   push @BAMALIGNSORTED, $SORTBAMALIGN.".bam";
   my @TMP = @SORT;
   push @TMP, $BAMALIGN;
   push @TMP, $SORTBAMALIGN;
   $COMM = "samtools sort @TMP";
   print "$COMM\n";
   system($COMM);
   unlink($BAMALIGN);
   print "Alignment Sorting Finished\n\n";
}


1;
