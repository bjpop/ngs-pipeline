package Sam2BamSorted;

use strict;
use warnings;
use File::Basename;

# Convert SAM to BAM
# Assumes the Binary directory has been created already.
sub sam2Bam 
{
   # XXX fix the VREGION argument
   my @VREGION = ();
   my ($VIEWFLAG, $REFERENCE, $SAMFILE) = @_;
   my ($name,$path,$suffix) = fileparse($SAMFILE, ".sam");
   my $BAMALIGN = $path."Binary/".$name.".bam";
   my $COMM = "samtools view $VIEWFLAG -t $REFERENCE.fai -o $BAMALIGN $SAMFILE @VREGION";

   print "Converting to BAM\n";
   print "$COMM\n";
   system($COMM);
   print "BAM Conversion Finished\n\n";

   $BAMALIGN;
}

# Sort a BAM file
# XXX fixme
sub sortBam
{
   # XXX fix the SORT argument
   my @SORT = ();
   my ($BAMALIGN) = @_;
   my $SORTBAMALIGN = $BAMALIGN;
   $SORTBAMALIGN =~ s/\.bam$/.sorted/;

   print "Sorting Bam alignments\n";
   my $COMM = "samtools sort @SORT $BAMALIGN $SORTBAMALIGN";
   print "$COMM\n";
   system($COMM);

   # I think we should leave the decision to remove the BAM file to the
   # caller of this function.
   #unlink($BAMALIGN);

   print "BAM alignment sorting finished\n\n";

   $SORTBAMALIGN.".bam";
}

1;
