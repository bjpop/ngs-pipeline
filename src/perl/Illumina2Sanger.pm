package Illumina2Sanger;

use strict;
use warnings;

# convert sequences to fastq format
sub convertManyFiles
{
   my @SEQFASTQ;
   my $SEQ;
   print "Converting Sequences From Illumina FASTQ to Standard/Sanger FASTQ Format\n";
   foreach $SEQ (@_) 
   {
       push @SEQFASTQ, convertOneFile($SEQ);
   }
   @SEQFASTQ;
}

sub convertOneFile
{
   my $ILL_FILE = $_[0];
   my $SANGER_FILE = $ILL_FILE;

   if(grep(/fastq$/, $ILL_FILE)) 
   {
      print "$ILL_FILE Already In Standard/Sanger FASTQ Format\n";
   } 
   elsif(grep(/txt$/, $ILL_FILE))
   {
      my $COMM = "maq ill2sanger $ILL_FILE $SANGER_FILE\n";
      print "$COMM";
      system($COMM);
      print "$ILL_FILE Converted\n";

      $SANGER_FILE=~s/\.txt$/\.fastq/ ;
   }
   else
   {
      die "sequence file not in illumina (.txt) or sanger (.fastq) format\n";
   }
   $SANGER_FILE;
}

# convert paired end sequences to fastq format
# my @PEFASTQ;
#if(defined $OPTIONS{p}) {
#	print "Converting Paired End Sequences to Fastq format\n";
#	foreach $PE (@PAIREDEND) {
#		if(!grep(/fastq$/, $PE)) {
#			$PE = convert($PE);
#			print "$PE Converted\n";
#		} else {
#			print "$PE Already In Standard/Sanger FASTQ Format\n";
#		}
#		push @PEFASTQ, $PE;
#	}
#	print "\n";
#}	

1;
