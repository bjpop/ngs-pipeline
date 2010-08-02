package Align;

use strict;
use warnings;
use File::Basename;

# Align a single sequence to the reference database 
sub align
{
   my ($ALG, $THREAD_FLAGS, $DIR, $REFERENCE, $SEQUENCE) = @_;
   my ($name,$path,$suffix) = fileparse($SEQUENCE, ".fastq");
   my $ALIGNMENT = $DIR.$name.".sai";		


   printf "alg = %s, thread_flags = %s, dir = %s, ref = %s, seq = %s", $ALG, $THREAD_FLAGS, $DIR, $REFERENCE, $SEQUENCE; 

   my $COMM = "bwa $ALG -t $THREAD_FLAGS $REFERENCE $SEQUENCE > $ALIGNMENT";

   print "Running Alignment\n";
   print "$COMM\n";
   system($COMM);
   $ALIGNMENT;
}

1;
