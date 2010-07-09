package Align2Sam;

use strict;
use warnings;
use File::Basename;

# Convert alignment to SAM format. 
sub align2Sam
{
   my ($SAMFLAG, $REFERENCE, $SEQUENCE, $ALIGNMENT) = @_;
   my ($name,$path,$suffix) = fileparse($ALIGNMENT, ".sai");
   my $SAMALIGN = $path.$name.".sam";
   my $COMM = "bwa $SAMFLAG $REFERENCE $ALIGNMENT $SEQUENCE > $SAMALIGN";


   print "Converting alignment to SAM format\n";
   print "$COMM\n";
   system($COMM);
   $SAMALIGN;
}

1;
