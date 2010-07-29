package IndexReference;

use strict;
use warnings;

# create reference database
sub mkRefDataBase
{
   my ($REFERENCE, @INDEX) = @_;
   my $REF_DATABASE = $REFERENCE.".bwt";

   if(-e $REF_DATABASE) {
      print "Reference database already exists: using $REF_DATABASE\n";
   } 
   else {
      print "Creating Reference Database\n";
      my $COMM = "bwa index @INDEX $REFERENCE";
      print "$COMM\n";
      system($COMM);
   }
}

# index the reference
sub indexReference
{
   my $REFERENCE = $_[0];
   my $REF_INDEX = $REFERENCE.".fai";

   if(-e $REF_INDEX) {
      print "Reflist already exists: using $REF_INDEX\n";
   }
   else {
      print "Creating Reflist\n";
      my $COMM = "samtools faidx $REFERENCE";
      print "$COMM\n";
      system($COMM);
   }
}

1;
