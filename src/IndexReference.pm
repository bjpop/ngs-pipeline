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
      my $COMM;
      print "Creating Reference Database\n";
      $COMM = "bwa index @INDEX $REFERENCE";
      print "$COMM\n";
      system($COMM);
   }
}

# index the reference
sub indexReference
{
   my $REFERENCE = @_[0];
   my $REF_INDEX = $REFERENCE.".fai";

   if(-e $REF_INDEX) {
      print "Reflist already exists: using $REF_INDEX\n";
   }
   else {
      my $COMM;
      print "Creating Reflist\n";
      $COMM = "samtools faidx $REFERENCE";
      print "$COMM\n";
      system($COMM);
   }
}

1;
