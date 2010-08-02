package OutputDir;

use strict;
use warnings;
use File::Basename;

# create output dirs
sub mk_output_dir 
{
   my ($SEQ, $OPTION_DIR) = @_;
   # ($name,$path,$suffix) = fileparse($SEQ, ".fastq");

   my ($name,$path,$suffix) = fileparse($SEQ, (".fastq", ".txt"));
   my $DIR = $path.$OPTION_DIR;
   unless(-d $DIR) 
   {
      print "Creating folder $DIR\n\n";
      mkdir($DIR, 0777) || print $!;
   }
   unless(-d "$DIR/Binary") 
   {
        print "Creating folder $DIR/Binary\n\n";
        mkdir("$DIR/Binary", 0777) || print $!;
   }
   $DIR."/";
}

1;
