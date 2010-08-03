package Logger;

use strict;
use warnings;

sub init_logger
{
   my ($LOG_FILE) = @_;

   # if log file specified
   if(defined $LOG_FILE)
   {
      open (FH, ">$LOG_FILE");
      close (FH);
      # Open the log file and redirect output to it
      # Not sure if redirecting STDERR is a good idea?
      open (STDERR, ">>$LOG_FILE");
      open (STDOUT, ">>$LOG_FILE");
      my $now = localtime time;
      print "Log File Created $now\n";
   }

   # this will be undef if no option was set.
   $LOG_FILE;
}

1;
