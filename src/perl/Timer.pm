package Timer;

use strict;
use warnings;

my $LAST_TIME;

sub init_timer 
{
   my ($TIMER_FILE) = @_;
   open (TIME, ">$TIMER_FILE");
   $LAST_TIME = time;
}

sub timer
{
   my ($MESSAGE) = @_;
   my $NOW = time;
   my $DELTA = $NOW - $LAST_TIME;
   $LAST_TIME = $NOW;
   printf TIME "%s: %s\n", $DELTA, $MESSAGE;
}

sub stop_timer
{
   close TIME;
}

1;
