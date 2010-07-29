#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Std;
use Align;
use File::Basename;
use IndexReference;
use Options;

$| = 1;

my %OPTIONS = Options::get_options();

# initialisation
my $REFERENCE = shift @ARGV;
my @ALN = ("-t $OPTIONS{t}");
my @INDEX = ("-a $OPTIONS{i}");
my @SAMSE = ();my @SAMPE = ();
my @VIEW = ("-b");
my @VREGION = ();
my @SORT = ();
my @PILEUP = ("-c");
my @VARFILTER = ("-p", "-D $OPTIONS{D}");
my $ALG = $OPTIONS{a};

# if log file specified
if(defined $OPTIONS{l}) {
   open (FH, ">$OPTIONS{l}");
   close (FH);
   # Open the log file and redirect output to it
   open (STDERR, ">>$OPTIONS{l}");
   open (STDOUT, ">>$OPTIONS{l}");
   my $now = localtime time;
   print "Log File Created $now\n";
}

IndexReference::mkRefDataBase($REFERENCE, @INDEX);
IndexReference::indexReference($REFERENCE);
