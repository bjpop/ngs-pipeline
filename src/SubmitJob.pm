package SubmitJob;

use strict;
use warnings;
use File::Temp qw(tempfile);

sub submitJob
{
    my ($WALLTIME, $PROG, @ARGS) = @_;

    my $PBS_SCRIPT = <<SCRIPT;
#!/bin/bash
#PBS -q smp
#PBS -l mem=24gb
#PBS -l walltime=$WALLTIME
#PBS -m abe
cd \$PBS_O_WORKDIR
$PROG @ARGS
SCRIPT

   my ($FH, $FN) = tempfile();
   print $FH "$PBS_SCRIPT\n";
   my $COMM = "qsub $FN";
   print "$COMM\n";
}

1;
