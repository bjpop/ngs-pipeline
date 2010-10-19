package SubmitJob;

use strict;
use warnings;
use File::Temp qw(tempfile);

# XXX remove bogus export of PERL5LIB
sub submitJob
{
    my ($JOB_NAME, $MODULE_NAMES, $JOBIDS, $WALLTIME, $PROG, @ARGS) = @_;

    my $DEPENDENCIES = ($JOBIDS eq '') ? '' : "#PBS -W depend=afterok$JOBIDS";
    my $MODULES = ($MODULE_NAMES eq '') ? '' : "module load $MODULE_NAMES";

    my $PBS_SCRIPT = <<SCRIPT;
#!/bin/bash
#PBS -q smp
#PBS -l mem=10gb
#PBS -l walltime=$WALLTIME
#PBS -m abe
#PBS -N $JOB_NAME
$DEPENDENCIES
cd \$PBS_O_WORKDIR
$MODULES
export PERL5LIB=/vlsci/VLSCI/bjpop/bin
$PROG @ARGS
SCRIPT

   print $PBS_SCRIPT;

   my ($FH, $FN) = tempfile();
   print $FH "$PBS_SCRIPT\n";
   my $COMM = "qsub $FN";
   my $RESULT = `$COMM`;
   #print $RESULT;
   $RESULT;
}

1;
