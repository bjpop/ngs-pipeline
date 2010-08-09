#!/usr/bin/env perl

use strict;
use warnings;
use Timer qw(init_timer timer stop_timer);
use Options;
use SubmitJob qw( submitJob );

$| = 1;

Timer::init_timer('Timer.log');
Timer::timer('Start of job');

my @OLDARGV = @ARGV;
my %OPTIONS = Options::get_options();
my ($REFERENCE, @SEQUENCES) = @ARGV;
my @ARGS = @OLDARGV[0..@OLDARGV - @ARGV];

my $makeDBJobID = SubmitJob::submitJob('1:00:00', 'MakeReferenceDatabase.pl', @ARGS);
foreach (@SEQUENCES)
{
    SubmitJob::submitJob('1:00:00', 'Sequence2BAM.pl', @ARGS, $_);
}
SubmitJob::submitJob('1:00:00', 'BAMs2Variations.pl', @ARGS);

Timer::timer('End of job');
Timer::stop_timer();
