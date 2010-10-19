#!/usr/bin/env perl

use strict;
use warnings;
use Timer;  #qw(init_timer timer stop_timer);
use Options;
use SubmitJob qw( submitJob );

$| = 1;

Timer::init_timer('Timer.log');
Timer::timer('Start of job');

my @OLDARGV = @ARGV;
my %OPTIONS = Options::get_options();
my ($REFERENCE, @SEQUENCES) = @ARGV;
my @ARGS = @OLDARGV[0..@OLDARGV - @ARGV];

my $makeDBJobID = SubmitJob::submitJob('MakeRefDatabase', 'bwa-gcc samtools-gcc', '', '2:00:00', 'MakeReferenceDatabase.pl', @OLDARGV);
chomp($makeDBJobID);

my $SEQ_JOB_IDS = '';
my $SEQ_JOB_ID;
foreach (@SEQUENCES)
{
    $SEQ_JOB_ID = SubmitJob::submitJob('Sequence2BAM', 'samtools-gcc', ':'.$makeDBJobID, '2:00:00', 'Sequence2BAM.pl', @ARGS, $_);
    chomp($SEQ_JOB_ID);
    $SEQ_JOB_IDS .= (':' . $SEQ_JOB_ID);
}
SubmitJob::submitJob('BAMs2Variations', 'samtools-gcc', $SEQ_JOB_IDS, '2:00:00', 'BAMs2Variations.pl', @ARGS);

Timer::timer('End of job');
Timer::stop_timer();
