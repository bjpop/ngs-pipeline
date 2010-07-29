package Options;

use strict;
use warnings;
use Getopt::Std;
use File::Basename;

sub get_options {

   my %OPTIONS = (a => "aln", d => "", D => 1000000, i => "bwtsw", t => 1);
   # XXX what about option D?
   getopts('a:b:c:d:i:l:pt:', \%OPTIONS);

   die qq(
Usage:   Pipeline.pl [OPTIONS] <ref.fa|bfa> <reads1.txt|fastq> [reads2.txt|fastq]

OPTIONS: -a STR    algorithm to use for alignment [null]
         -b STR    begin from here [null]
         -c STR    filename of a config file [null]
         -d STR    prefix for output directory [BWA]
         -D INT    maximum read depth for variant calling [$OPTIONS{D}]
         -i STR    Algorithm for constructing BWT index [$OPTIONS{i}]
         -l STR    filename of a log file [null]
         -p        mate-pair alignment (PE mode)
         -t INT    max number of threads to use [$OPTIONS{t}]
   ) if(@ARGV < 2);

   # determine output dir name
   if($OPTIONS{d} ne "") 
   {
      $OPTIONS{d} = "BWA_".$OPTIONS{d};
   } else 
   {
      $OPTIONS{d} = "BWA";
   }

   $OPTIONS
}
