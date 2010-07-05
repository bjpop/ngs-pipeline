use strict;
use warnings;
use Getopt::Std;
use threads;
use threads::shared;
use File::Basename;
$| = 1;

# Grab and set all options
my %OPTIONS = (a => "aln", d => "", D => 1000000, i => "bwtsw", t => 1);
getopts('a:b:c:d:i:l:pt:', \%OPTIONS);

# If the analysis dir or bed file have not been specified die
# otherwise store the results
die qq(
Usage:   RunBWA.pl [OPTIONS] <ref.fa|bfa> <reads1.txt|fastq> [reads2.txt|fastq]

OPTIONS: -a STR	   algorithm to use for alignment [null]
         -b STR	   begin from here [null]
         -c STR	   filename of a config file [null]
         -d STR	   prefix for output directory [BWA]
         -D INT    maximum read depth for variant calling [$OPTIONS{D}]
         -i STR    Algorithm for constructing BWT index [$OPTIONS{i}]
         -l STR    filename of a log file [null]
         -p 	   mate-pair alignment (PE mode)
         -t INT    max number of threads to use [$OPTIONS{t}]
) if(@ARGV < 2);

# initialisation
my $REFERENCE = shift @ARGV;
my @ALN = ("-t $OPTIONS{t}");
my @INDEX = ("-a $OPTIONS{i}"); 
my @SAMSE = ();
my @SAMPE = ();
my $COMM;
my @VIEW = ("-b");
my @VREGION = ();
my @SORT = ();
my @PILEUP = ("-c");
my @VARFILTER = ("-p", "-D $OPTIONS{D}");
my @SEQUENCE;
my @PAIREDEND;
my $ALG = $OPTIONS{a};
my $BAMALIGN;
my @BAMALIGNSORTED;
my $SORTBAMALIGN;
my $SEQIND;
my $SEQ;
my $PE;
my $file;
my $line;
my $command;
my $name;
my $path;
my $suffix;
my $BAMALIGN;

# grab the sequence files and make sure if using paired ends we have a matching sequence file
if(defined $OPTIONS{p}) {
	my $TMP;
	foreach $file (sort @ARGV) {
		($name,$path,$suffix) = fileparse($file, (".fastq", ".txt"));
		if(grep(/1$/, $name)) {
			push @SEQUENCE, $file;
			$TMP = $path.$name;
			$TMP =~ s/1$/2/;
			$TMP .= $suffix;
		} elsif($file eq $TMP) {
			push @PAIREDEND, $file;
		} else {
			die("** for PE mode, a matching paired end must be specified (i.e. sequence_1.txt sequence_2.txt or sequence_1.fastq sequence_2.fastq)\n");
		}
	}
} else {
	@SEQUENCE = @ARGV;
}


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

# grab any parameters needed from the config file
if(defined $OPTIONS{c}) {
	open(DAT, $OPTIONS{c}) || print "Could not open config file! Using Default Settings\n";
	while ($line = <DAT>) {
		chomp $line;
		if($line eq "INDEX" || $line eq "ALN" || $line eq "SAMSE" || $line eq "SAMPE" || $line eq "VIEW" || $line eq "VREGION" || $line eq "SORT" || $line eq "PILEUP" || $line eq "VARFILTER") {
			$command = $line;
		} else {
			if($command eq "INDEX") {
				push(@INDEX, $line);
			} elsif($command eq "ALN") {
				push(@ALN, $line);
			} elsif($command eq "SAMSE") {
				push(@SAMSE, $line);
			} elsif($command eq "SAMPE") {
				push(@SAMPE, $line);
			} elsif($command eq "VIEW") {
				push(@VIEW, $line);
			} elsif($command eq "VREGION") {
				push(@VREGION, $line);
			} elsif($command eq "SORT") {
				push(@SORT, $line);
			} elsif($command eq "PILEUP") {
				push(@PILEUP, $line);
			} elsif($command eq "VARFILTER") {
				if(grep(/^-D/, $command)) {
					$VARFILTER[1] = $command
				} else {
					push(@VARFILTER, $line);
				}
			}
		}
	}
	close(DAT);
}

# determine output dir name
if($OPTIONS{d} ne "") {
	$OPTIONS{d} = "BWA_".$OPTIONS{d};
} else {
	$OPTIONS{d} = "BWA";
}

# create reference database
if(-e $REFERENCE.".bwt") {
	print "Reference Database Already Exists\n";
} else {
	print "Creating Reference Database\n";
	$COMM = "bwa index @INDEX $REFERENCE";
	print "$COMM\n";
	system($COMM);
	print "Database Created\n\n";
}	

# convert sequences to fastq format
my @SEQFASTQ;
print "Converting Sequences to Fastq format\n";
foreach $SEQ (@SEQUENCE) {
	if(!grep(/fastq$/, $SEQ)) {
		$SEQ = convert($SEQ);
		print "$SEQ Converted\n";
	} else {
		print "$SEQ Already In Standard/Sanger FASTQ Format\n";
	}
	push @SEQFASTQ, $SEQ;
}
print "\n";

# convert paired end sequences to fastq format
my @PEFASTQ;
if(defined $OPTIONS{p}) {
	print "Converting Paired End Sequences to Fastq format\n";
	foreach $PE (@PAIREDEND) {
		if(!grep(/fastq$/, $PE)) {
			$PE = convert($PE);
			print "$PE Converted\n";
		} else {
			print "$PE Already In Standard/Sanger FASTQ Format\n";
		}
		push @PEFASTQ, $PE;
	}
	print "\n";
}	

# create output dirs
($name,$path,$suffix) = fileparse($SEQUENCE[0], ".fastq");
my $DIR = $path.$OPTIONS{d};
unless(-d $DIR) {
	print "Creating folder $DIR\n\n";
	mkdir($DIR, 0777) || print $!;
}
unless(-d "$DIR/Binary") {
	print "Creating folder $DIR/Binary\n\n";
	mkdir("$DIR/Binary", 0777) || print $!;
}

	
# index the reference and add it to SAMPE and SAMSE
if(!-e $REFERENCE.".fai") {
	$COMM = "samtools faidx $REFERENCE";
	print "Creating Reflist\n";
	print "$COMM\n\n";
	system($COMM);
} else {
	print "Reflist already exists: using $REFERENCE.fai\n";
}
push @SAMPE, $REFERENCE;
push @SAMSE, $REFERENCE;
$DIR = $DIR."/";

for $SEQIND (0..$#SEQFASTQ) {
	# Align each sequence and convert to SAMtools
	my $SAMALIGN = $SEQFASTQ[$SEQIND];
	($name,$path,$suffix) = fileparse($SAMALIGN, ".fastq");
	$SAMALIGN = $DIR.$name.".sam";

	print "Running Alignment\n";
	my $ALIGNMENT = $DIR.$name.".sai";		
	$COMM = "bwa $ALG @ALN $REFERENCE $SEQFASTQ[$SEQIND] > $ALIGNMENT";
	print "$COMM\n";
	system($COMM);

	if(defined $OPTIONS{p}) {
		print "Running Paired End Alignment\n";
		my $PAIRED = $PEFASTQ[$SEQIND];
		($name,$path,$suffix) = fileparse($PAIRED, ".fastq");
		$PAIRED = $DIR.$name.".sai";		
		$COMM = "bwa $ALG @ALN $REFERENCE $PEFASTQ[$SEQIND] > $PAIRED";
		print "$COMM\n";
		system($COMM);

		$SAMALIGN =~ s/\d\.(\w*)$/.$1/;
		$COMM = "bwa sampe @SAMPE $ALIGNMENT $PAIRED $SEQFASTQ[$SEQIND] $PEFASTQ[$SEQIND] > $SAMALIGN";
		print "Alignment Finished\n\n";
	} else {
		$COMM = "bwa samse @SAMSE $ALIGNMENT $SEQFASTQ[$SEQIND] > $SAMALIGN";
		print "Alignment Finished\n\n";
	}
	print "Running SAMtools Conversion\n";
	print "$COMM\n";
	system($COMM);
	print "Conversion Complete\n\n";
	
	# convert to BAM - a binary format
	print "Converting to BAM\n";
	($name,$path,$suffix) = fileparse($SAMALIGN, ".sam");
	$BAMALIGN = $path."Binary/".$name.".bam";
	$COMM = "samtools view @VIEW -t $REFERENCE.fai -o $BAMALIGN $SAMALIGN @VREGION";
	print "$COMM\n";
	system($COMM);
	print "BAM Conversion Finished\n\n";
	
	# sort the alignments
	print "Sorting Alignments\n";
	$SORTBAMALIGN = $BAMALIGN;
	$SORTBAMALIGN =~ s/\.bam$/.sorted/;
	push @BAMALIGNSORTED, $SORTBAMALIGN.".bam";
	my @TMP = @SORT;
	push @TMP, $BAMALIGN;
	push @TMP, $SORTBAMALIGN;
	$COMM = "samtools sort @TMP";
	print "$COMM\n";
	system($COMM);
	unlink($BAMALIGN);
	print "Alignment Sorting Finished\n\n";
}

# Merge and sort if we have multiple lanes
if($#SEQFASTQ > 0)	{
	print "Merging Bam files\n";
	$BAMALIGN = $DIR."Binary/all_reads_aligned.bam";
	$COMM = "samtools merge $BAMALIGN @BAMALIGNSORTED";
	print "$COMM\n";
	system($COMM);
	print "Bam files merged\n\n";
	print "Sorting Merged Alignments\n";
	$SORTBAMALIGN = $DIR."Binary/all_reads_aligned.sorted";
	my @TMP = @SORT;
	push @TMP, $BAMALIGN;
	push @TMP, $SORTBAMALIGN;
	$COMM = "samtools sort @TMP";
	print "$COMM\n";
	system($COMM);
	unlink($BAMALIGN);
	print "Merged Alignment Sorting Finished\n\n";
}

# index the alignments
$SORTBAMALIGN = $SORTBAMALIGN.".bam";
print "Indexing Sorted Alignment File\n";
$COMM = "samtools index $SORTBAMALIGN";
print "$COMM\n";
system($COMM);
print "Indexing Finished\n\n";

# construct a pileup file
print "Constructing Pileup File\n";
my $PFILE = $DIR."reads_aligned.pileup";
$COMM = "samtools pileup @PILEUP -f $REFERENCE $SORTBAMALIGN > $PFILE";
print "$COMM\n";
system($COMM);
print "Pileup Finished\n\n";

# run the variations filter
print "Running Variations Filter\n";
$COMM = "samtools.pl varFilter @VARFILTER $PFILE &> $PFILE.all.snps";
print "$COMM\n";
system($COMM);

shift @VARFILTER;
$COMM = "samtools.pl varFilter @VARFILTER $PFILE > $PFILE.snps";
print "$COMM\n";
system($COMM);
print "Variations Filter Finished\n";

if(defined $OPTIONS{l}) {
	my $now = localtime time;
	print "\nLog File Closed $now\n";
	close(STDERR);
	close(STDOUT);
}

sub convert{
	print "Converting Sequences From Illumina FASTQ to Standard/Sanger FASTQ Format\n";
	my $ILL = $_[0];
	$_[0] =~ s/\.txt$/\.fastq/; 
	print "maq ill2sanger $ILL $_[0]\n";
	system("maq", "ill2sanger", $ILL, $_[0]);
	$_[0]
}

