from os.path import (split, splitext, join, exists)
from os import (system, remove, mkdir)
from Logger import (logInfo, logWarn)
from sys import exit
from shutil import copy

def mkOutputDir(dir):
    def tryMakeDir(dir):
       if not exists(dir):
          logInfo('Creating folder %s' % dir)
          try:
              mkdir(dir, '0777')
          except IOError, e:
              exit('%s\nFailed to make directory %s' % (e, dir))
    tryMakeDir(dir)
    tryMakeDir(os.path.join(dir, 'Binary'))

def runCommand(msg, cmd):
    logInfo(msg)
    logInfo(command)
    status = system(command)
    if status != 0:
        logWarn("command '%s' returned non-zero status: %d'" % (command, status))
    return status

def align(algorithm, threadFlags, dir, reference, sequence):
    (root, ext) = splitext(sequence)
    if ext != ".fastq":
        exit("align: sequence file %s does not have .fastq extension" % sequence)
    alignmentFile = root + ".sai"
    command = "bwa %s -t %s %s %s > %s" % (algorithm, threadFlags, reference, sequence, alignmentFile)
    runCommand ("Running Alignment", command)
    return alignment

def align2Sam (flags, reference, sequence, alignment):
    (root, ext) = splitext(alignment)
    if ext != ".sai":
        exit("align2Sam: alignment file %s does not have .sai extension" % alignment)
    samfile = root + ".sam"
    command = "bwa %s %s %s %s > %s" % (flags, reference, alignment, sequence, samfile)
    runCommand("Converting alignment to SAM format", command)
    return samfile

def ills2FastQs(sequences):
    return map(illumina2FastQ, sequences)

def illumina2FastQ(file):
    (root, ext) = split(file)
    sangerFile = root + ".fastq"
    if ext == ".fastq":
        logInfo("illumina2FastQ: file %s already in Standard/Sanger FASTQ Format" % file)
    elif ext == ".txt":
        command = "maq ill2sanger %s %s" % (file, sangerFile)
        runCommand("Converting illumina file to Standard/Sanger FASTQ Format", command)
    else:
        exit("sequence file %s not in illumina (.txt) or sanger (.fastq) format" % file)
    return sangerFile

# create reference database
def mkRefDataBase(reference, indexFlags):
    refDatabase = referenece + ".bwt"
    if exists(refDatabase) 
        logInfo("Reference database already exists: using %s" % refDatabase)
    else:
        command = "bwa index %s %s" % (indexFlags, reference)
        runCommand("Creating Reference Database", command)
    return refDatabase 

# index the reference
def indexReference(reference):
    referenceIndex = reference + .".fai"
    if exists(referenceIndex):
        logInfo("Reflist already exists: using %s" % referenceIndex)
    else:
        command = "samtools faidx %s" % reference
        runCommand("Creating Reflist", command)
    return referenceIndex

def mergeBamsAndIndex(dir, bamFiles):
    bamAlignFile = dir + "Binary/all_reads_aligned.bam"
    numBams = len(bamFiles);
    # If there is more than one BAM file then merge them.
    if numBams > 1:
        command = "samtools merge %s %s" % (bamAlignFile, ' '.join(bamFiles))
        runCommand ("Merging Bam files", command)
    elif numBams == 1:
        copy(bamfiles[0], bamAlignFile)
    else:
        exit("Attempt to merge and index zero BAM files")
    command = "samtools index %s" % bamAlignFile
    runCommand("Indexing Alignment File", command)
    return bamAlignFile

def pileup(dir, pileupFlags, reference, bamAlignFile):
    pileupFile = dir + "reads_aligned.pileup"
    command = "samtools pileup %s -f %s %s > %s" % (pileupFlags, reference, bamAlignFile, pileupFile)
    runCommand("Constructing Pileup File", command)
    return pileupFile

# Convert SAM to BAM
# Assumes the Binary directory has been created already.
def sam2Bam(vregion, viewFlag, reference, samFile):
    (root, base) = split(samFile)
    (name, ext) = splitext(base)
    if ext != ".sam":
        exit("sam2Bam: file %s does not end in .sam extension" % samFile)
    bamAlignFile = os.path.join(root + "Binary", name + ".bam")
    command = "samtools view $VIEWFLAG -t $REFERENCE.fai -o $BAMALIGN $SAMFILE @VREGION";
    faiFile = reference + ".fai"
    command = "samtools view %s -t %s -o %s %s %s" % (viewFlag, faiFile, bamAlignFile, samFile, vregion) 
    runCommand("Converting to BAM", command)
    return bamAlignFile

# Sort a BAM file
def sortBam(bamAlignFile, sortFlags):
    (root, ext) = splitext(bamAlignFile)
    sortedBamFile = root + ".sorted"
    command = "samtools sort %s %s %s" % (sortFlags, bamAlignFile, sortedBamFile)
    runCommand("Sorting Bam alignments", command)
    # I think we should leave the decision to remove the BAM file to the
    # caller of this function.
    #remove(bamAlignFile)
    return sortedBamFile + ".bam"

# run the variations filter
def varFilter(varFilterFlags, pileupFile, outputExtension):
   snpsFile = pileupFile + outputExtension
   command = "samtools.pl varFilter $VARFILTERFLAGS $PILEUPFILE &> $SNPSFILE";
   command = "samtools.pl varFilter %s %s &> %s" % (varFilterFlags, pileupFile, snpsFile)
   runCommand("Running Variations Filter", command)
   return snpsFile 
