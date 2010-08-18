from os.path import (split, splitext, exists)
import os.path
import os
import sys
from shutil import copy
from yaml import load
from sys import argv
from tempfile import mkstemp
from logging import (info, warning, debug, DEBUG, WARNING, INFO, CRITICAL, ERROR, basicConfig)
import subprocess

################################################################################
#
# Processing stages.
#
################################################################################

logFile = None

def runCommand(message, command):
    info(message)
    info(command)
    process = subprocess.Popen(command, shell=True, stdout=logFile, stderr=logFile)
    status = os.waitpid(process.pid, 0)[1]
    if status != 0:
        warning("command '%s' returned non-zero status: %d'" % (command, status))
    return status

def align(algorithm, threadFlags, dir, reference, sequence):
    (root, ext) = splitext(sequence)
    if ext != ".fastq":
        sys.exit("align: sequence file %s does not have .fastq extension" % sequence)
    alignmentFile = root + ".sai"
    command = "bwa %s -t %s %s %s > %s" % (algorithm, threadFlags, reference, sequence, alignmentFile)
    runCommand ("Running Alignment", command)
    return alignment

def align2Sam (flags, reference, sequence, alignment):
    (root, ext) = splitext(alignment)
    if ext != ".sai":
        sys.exit("align2Sam: alignment file %s does not have .sai extension" % alignment)
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
        info("illumina2FastQ: file %s already in Standard/Sanger FASTQ Format" % file)
    elif ext == ".txt":
        command = "maq ill2sanger %s %s" % (file, sangerFile)
        runCommand("Converting illumina file to Standard/Sanger FASTQ Format", command)
    else:
        sys.exit("sequence file %s not in illumina (.txt) or sanger (.fastq) format" % file)
    return sangerFile

# create reference database
def mkRefDataBase(reference, indexFlags):
    refDatabase = reference + ".bwt"
    if exists(refDatabase):
        info("Reference database already exists: using %s" % refDatabase)
    else:
        command = "bwa index %s %s" % (indexFlags, reference)
        runCommand("Creating Reference Database", command)
    return refDatabase

# index the reference
def indexReference(reference):
    referenceIndex = reference + ".fai"
    if exists(referenceIndex):
        info("Reflist already exists: using %s" % referenceIndex)
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
        sys.exit("Attempt to merge and index zero BAM files")
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
        sys.exit("sam2Bam: file %s does not end in .sam extension" % samFile)
    bamAlignFile = os.path.join(root + "Binary", name + ".bam")
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
   command = "samtools.pl varFilter %s %s &> %s" % (varFilterFlags, pileupFile, snpsFile)
   runCommand("Running Variations Filter", command)
   return snpsFile

################################################################################
#
# Output directory
#
################################################################################

def mkOutputDir(dir):
    def tryMakeDir(dir):
       if not exists(dir):
          info('Creating folder %s' % dir)
          try:
              mkdir(dir, 0777)
          except IOError, e:
              sys.exit('%s\nFailed to make directory %s' % (e, dir))
    tryMakeDir(dir)
    tryMakeDir(os.path.join(dir, 'Binary'))

################################################################################
#
# Options processing
#
################################################################################

defaultBWA = {
   'align' : 'aln',
   'index' : 'bwtsw',
   'threads' : 1
}

defaultSamtools = {
   'maxReadDepth' : 1000000
}

defaultLogging = {
   'level': 'info',
   'file': 'pipeline.log'
}

defaultConfig = {
   'bwa' : defaultBWA,
   'samtools' : defaultSamtools,
   'logging': defaultLogging,
   'dir' : 'BWA'
}

defaultOptionsFile = 'ngs_pipeline.opt'

def getOptions():
    if len(argv) <= 1:
        filename = defaultOptionsFile
    else:
        filename = argv[1]
    try:
        file = open(filename)
        contents = file.read()
    except IOError, e:
        sys.exit('Could not open configuration file: %s' % filename)
    newConfig = load(contents)
    defaultConfig.update(newConfig)
    return defaultConfig

################################################################################
#
# Parallel job submission.
#
################################################################################

scriptTemplate = '''
#!/bin/bash
#PBS -q smp
#PBS -l mem=10gb
#PBS -l walltime=%s
#PBS -m abe
#PBS -N %s 
%s
cd $PBS_O_WORKDIR
%s
%s
'''

def makeScript(walltime, jobName, dependencies, modules, command): 
    return scriptTemplate % (walltime, jobName, dependencies, modules, command)

def submitJob(jobName, moduleList, dependList, wallTime, command): 
    if dependList == []:
        dependStr = ''
    else:
        dependStr = '#PBS -W depend=afterok:%s' % ':'.join(dependIds)
    if moduleList == []:
        moduleStr = ''
    else:
        modileStr = 'module load %s' % ' '.join(moduleList)
    script = makeScript(wallTime, jobName, dependStr, moduleStr, command)
    (handle, filepath) = mkstemp()
    handle.write(script)
    qsub = "qsub " + filepath
    status = system(qsub)
    handle.close()
    remove(filepath)
    return status

################################################################################
#
# Output logging.
#
################################################################################

def initLogger(options):
    file = options['file']
    levels = {
        'debug': DEBUG,
        'info': INFO,
        'warning': WARNING,
        'error': ERROR,
        'critical': CRITICAL }
    if file:
        try:
            handle = open(file, 'w')
        except IOError, e:
            sys.exit(str(e) + '\nCould not open logging file: %s' % file)
        global logFile
        logFile = handle
        basicConfig(stream=handle,level=levels[options['level']])

################################################################################
#
# Tie all the stages together.
#
################################################################################

def pipeline(options, reference, sequences):
    # Step 1. Create reference database.
    mkRefDataBase(reference, '-a ' + options['bwa']['index'])

'''
    # Step 2. Index the reference.
    indexReference(reference)

    sortedBams = []
    for seq in sequences:
        # Step 3. Convert sequence to fastq format.
        seqFastq = illumina2FastQ(seq)
        # Step 4. Align sequence to the reference database.
        seqAlign = align(options['bwa']['align'], options['bwa']['threads'], dir, reference, seqFastq)
        # Step 5. Convert alignment to SAM format.
        seqSam = align2Sam('samse', reference, seq, seqAlign)
        # Step 6. Convert SAM to BAM.
        seqBam = sam2Bam('-b', reference, seqSam)
        # Step 7. Sort the BAM file.
        seqBamSorted = sortBam(seqBam)
        sortedBams.append(seqBamSorted)

    # Steps 8-10. Merge BAM files, and align the result.
    bamAlign = mergeBamsAndIndex(dir, sortedBams)

    # Step 11. Construct pileup file.
    pileupFile = pileup(dir, '-c', reference, bamAlign)

    # Step 12-13. Run variation filter.
    maxReadDepth = options['samtools']['maxReadDepth']
    varFilter('-p -D %d' % maxReadDepth, pileupFile, '.all.snps')
    varFilter('-D %d' % maxReadDepth, pileupFile, '.snps')
'''

def main():
    options = getOptions()
    #print options
    reference = options['reference']
    sequences = options['sequences']
    if len(sequences) == 0:
        sys.exit('At least one sequence file must be specified')
    initLogger(options['logging'])
    # Create the output directory
    mkOutputDir(options['dir'])
    pipeline(options, reference, sequences)

main()
