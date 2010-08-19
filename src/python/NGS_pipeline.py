import os.path as path
import os
import sys
import shutil as sh
import yaml
import tempfile
import logging as log
import subprocess as proc

################################################################################
#
# Processing stages.
#
################################################################################

logFile = None

def fail(msg):
    log.error(msg)
    sys.exit()

def runCommand(message, command):
    log.info(message)
    log.info(command)
    process = proc.Popen(command, shell=True, stdout=logFile, stderr=logFile)
    status = os.waitpid(process.pid, 0)[1]
    if status != 0:
        log.warning("command '%s' returned non-zero status: %d'" % (command, status))
    return status

def align(algorithm, threadFlags, reference, sequence, outputDir):
    (path, name, ext) = splitPath(sequence)
    if ext != '.fastq':
        fail('align: sequence file %s does not have .fastq extension' % sequence)
    alignmentFile = os.path.join(outputDir, name + '.sai')
    command = 'bwa %s -t %s %s %s > %s' % (algorithm, threadFlags, reference, sequence, alignmentFile)
    runCommand ('Running Alignment', command)
    return alignmentFile

def align2Sam (flags, reference, sequence, alignment, outputDir):
    (path, name, ext) = splitPath(alignment)
    if ext != '.sai':
        fail('align2Sam: alignment file %s does not have .sai extension' % alignment)
    samfile = os.path.join(outputDir, name + '.sam')
    command = 'bwa %s %s %s %s > %s' % (flags, reference, alignment, sequence, samfile)
    runCommand('Converting alignment to SAM format', command)
    return samfile

def ills2FastQs(sequences):
    return map(illumina2FastQ, sequences)

def illumina2FastQ(file):
    (root, ext) = path.splitext(file)
    sangerFile = root + '.fastq'
    if ext == '.fastq':
        log.info('illumina2FastQ: file %s already in Standard/Sanger FASTQ Format' % file)
    elif ext == '.txt':
        command = 'maq ill2sanger %s %s' % (file, sangerFile)
        runCommand('Converting illumina file to Standard/Sanger FASTQ Format', command)
    else:
        fail('sequence file %s not in illumina (.txt) or sanger (.fastq) format' % file)
    return sangerFile

# create reference database
def mkRefDataBase(reference, indexFlags):
    refDatabase = reference + '.bwt'
    if path.exists(refDatabase):
        log.info('Reference database already exists: using %s' % refDatabase)
    else:
        command = 'bwa index %s %s' % (indexFlags, reference)
        runCommand('Creating Reference Database', command)
    return refDatabase

# index the reference
def indexReference(reference):
    referenceIndex = reference + '.fai'
    if path.exists(referenceIndex):
        log.info('Reflist already exists: using %s' % referenceIndex)
    else:
        command = 'samtools faidx %s' % reference
        runCommand('Creating Reflist', command)
    return referenceIndex

def mergeBamsAndIndex(bamFiles, outputDir):
    bamAlignFile = os.path.join(outputDir, 'Binary', 'all_reads_aligned.bam')
    numBams = len(bamFiles);
    # If there is more than one BAM file then merge them.
    if numBams > 1:
        command = 'samtools merge %s %s' % (bamAlignFile, ' '.join(bamFiles))
        runCommand ('Merging Bam files', command)
    elif numBams == 1:
        copy(bamfiles[0], bamAlignFile)
    else:
        fail('Attempt to merge and index zero BAM files')
    command = 'samtools index %s' % bamAlignFile
    runCommand('Indexing Alignment File', command)
    return bamAlignFile

def pileup(pileupFlags, reference, bamAlignFile, outputDir):
    pileupFile = os.path.join(outputDir, 'reads_aligned.pileup')
    command = 'samtools pileup %s -f %s %s > %s' % (pileupFlags, reference, bamAlignFile, pileupFile)
    runCommand('Constructing Pileup File', command)
    return pileupFile

# Convert SAM to BAM
# Assumes the Binary directory has been created already.
def sam2Bam(reference, samFile, viewRegions, outputDir):
    (root, base) = path.split(samFile)
    (name, ext) = path.splitext(base)
    if ext != '.sam':
        fail('sam2Bam: file %s does not end in .sam extension' % samFile)
    bamAlignFile = os.path.join(outputDir, 'Binary', name + '.bam')
    faiFile = reference + '.fai'
    command = 'samtools view -b -t %s -o %s %s %s' % (faiFile, bamAlignFile, samFile, ' '.join(viewRegions))
    runCommand('Converting to BAM', command)
    return bamAlignFile

# Sort a BAM file
def sortBam(bamAlignFile, sortFlags, outputDir):
    (path, name, ext) = splitPath(bamAlignFile)
    sortedBamFilePrefix = os.path.join(outputDir, 'Binary', name + '.sorted')
    command = 'samtools sort %s %s %s' % (sortFlags, bamAlignFile, sortedBamFilePrefix)
    runCommand('Sorting Bam alignments', command)
    # I think we should leave the decision to remove the BAM file to the
    # caller of this function.
    #remove(bamAlignFile)
    return sortedBamFilePrefix + '.bam'

# run the variations filter
def varFilter(varFilterFlags, pileupFile, outputExtension, outputDir):
   (path, name, ext) = splitPath(pileupFile)
   snpsFile = os.path.join(outputDir, name + outputExtension)
   command = 'samtools.pl varFilter %s %s &> %s' % (varFilterFlags, pileupFile, snpsFile)
   runCommand('Running Variations Filter', command)
   return snpsFile

################################################################################
#
# Output directory
#
################################################################################

def mkOutputDir(dir):
    def tryMakeDir(dir):
       if not path.exists(dir):
          log.info('Creating folder %s' % dir)
          try:
              os.mkdir(dir, 0777)
          except IOError, e:
              fail('%s\nFailed to make directory %s' % (e, dir))
    tryMakeDir(dir)
    tryMakeDir(os.path.join(dir, 'Binary'))

################################################################################
#
# Filepath splitting
#
################################################################################

def splitPath(path):
    (prefix, base) = os.path.split(path)
    (name, ext) = os.path.splitext(base)
    return (prefix, name, ext)

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
   'maxReadDepth' : 1000000,
   'viewRegions' : [],
   'sort' : '',
   'pileup' : '-c'
}

defaultLogging = {
   'level': 'info',
   'file': 'pipeline.log'
}

defaultConfig = {
   'bwa' : defaultBWA,
   'samtools' : defaultSamtools,
   'logging': defaultLogging,
   'dir' : 'BWA',
   'reference' : None,
   'sequences' : []
}

defaultOptionsFile = 'ngs_pipeline.opt'

def getOptions():
    if len(sys.argv) <= 1:
        filename = defaultOptionsFile
    else:
        filename = sys.argv[1]
    try:
        file = open(filename)
        contents = file.read()
    except IOError, e:
        fail('Could not open configuration file: %s' % filename)
    newConfig = yaml.load(contents)
    #defaultConfig.update(newConfig)
    updateDict(defaultConfig, newConfig)
    return defaultConfig

def updateDict(old, new):
    for k in new:
        newVal = new[k]
        if k in old:
            oldVal = old[k]
            if type(oldVal) == dict:
                if type(newVal) == dict:
                    updateDict(oldVal, newVal)
                else:
                    fail('incorrect options for %s, was expecting multiple items but found only one' % k)
            else:
                old[k] = newVal

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
    (handle, filepath) = tempfile.mkstemp()
    handle.write(script)
    qsub = 'qsub ' + filepath
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
    global logFile
    file = options['file']
    levels = {
        'debug': log.DEBUG,
        'info': log.INFO,
        'warning': log.WARNING,
        'error': log.ERROR,
        'critical': log.CRITICAL }
    if file == 'stdout':
        logfile = sys.stdout
    else:
        try:
            logFile = open(file, 'w')
        except IOError, e:
            fail(str(e) + '\nCould not open logging file: %s' % file)
    log.basicConfig(stream=logFile,level=levels[options['level']])

################################################################################
#
# Tie all the stages together.
#
################################################################################

def pipeline(options, outputDir, reference, sequences):
    # Step 1. Create reference database.
    mkRefDataBase(reference, '-a ' + options['bwa']['index'])

    # Step 2. Index the reference.
    indexReference(reference)

    sortedBams = []
    for seq in sequences:
        # Step 3. Convert sequence to fastq format.
        seqFastq = illumina2FastQ(seq)
        # Step 4. Align sequence to the reference database.
        seqAlign = align(options['bwa']['align'], options['bwa']['threads'], reference, seqFastq, outputDir)
        # Step 5. Convert alignment to SAM format.
        seqSam = align2Sam('samse', reference, seq, seqAlign, outputDir)
        # Step 6. Convert SAM to BAM.
        seqBam = sam2Bam(reference, seqSam, options['samtools']['viewRegions'], outputDir)
        # Step 7. Sort the BAM file.
        seqBamSorted = sortBam(seqBam, options['samtools']['sort'], outputDir)
        sortedBams.append(seqBamSorted)

    # Steps 8-10. Merge BAM files, and align the result.
    bamAlign = mergeBamsAndIndex(sortedBams, outputDir)

    # Step 11. Construct pileup file.
    pileupFile = pileup(options['samtools']['pileup'], reference, bamAlign, outputDir)

    # Step 12-13. Run variation filter.
    maxReadDepth = options['samtools']['maxReadDepth']
    varFilter('-p -D %d' % maxReadDepth, pileupFile, '.all.snps', outputDir)
    varFilter('-D %d' % maxReadDepth, pileupFile, '.snps', outputDir)

def main():
    options = getOptions()
    reference = options['reference']
    if not reference:
        fail('One reference file must be specified')
    sequences = options['sequences']
    if len(sequences) == 0:
        fail('At least one sequence file must be specified')
    initLogger(options['logging'])
    outputDir = options['dir']
    mkOutputDir(outputDir)
    pipeline(options, outputDir, reference, sequences)

main()
