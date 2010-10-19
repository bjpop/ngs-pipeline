from Logger import initLogger
from Pipeline import *
from OutputDir import mkOutputDir
from Options import getConfig
from sys import exit

options = getOptions()

initLogger(options['logging'])
reference = options['reference']
sequences = options['sequences']

if len(sequences == 0):
    exit('At least one sequence file must be specified')

# Create the output directory
mkOutputDir(options['dir'])

# Step 1. Create reference database.
mkRefDataBase(reference, '-a ' + options['bwa']['index'])

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
