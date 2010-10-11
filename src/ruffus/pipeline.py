from ruffus import *
import os.path
import utils

options = utils.getOptions()
reference = options['reference']
sequences = options['sequences']
logger = utils.initLog(options)

@files(reference, reference + '.bwt', logger)
def mkRefDataBase(reference, output, logger):
    utils.runCommand('bwa index %s -a bwtsw' % reference, logger)

@follows(mkRefDataBase)
@files(reference, reference + '.fai', logger)
def indexReference(reference, output, logger):
    utils.runCommand('samtools faidx %s' % reference, logger)

@follows(indexReference)
@transform(sequences, regex(r"(.*).fastq"), r"\1.sai", logger)
def alignSequence(sequence, output, logger):
    algorithm = options['bwa']['align']
    threads = options['bwa']['threads']
    utils.runCommand('bwa %s -t %s %s %s > %s' % (algorithm, threads, reference, sequence, output),
                     logger)

@transform(alignSequence, regex(r"(.*).sai"), r"\1.sam", r"\1.fastq", logger)
def alignToSam(alignment, output, sequence, logger):
    utils.runCommand('bwa samse %s %s %s > %s' % (reference, alignment, sequence, output),
                     logger)

@transform(alignToSam, regex(r"(.*).sam"), r"\1.bam", logger)
def samToBam(samFile, output, logger):
    regions = ' '.join(options['samtools']['viewRegions'])
    indexedReference = reference + '.fai'
    utils.runCommand('samtools view -b -t %s -o %s %s %s' % (indexedReference, output, samFile, regions),
                     logger)

@transform(samToBam, regex(r"(.*).bam"), r"\1.sorted.bam", logger)
def sortBam(bamFile, output, logger):
    sortFlags = options['samtools']['sort']
    (prefix, name, ext) = utils.splitPath(output)
    utils.runCommand('samtools sort %s %s %s' % (sortFlags, bamFile, os.path.join(prefix,name)),
                     logger)

@collate(sortBam, regex(r"(.*/|^).*.bam"), r"\1/all_reads_aligned.sorted.bam", logger)
def mergeBams(sortedBams, output, logger):
    bams = ' '.join(sortedBams)
    utils.runCommand('samtools merge %s %s' % (output, bams), logger)

@transform(mergeBams, regex(r"(.*.bam)"), r"\1.bai", logger)
def indexMergedBams(mergedBamsFile, output, logger):
    utils.runCommand('samtools index %s' % mergedBamsFile, logger)

@transform(indexMergedBams, regex(r"(.*).bam.bai"), r"\1.pileup", logger)
def pileup(baiFile, output, logger):
    (prefix, name, ext) = utils.splitPath(baiFile)
    bamAlignFile = os.path.join(prefix, name)
    pileupFlags = options['samtools']['pileup']
    utils.runCommand('samtools pileup %s -f %s %s > %s' % (pileupFlags, reference, bamAlignFile, output),
                     logger)

@transform(pileup, regex(r"(.*.pileup)"), r"\1.all.snps", logger)
def variationAll(pileupFile, output, logger):
    depth = options['samtools']['maxReadDepth']
    utils.runCommand('samtools.pl varFilter -p -D %d %s &> %s' % (depth, pileupFile, output),
                     logger)

@transform(pileup, regex(r"(.*.pileup)"), r"\1.snps", logger)
def variation(pileupFile, output, logger):
    depth = options['samtools']['maxReadDepth']
    utils.runCommand('samtools.pl varFilter -D %d %s &> %s' % (depth, pileupFile, output),
                     logger)

pipeline = [variationAll, variation]

pipelineOptions = options['pipeline']
if pipelineOptions['style'] == 'run':
    pipeline_run(pipeline, multiprocess = pipelineOptions['procs'], logger = logger['proxy'])
elif pipelineOptions['style'] == 'flowchart':
    pipeline_printout_graph ('flowchart.svg', 'svg', pipeline, no_key_legend = False)
