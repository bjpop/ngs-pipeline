from ruffus import *
import os.path
from utils import (runStage, splitPath, getOptions, initLog, getCommand)

options = getOptions()
reference = options['reference']
sequences = options['sequences']
logger = initLog(options)

@files(reference, reference + '.bwt', logger)
def mkRefDataBase(reference, output, logger):
    runStage('mkRefDataBase', logger, options, reference, output)

@follows(mkRefDataBase)
@files(reference, reference + '.fai', logger)
def indexReference(reference, output, logger):
    runStage('indexReference', logger, options, reference, output)

@follows(indexReference)
@transform(sequences, regex(r"(.*).fastq"), r"\1.sai", logger)
def alignSequence(sequence, output, logger):
    runStage('alignSequence', logger, options, reference, sequence, output)

@transform(alignSequence, regex(r"(.*).sai"), r"\1.sam", r"\1.fastq", logger)
def alignToSam(alignment, output, sequence, logger):
    runStage('alignToSam', logger, options, reference, alignment, sequence, output)

@transform(alignToSam, regex(r"(.*).sam"), r"\1.bam", logger)
def samToBam(samFile, output, logger):
    indexedReference = reference + '.fai'
    runStage('samToBam', logger, options, indexedReference, samFile, output)

@transform(samToBam, regex(r"(.*).bam"), r"\1.sorted.bam", logger)
def sortBam(bamFile, output, logger):
    (prefix, name, ext) = splitPath(output)
    outFile = os.path.join(prefix,name)
    runStage('sortBam', logger, options, bamFile, outFile)

@collate(sortBam, regex(r"(.*/|^).*.bam"), r"\1/all_reads_aligned.sorted.bam", logger)
def mergeBams(sortedBams, output, logger):
    bams = ' '.join(sortedBams)
    runStage('mergeBams', logger, options, bams, output)

@transform(mergeBams, regex(r"(.*.bam)"), r"\1.bai", logger)
def indexMergedBams(mergedBamsFile, output, logger):
    runStage('indexMergedBams', logger, options, mergedBamsFile, output)

@transform(indexMergedBams, regex(r"(.*).bam.bai"), r"\1.pileup", logger)
def pileup(baiFile, output, logger):
    (prefix, name, ext) = splitPath(baiFile)
    bamAlignFile = os.path.join(prefix, name)
    runStage('pileup', logger, options, reference, bamAlignFile, output)

@transform(pileup, regex(r"(.*.pileup)"), r"\1.all.snps", logger)
def variationAll(pileupFile, output, logger):
    runStage('variationAll', logger, options, pileupFile, output)

@transform(pileup, regex(r"(.*.pileup)"), r"\1.snps", logger)
def variation(pileupFile, output, logger):
    runStage('variation', logger, options, pileupFile, output)

pipeline = [variationAll, variation]

pipelineOptions = options['pipeline']
if pipelineOptions['style'] == 'run':
    pipeline_run(pipeline, multiprocess = pipelineOptions['procs'], logger = logger['proxy'])
elif pipelineOptions['style'] == 'flowchart':
    pipeline_printout_graph ('flowchart.svg', 'svg', pipeline, no_key_legend = False)
