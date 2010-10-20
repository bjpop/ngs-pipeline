from ruffus import *
import os.path
from utils import (runCommand, splitPath, getOptions, initLog, getCommand)

options = getOptions()
reference = options['reference']
sequences = options['sequences']
logger = initLog(options)

@files(reference, reference + '.bwt', logger)
def mkRefDataBase(reference, output, logger):
    command = getCommand('mkRefDataBase', options)
    runCommand(command(reference,output), logger)

@follows(mkRefDataBase)
@files(reference, reference + '.fai', logger)
def indexReference(reference, output, logger):
    command = getCommand('indexReference', options)
    runCommand(command(reference, output), logger)

@follows(indexReference)
@transform(sequences, regex(r"(.*).fastq"), r"\1.sai", logger)
def alignSequence(sequence, output, logger):
    command = getCommand('alignSequence', options)
    runCommand(command(reference, sequence, output), logger)

@transform(alignSequence, regex(r"(.*).sai"), r"\1.sam", r"\1.fastq", logger)
def alignToSam(alignment, output, sequence, logger):
    command = getCommand('alignToSam', options)
    runCommand(command(reference, alignment, sequence, output), logger)

@transform(alignToSam, regex(r"(.*).sam"), r"\1.bam", logger)
def samToBam(samFile, output, logger):
    indexedReference = reference + '.fai'
    command = getCommand('samToBam', options)
    runCommand(command(indexedReference, samFile, output), logger)

@transform(samToBam, regex(r"(.*).bam"), r"\1.sorted.bam", logger)
def sortBam(bamFile, output, logger):
    (prefix, name, ext) = splitPath(output)
    outFile = os.path.join(prefix,name)
    command = getCommand('sortBam', options)
    runCommand(command(bamFile, outFile), logger)

@collate(sortBam, regex(r"(.*/|^).*.bam"), r"\1/all_reads_aligned.sorted.bam", logger)
def mergeBams(sortedBams, output, logger):
    bams = ' '.join(sortedBams)
    command = getCommand('mergeBams', options)
    runCommand(command(bams, output), logger)

@transform(mergeBams, regex(r"(.*.bam)"), r"\1.bai", logger)
def indexMergedBams(mergedBamsFile, output, logger):
    command = getCommand('indexMergedBams', options)
    runCommand(command(mergedBamsFile, output), logger)

@transform(indexMergedBams, regex(r"(.*).bam.bai"), r"\1.pileup", logger)
def pileup(baiFile, output, logger):
    (prefix, name, ext) = splitPath(baiFile)
    bamAlignFile = os.path.join(prefix, name)
    command = getCommand('pileup', options)
    runCommand(command(reference, bamAlignFile, output), logger)

@transform(pileup, regex(r"(.*.pileup)"), r"\1.all.snps", logger)
def variationAll(pileupFile, output, logger):
    command = getCommand('variationAll', options)
    runCommand(command(pileupFile, output), logger)

@transform(pileup, regex(r"(.*.pileup)"), r"\1.snps", logger)
def variation(pileupFile, output, logger):
    command = getCommand('variation', options)
    runCommand(command(pileupFile, output), logger)

pipeline = [variationAll, variation]

pipelineOptions = options['pipeline']
if pipelineOptions['style'] == 'run':
    pipeline_run(pipeline, multiprocess = pipelineOptions['procs'], logger = logger['proxy'])
elif pipelineOptions['style'] == 'flowchart':
    pipeline_printout_graph ('flowchart.svg', 'svg', pipeline, no_key_legend = False)
