#!/bin/env python

'''
Next generation sequencing pipeline.

Authors: Bernie Pope, Maria Doyle, Jason Ellul and Jason Li.

Description:

This program implements a workflow pipeline for next generation
sequencing, predominantly the read mapping task, up until variant
detection.

It uses the Ruffus library to make the description of the pipeline
more declarative.

It supports parallel evaluation of independent pipeline stages,
and can run stages on a cluster environment.

The pipeline is configured by an options file in YAML format,
including the actual commands which are run at each stage.
'''

from ruffus import *
import os.path
import shutil
from utils import (runStage, splitPath, getOptions, initLog, getCommand)

# Read the configuation options from file, determine the reference file
# and list of sequence files.
options = getOptions()
reference = options['reference']
sequences = options['sequences']
isPairedEnd = options['paired']
# Start the logging process.
logger = initLog(options)

# Index the reference file.
@files(reference, reference + '.bwt', logger)
def mkRefDataBase(reference, output, logger):
    runStage('mkRefDataBase', logger, options, reference, output)

# Index the reference file.
# XXX not sure why we need to do both mkRefDataBase and indexReference.
@follows(mkRefDataBase)
@files(reference, reference + '.fai', logger)
def indexReference(reference, output, logger):
    runStage('indexReference', logger, options, reference, output)

# Convert illumiar data (.txt) to sanger format (.fastq) if necessary
# XXX the input regex should only match .txt or .fastq.
@transform(sequences, regex('(.*)\..*$'), r'\1.fastq', logger)
def illToSanger(sequence, output, logger):
    if sequence.endswith('.txt'):
        runStage('illToSanger', logger, options, sequence, output)

# Align sequence reads to the reference genome.
@follows(indexReference)
@transform(illToSanger, suffix('.fastq'), '.sai', logger)
def alignSequence(sequence, output, logger):
    runStage('alignSequence', logger, options, reference, sequence, output)

if isPairedEnd:
   input = r'(.+)_1\.sai'
   extraInputs = [r'\1_2.sai', r'\1_1.fastq', r'\1_2.fastq']
else:
   input = r'(.+)\.sai'
   extraInputs = [r'\1.fastq']

# Convert alignments to SAM format.
@transform(alignSequence, regex(input), add_inputs(extraInputs), r'\1.sam', logger)
def alignToSam(inputs, output, logger):
    if isPairedEnd:
        align1, [align2, seq1, seq2] = inputs
        runStage('alignToSamPE', logger, options, reference, align2, align2, seq1, seq2, output)
    else:
        align, seq = inputs
        runStage('alignToSamSE', logger, options, reference, align, seq, output)

#@transform(alignSequence, suffix('.sai'), '.sam', '.fastq', logger)
#def alignToSam(alignment, output, sequence, logger):
#    runStage('alignToSam', logger, options, reference, alignment, sequence, output)

# Convert SAM alignments to BAM format.
@transform(alignToSam, suffix('.sam'), '.bam', logger)
def samToBam(samFile, output, logger):
    indexedReference = reference + '.fai'
    runStage('samToBam', logger, options, indexedReference, samFile, output)

# Sort BAM alignments by (leftmost?) coordinates.
@transform(samToBam, suffix('.bam'), '.sorted.bam', logger)
def sortBam(bamFile, output, logger):
    (prefix, name, ext) = splitPath(output)
    outFile = os.path.join(prefix,name)
    runStage('sortBam', logger, options, bamFile, outFile)

# Merge all the sorted BAM alignments.
@collate(sortBam, regex(r"(.*/|^).*.bam"), r"\1/all_reads_aligned.sorted.bam", logger)
def mergeBams(sortedBams, output, logger):
    if len(sortedBams) == 1:
        shutil.copyfile(sortedBams[0], output)
    else: 
        bams = ' '.join(sortedBams)
        runStage('mergeBams', logger, options, bams, output)

# Index sorted (merged) BAM alignment for fast access.
@transform(mergeBams, suffix('.bam'), '.bam.bai', logger)
def indexMergedBams(mergedBamsFile, output, logger):
    runStage('indexMergedBams', logger, options, mergedBamsFile, output)

# Convert (sorted merged) BAM alignment to pileup format.
@transform(indexMergedBams, suffix('.bam.bai'), '.pileup', logger)
def pileup(baiFile, output, logger):
    (prefix, name, ext) = splitPath(baiFile)
    bamAlignFile = os.path.join(prefix, name)
    runStage('pileup', logger, options, reference, bamAlignFile, output)

# Determine variants.
@transform(pileup, suffix('.pileup'), '.all.snps', logger)
def variationAll(pileupFile, output, logger):
    runStage('variationAll', logger, options, pileupFile, output)

# Determine variants.
@transform(pileup, suffix('.pileup'), '.snps', logger)
def variation(pileupFile, output, logger):
    runStage('variation', logger, options, pileupFile, output)

# Define the end-points of the pipeline.
pipeline = [variationAll, variation]

# Invoke the pipeline.
pipelineOptions = options['pipeline']
if pipelineOptions['style'] == 'run':
    # Perform the pipeline steps.
    pipeline_run(pipeline, multiprocess = pipelineOptions['procs'], logger = logger['proxy'])
elif pipelineOptions['style'] == 'flowchart':
    # Draw the pipeline as a diagram.
    pipeline_printout_graph ('flowchart.svg', 'svg', pipeline, no_key_legend = False)
