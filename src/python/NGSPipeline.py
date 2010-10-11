#!/usr/bin/env python

from tempfile import (NamedTemporaryFile, mkstemp)
import logging as log
import os

from NGSPipelineParts import (
    getOptions, referenceDatabase, seqs2Bams,
    mergeAndPileup, mkOutputDir, initLogger,
    defaultOptionsFile, runCommand)

def sequential(options):
    referenceDatabase(options)
    bams = seqs2Bams(options)
    mergeAndPileup(options, bams)

def referenceDBScript(options):
    return PBS_Script (
        #command = 'ReferenceDatabase.py -c %s' % options['optionsFile'],
        command = '/vlsci/VLSCI/bjpop/code/ngs_pipeline/src/python/ReferenceDatabase.py -c %s' % options['optionsFile'],
        #name = 'ReferenceDatabaseRun',
        memInGB = '10',
        walltime = '5:00:00',
        moduleList = ['python-gcc/2.6.4', 'bwa-gcc', 'samtools-gcc'],
    )

class PBS_Script(object):
    def __init__(self, command, queue='batch', name=None, memInGB=None, walltime=None, moduleList=None):
        self.command = command
        if queue in ['batch', 'smp']:
            self.queue = queue
        else:
            self.queue = 'batch'
        self.name = name
        self.memInGB = memInGB
        self.walltime = walltime
        self.moduleList = moduleList

    def __str__(self):
        script = ['#!/bin/bash']
        script.append('#PBS -q %s' % self.queue)
        if self.name:
            script.append('#PBS -N %s' % self.name)
        if self.memInGB:
            if self.queue == 'smp':
                script.append('#PBS -l mem=%sgb' % self.memInGB)
            else:
                script.append('#PBS -l pvmem=%sgb' % self.memInGB)
        if self.walltime:
            script.append('#PBS -l walltime=%s' % self.walltime)
        if type(self.moduleList) == list and len(self.moduleList) > 0:
            script.append('module load %s' % ' '.join(self.moduleList))
        script.append('cd $PBS_O_WORKDIR')
        script.append(self.command)
        return '\n'.join(script) + '\n'

    def launch(self):
        file = NamedTemporaryFile()
        file.write(str(self))
        file.flush()
        command = 'qsub ' + file.name
        log.info(command)
        commandHandle = os.popen(command)
        jobID = commandHandle.read()
        commandHandle.close()
        file.close()
        return jobID

def parallel(options):
    script = referenceDBScript(options)
    print script
    jobID = script.launch()
    print jobID
    #sequences = options['sequences']

def main():
    options = getOptions()
    outputDir = options['dir']
    mkOutputDir(outputDir)
    initLogger(options['logging'])
    control = options['control']
    if control == 'sequential':
        sequential(options)
    elif control == 'parallel':
        parallel(options)

if __name__ == '__main__':
    main()
