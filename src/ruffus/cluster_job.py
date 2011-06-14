from shell_command import shellCommand
import sys
from time import sleep
from tempfile import NamedTemporaryFile

def isJobCompleted(jobID):
    (stdout, stderr, returnCode) = shellCommand("qstat %s" % jobID)
    if returnCode != 0:
        return True
    else:
        try:
           lines = stdout.split('\n')
           statusLine = lines[2]
           statusVal = statusLine.split()[4]
        except:
           return "bad result from qstat"
        return statusVal is 'C'

def waitForJobCompletion(jobID):
    while(not isJobCompleted(jobID)):
        sleep(10)

def runJobAndWait(script):
    #print jobScript
    jobID = script.launch()
    #print jobID
    waitForJobCompletion(jobID)

class PBS_Script(object):
    def __init__(self, command, walltime=None, name=None, memInGB=None, queue='batch', moduleList=None, logDir=None):
        self.command = command
        if queue in ['batch', 'smp']:
            self.queue = queue
        else:
            self.queue = 'batch'
        self.name = name
        self.memInGB = memInGB
        self.walltime = walltime
        self.moduleList = moduleList
        self.logDir = logDir

    def __str__(self):
        script = ['#!/bin/bash']
        # XXX fixme 
        # should include job id in the output name.
        # should use the proper log directory.
        script.append('#PBS -q %s' % self.queue)
        if self.logDir:
           script.append('#PBS -o %s' % self.logDir)
           script.append('#PBS -e %s' % self.logDir)
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
        (stdout, stderr, returnCode) = shellCommand(command)
        file.close()
        if returnCode == 0:
            return stdout
        else:
            raise(Exception('qsub command failed with exit status: ' + str(returnCode)))
