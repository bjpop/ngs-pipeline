import os.path
import sys
import yaml
import subprocess
from ruffus.proxy_logger import *
import logging
import os
from shell_command import shellCommand
from cluster_job import (PBS_Script, runJobAndWait)

defaultOptionsFile = 'ngs_pipeline.opt'
defaultWalltime = None
defaultModules = []
defaultQueue = 'batch'
defaultMemInGB = None

def mkDir(dir):
    if not os.path.exists(dir):
        try:
           os.mkdir(dir, 0777)
        except IOError, e:
           sys.exit('%s\nFailed to make directory %s' % (e, dir))

def initLog(options):
    logDir = options['logging']['dir']
    logFile = os.path.join(logDir, options['logging']['file'])
    mkDir(logDir)

    loggerArgs={}
    loggerArgs["file_name"] = logFile
    loggerArgs["level"] = logging.DEBUG
    loggerArgs["rotating"] = True
    loggerArgs["maxBytes"]=20000
    loggerArgs["backupCount"]=10
    loggerArgs["formatter"]="%(asctime)s - %(message)s"

    (proxy, mutex) = make_shared_logger_and_proxy (setup_std_shared_logger, "NGS_pipeline", loggerArgs)
    return { 'proxy': proxy, 'mutex': mutex }

def distributedCommand(stage, comm, options):
    stageOptions = options['stages'][stage]
    time = stageOptions.get('walltime', defaultWalltime)
    mods = stageOptions.get('modules', defaultModules)
    q = stageOptions.get('queue', defaultQueue)
    mem = stageOptions.get('memInGB', defaultMemInGB)
    jobScript = PBS_Script(command=comm, walltime=time, name=stage, memInGB=mem, queue=q, moduleList=mods)
    #print "running job:"
    #print jobScript
    runJobAndWait(jobScript)

def runStage(stage, logger, options, *args):
    command = getCommand(stage, options)
    commandStr = command(*args)
    logInfo(stage + ': ' + commandStr, logger)
    if options['pipeline']['distributed']:
        distributedCommand(stage, commandStr, options)
    else:
        (stdoutStr, stderrStr, returncode) = shellCommand(commandStr)
        if returncode != 0:
            msg = ("Failed to run '%s'\n%s%sNon-zero exit status %s" %
                   (commandStr, stdoutStr, stderrStr, returncode))
            logInfo(msg, logger)


def getCommand(name, options):
    funStr = options['stages'][name]['command']
    return eval(funStr)

def logInfo(msg, logger):
    with logger['mutex']: logger['proxy'].info(msg)

def splitPath(path):
    (prefix, base) = os.path.split(path)
    (name, ext) = os.path.splitext(base)
    return (prefix, name, ext)

defaultLogging = {
   'dir': 'log',
   'file': 'NGS_pipeline.log',
}

defaultConfig = {
   'logging': defaultLogging,
#   'dir' : 'BWA',
   'reference' : None,
   'sequences' : [],
#   'control' : 'sequential',
   'optionsFile' : defaultOptionsFile
}

def validateOptions(options):
    reference = options['reference']
    if not reference:
        fail('One reference file must be specified')
    sequences = options['sequences']
    if len(sequences) == 0:
        fail('At least one sequence file must be specified')
    return options

def getOptionsFile():
    if len(sys.argv) > 1:
        return sys.argv[1]
    else:
        return defaultOptionsFile

def getOptions():
    try:
        file = open(getOptionsFile())
        contents = file.read()
        file.close()
    except IOError, e:
        fail('Could not open configuration file: %s' % optionsFile)
    newConfig = yaml.load(contents)
    updateDict(defaultConfig, newConfig)
    validateOptions(defaultConfig)
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
        else:
            old[k] = newVal
