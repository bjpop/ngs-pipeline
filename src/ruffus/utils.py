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
defaultWalltime = None # use the default walltime of the scheduler
defaultModules = []
defaultQueue = 'batch'
defaultMemInGB = None # use the default mem of the scheduler
defaultDistributed = False
defaultLogDir = 'log'
defaultLogFile = 'pipeline.log'
defaultStyle = 'run'
defaultProcs = 4
defaultPaired = False

stageDefaults = {
   'distributed': defaultDistributed,
   'walltime': defaultWalltime,
   'memInGB': defaultMemInGB,
   'modules': defaultModules,
   'queue': defaultQueue
}

pipeline = {
   'logDir': defaultLogDir,
   'logFile': defaultLogFile,
   'style': defaultStyle,
   'procs': defaultProcs,
   'paired': defaultPaired
}

defaultConfig = {
   'reference': None,
   'sequences': [],
   'optionsFile': defaultOptionsFile,
   'stageDefaults': stageDefaults,
   'pipeline': pipeline
}

def mkDir(dir):
    if not os.path.exists(dir):
        try:
           os.mkdir(dir, 0777)
        except IOError, e:
           sys.exit('%s\nFailed to make directory %s' % (e, dir))

def initLog(options):
    logDir = options['pipeline']['logDir']
    logFile = os.path.join(logDir, options['pipeline']['logFile'])
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

# Look for a particular option in a stage, otherwise return the result
def getStageOptions(options, stage, optionName):
    try:
        return options['stages'][stage][optionName]
    except KeyError:
        return options['stageDefaults'][optionName]

def distributedCommand(stage, comm, options):
    #stageOptions = options['stages'][stage]
    #time = stageOptions.get('walltime', defaultWalltime)
    #mods = stageOptions.get('modules', defaultModules)
    #queue = stageOptions.get('queue', defaultQueue)
    #mem = stageOptions.get('memInGB', defaultMemInGB)
    time = getStageOptions(options, stage, 'walltime')
    mods = getStageOptions(options, stage, 'modules')
    queue = getStageOptions(options, stage, 'queue')
    mem = getStageOptions(options, stage, 'memInGB')
    script = PBS_Script(command=comm, walltime=time, name=stage, memInGB=mem, queue=queue, moduleList=mods)
    print script
    runJobAndWait(script)

def runStage(stage, logger, options, *args):
    command = getCommand(stage, options)
    commandStr = command(*args)
    logInfo(stage + ': ' + commandStr, logger)
    if getStageOptions(options, stage, 'distributed'):
        distributedCommand(stage, commandStr, options)
    else:
        (stdoutStr, stderrStr, returncode) = shellCommand(commandStr)
        if returncode != 0:
            msg = ("Failed to run '%s'\n%s%sNon-zero exit status %s" %
                   (commandStr, stdoutStr, stderrStr, returncode))
            logInfo(msg, logger)


def getCommand(name, options):
    try:
        funStr = options['stages'][name]['command']
        return eval(funStr)
    except KeyError:
        fail("command: %s, is not defined in the options file" % name)

def logInfo(msg, logger):
    with logger['mutex']: logger['proxy'].info(msg)

def splitPath(path):
    (prefix, base) = os.path.split(path)
    (name, ext) = os.path.splitext(base)
    return (prefix, name, ext)

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
