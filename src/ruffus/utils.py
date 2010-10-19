import os.path
import sys
import yaml
import subprocess
from ruffus.proxy_logger import *
import logging
import os

defaultOptionsFile = 'ngs_pipeline.opt'

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

    (proxy, mutex) = make_shared_logger_and_proxy (setup_std_shared_logger,
                                                "NGS_pipeline", loggerArgs)
    return { 'proxy': proxy, 'mutex': mutex }

def shellCommand(command):
    process = subprocess.Popen(command, stdout = subprocess.PIPE,
                                stderr = subprocess.PIPE, shell = True)
    stdoutStr, stderrStr = process.communicate()
    return(stdoutStr, stderrStr, process.returncode)

def runCommand(command, logger):
    (stdoutStr, stderrStr, returncode) = shellCommand(command)
    if returncode != 0:
        msg = ("Failed to run '%s'\n%s%sNon-zero exit status %s" %
               (command, stdoutStr, stderrStr, process.returncode))
        logInfo(msg, logger)
    logInfo(command, logger)

def logInfo(msg, logger):
    with logger['mutex']: logger['proxy'].info(msg)

def splitPath(path):
    (prefix, base) = os.path.split(path)
    (name, ext) = os.path.splitext(base)
    return (prefix, name, ext)

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
   'dir': 'log',
   'file': 'NGS_pipeline.log',
}

defaultConfig = {
   'bwa' : defaultBWA,
   'samtools' : defaultSamtools,
   'logging': defaultLogging,
   'dir' : 'BWA',
   'reference' : None,
   'sequences' : [],
   'control' : 'sequential',
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
