from yaml import load
from sys import argv

defaultBWA = {
   'align' : 'aln',
   'index' : 'bwtsw',
   'threads' : 1
}

defaultSamtools = {
   'maxReadDepth' : 1000000
}

defaultLogging = {
   'level': 'info',
   'file': 'pipeline.log'
}

defaultConfig = {
   'bwa' : defaultBWA,
   'samtools' : defaultSamtools,
   'logging': defaultLogging,
   'dir' : 'BWA'
}

def getConfig():
    if len(argv) <= 1:
        filename = 'ngs_pipeline.opt'
    else:
        filename = argv[1]
    try:
        file = open(filename)
        contents = file.read()
    except IOError, e:
        exit('Could not open configuration file: %s' % filename)
    newConfig = load(contents)
    defaultConfig.update(newConfig)
    return defaultConfig 
