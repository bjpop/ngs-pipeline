import logging

def initLogger(options):
    logFile = options['file']
    logLevel = options['level'] 
    levels = {
        'debug': logging.DEBUG,
        'info': logging.INFO,
        'warning': logging.WARNING,
        'error': logging.ERROR,
        'critical': logging.CRITICAL }
   if logFile:
        logging.basicConfig(filemode='w',filename=logFile,level=levels[logLevel])

def logInfo(msg):
    logging.info(msg)

def logWarn(msg):
    logging.warning(msg)
