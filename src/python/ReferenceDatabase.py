#!/usr/bin/env python

from NGSPipelineParts import (referenceDatabase, getOptions, initLogger)

def main():
    options = getOptions()
    initLogger(options['logging'])
    referenceDatabase(options)

if __name__ == '__main__':
    main()
