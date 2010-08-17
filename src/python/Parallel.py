from tempfile import mkstemp
from os import system

scriptTemplate = '''
#!/bin/bash
#PBS -q smp
#PBS -l mem=10gb
#PBS -l walltime=%s
#PBS -m abe
#PBS -N %s 
%s
cd $PBS_O_WORKDIR
%s
%s
'''

def makeScript(walltime, jobName, dependencies, modules, command): 
    return scriptTemplate % (walltime, jobName, dependencies, modules, command)

def submitJob(jobName, moduleList, dependList, wallTime, command): 
    if dependList = []:
        dependStr = ''
    else:
        dependStr = '#PBS -W depend=afterok:%s' % ':'.join(dependIds)
    if moduleList = []:
        moduleStr = ''
    else:
        modileStr = 'module load %s' % ' '.join(moduleList)
   script = makeScript(wallTime, jobName, dependStr, moduleStr, command) 
   (handle, filename) = mkstemp()
   handle.write(script)
   qsub = "qsub " + filename 
   status = system(qsub) 
   return status
