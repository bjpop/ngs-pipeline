import subprocess

def shellCommand(command):
    process = subprocess.Popen(command, stdout = subprocess.PIPE,
                                stderr = subprocess.PIPE, shell = True)
    stdoutStr, stderrStr = process.communicate()
    return(stdoutStr, stderrStr, process.returncode)
