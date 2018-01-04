import commands
import os
import tarfile
import subprocess


PCHECK='voms-proxy-info -timeleft'

try:
	x = subprocess.check_output(PCHECK, shell=True)
except subprocess.CalledProcessError as xExc:
	print "Error Code ",xExc.returncode, grepexc.output

if int(x) == 0:
	print "Need to set proxy: voms-proxy-init -voms cms --valid 168:00"
	print "Please enter the voms-proxy-init password..."
	os.system("voms-proxy-init -voms cms --valid 168:00")

# Function to create a tar file
def make_tarfile(output_filename, source_dir):
    with tarfile.open(output_filename, "w:gz") as tar:
                tar.add(source_dir, arcname=os.path.basename(source_dir))

# Get CMSSW directory path and name
cmsswDirPath = commands.getstatusoutput('echo ${CMSSW_BASE}')
CMSSWRel = os.path.basename(cmsswDirPath[1])

print "CMSSW release used : ",CMSSWRel

# create tarball of present working CMSSW base directory
os.system('rm CMSSW*.tgz')
make_tarfile(CMSSWRel+".tgz", cmsswDirPath[1])
