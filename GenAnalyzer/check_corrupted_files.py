import os
from optparse import OptionParser
import ROOT
from ROOT import TFile

DEBUG = 1

def files_to_remove(files,dir):
  filelist_to_remove = []
  for file in files:
    try:
      tfile = TFile.Open(file);
    except:
      pass
    if tfile:
      if (tfile.IsZombie()):
        filelist_to_remove.append(file)
    else:
      print('File could not be opened, adding it to missing files')
      filelist_to_remove.append(file)
      
  if DEBUG: print(filelist_to_remove)
  return filelist_to_remove

def list_root(directory):
  flist = []
  flistWithPath = []
  for root, directories, filenames in os.walk(directory): 
    for filename in filenames: 
      if filename.endswith(".root"):
        if "inLHE" in filename: flist.append(filename)
        fileWithPath = os.path.join(root,filename)  # Get file name with path
        if "inLHE" in filename: flistWithPath.append(fileWithPath)
  return flist,flistWithPath


stageDir = "/eos/uscms/store/user/lnujj/rasharma/DoubleHiggs_NonResonant/"
present_output, present_output_WithPath =  list_root(stageDir)

print present_output
print "\n\n=====\n\n"
print present_output_WithPath
corrupted_files = files_to_remove(present_output_WithPath,stageDir)
print "\n\n=====\n\n"
print "corrupted_files:\n"
print corrupted_files
