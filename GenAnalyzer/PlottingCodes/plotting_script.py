import os, sys
from runnning_plots import variables
import runnning_plots2


# Group of Different functions for different styles
os.system("")
class style():
    BLACK = '\033[30m'
    RED = '\033[31m'
    GREEN = '\033[32m'
    YELLOW = '\033[33m'
    BLUE = '\033[34m'
    MAGENTA = '\033[35m'
    CYAN = '\033[36m'
    WHITE = '\033[37m'
    UNDERLINE = '\033[4m'
    RESET = '\033[0m'

# print(style.YELLOW + "Hello, World!")# input_root_file = '/eos/uscms/store/user/rasharma/double-higgs/SignalSample/GEN_reducedNtuples/'

# InputFileList = ['GF_HH_Benchmark1', 'GF_HH_Benchmark2', 'GF_HH_Benchmark3', 'GF_HH_Benchmark4', 'GF_HH_Benchmark5', 'GF_HH_Benchmark6', 'GF_HH_Benchmark7', 'GF_HH_Benchmark8', 'GF_HH_Benchmark9', 'GF_HH_Benchmark10', 'GF_HH_Benchmark11', 'GF_HH_Benchmark12', 'GF_HH_BenchmarkSM']
# InputFileList = ['GF_HH_Benchmark1', 'GF_HH_Benchmark7', 'GF_HH_Benchmark11', 'GF_HH_BenchmarkSM']

input_root_file = '../'
InputFileList = ['GF_HH_Benchmark3']

singleVariable = True
multipleVariable = False
multipleVariableDifferentFiles = False

if (singleVariable):
  for files in InputFileList:
    for key in variables:
      #print variables[key]
      #python make_aQGC_plots.py -n ../GF_HH_Benchmark222.root --grid -v gen_leading_photon_Pt gen_Subleading_photon_Pt --leg "Leading" "Sub Leading" --xmin 0 --xmax 500 --xlabel "Photon p_{T}" -o photon_pt.png --legPos "tr" --nbin 51
        print(style.RED+ '\n\n==> python make_aQGC_plots.py -n '+input_root_file +os.sep+files+'.root'+ ' --grid  -v "' + key + '" --leg "' + variables[key][1] +'" --xmin ' + str(variables[key][4]) + ' --xmax ' + str(variables[key][5]) + ' --xlabel "' + variables[key][1] + '" -o ' +'plots/'+files+'_'+ key + '.png  --legPos "tr"  --nbin ' + str(variables[key][3])+ style.RESET)
        toPlot = 'python make_aQGC_plots.py -n '+input_root_file +os.sep+files+'.root'+ ' --grid  -v "' + key + '" --leg "' + variables[key][1] +'" --xmin ' + str(variables[key][4]) + ' --xmax ' + str(variables[key][5]) + ' --xlabel "' + variables[key][1] + '" -o ' +'plots/'+files+'_'+ key + '.png  --legPos "tr"  --nbin ' + str(variables[key][3])
        os.system(toPlot)

if (multipleVariable):
  for var_arr in runnning_plots2.variables:
    varToAdd = ""
    legendToAdd = ""
    key = []
    for local_var in var_arr:
      key = local_var
      varToAdd += '"' + local_var + '" '
      legendToAdd += '"' + variables[key][1] + '" '
    #print varToAdd
    toPlot = 'python make_aQGC_plots.py -n '+input_root_file + ' --grid  -v ' + varToAdd + ' --leg ' + legendToAdd +' --xmin ' + str(variables[key][4]) + ' --xmax ' + str(variables[key][5]) + ' --xlabel "' + (variables[key][1]).replace('Leading','').replace('SubLeading','').replace('Sub','') + '" -o ' + key + '_multi.png  --legPos "tr"  --nbin ' + str(variables[key][3])
    print "\n\n","*"*51
    print toPlot
    os.system(toPlot)

if (multipleVariableDifferentFiles):
  fileListsToPlot = ""
  legendToAdd = ""
  for files in InputFileList:
    fileListsToPlot +=input_root_file + os.sep + files +'.root  '
    legendToAdd += files.replace('GF_HH_','') + ' '
  print fileListsToPlot
  print legendToAdd
  for key in variables:
    toPlot = 'python make_plots_diffrootfiles.py -n '+fileListsToPlot+ ' --grid  -v ' + key + ' --leg ' + legendToAdd +' --xmin ' + str(variables[key][4]) + ' --xmax ' + str(variables[key][5]) + ' --xlabel "' + (variables[key][1]).replace('Leading','').replace('SubLeading','').replace('Sub','') + '" -o ' + key + '_DiffBenchMark.png  --legPos "tr"  --nbin ' + str(variables[key][3])
    print "\n\n","*"*51
    print toPlot
    os.system(toPlot)
#if (twoDPlots):
#  
