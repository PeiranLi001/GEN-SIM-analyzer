import os, sys
from runnning_plots import variables
import runnning_plots2
from datetime import datetime
current_datetime = datetime.now()


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

dirName =(str(current_datetime.year)[-2:]
         +str(format(current_datetime.month,'02d'))
         +str(format(current_datetime.day,'02d'))
         +"_"
         +str(format(current_datetime.hour,'02d'))
         +str(format(current_datetime.minute,'02d'))
         +str(format(current_datetime.second,'02d'))
         )
print dirName

# print(style.YELLOW + "Hello, World!")# input_root_file = '/eos/uscms/store/user/rasharma/double-higgs/SignalSample/GEN_reducedNtuples/'

# InputFileList = ['GF_HH_Benchmark1', 'GF_HH_Benchmark2', 'GF_HH_Benchmark3', 'GF_HH_Benchmark4', 'GF_HH_Benchmark5', 'GF_HH_Benchmark6', 'GF_HH_Benchmark7', 'GF_HH_Benchmark8', 'GF_HH_Benchmark9', 'GF_HH_Benchmark10', 'GF_HH_Benchmark11', 'GF_HH_Benchmark12', 'GF_HH_BenchmarkSM']
# InputFileList = ['GF_HH_Benchmark1', 'GF_HH_Benchmark7', 'GF_HH_Benchmark11', 'GF_HH_BenchmarkSM']

input_root_file = '../'
InputFileList = ['Radion_hh_narrow_M500_bqrk']
# InputFileList = ['Radion_hh_narrow_M500']

singleVariable = True
multipleVariable = False
multipleVariableDifferentFiles = False

# CutToApply = 'AK8Gen_HiggsJet_MaxPt_Pt>400 && AK8Gen_HiggsJet_MaxPt_M>100 && AK8Gen_HiggsJet_MaxPt_M<160 && AK8Gen_HiggsJet_MaxPt_deltaR_H1<1.0'
# CutToApply = 'AK8Gen_HiggsJet_minDMass_Pt>400 && AK8Gen_HiggsJet_minDMass_M>100 && AK8Gen_HiggsJet_minDMass_M<160 && AK8Gen_HiggsJet_minDMass_deltaR_H1<1.0'
# CutToApply = 'AK8Gen_MergedWjets_MaxPt_Higgs_M>100 && AK8Gen_MergedWjets_MaxPt_Higgs_M<160'
# CutToApply = 'AK8Gen_MergedWjets_MaxPt_Leading_Pt>200 && AK8Gen_MergedWjets_MaxPt_SubLeading_Pt>200'
# CutToApply = 'AK8Gen_MergedWjets_MaxPt_Leading_M<105 && AK8Gen_MergedWjets_MaxPt_SubLeading_M < 105 && AK8Gen_MergedWjets_MaxPt_Higgs_M>100 && AK8Gen_MergedWjets_MaxPt_Higgs_M<160'
# CutToApply = 'AK4GEN_AllResolved_Higgs_M>100 && AK4GEN_AllResolved_Higgs_M<160 && AK4GEN_AllResolved_onShellWboson_M > 60 && AK4GEN_AllResolved_onShellWboson_M<105'
# CutToApply = 'OneAK8TwoAK4_pTMax_ReconsH_M>100 && OneAK8TwoAK4_pTMax_ReconsH_M<160 && OneAK8TwoAK4_pTMax_AK8_Pt>100'
# CutToApply = 'OneAK8TwoAK4_pTMax_AK8_Pt>100'
# CutToApply = 'OneAK8TwoAK4_minMass_ReconsH_M>100 && OneAK8TwoAK4_minMass_ReconsH_M<160 && OneAK8TwoAK4_minMass_ReconsW_AK4_M<100 && OneAK8TwoAK4_minMass_AK8_M<100 && OneAK8TwoAK4_minMass_AK8_Pt>40 && OneAK8TwoAK4_minMass_leadingAK4_Pt>30 && OneAK8TwoAK4_minMass_subleadingAK4_Pt>20'
# CutToApply = 'OneAK8TwoAK4_minMass_AK8_M<120 && OneAK8TwoAK4_minMass_ReconsW_AK4_M<120 && OneAK8TwoAK4_minMass_AK8_Pt > 200 '
# CutToApply = 'OneAK8TwoAK4_pTMax_AK8_M<120 && OneAK8TwoAK4_pTMax_ReconsW_AK4_M<120  && OneAK8TwoAK4_pTMax_AK8_Pt<300 && OneAK8TwoAK4_pTMax_ReconsH_M>100'
CutToApply = ''
keyToUse = 'all'
# keyToUse = 'ak8higgs_minDmass'
# keyToUse = 'oneAK8twoAK4_minMass'
# keyToUse = 'oneAK8twoAK4_pTmax'


if (singleVariable):
  for files in InputFileList:
    for key2 in variables:
      # print variables[key]
      if key2 == keyToUse:
        print(style.BLUE+ '='*51 + '\n=' +key2 + '\n'+ '='*51 +style.RESET)
        os.system('mkdir -p '+key2 + os.sep + dirName)
        newDict = variables[key2]
        for key in newDict:
          # print key
          # print newDict[key]
          toPlot = 'python make_aQGC_plots.py --logy -n '+input_root_file +os.sep+files+'.root'+ ' --grid  -v "' + key + '" --leg "' + newDict[key][1] +'" --xmin ' + str(newDict[key][4]) + ' --xmax ' + str(newDict[key][5]) + ' --xlabel "' + newDict[key][1] + '" -o ' +key2 + os.sep + dirName+ os.sep +files+'_'+ key + '.png  --legPos "tr"  --nbin ' + str(newDict[key][3]) + '  --cut "' + CutToApply +'"'
          print(style.RED + '\n\n==> ' + toPlot + style.RESET)
          os.system(toPlot)

if (multipleVariable):
  newDict = variables[keyToUse]
  for var_arr in runnning_plots2.variables:
    varToAdd = ""
    legendToAdd = ""
    key = []
    for local_var in var_arr:
      key = local_var
      varToAdd += '"' + local_var + '" '
      legendToAdd += '"' + newDict[key][1] + '" '
    # print varToAdd
    # twoDPlots' + os.sep +keyToUse + os.sep + dirName + os.sep+ key + '_multi.png
    os.system('mkdir -p twoDPlots'+os.sep+keyToUse+os.sep+dirName)
    toPlot = 'python make_aQGC_plots.py -n '+input_root_file +InputFileList[0]+'.root'+ ' --grid  -v ' + varToAdd + ' --leg ' + legendToAdd +' --xmin ' + str(newDict[key][4]) + ' --xmax ' + str(newDict[key][5]) + ' --xlabel "' + (newDict[key][1]).replace('Leading','').replace('SubLeading','').replace('Sub','') + '" -o  twoDPlots' + os.sep +keyToUse + os.sep + dirName + os.sep+ key + '_multi.png  --legPos "tr"  --nbin ' + str(newDict[key][3]) + '  --cut "' + CutToApply  + '"'
    # +'" --logy'
    print "\n\n","*"*51
    print var_arr
    print(style.RED + toPlot + style.RESET)
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
