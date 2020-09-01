import os, sys
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(0)

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

input_root_file = '/eos/uscms/store/user/rasharma/double-higgs/SignalSample/GEN_reducedNtuples/'
InputFileList = ['GF_HH_Benchmark1', 'GF_HH_Benchmark2', 'GF_HH_Benchmark3', 'GF_HH_Benchmark4', 'GF_HH_Benchmark5', 'GF_HH_Benchmark6', 'GF_HH_Benchmark7', 'GF_HH_Benchmark8', 'GF_HH_Benchmark9', 'GF_HH_Benchmark10', 'GF_HH_Benchmark11', 'GF_HH_Benchmark12', 'GF_HH_BenchmarkSM']

input_root_file = '../'
InputFileList = ['Radion_hh_narrow_M1500']


# file = ROOT.TFile(input_root_file + os.sep + InputFileList[0]+'.root','READ')
# file = ROOT.TFile('../Radion_hh_narrow_M1500.root','READ')
file = ROOT.TFile('../Radion_hh_narrow_M500.root','READ')

if not file:
  print 'Failed to open %s' % (input_root_file + os.sep + InputFileList[0]+'.root')
  exit(0)

tree = file.Get("otree")

from twoDplots_info import varlists
from runnning_plots import variables
from plot_functions import getCanvas

key1 = 'ak8higgs_minDmass'
# key1 = 'oneAK8twoAK4_minMass'
# key1 = 'oneAK8twoAK4_pTmax'
key2 = 'all'
# CutToApply = 'AK8Gen_HiggsJet_minDMass_M>100 && '\
#               'AK8Gen_HiggsJet_minDMass_M<160 && '\
#               'AK8Gen_HiggsJet_minDMass_Pt>400'
#               # 'AK8Gen_HiggsJet_minDMass_deltaR_H1<1.0 && '\
#               # 'gen_leading_Higgs_Pt>400 '

# CutToApply = 'AK8Gen_HiggsJet_MaxPt_M>100 && '\
#               'AK8Gen_HiggsJet_MaxPt_M<160 && '\
#               'AK8Gen_HiggsJet_MaxPt_deltaR_H1<1.0 && '\
#               'AK8Gen_HiggsJet_MaxPt_Pt>400 && '\
#               'gen_leading_Higgs_Pt>400 '

# CutToApply = 'AK8Gen_MergedWjets_MaxPt_Leading_Pt>200 &&  AK8Gen_MergedWjets_MaxPt_SubLeading_Pt>200'
# CutToApply =  'AK4GEN_AllResolved_Higgs_M>100 && AK4GEN_AllResolved_Higgs_M<160 '\
#               ' && AK4GEN_AllResolved_onShellWboson_M > 60 && AK4GEN_AllResolved_onShellWboson_M<105 '\
#               ' && AK4GEN_AllResolved_offShellWboson_M < 60 '\
#               ' && gen_leading_WBoson_M > 60 && gen_leading_WBoson_M < 100 '\
#               ' && gen_Subleading_WBoson_M < 60 '\
#               ' && AK4GEN_AllResolved_onShellWboson_Pt < 200 '\
#               ' && AK4GEN_AllResolved_offShellWboson_Pt < 200 '\
#               ' && gen_Subleading_WBoson_Pt<200 '\
#               ' && gen_leading_WBoson_Pt < 200 '

CutToApply = ''
# CutToApply = 'OneAK8TwoAK4_pTMax_AK8_M<120 && OneAK8TwoAK4_pTMax_ReconsW_AK4_M<120'
# CutToApply = 'OneAK8TwoAK4_pTMax_AK8_M<120 && OneAK8TwoAK4_pTMax_ReconsW_AK4_M<120  && OneAK8TwoAK4_pTMax_AK8_Pt<200'
# CutToApply = 'OneAK8TwoAK4_pTMax_AK8_M<120 && OneAK8TwoAK4_pTMax_ReconsW_AK4_M<120  && OneAK8TwoAK4_pTMax_AK8_Pt<300 && OneAK8TwoAK4_pTMax_ReconsH_M>100'
# CutToApply = 'OneAK8TwoAK4_pTMax_AK8_Pt>200'
# CutToApply = 'OneAK8TwoAK4_minMass_ReconsH_M>100 && OneAK8TwoAK4_minMass_ReconsH_M<160 && OneAK8TwoAK4_minMass_ReconsW_AK4_M<100 && OneAK8TwoAK4_minMass_AK8_M<100 && OneAK8TwoAK4_minMass_AK8_Pt>40 && OneAK8TwoAK4_minMass_leadingAK4_Pt>30 && OneAK8TwoAK4_minMass_subleadingAK4_Pt>20'

for var in varlists:
  print(var)
  title = variables[key1][var[0]][1]
  nbin1 = variables[key1][var[0]][3]
  nbin2 = variables[key2][var[1]][3]
  xmin = variables[key1][var[0]][4]
  xmax = variables[key1][var[0]][5]
  ymin = variables[key2][var[1]][4]
  ymax = variables[key2][var[1]][5]
  print(title,nbin1,xmin,xmax, nbin2,ymin,ymax)

  c1 = getCanvas()
  hist = ROOT.TH2F("hist", ';'+var[0]+';'+var[1], nbin1, xmin, xmax, nbin2, ymin, ymax)
  print(style.RED+'tree.Draw('+var[0]+':'+var[1]+'">>hist",'+CutToApply+' , "+colz+"'+style.RESET)
  tree.Draw(var[0]+":"+var[1]+">>hist",CutToApply,'colz')
  print('Entries = ',hist.GetEntries())
  c1.SaveAs('twoD_'+key1+'_'+var[0]+'_'+var[1]+'.png')



