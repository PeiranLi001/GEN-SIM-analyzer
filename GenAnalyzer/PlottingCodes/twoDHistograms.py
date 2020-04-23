import os, sys
import ROOT

def getCanvas():
    H_ref = 600; 
    W_ref = 600; 
    W = W_ref
    H = H_ref

    T = 0.08*H_ref
    B = 0.12*H_ref 
    L = 0.18*W_ref
    R = 0.01*W_ref

    canvas = ROOT.TCanvas("c2","c2",50,50,W,H)
    canvas.SetFillColor(0)
    canvas.SetBorderMode(0)
    canvas.SetFrameFillStyle(0)
    canvas.SetFrameBorderMode(0)
    canvas.SetLeftMargin( L/W *0.6 )
    canvas.SetRightMargin( R/W*5 )
    canvas.SetTopMargin( T/H )
    canvas.SetBottomMargin( B/H )
    canvas.SetTickx(0)
    canvas.SetTicky(0) 
    return canvas

input_root_file = '/eos/uscms/store/user/rasharma/double-higgs/SignalSample/GEN_reducedNtuples/'
InputFileList = ['GF_HH_Benchmark1', 'GF_HH_Benchmark2', 'GF_HH_Benchmark3', 'GF_HH_Benchmark4', 'GF_HH_Benchmark5', 'GF_HH_Benchmark6', 'GF_HH_Benchmark7', 'GF_HH_Benchmark8', 'GF_HH_Benchmark9', 'GF_HH_Benchmark10', 'GF_HH_Benchmark11', 'GF_HH_Benchmark12', 'GF_HH_BenchmarkSM']

# file = ROOT.TFile(input_root_file + os.sep + InputFileList[0]+'.root','READ')
file = ROOT.TFile('../GF_HH_Benchmark11.root','READ')

if not file:
  print 'Failed to open %s' % (input_root_file + os.sep + InputFileList[0]+'.root')
  exit(0)

tree = file.Get("otree")

from twoDplots_info import varlists
from runnning_plots import variables
# import plot_functions as plotter
# from plot_functions import getCanvas

for var in varlists:
  # print(var)
  title = variables[var[0]][1]
  nbin1 = variables[var[0]][3]
  nbin2 = variables[var[1]][3]
  xmin = variables[var[0]][4]
  xmax = variables[var[0]][5]
  ymin = variables[var[1]][4]
  ymax = variables[var[1]][5]
  print(title,nbin1,xmin,xmax, nbin2,ymin,ymax)
  
  # c1 = getCanvas()
  # tree.Print()
  # hist = ROOT.TH2F("hist", "", nbin1, xmin, xmax, nbin2, ymin, ymax)
  # genJetAK8_MaxPt_deltaR_H1
  # genJetAK8_minDMass_deltaR_H1
  print('otree->Draw("'+var[0]+':'+var[1]+'","genJetAK8_MaxPt_deltaR_H1<0.8 && genJetAK8_MaxPt_M>100.","colz")')
  # tree.Draw('"'+var[0]+':'+var[1]+'>>hist")')
  # c1.SaveAs(var[0]+'_'+var[1]+'.png')
  


