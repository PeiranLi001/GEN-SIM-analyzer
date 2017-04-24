
# coding: utf-8

# # Make $m_{TT}$ plot for LL+LT+TT normalized to cross-section

# In[1]:

from ConfigParser import RawConfigParser
config = RawConfigParser()   
config.optionxform = str       # Last two lines are done because ConfigParser will not preserve case
config.read("DataMCInfo.ini")


# In[2]:

crossSections = dict([sample, float(xsec)] for sample, xsec in config.items('CrossSection'))
nProcessed    = dict([sample, int(nPro)] for sample, nPro in config.items('nEventsProcessed'))

#from pprint import pprint
#print "cross sections:" 
#pprint(crossSections)
#print "number of events processed:"
#pprint(nProcessed)


# In[3]:

weights = {}
luminosity = 36000.0    # This is just an example value
for sample in crossSections.keys():
    weights[sample] = luminosity * crossSections[sample]/nProcessed[sample]
    
#pprint(weights)    
#print "Weight for LL = ",weights["Signal_LL"]


# In[8]:

import ROOT as r
r.gROOT.SetBatch(True)

tfile_LL = r.TFile('root:://cmseos.fnal.gov//eos/uscms/store/user/rasharma/LHE_GEN_Analyzer_Output_13April/Signal_LL.root')
tfile_LT = r.TFile('root:://cmseos.fnal.gov//eos/uscms/store/user/rasharma/LHE_GEN_Analyzer_Output_13April/Signal_LT.root')
tfile_TT = r.TFile('root:://cmseos.fnal.gov//eos/uscms/store/user/rasharma/LHE_GEN_Analyzer_Output_13April/Signal_TT.root')
tfile_SM = r.TFile('root:://cmseos.fnal.gov//eos/uscms/store/user/rasharma/LHE_GEN_Analyzer_Output/aQGC_WPlepWMhadJJ_EWK_SM_mjj200_Pythia8CUEP8M1_13TeV_Madgraph.root')
tree_LL = tfile_LL.Get("otree")
tree_LT = tfile_LT.Get("otree")
tree_TT = tfile_TT.Get("otree")
tree_SM = tfile_TT.Get("otree")

print tree_LL.GetEntries()  


# In[ ]:

clock = r.TStopwatch()
clock.Start()

#for branch in tree_LL.GetListOfBranches():
#  print branch.GetName()

  
r.gStyle.SetOptStat(0)
#
canvas = r.TCanvas() 
canvas.SetLogy(1)
canvas.SetTickx(1)
canvas.SetTicky(1) 
canvas.SetGrid()
#
leg = r.TLegend(0.7,0.6,0.9,0.9)
#
h_LL_mtWW = r.TH1F("h_LL_mtWW","",30,0,6000)
h_LL_mtWW.SetLineColor(1)
h_LL_mtWW.SetMarkerColor(1)
#leg.AddEntry(h_LL_mtWW,"LL","lpe")
h_LT_mtWW = r.TH1F("h_LT_mtWW","",30,0,6000)
h_LT_mtWW.SetLineColor(2)
h_LT_mtWW.SetMarkerColor(2)
#leg.AddEntry(h_LT_mtWW,"LT","lpe")
h_TT_mtWW = r.TH1F("h_TT_mtWW","",30,0,6000)
h_TT_mtWW.SetLineColor(3)
h_TT_mtWW.SetMarkerColor(3)
#leg.AddEntry(h_TT_mtWW,"TT","lpe")
#h_LL_mtWW.SetMaximum(1.0)
#h_LT_mtWW.SetMaximum(1.0)
#h_TT_mtWW.SetMaximum(1.0)
#h_LL_mtWW.SetMinimum(0.0)

tree_LL.Draw("LHE_MothInfo_mtWW>>h_LL_mtWW","")
tree_LT.Draw("LHE_MothInfo_mtWW>>h_LT_mtWW","")
tree_TT.Draw("LHE_MothInfo_mtWW>>h_TT_mtWW","")


h_LL_mtWW.Scale(weights["Signal_LL"])
h_LT_mtWW.Scale(weights["Signal_LT"])
h_TT_mtWW.Scale(weights["Signal_TT"])
#h_LL_mtWW.GetYaxis().SetRangeUser(0.0,1.0)
#h_LT_mtWW.GetYaxis().SetRangeUser(0.0,1.0)
#h_TT_mtWW.GetYaxis().SetRangeUser(0.0,1.0)
#h_LL_mtWW.Draw()
#h_LT_mtWW.Draw("same")
#h_TT_mtWW.Draw("same")

h_add = r.TH1F("h_add","",30,0,6000)
h_add.Sumw2()
h_add.GetXaxis().SetTitle("mt (WW) GeV")
h_add.Add(h_LL_mtWW)
h_add.Add(h_LT_mtWW)
h_add.Add(h_TT_mtWW)
h_add.SetLineColor(4)
h_add.SetMarkerColor(4)
leg.AddEntry(h_add,"LL+LT+TT","lpe")
h_add.Draw()
print "total integral = ",h_add.Integral()
print "Integral (2000,6000) = ",h_add.Integral(h_add.GetXaxis().FindBin(2000),h_add.GetXaxis().FindBin(6000))
#h_add.Draw("same")

leg.Draw()

canvas.Draw()
canvas.SaveAs("mt_plot_log.png")
canvas.SaveAs("mt_plot_log.pdf")

clock.Stop()


# In[ ]:

print 'Tree loop profiling stats:'
#print '# entries examined:', evaluated
#print 'Real Time used:', clock.RealTime(), '(per event: %.5f)' % (clock.RealTime()/evaluated)
print 'Real Time used:', clock.RealTime()
print 'CPU Time used:', clock.CpuTime()
