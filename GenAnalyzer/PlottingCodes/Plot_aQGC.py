import ROOT as r
r.gROOT.SetBatch(True)
clock = r.TStopwatch()
clock.Start()

r.gStyle.SetOptStat(0)
import CMS_lumi, tdrstyle
tdrstyle.setTDRStyle() 


def createRatio(h1, h2, xlabel, col):
    h3 = h1.Clone("h3")
    h3.SetLineColor(col)
    h3.SetMarkerColor(col)
    h3.SetMarkerStyle(21)
    h3.SetTitle("")
    h3.SetMinimum(0.5)
    h3.SetMaximum(1.5)
    # Set up plot for markers and errors
    #h3.Sumw2()
    h3.SetStats(0)
    h3.Divide(h2)

    # Adjust y-axis settings
    y = h3.GetYaxis()
    y.SetTitle("ratio (SM/aQGC)")
    y.SetNdivisions(505)
    y.SetTitleSize(20)
    y.SetTitleFont(43)
    y.SetTitleOffset(1.55)
    y.SetLabelFont(43)
    y.SetLabelSize(15)

    # Adjust x-axis settings
    x = h3.GetXaxis()
    x.SetTitle(xlabel)
    x.SetTitleSize(20)
    x.SetTitleFont(43)
    x.SetTitleOffset(4.0)
    x.SetLabelFont(43)
    x.SetLabelSize(15)

    return h3

def createCanvasPads():
    c = r.TCanvas("c", "canvas", 600, 700)
    # Upper histogram plot is pad1
    pad1 = r.TPad("pad1", "pad1", 0.0, 0.3, 1.0, 1.0)
    pad1.SetBottomMargin(0)  # joins upper and lower plot
    pad1.SetLeftMargin(0.1)
    pad1.SetRightMargin(0.03)
    pad1.SetLogy(1)
    pad1.SetGridx()
    pad1.SetTickx(1)
    pad1.SetTicky(1)
    pad1.Draw()
    # Lower ratio plot is pad2
    c.cd()  # returns to main canvas before defining pad2
    pad2 = r.TPad("pad2", "pad2", 0.0, 0.00, 1.0, 0.3)
    #pad2.SetLogy(1)
    pad2.SetTopMargin(0)  # joins upper and lower plot
    pad2.SetBottomMargin(0.25)
    pad2.SetLeftMargin(0.1)
    pad2.SetRightMargin(0.03)
    pad2.SetGridx()
    pad2.SetTickx(1)
    pad2.SetTicky(1)
    pad2.Draw()

    return c, pad1, pad2

def getHistFromFile (tree, var, weight, nbin, xmin, xmax, xlabel, col):  
    hist = r.TH1F("hist", xlabel, nbin, xmin, xmax)
    tree.Draw(var+ ">>hist","LHEWeights["+str(weight)+"]","")
    #tree.Draw(var+ ">>hist","LHEWeights["+str(weight)+"]","",1000)
    if not hist:
        print 'Failed to get hist from file'
        exit(0)
    hist.GetXaxis().SetTitle(xlabel)
    ylabel = "Events / %s GeV" % int(hist.GetBinWidth(1))
    hist.GetYaxis().SetTitle(ylabel)
    hist.SetLineColor(col)
    hist.SetMarkerColor(col)
    hist.SetDirectory(r.gROOT) # detach "hist" from the file
    return hist

def CreatePlots(ttree, variables, outputFileName, xlabel="x-axis", nbin=20, minx=0, maxx=5000, nEvents=-1):
	colorArray = [1, 632, 416, 600, 400, 616, 432, 920, 800, 880, 820, 840, 860, 900, 940, 960, 980, 640, 660, 680, 700, 720, 740, 760]
	## Define 3 histogram for three different samples
	aQGC_pars = { "FS0": [(491, 490, 489, 488), ("SM", -20, -40, -60)],	
		      "FS1": [(570, 569, 568, 567), ("SM", -10, -20, -30)],
		      "FT0": [(995, 994, 993, 992), ("SM", -0.2, -0.4, -0.6)],	
		      "FT1": [(1055, 1054, 1053, 1052), ("SM", -0.5, -1.0, -1.5)],
		      "FT2": [(1122, 1121, 1120, 1119), ("SM", -0.5, -1.0, -1.5)],	
		      "FM0": [(646, 645, 644, 643), ("SM", -1, -2, -3)],
		      "FM1": [(722, 721, 720, 719), ("SM", -5, -10, -15)],	
		      "FM6": [(798, 797, 796, 795), ("SM", -2, -4, -6)],
		      "FM7": [(900, 899, 898, 897), ("SM", -5, -10, -15)]
	    }	
	for i,par1 in enumerate(aQGC_pars):
	   #leg = r.TLegend(0.6,0.70,0.9,0.95)
	   leg = r.TLegend(0.10,0.00,0.50,0.35)
	   hist1 = []
	   ratio = []
	   colorPick = 0
	   print i, par1, aQGC_pars[par1][0]
	   for j,par2 in enumerate(aQGC_pars[par1][0]):
	   	hist1.append(getHistFromFile(ttree, variables, par2, nbin, minx, maxx, xlabel, colorArray[colorPick]))
	   	#hist1[colorPick].Sumw2()
	   	hist1[colorPick].Scale(1/hist1[colorPick].Integral())
		if j == 0:
	   		leg.AddEntry(hist1[colorPick], "SM","lpe")
		else:
	   		leg.AddEntry(hist1[colorPick], par1+"= "+str(aQGC_pars[par1][1][j])+"x10^{-12} TeV^{-4}","lpe")
	   	ratio.append(createRatio(hist1[0],hist1[colorPick], xlabel, colorArray[colorPick]))
		colorPick+=1
	
	   canvas, pad1, pad2 = createCanvasPads()
	   
	   # Draw all hitogram using sames.
	   pad1.cd()
	   for a,list1 in enumerate(hist1):
	      if a == 0:
	         list1.Draw()
	      else:
	         list1.Draw("same")
	   leg.Draw()
	   
	   pad2.cd()

	   for a,list1 in enumerate(ratio):
	      if a != 0:
	         if a == 1:
	            list1.Draw("ep")
	         else:
	            list1.Draw("ep same")
	   
	   l = r.TLine( minx, 1, maxx, 1 )
	   l.SetLineStyle(3)
	   l.Draw("same")
	   
	   #canvas.Draw()
	   canvas.SaveAs(outputFileName+"_"+par1+".png")
	   canvas.SaveAs(outputFileName+"_"+par1+".pdf")



# Read three fils and their tree
tfile = r.TFile("root:://cmseos.fnal.gov//eos/uscms/store/user/rasharma/LHE_GEN_Analyzer_Output/WPlepWMhad_aQGC.root")

ttree = tfile.Get("otree")

# See which weights we are going to apply
for n, event in enumerate(ttree):
	if n>0:
		break
	print "Weight Details for OnlyFT1_initialize: \tWeight = ",event.LHEWeights[491]


#CreatePlots(ttree, variables, outputFileName, xlabel="x-axis", nbin=20, minx=0, maxx=5000, nEvents=5000):
CreatePlots(ttree, "LHE_MothInfo_mtWW", "mT_WW", "m_{T} (WW)", 20, 0, 5000, -1)
CreatePlots(ttree, "LHE_MothInfo_mjj", "LHE_MothInfo_mjj", "VBF di-jet invariant mass", 20, 0, 5000, -1)
CreatePlots(ttree, "LHE_MothInfo_dEtajj", "LHE_MothInfo_dEtajj", "#Delta #eta_{jj}", 20, 0, 10, -1)
CreatePlots(ttree, "LHE_MothInfo_mWW", "LHE_MothInfo_mWW", "Invariant mass of WW system", 20, 0, 5000, -1)

# Total 9 parameters:
# For each parameters 4 histograms

clock.Stop()

print 'Tree loop profiling stats:'
print 'Real Time used:', clock.RealTime(),"seconds"
print 'CPU Time used:', clock.CpuTime()
