#!/usr/bin/env python
import ROOT as ROOT
import CMS_lumi, tdrstyle
import argparse

def PrintaQGC_parameters(plot_info, aQGC_par):
    	file = ROOT.TFile(plot_info["file_name"][0])
    	if not file:
    	    print 'Failed to open %s' % plot_info["file_name"][0]
    	    exit(0)
    	tree = file.Get(plot_info["tree_folder"] + plot_info["tree_name"])
	FT0_key = []
	FT0_val = []
	for n, event in enumerate(tree):
		if n>0:
			break
		indices1 = [i for i, s in enumerate(event.LHEWeightIDs) if aQGC_par in s]
	
		#for i in range(0,len(indices1), len(indices1)/6):
		for i in range(0,len(indices1)):
			print indices1[i],"\t",event.LHEWeightIDs[indices1[i]]

def getaQGC_parameters(plot_info, aQGC_par):
    	file = ROOT.TFile(plot_info["file_name"][0])
    	if not file:
    	    print 'Failed to open %s' % plot_info["file_name"][0]
    	    exit(0)
    	tree = file.Get(plot_info["tree_folder"] + plot_info["tree_name"])
	FT0_key = []
	FT0_val = []
	for n, event in enumerate(tree):
		if n>0:
			break
		indices1 = [i for i, s in enumerate(event.LHEWeightIDs) if aQGC_par in s]
	
		#for i in range(0,len(indices1), len(indices1)/6):
		for i in range(len(indices1)/2,len(indices1)):
			FT0_key.append(indices1[i])
			FT0_val.append(event.LHEWeightIDs[indices1[i]])
		indices1s = [i for i, s in enumerate(event.LHEWeightIDs) if aQGC_par.replace("_m","_0p0") in s]
		FT0_key.append(indices1s[0])
		FT0_val.append(event.LHEWeightIDs[indices1s[0]]+" (SM)")
	return FT0_key,FT0_val

	
def aQGC_GetHisto(plot_info, key1, ColNum):
    file1 = ROOT.TFile(plot_info["file_name"][0])
    if not file1:
        print 'Failed to open %s' % plot_info["file_name"][0]
        exit(0)
    tree1 = file1.Get(plot_info["tree_folder"] + plot_info["tree_name"])
    hist1 = ROOT.TH1F("hist1", "Test", plot_info["nbin"], plot_info["xmin"], plot_info["xmax"])    
    tree1.Draw(plot_info["tree_var"][0] + ">>hist1",plot_info["weight1"],"")
    hist1.SetDirectory(0)
    hist1.Scale(1/hist1.Integral())
    #colorArray = ["ROOT.kBlack", "ROOT.kGray", "ROOT.KRed", "ROOT.kGreen", "ROOT.kBlue", "ROOT.kYellow", "ROOT.kMagenta", "ROOT.kCyan", "ROOT.Orange", "ROOT.kViolet", "ROOT.kSpring", "ROOT.kTeal", "ROOT.kAzure", "ROOT.kPink"]
    #colorArray = [ROOT.kBlack, ROOT.kGray, ROOT.KRed, ROOT.kGreen, ROOT.kBlue, ROOT.kYellow, ROOT.kMagenta, ROOT.kCyan, ROOT.Orange, ROOT.kViolet, ROOT.kSpring, ROOT.kTeal, ROOT.kAzure, ROOT.kPink]
    colorArray = [1, 920, 632, 416, 600, 400, 616, 432, 800, 880, 820, 840, 860, 900, 940, 960, 980, 640, 660, 680, 700, 720, 740, 760]
    setHistAttributes(hist1, plot_info, colorArray[ColNum],0)
    return hist1;
	
def aQGC_plotting (plot_info, aQGC_key, aQGC_val, outputNameString, skip):
    print "size of aQGC parameter: ",len(aQGC_key)
    print "==> Key: ",aQGC_key
    print "==> Val: ",aQGC_val
    if len(aQGC_key) == 0:
        print 'Error: No aQGC parameter for this input root file'
        exit(0)
    	
    c2 = getCanvas()
    legend = ROOT.TLegend(.83 ,.20 ,1.0 ,.930)
    legend.SetFillColor(ROOT.kWhite)
    legend.SetTextSize(0.04)

    hist1 = []
    new_aQGC_val = []
    for i in range(0,len(aQGC_key),skip):
	print "Working on hist number : ",i
    	hist1.append(aQGC_GetHisto(plot_info,aQGC_key[i],i/skip))
	new_aQGC_val.append(aQGC_val[i])
    for a,list1 in enumerate(hist1):
    	#print a,list1
	if new_aQGC_val[a].find("SM") == -1:
    		legend.AddEntry(list1, new_aQGC_val[a].replace("_m"," = -").replace("p",".").replace("_0"," = 0"),"lpe")
	else:
    		legend.AddEntry(list1, "SM", "lpe")
	if i==0:
		list1.Draw()
	else:
		list1.Draw("sames")

    legend.Draw("same")
    if plot_info["logy"]:
        c2.SetLogy()
    if plot_info["logx"]:
        c2.SetLogx()
    if plot_info["grid"]:
    	c2.SetGrid()
    #setTDRStyle(c2, 1, 13, plot_info["printCMS"]) 
    setTDRStyle(c2, 1, 13, "No") 
    hist1[0].GetXaxis().SetTitle(plot_info["xlabel"])
    if plot_info["ylabel"] == "":
        plot_info["ylabel"] = "Events / %s GeV" % int(hist1[0].GetBinWidth(1))
    hist1[0].GetYaxis().SetTitle(plot_info["ylabel"])
    c2.SaveAs(str(outputNameString)+plot_info["output_file"])
    #c2.SaveAs("test.pdf");

def CompHistFromTwoBranchSameFile (plot_info):
    c1 = ROOT.TCanvas()
    legend = ROOT.TLegend(.7 ,.70 ,.985 ,.950)
    legend.SetFillColor(ROOT.kWhite)
    file1 = ROOT.TFile(plot_info["file_name"][0])
    if not file1:
        print 'Failed to open %s' % plot_info["file_name"][0]
        exit(0)
    tree1 = file1.Get(plot_info["tree_folder"] + plot_info["tree_name"])
    hist1 = []
    new_aQGC_val = []
    for i in range(0,len(plot_info["tree_var"])):
	print "Working on hist number : ",i
    	hist1.append(getHistFromFile(plot_info,0,i,i))
	new_aQGC_val.append(plot_info["tree_var"][i])
    for a,list1 in enumerate(hist1):
    	list1.GetXaxis().SetTitle(plot_info["xlabel"])
    	if plot_info["ylabel"] == "":
    	    plot_info["ylabel"] = "Events / %s GeV" % int(list1.GetBinWidth(1))
    	list1.GetYaxis().SetTitle(plot_info["ylabel"])
    	list1.Scale(1/list1.Integral())
    	#print a,list1
	#if new_aQGC_val[a].find("SM") == -1:
    	#legend.AddEntry(list1, new_aQGC_val[a].replace("_m"," = -").replace("p",".").replace("_0"," = 0"),"lpe")
    	legend.AddEntry(list1, new_aQGC_val[a],"lpe")
	#else:
    	#	legend.AddEntry(list1, "SM", "lpe")
	if i==0:
		list1.Draw()
	else:
		list1.Draw("sames")
	
    legend.Draw()
    if plot_info["logy"]:
        c1.SetLogy()
    if plot_info["logx"]:
        c1.SetLogx()
    if plot_info["grid"]:
    	c1.SetGrid()
    c1.SaveAs(plot_info["output_file"])
    #if not hist1:
    #    print 'Failed to get hist from file'
    #    exit(0)
    #hist1.SetDirectory(ROOT.gROOT) # detach "hist" from the file
    #return hist1
def CompHistFromTwoFile (plot_info):
    c1 = ROOT.TCanvas()
    legend = ROOT.TLegend(.6 ,.70 ,.885 ,.875)
    legend.SetFillColor(ROOT.kWhite)
    file1 = ROOT.TFile(plot_info["file_name"][0])
    file2 = ROOT.TFile(plot_info["file_name"][1])
    #print "=======================\n\n"
    #print plot_info["file_name"][0],"\t",plot_info["leg1"]
    #print plot_info["file_name"][1],"\t",plot_info["leg2"]
    #print "\n\n======================="
    if not file1:
        print 'Failed to open %s' % plot_info["file_name"][0]
        exit(0)
    if not file2:
        print 'Failed to open %s' % plot_info["file_name"][1]
        exit(0)
    tree1 = file1.Get(plot_info["tree_folder"] + plot_info["tree_name"])
    tree2 = file2.Get(plot_info["tree_folder"] + plot_info["tree_name"])
    hist1 = ROOT.TH1F("hist1", "Test", plot_info["nbin"], plot_info["xmin"], plot_info["xmax"])    
    hist2 = ROOT.TH1F("hist2", "Test", plot_info["nbin"], plot_info["xmin"], plot_info["xmax"])    
    tree1.Draw(plot_info["tree_var"][0] + ">>hist1",plot_info["weight1"])
    #tree1.Draw(plot_info["tree_var"][0] + ">>hist1","")
    legend.AddEntry(hist1, plot_info["leg1"])

    tree2.Draw(plot_info["tree_var"][0] + ">>hist2",plot_info["weight2"])
    legend.AddEntry(hist2, plot_info["leg2"])

    hist1.Sumw2()
    hist2.Sumw2()
    hist1.Scale(1/hist1.Integral())
    hist2.Scale(1/hist2.Integral())
    setHistAttributes(hist1, plot_info, ROOT.kBlack,0)
    setHistAttributes(hist2, plot_info, ROOT.kRed,0)
    hist1.Draw()
    hist2.Draw("sames")
    legend.Draw("same")
    if plot_info["logy"]:
        c1.SetLogy()
    if plot_info["logx"]:
        c1.SetLogx()
    if plot_info["grid"]:
    	c1.SetGrid()
    #setTDRStyle(c1, 1, 13, plot_info["printCMS"]) 
    setTDRStyle(c1, 1, 13, "No") 
    hist1.GetXaxis().SetTitle(plot_info["xlabel"])
    if plot_info["ylabel"] == "":
        plot_info["ylabel"] = "Events / %s GeV" % int(hist1.GetBinWidth(1))
    hist1.GetYaxis().SetTitle(plot_info["ylabel"])
    c1.SaveAs(plot_info["output_file"])
    #if not hist1:
    #    print 'Failed to get hist from file'
    #    exit(0)
    #hist1.SetDirectory(ROOT.gROOT) # detach "hist" from the file
    #return hist1
def getHistFromFile (plot_info,i,j,k):
    file = ROOT.TFile(plot_info["file_name"][i])
    print "File Name: ",plot_info["file_name"][i]
    print "Tree Folder: ",plot_info["tree_folder"]
    if not file:
        print 'Failed to open %s' % plot_info["file_name"][i]
        exit(0)
    tree = file.Get(plot_info["tree_folder"] + plot_info["tree_name"])
    hist = ROOT.TH1F("hist", "Test", plot_info["nbin"], plot_info["xmin"], plot_info["xmax"])    
    tree.Draw(plot_info["tree_var"][j] + ">>hist")
    if not hist:
        print 'Failed to get hist from file'
        exit(0)
    colorArray = [1, 920, 632, 416, 600, 400, 616, 432, 800, 880, 820, 840, 860, 900, 940, 960, 980, 640, 660, 680, 700, 720, 740, 760]
    setHistAttributes(hist, plot_info, colorArray[k],0)
    hist.SetDirectory(ROOT.gROOT) # detach "hist" from the file
    return hist
def setHistAttributes (hist, plot_info, line_color, fill_color):
    #hist.SetFillColor(fill_color)
    hist.SetLineColor(line_color)
    hist.SetMarkerColor(line_color)
    #hist.SetLineWidth(2)
    if plot_info["rebin"] != 0:
        if type(hist) != "<class '__main__.TH2F'>":
            hist.Rebin(plot_info["rebin"])
        else:
            print 'Rebin only defined for 1D hist. Use --rebin2D instead.'
    #hist.Draw()
    if plot_info["xmin"] < plot_info["xmax"]:
        hist.GetXaxis().SetRangeUser(plot_info["xmin"], plot_info["xmax"])
    if plot_info["ymin"] < plot_info["ymax"]:
        hist.GetYaxis().SetRangeUser(plot_info["ymin"], plot_info["ymax"])
    #hist.GetYaxis().SetTitle(plot_info["ylabel"])
def addHistToStack (hist_stack, plot_info, hist_opts, line_color, fill_color):
    hist = getHistFromFile(plot_info)
    setHistAttributes(hist, plot_info, line_color, fill_color)
    hist_stack.Add(hist, hist_opts)
def makePlot (hist, hist_opts, plot_info):
    #legend = ROOT.TLegend(.5 ,.65 ,.885 ,.875)
    canvas = getCanvas()
    setTDRStyle(canvas, 1, 13, plot_info["printCMS"]) 
    if plot_info["logy"]:
        canvas.SetLogy()
    if plot_info["logx"]:
        canvas.SetLogx()
    #draw the lumi text on the canvas
    hist.Draw(hist_opts)
    setTDRStyle(canvas, 1, 13, plot_info["printCMS"]) 
    hist.GetXaxis().SetTitle(plot_info["xlabel"])
    if plot_info["ylabel"] == "":
        plot_info["ylabel"] = "Events / %s GeV" % int(hist.GetBinWidth(1))
    #hist.GetYaxis().SetTitle(plot_info["ylabel"])
    #hist.SetTitleOffset(1.3, "y")
    #hist.SetTitleOffset(1.1, "x")
    canvas.cd()
    canvas.Update()
    canvas.RedrawAxis()
    #frame = canvas.GetFrame()
    #frame.Draw()
    #legend.SetFillColor(ROOT.kWhite)
    #legend.AddEntry(hist, legendName)

    #legend.Draw("same")
    canvas.Print(plot_info["output_file"]) 
def makeStackPlots (stacked, unstacked, hist_opts, plot_info):
    #legend = ROOT.TLegend(.5 ,.65 ,.885 ,.875)
    canvas = getCanvas()
    if plot_info["logy"]:
        canvas.SetLogy()
    if plot_info["logx"]:
        canvas.SetLogx()
    #draw the lumi text on the canvas
    hist1.Draw(hist_opts)
    hist1.GetXaxis().SetTitle(plot_info["xlabel"])
    if plot_info["ylabel"] == "":
        plot_info["ylabel"] = "Events / %s GeV" % int(hist.GetBinWidth(1))
    hist1.GetYaxis().SetTitle(plot_info["ylabel"])
    hist1.SetTitleOffset(1.3, "y")
    hist.SetTitleOffset(1.1, "x")
    setTDRStyle(canvas, 1, 13, plot_info["printCMS"]) 
    canvas.cd()
    canvas.Update()
    canvas.RedrawAxis()
    #frame = canvas.GetFrame()
    #frame.Draw()
    #legend.SetFillColor(ROOT.kWhite)
    #legend.AddEntry(hist, legendName)

    #legend.Draw("same")
    canvas.Print(plot_info["output_file"]) 

def setTDRStyle(canvas, luminosity, energy, printCMS):
    tdrstyle.setTDRStyle() 
    if printCMS == "right" or printCMS == "left":
        if energy == 13:
            CMS_lumi.lumi_13TeV = "%s fb^{-1}" % luminosity
            if printCMS == "left":
                iPos = 11
            else:
                iPos = 13
            CMS_lumi.writeExtraText = 1
            CMS_lumi.extraText = "Preliminary"
            CMS_lumi.CMS_lumi(canvas, 4, iPos)
def getCanvas():
    H_ref = 600; 
    W_ref = 900; 
    W = W_ref
    H  = H_ref

    T = 0.08*H_ref
    B = 0.12*H_ref 
    L = 0.12*W_ref
    R = 0.04*W_ref

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
def getBasicParser():
    parser = argparse.ArgumentParser()
    #parser.add_argument('-n', '--file_name', type=str, required=True, 
    #                    help="Name of root file in which histogram is stored.")
    parser.add_argument('-o', '--output_file', type=str, required=True,
                        help="Name of produced plot file (type pdf/png/jpg etc.).")
    parser.add_argument('-t', '--tree_name', type=str, required=False, default="otree",
                        help="Name of root tree")  
    parser.add_argument('-f', '--tree_folder', type=str, required=False, default="",
                        help="Folder where tree is stored")  
    parser.add_argument('-v', '--tree_var', nargs='+', type=str, required=False, default="LHE_mWW",
                        help="Variable name in root tree")  
    parser.add_argument('--xlabel', type=str, required=False, default="Invariant mass of WW", 
                        help="x axis label")
    parser.add_argument('--ylabel', type=str, required=False, default="", 
                        help="y axis label")
    parser.add_argument('--leg1', type=str, required=False, default="Legend 1", 
                        help="First Histo/Rootfile legend")
    parser.add_argument('--leg2', type=str, required=False, default="Legend 2", 
                        help="Second Histo/Rootfile legend")
    parser.add_argument('--weight1', type=str, required=False, default="", 
                        help="Weight corresponding to variable")
    parser.add_argument('--weight2', type=str, required=False, default="", 
                        help="Weight corresponding to variable")
    parser.add_argument('--xmin', type=float, required=False, default=0, 
                        help="minimum x value")
    parser.add_argument('--xmax', type=float, required=False, default=600, 
                        help="maximum x value")   
    parser.add_argument('--ymin', type=float, required=False, default=0, 
                        help="minimum y value")
    parser.add_argument('--ymax', type=float, required=False, default=0, 
                        help="maximum y value")
    parser.add_argument('--nbin', type=int, required=False, default=70, 
                        help="Number of bins to group together (1D only)")
    parser.add_argument('--rebin', type=int, required=False, default=0, 
                        help="Number of bins to group together (1D only)")
    parser.add_argument('--logy', action='store_true',
                        help="Set y axis to logarithmic scale")
    parser.add_argument('--logx', action='store_true', 
                        help="Set x axis to logarithmic scale")
    parser.add_argument('--grid', action='store_true', 
                        help="Set grid for x-y axis")
    parser.add_argument('--printCMS', type=str, default="left",required=False,
                        choices=["left","right"], help="""print 'CMS preliminary' 
                        in left (or right) upper corner""")
    return parser
