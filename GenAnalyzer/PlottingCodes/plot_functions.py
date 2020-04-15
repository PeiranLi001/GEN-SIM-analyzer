#!/usr/bin/env python
import ROOT as ROOT
import CMS_lumi, tdrstyle
import argparse


#############	Define Legend, Canvas, 	##########################

def getLegends(pos,ncol,nvar,fontsize):
    if pos == "tr":
    	legend = ROOT.TLegend(.70-(0.15*(ncol-1)), 0.80-(0.02*(nvar/ncol-1)) ,.95 ,.920)
    elif pos == "tl":
    	legend = ROOT.TLegend(0.11, 0.85-(0.02*(nvar/ncol-1)) ,.80+(0.15*(ncol-1)) ,.920)
    elif pos == "tc":
    	legend = ROOT.TLegend(0.30, 0.80-(0.02*(nvar/ncol-1)) ,.80 ,.920)
    elif pos == "bl":
    	legend = ROOT.TLegend(0.11, 0.25+(0.02*(nvar/ncol-1)) ,.60+(0.15*(ncol-1)) ,.130)
    elif pos == "bc":
    	legend = ROOT.TLegend(0.30, 0.25+(0.02*(nvar/ncol-1)) ,.80 ,.130)
    elif pos == "br":
    	legend = ROOT.TLegend(.45-(0.15*(ncol-1)), 0.25+(0.02*(nvar/ncol-1)) ,.95 ,.130)
    else:
    	print "Invalid default position: Switching to default legend position top-right"
    	legend = ROOT.TLegend(.85-(0.15*(ncol-1)), 0.88-(0.02*(nvar/ncol-1)) ,.50 ,.950)
    legend.SetFillColor(ROOT.kWhite)
    legend.SetFillStyle(0)
    legend.SetTextSize(fontsize)
    legend.SetNColumns(ncol)
    return legend

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

def setHistAttributes (hist, plot_info, line_color, fill_color):
    #hist.SetFillColor(fill_color)
    hist.SetLineColor(line_color)
    hist.SetMarkerColor(line_color)
    #colorArray = [1, 632, 416, 600, 400, 616, 432, 920, 800, 880, 820, 840, 860, 900, 940, 960, 980, 640, 660, 680, 700, 720, 740, 760]
    #stl=1
    #if(line_color==1):
    #	stl=1
    #if(line_color==632):
    #	stl=2
    #if(line_color==416):
    #	stl=3
    #if(line_color==600):
    #	stl=4
    #hist.SetMarkerStyle(stl)
    #hist.SetMarkerSize(1.5)
    #hist.SetLineWidth(2)
    hist.SetTitleOffset(1.3, "x")
    hist.SetTitleOffset(1.3, "y")
    if plot_info["rebin"] != 0:
        if type(hist) != "<class '__main__.TH2F'>":
            hist.Rebin(plot_info["rebin"])
        else:
            print 'Rebin only defined for 1D hist. Use --rebin2D instead.'
    if plot_info["xmin"] < plot_info["xmax"]:
        hist.GetXaxis().SetRangeUser(plot_info["xmin"], plot_info["xmax"])
    if plot_info["ymin"] < plot_info["ymax"]:
        #hist.GetYaxis().SetRangeUser(plot_info["ymin"], plot_info["ymax"])
	#hist.SetMinimum(plot_info["ymin"])
	hist.SetMaximum(plot_info["ymax"])
	#hist.SetMaximum(30)

    #hist.GetYaxis().SetTitle(plot_info["ylabel"])

def createRatio(h1, h2,col):	# h1/h2
    h3 = h1.Clone("h3")
    h3.SetLineColor(col)
    h3.SetMarkerColor(col)
    h3.SetMarkerStyle(21)
    h3.SetTitle("")
    h3.SetMinimum(0.0)
    h3.SetMaximum(1.10)
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
    x.SetTitle("mT (WW) GeV")
    x.SetTitleSize(20)
    x.SetTitleFont(43)
    x.SetTitleOffset(4.0)
    x.SetLabelFont(43)
    x.SetLabelSize(15)

    return h3

def createCanvasPads():	# Create Canvas having two pads
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
# i == File Number, j == variable number, k == color number
def getHistFromFile (plot_info,i,j,k):	
    print "\n===> weight = ",plot_info["weight"],"\n\n"
    file = ROOT.TFile(plot_info["file_name"][i])
    print "File Name: ",plot_info["file_name"][i]
    if not file:
        print 'Failed to open %s' % plot_info["file_name"][i]
        exit(0)
    if plot_info["tree_var"][j] == "":
    	print 'Branch name is blank. Please check!!!'
	exit(0)
    tree = file.Get(plot_info["tree_folder"] + plot_info["tree_name"])
    hist = ROOT.TH1F("hist", "Test", plot_info["nbin"], plot_info["xmin"], plot_info["xmax"])    
    #tree.Draw(plot_info["tree_var"][j] + ">>hist",plot_info["weight1"])
    if plot_info["weight"] == "":
    	tree.Draw(plot_info["tree_var"][j] + ">>hist",plot_info["cut"],"")
	#if j == 0:
	#	cut = "AK4_jj_DeltaEta_gen"
	#	#cut = ""
	#	tree.Draw(plot_info["tree_var"][j] + ">>hist",cut)
	#if j == 1:
	#	cut = "vbf_maxpt_jj_Deta<1.0"
	#	tree.Draw(plot_info["tree_var"][j] + ">>hist",cut)
    else:
	for n, event in enumerate(tree):
		if n>0:
			break
		#indices1 = [i for i, s in enumerate(event.LHEWeightIDs) if aQGC_par in s]
		print "Weight Details: ID = ",event.LHEWeightIDs[491],"\tWeight = ",event.LHEWeights[491]
    	tree.Draw(plot_info["tree_var"][j] + ">>hist",plot_info["weight"][0]+"*("+plot_info["cut"]+")","")
	print plot_info["weight"][0]+"*("+plot_info["cut"]+")"
    if not hist:
        print 'Failed to get hist from file'
        exit(0)
    colorArray = [1, 632, 416, 600, 400, 616, 432, 920, 800, 880, 820, 840, 860, 900, 940, 960, 980, 640, 660, 680, 700, 720, 740, 760]
    setHistAttributes(hist, plot_info, colorArray[k],0)
    hist.SetDirectory(ROOT.gROOT) # detach "hist" from the file
    return hist

def getHistFromFileaQGC(plot_info, key1, ColNum):
    file1 = ROOT.TFile(plot_info["file_name"][0])
    if not file1:
        print 'Failed to open %s' % plot_info["file_name"][0]
        exit(0)
    tree1 = file1.Get(plot_info["tree_folder"] + plot_info["tree_name"])
    hist1 = ROOT.TH1F("hist1", "Test", plot_info["nbin"], plot_info["xmin"], plot_info["xmax"])    
    if plot_info["cut"] == "":
    	tree1.Draw(plot_info["tree_var"][0] + ">>hist1","LHEWeights["+str(key1)+"]","",50000)
    else:
    	tree1.Draw(plot_info["tree_var"][0] + ">>hist1","LHEWeights["+str(key1)+"]"+"*("+plot_info["cut"]+")","",50000)
    print "Number of events = ",hist1.GetEntries()
    print "Integral of hist = ",hist1.Integral()
    #hist1.Scale(1/hist1.Integral())
    #colorArray = ["ROOT.kBlack", "ROOT.kGray", "ROOT.KRed", "ROOT.kGreen", "ROOT.kBlue", "ROOT.kYellow", "ROOT.kMagenta", "ROOT.kCyan", "ROOT.Orange", "ROOT.kViolet", "ROOT.kSpring", "ROOT.kTeal", "ROOT.kAzure", "ROOT.kPink"]
    #colorArray = [ROOT.kBlack, ROOT.kGray, ROOT.KRed, ROOT.kGreen, ROOT.kBlue, ROOT.kYellow, ROOT.kMagenta, ROOT.kCyan, ROOT.Orange, ROOT.kViolet, ROOT.kSpring, ROOT.kTeal, ROOT.kAzure, ROOT.kPink]
    colorArray = [1, 920, 632, 416, 600, 400, 616, 432, 800, 880, 820, 840, 860, 900, 940, 960]
    setHistAttributes(hist1, plot_info, colorArray[ColNum],0)
    hist1.SetDirectory(0)	# detatch hist1 from the file
    return hist1;

##################################################################

################	aQGC Specific Function	##################
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

def aQGC_plotting (plot_info, aQGC_key, aQGC_val, outputNameString, skip):
    print "size of aQGC parameter: ",len(aQGC_key)
    print "==> Key: ",aQGC_key
    print "==> Val: ",aQGC_val
    if len(aQGC_key) == 0:
        print 'Error: No aQGC parameter for this input root file'
        exit(0)
    	
    c2 = getCanvas()
    setTDRStyle(c2, 1, 13, "No") 
    legend = getLegends(plot_info["legPos"],2,len(aQGC_key)/skip,plot_info["legFontSize"])

    hist1 = []
    new_aQGC_val = []
    for i in range(0,len(aQGC_key),skip):
	print "Working on hist number : ",i
    	hist1.append(getHistFromFileaQGC(plot_info, aQGC_key[i] ,i/skip))   
	new_aQGC_val.append(aQGC_val[i])
    for a,list1 in enumerate(hist1):
    	#print a,list1
    	list1.GetXaxis().SetTitle(plot_info["xlabel"])
    	if plot_info["ylabel"] == "":
    	    plot_info["ylabel"] = "Events / %s GeV" % int(list1.GetBinWidth(1))
    	list1.GetYaxis().SetTitle(plot_info["ylabel"])
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


##################################################################
	
	

def CompHistFromTwoBranchSameFile (plot_info):
    c1 = getCanvas()
    setTDRStyle(c1, 1, 13, "No") 
    legend = getLegends(plot_info["legPos"],2,len(plot_info["tree_var"]),plot_info["legFontSize"])

    hist1 = []
    new_aQGC_val = []
    for i in range(0,len(plot_info["tree_var"])):
	print "Working on hist number : ",i
    	hist1.append(getHistFromFile(plot_info,0,i,i))
	print "=======\n"
	# below pattern is for the legend...
	if (plot_info["leg"][i] != ""):
		new_aQGC_val.append(plot_info["leg"][i])
	else:
		new_aQGC_val.append(plot_info["tree_var"][i])
    maximumY = 0.0
    for a,list1 in enumerate(hist1):
    	list1.Sumw2()
    	list1.GetXaxis().SetTitle(plot_info["xlabel"])
    	if plot_info["ylabel"] == "":
    	    plot_info["ylabel"] = "Events / %s GeV" % int(list1.GetBinWidth(1))
    	list1.GetYaxis().SetTitle(plot_info["ylabel"])
    	list1.Scale(1/list1.Integral())
    	legend.AddEntry(list1, new_aQGC_val[a],"lpe")
	if (list1.GetMaximum() > maximumY ):
		maximumY=list1.GetMaximum()*1.5
	list1.SetMaximum(maximumY)
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

def CompHistFromTwoFile (plot_info):
    c1 = getCanvas()
    setTDRStyle(c1, 1, 13, "No") 
    legend = getLegends(plot_info["legPos"],1,len(plot_info["file_name"]),plot_info["legFontSize"])
    hist = []
    for i in range(0,len(plot_info["file_name"])):
    	print "working on file number : ",i
	# Below if condition is for specific perpose
	if i == 1:
		plot_info["weight"] = ""
	print "weight = ",plot_info["weight"],"\n\n"
	hist.append(getHistFromFile(plot_info,i,0,i))

    for a,list1 in enumerate(hist):
    	list1.Sumw2()
	list1.GetXaxis().SetTitle(plot_info["xlabel"])
	if plot_info["ylabel"] == "":
		plot_info["ylabel"] = "Events / %s GeV" % int(list1.GetBinWidth(1))
	list1.GetYaxis().SetTitle(plot_info["ylabel"])
	list1.Scale(1.0/list1.Integral())
	#list1.Scale((36.0*174.0)/list1.Integral())
	#print "Print Bin Number = ",list1.GetXaxis().FindBin(1100)
	#print "HISTO: ",list1
	print "Total Integral = ",list1.Integral()
	#list1.Fit("gaus")
	#print "Integral (2000,6000) = ",list1.Integral(list1.GetXaxis().FindBin(2000),list1.GetXaxis().FindBin(6000))
	legend.AddEntry(list1,plot_info["leg"][a],"lpe")
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
    parser.add_argument('--cut',  type=str, required=False, default="",
                        help="cut to be applied")  
    parser.add_argument('--xlabel', type=str, required=False, default="Invariant mass of WW", 
                        help="x axis label")
    parser.add_argument('--ylabel', type=str, required=False, default="", 
                        help="y axis label")
    parser.add_argument('--leg', nargs='+', type=str, required=False, default="", 
                        help="Histo/Rootfile legend name")
    parser.add_argument('--legPos', type=str, required=False, default="bc", 
                        help="Legend Position: tr, tl, tc, br, bc, bl")
    parser.add_argument('--legFontSize', type=float, required=False, default=0.04, 
                        help="Legend font size")
    parser.add_argument('--weight', nargs='+', type=str, required=False, default="", 
                        help="Weight corresponding to variable")
    parser.add_argument('--xmin', type=float, required=False, default=0, 
                        help="minimum x value")
    parser.add_argument('--xmax', type=float, required=False, default=600, 
                        help="maximum x value")   
    parser.add_argument('--ymin', type=float, required=False, default=0, 
                        help="minimum y value")
    parser.add_argument('--ymax', type=float, required=False, default=0, 
                        help="maximum y value")
    parser.add_argument('--nbin', type=int, required=False, default=15, 
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
