#! /usr/bin/env python2.6

import sys
import ROOT as rt

from Template import *

######################
#  Print C++ header  #
######################


sys.stdout = open(OutPutCodeName+".C",'w')

ColorChoice = ["kBlue", "kSpring", "kTeal", "kGray", "kMagenta", "kRed", "kCyan+2", "kViolet+2", "kGreen+2", "kOrange", "kPink", "kBlack", "kYellow"]

print "// ############################################################"
print "// # Usage                                                    #"
print "// #                                                          #"
print "// ############################################################"
print ''
print '#include "TGraphErrors.h"'
print '#include "setTDRStyle.C"'
print ''
print ''

print 'void %s(){'%OutPutCodeName
print '\tifstream in;'
if iffit == 1:
    print '\tofstream outfile;'
    print '\toutfile.open ("%s.txt");'%(OutPutCodeName)
print ''
print '\t// Set TDR Style'
print '\tsetTDRStyle();'
print ''
print '\tTGaxis::SetMaxDigits(%i);'%setMaxdigit
print '\tgStyle->SetOptStat(%i);'%getstat
print '\tgStyle->SetOptFit(%i);'%getstat
print ''
print '\t//Define the variables'
print '\tFloat_t ',', '.join(VarInTextFile),' ;'
pre = "v_"
Pre_VarInTextFile = [pre + x for x in VarInTextFile ]
#print Pre_VarInTextFile
if len(VarToPlot) == 0:
	VarToPlot = VarInTextFile
Pre_VarToPlot = [pre + x for x in VarToPlot ]
#print "\n\n",Pre_VarToPlot,"\n\n"
print '\tvector<Float_t> ',', '.join(Pre_VarInTextFile),' ;'
print ''
print ''
print '\tInt_t nlines = 0;'
print '\tTFile *f = new TFile("%s.root","RECREATE");'%OutPutCodeName
print '\tTNtuple *ntuple = new TNtuple("ntuple","data from ascii file","',':'.join(VarInTextFile),'");'
print ''
print ''
print '\tTCanvas* c1 = new TCanvas("c1","",1);'
print '\tc1->Range(0,0,1,1);'
print '\t//c1->SetLogy(1);'
print '\t//c1->SetLogx(1);'
print '\t//TPad *pad = new TPad("pad","",0,0,1,1);'
print '\t//pad->SetGrid();'
print '\t//pad->Draw();'
print '\t//pad->cd();'
print '\tTLegend * leg1 = new TLegend(0.530000,0.900000,0.630000,0.800000);'
print ''
print ''

print '\tstring line;'
for f in range(0,len(InputData)):
	#print InputData[f]
	print '\tin.open("%s%s");'%(datapath,InputData[f])
	print '\twhile(getline(in,line))'
	print '\t{'
	print "\t\tif(line[0] == '#') continue;"
	print '\t\t'
	print '\t\tstringstream(line) >>',' >> '.join(VarInTextFile),';'
	print '\t\t'
	print '\t\t//cout<<"===> "<<',' <<"\\t" <<  '.join(VarInTextFile),'<<endl;'
	print ''
	#print '\t\tif(%s < 11.0){'%VarInTextFile[1]
	for var in range(0,len(VarInTextFile)):
		if (xscaleFactor != 0.0 and yscaleFactor != 0.0):
			if (var == 0):
				print '\t\t%s.push_back(%s/%0.15f);'%(Pre_VarInTextFile[var],VarInTextFile[var],xscaleFactor)
			else:
				print '\t\t%s.push_back(%s/%0.15f);'%(Pre_VarInTextFile[var],VarInTextFile[var],yscaleFactor)
		else:
			if (xscaleFactor != 0.0 ):
				if (var == 0):
					print '\t\t%s.push_back(%s/%0.15f);'%(Pre_VarInTextFile[var],VarInTextFile[var],xscaleFactor)
				else:
					print '\t\t%s.push_back(%s);'%(Pre_VarInTextFile[var],VarInTextFile[var])
			else:
				if (yscaleFactor != 0.0 ):
					if (var == 0):
						print '\t\t%s.push_back(%s);'%(Pre_VarInTextFile[var],VarInTextFile[var])
					else:
						print '\t\t%s.push_back(%s/%0.15f);'%(Pre_VarInTextFile[var],VarInTextFile[var],yscaleFactor)
				else:
					print '\t\t%s.push_back(%s);'%(Pre_VarInTextFile[var],VarInTextFile[var])
	print '\t\tntuple->Fill(',' , '.join(VarInTextFile),');'
	#print '\t\t}'
	print '\t}'
	
	print '\tin.close();'
	print '\t'
	print '\tTGraphErrors * gr%i = new TGraphErrors(%s.size()); '%(f,Pre_VarInTextFile[0])
	print '\t    '
	print '\tfor (unsigned int i = 0; i<%s.size();i++)'%Pre_VarInTextFile[0]
	print '\t{'
	#print '\t       gr%i->SetPoint(i,%s[i],%s[i]);'%(f,Pre_VarInTextFile[0],Pre_VarInTextFile[1])
	print '\t        gr%i->SetPoint(i,%s[i],%s[i]);'%(f,Pre_VarToPlot[0],Pre_VarToPlot[1])
	if (len(Pre_VarToPlot)==3):
		print '\t        gr%i->SetPointError(i,0,%s[i]);'%(f,Pre_VarToPlot[2])
	else:
		print '\t        gr%i->SetPointError(i,0,0);'%f
	print '\t}'
	print '\t'
	print '\t'
	print '\tgr%i->SetTitle("");'%f
	print '\t//gr%i->SetTitle("%s vs %s");'%(f,VarInTextFile[0],VarInTextFile[1])
	print '\tgr%i->GetXaxis()->SetTitle("%s");'%(f,xlabel)
	print '\tgr%i->GetYaxis()->SetTitle("%s");'%(f,ylabel)
	print '\tgr%i->GetYaxis()->SetTitleOffset(%f);'%(f,yoffset)
	print '\tgr%i->GetXaxis()->SetTitleOffset(%f);'%(f,xoffset)
	print '\tgr%i->GetYaxis()->SetTitleSize(0.05);'%f
	print '\tgr%i->GetXaxis()->SetTitleSize(0.05);'%f
	print '\t//gr%i->GetXaxis()->SetLabelSize(0.05);'%f
	if yrange[0]!=yrange[1]:
		print '\tgr%i->GetYaxis()->SetRangeUser(%f,%f);'%(f,yrange[0],yrange[1])
	print '\tgr%i->SetMarkerSize(1);'%f
	print '\tgr%i->SetMarkerColor(%s);'%(f,ColorChoice[f])
	print '\tgr%i->SetLineColor(%s);'%(f,ColorChoice[f])
	print '\tgr%i->SetMarkerStyle(21);'%f
	print '\t'
	print '\tgr%i->Draw("AP");'%f
	print '\t//gr%i->Draw("ACP");'%f
        print '\tleg1->AddEntry(gr%i, "%s","lep");'%(f,legends[f])
        print '\tleg1->SetFillColor(0);'
        print '\tleg1->SetFillStyle(0);'
        print '\tleg1->SetBorderSize(0);'
        print '\tleg1->SetTextFont(42);'
        print '\tleg1->SetTextSize(0.05);'
        print '\tleg1->Draw();'
	if iffit == 1:
		try:
		    fitXrange
		except :
		    print '\tgr%i->Fit("%s","","");'%(f,fitfunction)
		else:
		    print '\tgr%i->Fit("%s","","",%f,%f);'%(f,fitfunction,fitXrange[0],fitXrange[1])
		print '\tgr%i->GetFunction("%s")->SetLineColor(%s);'%(f,fitfunction,ColorChoice[f])
		print '\n'
		print '\tTF1 *fit%i = gr%i->GetFunction("%s");'%(f,f,fitfunction)
		#print '\tcout<<"Fit value at 0 = "<<fit%i->Eval(-3.0)<<"\t"<<fit%i->Eval(3.0)<<endl;'%(f,f)
		print '\tcout<<"aqgc_limits_For %s:  "<<fit%i->GetX(3,fit%i->GetXmin(),0.0)<<"\t"<<fit%i->GetX(3,0.0,fit%i->GetXmax())<<endl;'%(legends[f],f,f,f,f)
		print '\toutfile<<"aqgc_limits_For %s  "<<fit%i->GetX(3,fit%i->GetXmin(),0.0)<<"\t"<<fit%i->GetX(3,0.0,fit%i->GetXmax())<<endl;'%(legends[f],f,f,f,f)
	print '\n'
	
	print '\tc1->SetName("%s");'%legends[f]
        print '\tc1->Write();'
        print '\tc1->SaveAs("%s_%s.pdf");'%(OutPutCodeName,legends[f])
        print '\tc1->SaveAs("%s_%s.png");'%(OutPutCodeName,legends[f])
        print '\tc1->Clear();'
        print '\tleg1->Clear();'

	#print '\tgr->Draw("ALP");'
	print '\n'
	for size in range(0,len(Pre_VarInTextFile)):
		print '\t%s.clear();'%Pre_VarInTextFile[size]
print '\t'

print '\tTMultiGraph *mg = new TMultiGraph("mg",";%s;%s");'%(xlabel,ylabel)
for f in range(0,len(InputData)):
	print '\tmg->Add(gr%i);'%f
if yrange[0]!=yrange[1]:
	print '\tmg->SetMaximum(%f);'%yrange[1]
	print '\tmg->SetMinimum(%f);'%yrange[0]
print '\t//mg->GetXaxis()->SetLimits(750,815);'
print '\tmg->Draw("AP");'

LegY2Pos = 0.0

if len(tlatexx) > 0:
	for text in range(0,len(tlatexx)):
		print '\tTLatex *text%i = new TLatex(%f,%f,"%s");'%(text,pos1,pos2-((pos2/100.)*sep*text),tlatexx[text])
		LegY2Pos = pos2-((pos2/100.)*sep*text)
	print ''
	for text in range(0,len(tlatexx)):
		print '\ttext%i->SetNDC();'%text
	print ''
	for text in range(0,len(tlatexx)):
		print '\ttext%i->SetTextFont(42);'%text
	print ''
	for text in range(0,len(tlatexx)):
		print '\ttext%i->SetTextSize(0.05);'%text
	print ''
	for text in range(0,len(tlatexx)):
		print '\ttext%i->Draw("same");'%text
	print ''

if len(legends) > 0:
    if len(tlatexx) > 0:
	print '\tTLegend * leg = new TLegend(%f,%f,%f,%f);'%(pos1,LegY2Pos-0.06*len(tlatexx),pos1+0.20,LegY2Pos-0.02)
    else:
	print '\tTLegend * leg = new TLegend(%f,%f,%f,%f);'%(pos1,0.90,pos1+0.30,0.50)
    print '\tleg-> SetNColumns(3);'
    print '\t//leg->SetHeader("The Legend Title","C"); // option "C" allows to center the header'
    for size in range(0,len(legends)):
    	print '\tleg->AddEntry(gr%i, "%s","lep");'%(size,legends[size])
    print '\tleg->SetFillColor(0);'
    print '\tleg->SetFillStyle(0);'
    print '\tleg->SetBorderSize(0);'
    print '\tleg->SetTextFont(42);'
    print '\tleg->SetTextSize(0.05);'
    print '\tleg->Draw();'

print ''
print ''
print '\tTLatex *cmsprem = new TLatex(%f,%f,"#it{CMS Preliminary}");'%(cpos1,cpos2)
print '\tcmsprem->SetNDC();'
print ''
print '\tTLatex *gen = new TLatex(%f,%f,"%s");'%(cpos3,cpos2,DetInfo)
print '\tgen->SetNDC();'
print ''
print '\tcmsprem->Draw("same");'
print '\tgen->Draw("same");'
print ''
print ''
print ''
print '\tc1->Write();'
print '\tc1->SaveAs("%s.pdf");'%OutPutCodeName
print '\tc1->SaveAs("%s.png");'%OutPutCodeName
print '\tf->Write();'
if iffit == 1:
    print '\toutfile.close();'

print '}'


################
#  Close file  #
################
