#!/usr/bin/env python
import ROOT
import CMS_lumi, tdrstyle
import argparse
import plot_functions as plotter

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(0)

def main():
    plot_info = getPlotArgs()
    print plot_info
    c1 = ROOT.TCanvas()
    plotter.setTDRStyle(c1, 1, 13, "No")
    print "\n======= Start Makeing plots =================\n"
    
    plotter.CompHistFromTwoBranchSameFile(plot_info)



def getPlotArgs():
    parser = plotter.getBasicParser()
    parser.add_argument("-n", "--file_name", nargs='+', type=str, required=False,
                        help="Name of root file where plots are stored")
    parser.add_argument("--is_data", action='store_true',
                        help="Plot histogram with data points")
    return vars(parser.parse_args())
    
if __name__ == "__main__":
    main()
