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
    print "\n======= Start Makeing plots =================\n"

    #canvas = ROOT.TCanvas("c2","c2",50,50,800,800)
    #plotter.getHistFromFile(plot_info, 0, 1, 1)
    #canvas.Print(plot_info["output_file"]) 
    #plotter.CompHistFromTwoFile(plot_info)
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
