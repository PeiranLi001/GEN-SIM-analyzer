#include "TH1F.h"
void makeplot(){

  TFile *f = TFile::Open("../Radion_hh_narrow_M270_LHEBqrk.root");
  TTree *t1 = (TTree*)f->Get("otree");

  

  //create histograms
  TH1F *leading_photon_Pt = new TH1F("leading_photon_Pt","leading photon Pt",100,-10,1000);
  TH1F *leading_photon_M = new TH1F("leading_photon_M","leading photon M",100,-100,100);
  TH1F *njets = new TH1F("njets","njets",10,0,10);
  TH1F *deltaR_H1_H2 = new TH1F("deltaR_H1_H2","dR(H1,H2)",100,-15,15);



  double gen_leading_photon_Pt;
  double gen_leading_photon_M;
  double genJetAK4_njets;
  double gen_deltaR_H1_H2;
  //double weight;

  t1->SetBranchAddress("gen_leading_photon_Pt",& gen_leading_photon_Pt);
  t1->SetBranchAddress("gen_leading_photon_M",& gen_leading_photon_M);
  t1->SetBranchAddress("genJetAK4_njets",& genJetAK4_njets);
  t1->SetBranchAddress("gen_deltaR_H1_H2",& gen_deltaR_H1_H2);

 //t1->SetBranchAddress("weight",& weight);

  int entries = t1->GetEntries();
 
  std::cout << "entries:" << entries << std::endl;

  for (int i=0; i<entries; i++){
    
    t1->GetEntry(i);

    leading_photon_Pt->Fill(gen_leading_photon_Pt);
    leading_photon_M->Fill(gen_leading_photon_M);
    njets->Fill(genJetAK4_njets);
    deltaR_H1_H2->Fill(gen_deltaR_H1_H2);

  }
  
  TCanvas *c1 = new TCanvas();
  leading_photon_Pt->Draw();
  c1->SaveAs("~/private/HHtoWWgg/CMSSW_10_2_22/src/GEN-SIM-analyzer/GenAnalyzer/PlottingCodes/all/leading_photon_Pt.pdf");

  TCanvas *c2 = new TCanvas();
  leading_photon_M->Draw();
  c2->SaveAs("~/private/HHtoWWgg/CMSSW_10_2_22/src/GEN-SIM-analyzer/GenAnalyzer/PlottingCodes/all/leading_photon_M.pdf");

  TCanvas *c3 = new TCanvas();
  njets->Draw();
  c3->SaveAs("~/private/HHtoWWgg/CMSSW_10_2_22/src/GEN-SIM-analyzer/GenAnalyzer/PlottingCodes/all/njets.pdf");

  TCanvas *c4 = new TCanvas();
  deltaR_H1_H2->Draw();
  c4->SaveAs("~/private/HHtoWWgg/CMSSW_10_2_22/src/GEN-SIM-analyzer/GenAnalyzer/PlottingCodes/all/deltaR_H1_H2.pdf");


}
  
  
