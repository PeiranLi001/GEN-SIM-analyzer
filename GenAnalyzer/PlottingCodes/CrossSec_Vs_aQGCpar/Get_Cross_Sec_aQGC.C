#define Get_Cross_Sec_aQGC_cxx
#include "Get_Cross_Sec_aQGC.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <fstream>

void Get_Cross_Sec_aQGC::Loop(std::string outFileName, float Total_Cross_Section)
{
//   In a ROOT session, you can do:
//      root> .L Get_Cross_Sec_aQGC.C
//      root> Get_Cross_Sec_aQGC t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   ofstream file;
   if (outFileName.c_str() == "")
   	outFileName = "temp_data.txt";
   file.open (outFileName.c_str());

   fChain->SetBranchStatus("*",0);
   fChain->SetBranchStatus("LHEWeights",1);
   fChain->SetBranchStatus("LHEWeightIDs",1);

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   double SumWeight[718] = {0};
   double AvgWeight[718] = {0};
   double DevWeight[718] = {0};
   int countWgt ;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      if (jentry%10000 == 0)
      	cout<<"entry = "<<jentry<<endl;
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      countWgt = 0;
      for (int i =446; i<LHEWeights->size();i++)
      {
      	//WeightIDs.push_back(LHEWeightIDs->at(i));
	//Weights.push_back(LHEWeights->at(i));
	//cout<<LHEWeightIDs->at(i)<<endl;
	SumWeight[countWgt] += (LHEWeights->at(i))/(LHEWeights->at(0));
	countWgt += 1;
      }
   }

   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      if (jentry%10000 == 0)
      	cout<<"entry = "<<jentry<<endl;
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      countWgt = 0;
      for (int i =446; i<LHEWeights->size();i++)
      {
	AvgWeight[countWgt] = SumWeight[countWgt]/nentries;
	DevWeight[countWgt] = ( (LHEWeights->at(i))/(LHEWeights->at(0)) - AvgWeight[countWgt]) * ( (LHEWeights->at(i))/(LHEWeights->at(0)) - AvgWeight[countWgt]);
	countWgt += 1;
      }
   }

   cout<<"\n\n============================================================\n\n"<<endl;
   cout<<"\n\n\t\t\t Print Cross-Section	"<<endl;
   cout<<"============================================================\n\n"<<endl;
   cout<<"aQGC Par Name \t Cross-Section \t +/- \t ErrorOnCross-section"<<endl;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      if (jentry > 0) break;
      countWgt = 0;
      for (int i =446; i<LHEWeights->size();i++)
      {
      	//cout<<LHEWeightIDs->at(i)<<"\t"<<SumWeight[countWgt]*Total_Cross_Section<<" +/- "<<DevWeight[countWgt]/sizeof(DevWeight)<<endl;
      	cout<<LHEWeightIDs->at(i)<<endl;
      	cout<<"\t\tSum Wgt   = "<<SumWeight[countWgt]*Total_Cross_Section<<endl;
      	cout<<"\t\tSum Wgt/N = "<<(SumWeight[countWgt]*Total_Cross_Section)/nentries<<endl;
      	cout<<"\t\tError     = "<<DevWeight[countWgt]/sizeof(DevWeight)<<endl;
	file<<LHEWeightIDs->at(i)<<"\t"<<(SumWeight[countWgt]*Total_Cross_Section)/nentries<<"\t"<<DevWeight[countWgt]/sizeof(DevWeight)<<endl;
	countWgt++;
      }
   }
   file.close();
}
