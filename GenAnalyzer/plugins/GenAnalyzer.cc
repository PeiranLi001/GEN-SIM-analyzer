// -*- C++ -*-
//
// Package:    GenAnalyzerNew/GenAnalyzer
// Class:      GenAnalyzer
// 
/**\class GenAnalyzer GenAnalyzer.cc GenAnalyzerNew/GenAnalyzer/plugins/GenAnalyzer.cc

 Description: Read the GENSIM file

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Ram Krishna Sharma
//         Created:  Wed, 21 Sep 2016 22:13:15 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "GenAnalyzer.h"

#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/METReco/interface/GenMETCollection.h"
#include "DataFormats/METReco/interface/GenMET.h"

#include "DataFormats/Math/interface/deltaR.h"

#include <vector>
#include "TMath.h"
#include "TLorentzVector.h"
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

using namespace std;
using namespace reco;
using namespace edm;
//

// ------------ method called for each event  ------------
void
GenAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

  edm::Handle<LHEEventProduct> LHEEventHandle;
  iEvent.getByToken(LHEEventToken, LHEEventHandle);
  const LHEEventProduct* LHE = 0;

  if(LHEEventHandle.isValid()){
  	LHE = LHEEventHandle.product();

	for(const auto& weight : LHE->weights()) {
		LHEWeightIDs_.push_back(weight.id);
		LHEWeights_.push_back(weight.wgt);
	}
	
        std::vector<int> leptons ;      
        std::vector<int> finalQuarks ;      
        std::vector<int> intermediates ;
        std::vector<int> tops ;        
	TLorentzVector Is_Iqrk1,Is_Iqrk0;
        //PG loop over particles in the event
	int incomingPart = 0;
        for (int iPart = 0 ; iPart < LHE->hepeup().NUP; ++iPart){
            
            //PG incoming particle          
            if (LHE->hepeup().ISTUP.at (iPart) == -1){
		incomingPart++;
		if (incomingPart == 1)
		{
			Is_Iqrk0.SetPxPyPzE
			(
                 	LHE->hepeup().PUP[incomingPart][0], //PG px
                 	LHE->hepeup().PUP[incomingPart][1], //PG py
                 	LHE->hepeup().PUP[incomingPart][2], //PG pz
                 	LHE->hepeup().PUP[incomingPart][3] //PG E
			);
		}
		if (incomingPart == 2)
		{
			Is_Iqrk1.SetPxPyPzE
			(
                 	LHE->hepeup().PUP[incomingPart][0], //PG px
                 	LHE->hepeup().PUP[incomingPart][1], //PG py
                 	LHE->hepeup().PUP[incomingPart][2], //PG pz
                 	LHE->hepeup().PUP[incomingPart][3] //PG E
			);
		}
            }

            //PG outgoing particles          
            if (LHE->hepeup().ISTUP.at (iPart) == 1){
                //PG leptons
                if (abs (LHE->hepeup().IDUP.at (iPart)) == 11 ||   //PG electron
                    abs (LHE->hepeup().IDUP.at (iPart)) == 13 ||   //PG muon
                    abs (LHE->hepeup().IDUP.at (iPart)) == 15 ||   //PG tau
                    abs (LHE->hepeup().IDUP.at (iPart)) == 12 ||   //PG neutrino
                    abs (LHE->hepeup().IDUP.at (iPart)) == 14 ||   //PG neutrino
                    abs (LHE->hepeup().IDUP.at (iPart)) == 16)     //PG neutrino                    
                    {
                    leptons.push_back (iPart) ;
                    } //PG leptons
                else
                    {
                    finalQuarks.push_back (iPart) ;
                    }
                
            } 
            
            //PG intermediates
            if (LHE->hepeup().ISTUP.at(iPart) == 2){
                intermediates.push_back (iPart) ;
            }
            
            //PG tops
            if (abs(LHE->hepeup().IDUP.at(iPart)) == 6){
                tops.push_back (iPart) ;
            }
        } //PG loop over particles in the event
    

        // --------------- Indices for final state particles  -------------------
        int i_olep_part = -1;
        int i_olep_anti = -1;        
        int i_wqrk_1 = -1;
        int i_wqrk_2 = -1;        
        int i_iqrk_1 = -1;
        int i_iqrk_2 = -1;        

        int tmpfill1 = 0;
        int signalFlag = 0;
        int signalWCtr = 0;
        
        if (leptons.size() == 2){
            if (LHE->hepeup().IDUP.at(leptons.at(0)) > 0){
                i_olep_part = leptons.at(0);
                i_olep_anti = leptons.at(1);
                if (LHE->hepeup().IDUP.at(leptons.at(0)) == 11 ||   //PG electron
                    LHE->hepeup().IDUP.at(leptons.at(0)) == 13 ||   //PG muon
                    LHE->hepeup().IDUP.at(leptons.at(0)) == 15) isMuMinus_ = 1;           
            }
            else{
                i_olep_part = leptons.at(1);
                i_olep_anti = leptons.at(0);    
                if (LHE->hepeup().IDUP.at(leptons.at(1)) == 11 ||   //PG electron
                    LHE->hepeup().IDUP.at(leptons.at(1)) == 13 ||   //PG muon
                    LHE->hepeup().IDUP.at(leptons.at(1)) == 15) isMuMinus_ = 1;           
                            
            }
	    //std::cout << "lep1 = " << LHE->hepeup().IDUP.at(i_olep_part) << std::endl;
	    //std::cout << "lep1 = " << LHE->hepeup().IDUP.at(i_olep_anti) << std::endl;
        }
        else{
            std::cout << "Problem!" << std::endl;
        }
        
        // --------------- If signal, find the quarks from the W  -------------------
        if (finalQuarks.size() == 4){
            for (unsigned int a = 0; a < finalQuarks.size(); ++a ){
                int tmpMoth = LHE->hepeup().MOTHUP.at(finalQuarks.at(a)).first - 1;

                if ( abs(LHE->hepeup().IDUP.at(tmpMoth)) == 24 && LHE->hepeup().IDUP.at(finalQuarks.at(a)) > 0 ){
                    signalWCtr++;
                }
                if ( abs(LHE->hepeup().IDUP.at(tmpMoth)) == 24 && LHE->hepeup().IDUP.at(finalQuarks.at(a)) < 0 ){
                    signalWCtr++;
                }
                if ( abs(LHE->hepeup().IDUP.at(tmpMoth)) != 24 && tmpfill1){
                    signalWCtr++;
                }
                if ( abs(LHE->hepeup().IDUP.at(tmpMoth)) != 24 && !tmpfill1){
                    signalWCtr++;
                    tmpfill1++;
                }
            }
        }
        else{
            std::cout << "Problem!" << std::endl;
        }
        
        //if (signalWCtr == 4 && tops.size() == 0){ signalFlag = 1; NSignal++; }
        
        // --------------- Assign quarks based on W invariant mass -------------------
        float distanceToWMass = 9999.;
        if (finalQuarks.size() == 4){
        for (unsigned int a = 0; a < finalQuarks.size(); ++a ){
            for (unsigned int b = a+1; b < finalQuarks.size(); ++b ){
                TLorentzVector qrk0
                (
                 LHE->hepeup().PUP[a][0], //PG px
                 LHE->hepeup().PUP[a][1], //PG py
                 LHE->hepeup().PUP[a][2], //PG pz
                 LHE->hepeup().PUP[a][3] //PG E
                 ) ;
                TLorentzVector qrk1
                (
                 LHE->hepeup().PUP[b][0], //PG px
                 LHE->hepeup().PUP[b][1], //PG py
                 LHE->hepeup().PUP[b][2], //PG pz
                 LHE->hepeup().PUP[b][3] //PG E
                 ) ;

                float tmpDistance = std::abs( (qrk0+qrk1).M() - 80.2 );
                if (tmpDistance < distanceToWMass){
                    i_wqrk_1 = finalQuarks.at(a);
                    i_wqrk_2 = finalQuarks.at(b);                    
		    //cout<<"i_wqrk_1 = "<<i_wqrk_1<<"\ti_wqrk_2 = "<<i_wqrk_2<<endl;
                    distanceToWMass = tmpDistance;
                }
            }
        }
	}
        
        // quarks based on particle-antiparticle
        if (LHE->hepeup().IDUP.at(i_wqrk_1) < 0){
            int tmpint = i_wqrk_1;
            i_wqrk_1 = i_wqrk_2;
            i_wqrk_2 = tmpint;
        }
        
        // assign the other quarks
        std::vector<int> finalQuarksNotW;
        if (finalQuarks.size() == 4){
        for (unsigned int a = 0; a < finalQuarks.size(); ++a ){
            if (finalQuarks.at(a) != i_wqrk_1 && finalQuarks.at(a) != i_wqrk_2){
                finalQuarksNotW.push_back( finalQuarks.at(a) );
            }
        }
        if (finalQuarksNotW.size() == 2){
            i_iqrk_1 = finalQuarksNotW.at(0);
            i_iqrk_2 = finalQuarksNotW.at(1);                
        }
	}

        //std::cout << "signalFlag = " << signalFlag << std::endl;
        //std::cout << "isMuMinus = " << isMuMinus << std::endl;        
        //std::cout << "lep1 = " << LHE->hepeup().IDUP.at(i_olep_part) << std::endl;
        //std::cout << "lep2 = " << LHE->hepeup().IDUP.at(i_olep_anti) << std::endl;        
        //std::cout << "qrk1 = " << LHE->hepeup().IDUP.at(i_wqrk_1) << std::endl;
        //std::cout << "qrk2 = " << LHE->hepeup().IDUP.at(i_wqrk_2) << std::endl;        
        //std::cout << "i-qrk1 = " << LHE->hepeup().IDUP.at(i_iqrk_1) << std::endl;
        //std::cout << "i-qrk2 = " << LHE->hepeup().IDUP.at(i_iqrk_2) << std::endl;        

        if (finalQuarks.size() == 4){
        if (leptons.size() == 2){
        TLorentzVector fs_lep0
        (
         LHE->hepeup().PUP[i_olep_part][0], //PG px
         LHE->hepeup().PUP[i_olep_part][1], //PG py
         LHE->hepeup().PUP[i_olep_part][2], //PG pz
         LHE->hepeup().PUP[i_olep_part][3] //PG E
         ) ;
		LHELeptPt_	= fs_lep0.Pt();
		LHELeptEta_	= fs_lep0.Eta();
		LHELeptPhi_	= fs_lep0.Phi();
		LHELeptM_	= fs_lep0.M();
		LHELeptE_	= fs_lep0.E();
	//cout<<LHELeptPt_<<LHELeptEta_<<LHELeptPhi_<<LHELeptM_<<LHELeptE_<<end;

        TLorentzVector fs_lep1
        (
         LHE->hepeup().PUP[i_olep_anti][0], //PG px
         LHE->hepeup().PUP[i_olep_anti][1], //PG py
         LHE->hepeup().PUP[i_olep_anti][2], //PG pz
         LHE->hepeup().PUP[i_olep_anti][3] //PG E
         ) ;
		LHENuPt_  = fs_lep1.Pt();
		LHENuEta_ = fs_lep1.Eta();
		LHENuPhi_ = fs_lep1.Phi();
		LHENuM_   = fs_lep1.M();
		LHENuE_   = fs_lep1.E();
    
        TLorentzVector fs_Wqrk0
        (
         LHE->hepeup().PUP[i_wqrk_1][0], //PG px
         LHE->hepeup().PUP[i_wqrk_1][1], //PG py
         LHE->hepeup().PUP[i_wqrk_1][2], //PG pz
         LHE->hepeup().PUP[i_wqrk_1][3] //PG E
         ) ;
	 	LHE_Wqrk0_E_  = LHE->hepeup().PUP[i_wqrk_1][3];
		LHE_Wqrk0_pt_ = fs_Wqrk0.Pt();
		LHE_Wqrk0_eta_= fs_Wqrk0.Eta();
		LHE_Wqrk0_phi_= fs_Wqrk0.Phi();
		LHE_Wqrk0_M_  = fs_Wqrk0.M();
		LHE_Wqrk0_Mt_  = fs_Wqrk0.Mt();

        TLorentzVector fs_Wqrk1
        (
         LHE->hepeup().PUP[i_wqrk_2][0], //PG px
         LHE->hepeup().PUP[i_wqrk_2][1], //PG py
         LHE->hepeup().PUP[i_wqrk_2][2], //PG pz
         LHE->hepeup().PUP[i_wqrk_2][3] //PG E
         ) ;
	 	LHE_Wqrk1_E_   = LHE->hepeup().PUP[i_wqrk_2][3];
		LHE_Wqrk1_pt_  = fs_Wqrk1.Pt();
		LHE_Wqrk1_eta_ = fs_Wqrk1.Eta();
		LHE_Wqrk1_phi_ = fs_Wqrk1.Phi();
		LHE_Wqrk1_M_   = fs_Wqrk1.M();
		LHE_Wqrk1_Mt_  = fs_Wqrk1.Mt();

        TLorentzVector fs_Iqrk0
        (
         LHE->hepeup().PUP[i_iqrk_1][0], //PG px
         LHE->hepeup().PUP[i_iqrk_1][1], //PG py
         LHE->hepeup().PUP[i_iqrk_1][2], //PG pz
         LHE->hepeup().PUP[i_iqrk_1][3] //PG E
         ) ;
	 	LHE_Iqrk0_E_   = LHE->hepeup().PUP[i_iqrk_1][3];
		LHE_Iqrk0_pt_  = fs_Iqrk0.Pt();
		LHE_Iqrk0_eta_ = fs_Iqrk0.Eta();
		LHE_Iqrk0_phi_ = fs_Iqrk0.Phi();
		LHE_Iqrk0_M_   = fs_Iqrk0.M();
		LHE_Iqrk0_Mt_  = fs_Iqrk0.Mt();


        TLorentzVector fs_Iqrk1
        (
         LHE->hepeup().PUP[i_iqrk_2][0], //PG px
         LHE->hepeup().PUP[i_iqrk_2][1], //PG py
         LHE->hepeup().PUP[i_iqrk_2][2], //PG pz
         LHE->hepeup().PUP[i_iqrk_2][3] //PG E
         ) ;
	 	LHE_Iqrk1_E_   = LHE->hepeup().PUP[i_iqrk_2][3];
		LHE_Iqrk1_pt_  = fs_Iqrk1.Pt();
		LHE_Iqrk1_eta_ = fs_Iqrk1.Eta();
		LHE_Iqrk1_phi_ = fs_Iqrk1.Phi();
		LHE_Iqrk1_M_   = fs_Iqrk1.M();
		LHE_Iqrk1_Mt_  = fs_Iqrk1.Mt();

        TLorentzVector p4_WHad = fs_Wqrk0 + fs_Wqrk1;
        TLorentzVector p4_WLep = fs_lep0 + fs_lep1;        
        TLorentzVector p4_WW = p4_WHad + p4_WLep;
        
        double a_costheta1, a_costheta2, a_costhetastar, a_Phi, a_Phi1;
        computeAngles( p4_WW, p4_WLep, fs_lep0, fs_lep1, p4_WHad, fs_Wqrk0, fs_Wqrk1, 
                      a_costheta1, a_costheta2, a_Phi, a_costhetastar, a_Phi1);
        
        LHE_mWW_ = (float) p4_WW.M();
        LHE_mtWW_ = (float) p4_WW.Mt();
        LHE_mWLep_ = (float) p4_WLep.M();
        LHE_mtWLep_ = (float) p4_WLep.Mt();
        LHE_mWHad_ = (float) p4_WHad.M();        
        LHE_mtWHad_ = (float) p4_WHad.Mt();
        LHE_costheta1_ = (float) a_costheta1;                
        LHE_costheta2_ = (float) a_costheta2;
        LHE_phi_ = (float) a_Phi;
        LHE_costhetastar_ = (float) a_costhetastar;
        LHE_phi1_ = (float) a_Phi1;

        LHE_dEtajj_ = (float) fabs( fs_Iqrk0.Eta() - fs_Iqrk1.Eta() );
        LHE_dPhijj_ = (float) deltaPhi(fs_Iqrk0.Phi(),fs_Iqrk1.Phi());     
        LHE_mjj_ = (float) (fs_Iqrk0 + fs_Iqrk1).M();

	if (fabs(fs_Iqrk0.Eta()-fs_Iqrk1.Eta()) == 0.0)
		LHE_VBSCentrality_ = -999.0;
	else
		LHE_VBSCentrality_ = fabs(fs_Iqrk0.Eta()-(fs_Wqrk0.Rapidity()+fs_Wqrk1.Rapidity())-fs_Iqrk1.Eta())/fabs(fs_Iqrk0.Eta()-fs_Iqrk1.Eta());

        //isSignal = signalFlag;

	//initialQuarks_.clear();
		
	}
  }
  }



tree->Fill();
Clear();   
}

// ------------ method called once each job just before starting event loop  ------------
void 
GenAnalyzer::beginJob()
{
    std::cout<<"Inside beginJob()"<<std::endl;
    //outputFile_ = new TFile("aQGC_WPhadWMlepJJ_Madspin_EWK_MadDefCard_Pythia8_CUEP8M1_13TeV_Madgraph_NoMatching_New.root","RECREATE"); 
    outputFile_ = new TFile("LHEinfo.root","RECREATE"); 
    outputFile_->SetCompressionLevel(2);
    tree = new TTree("otree","GenParticles Basic Info"); 

    SetBranches();
}

// ------------ method called once each job just after ending the event loop  ------------
void 
GenAnalyzer::endJob() 
{
    outputFile_->Write();
    outputFile_->Close();
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
GenAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(GenAnalyzer);
