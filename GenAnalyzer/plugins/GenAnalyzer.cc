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
   
   edm::Handle<reco::GenParticleCollection> genpsHandle;
   iEvent.getByToken(genParticlesToken, genpsHandle);
   
   const vector<reco::GenParticle>* genps_coll = genpsHandle.product();
   
   edm::Handle<LHEEventProduct> LHEEventHandle;
   iEvent.getByToken(LHEEventToken, LHEEventHandle);
   const LHEEventProduct* LHE = 0;
   
   std::vector<int> leptons ;      
   std::vector<int> finalQuarks ;      
   std::vector<int> intermediates ;
   std::vector<int> tops ;        
   
   if (Verbose_)  std::cout<<"===================\n\n"<<endl;
   
   if(LHEEventHandle.isValid()){
   	// clear the defined vectors before start
  	leptons.clear();
	finalQuarks.clear();
	intermediates.clear();
	tops.clear();
	int WDaughter_MothInfo1 = 0;
	int WDaughter_MothInfo2 = 0;
	int FQuark_MothInfo1 = 0;
	int FQuark_MothInfo2 = 0;

  	LHE = LHEEventHandle.product();

	for(const auto& weight : LHE->weights()) {
		//LHEWeightIDs_.push_back(weight.id);
		LHEWeights_.push_back(weight.wgt);
	}
	//std::cout<<"size of LHEWeightIDS:\t"<<LHEWeightIDs_.size()<<std::endl;
	//std::cout<<"size of LHEWeight: \t"<<LHEWeights_.size()<<std::endl;
	
	TLorentzVector Is_Iqrk1,Is_Iqrk0;
        //PG loop over particles in the event
	int incomingPart = 0;
	if (Verbose_) std::cout<<"Total No. of particles = "<< LHE->hepeup().NUP <<std::endl;
        for (int iPart = 0 ; iPart < LHE->hepeup().NUP; ++iPart){

		int mother1 = LHE->hepeup().MOTHUP[iPart].first;
		int mother2 = LHE->hepeup().MOTHUP[iPart].second;
		if (Verbose_)
		if (LHE->hepeup().ISTUP.at (iPart) != -1)
		std::cout<<"PDGID = "<<LHE->hepeup().IDUP.at(iPart)<<"\tStatus = "<< LHE->hepeup().ISTUP.at(iPart)<<"\tMother1 pos = "<<mother1<<"\t"<<mother2<<"\tPDGID = "<< LHE->hepeup().IDUP.at(mother1-1) <<"\t"<< LHE->hepeup().IDUP.at(mother2-1) <<std::endl;
		else
		std::cout<<"PDGID = "<<LHE->hepeup().IDUP.at(iPart)<<"\tStatus = "<< LHE->hepeup().ISTUP.at(iPart)<<std::endl;
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
            std::cout << "Warning!!!  Lepton size > 2" << std::endl;
        }
        
        // --------------- If signal, find the quarks from the W  -------------------
        if (finalQuarks.size() == 4){
            for (unsigned int a = 0; a < finalQuarks.size(); ++a ){
                int tmpMoth = LHE->hepeup().MOTHUP.at(finalQuarks.at(a)).first - 1;
		if (Verbose_)
			std::cout<<LHE->hepeup().IDUP.at(finalQuarks.at(a))<<"\tMOthere = "<<LHE->hepeup().IDUP.at(tmpMoth)<<std::endl;

                if ( abs(LHE->hepeup().IDUP.at(tmpMoth)) == 24 && LHE->hepeup().IDUP.at(finalQuarks.at(a)) > 0 ){
                    signalWCtr++;
		    WDaughter_MothInfo1 = finalQuarks.at(a);
		    if (Verbose_)	std::cout<<"First"<< WDaughter_MothInfo1 <<std::endl;
                }
                if ( abs(LHE->hepeup().IDUP.at(tmpMoth)) == 24 && LHE->hepeup().IDUP.at(finalQuarks.at(a)) < 0 ){
                    signalWCtr++;
		    WDaughter_MothInfo2 = finalQuarks.at(a);
		    if (Verbose_)	std::cout<<"Second"<<WDaughter_MothInfo2 <<std::endl;
                }
                if ( abs(LHE->hepeup().IDUP.at(tmpMoth)) != 24 && tmpfill1){
                    signalWCtr++;
		    FQuark_MothInfo1 = finalQuarks.at(a);
		    if (Verbose_)	std::cout<<"Third"<< FQuark_MothInfo1 <<std::endl;
                }
                if ( abs(LHE->hepeup().IDUP.at(tmpMoth)) != 24 && !tmpfill1){
                    signalWCtr++;
		    FQuark_MothInfo2 = finalQuarks.at(a);
		    if (Verbose_)	std::cout<<"Forth"<< FQuark_MothInfo2 <<std::endl;
                    tmpfill1++;
                }
            }
        }
        else{
            std::cout << "Warning!!! Final quark size > 4" << std::endl;
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
	finalQuarksNotW.clear();
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
	 	LHE_DeltaM_Wqrk0_E_  = LHE->hepeup().PUP[i_wqrk_1][3];
		LHE_DeltaM_Wqrk0_pt_ = fs_Wqrk0.Pt();
		LHE_DeltaM_Wqrk0_eta_= fs_Wqrk0.Eta();
		LHE_DeltaM_Wqrk0_phi_= fs_Wqrk0.Phi();
		LHE_DeltaM_Wqrk0_M_  = fs_Wqrk0.M();
		LHE_DeltaM_Wqrk0_Mt_  = fs_Wqrk0.Mt();

        TLorentzVector fs_Wqrk1
        (
         LHE->hepeup().PUP[i_wqrk_2][0], //PG px
         LHE->hepeup().PUP[i_wqrk_2][1], //PG py
         LHE->hepeup().PUP[i_wqrk_2][2], //PG pz
         LHE->hepeup().PUP[i_wqrk_2][3] //PG E
         ) ;
	 	LHE_DeltaM_Wqrk1_E_   = LHE->hepeup().PUP[i_wqrk_2][3];
		LHE_DeltaM_Wqrk1_pt_  = fs_Wqrk1.Pt();
		LHE_DeltaM_Wqrk1_eta_ = fs_Wqrk1.Eta();
		LHE_DeltaM_Wqrk1_phi_ = fs_Wqrk1.Phi();
		LHE_DeltaM_Wqrk1_M_   = fs_Wqrk1.M();
		LHE_DeltaM_Wqrk1_Mt_  = fs_Wqrk1.Mt();

        TLorentzVector fs_Iqrk0
        (
         LHE->hepeup().PUP[i_iqrk_1][0], //PG px
         LHE->hepeup().PUP[i_iqrk_1][1], //PG py
         LHE->hepeup().PUP[i_iqrk_1][2], //PG pz
         LHE->hepeup().PUP[i_iqrk_1][3] //PG E
         ) ;
	 	LHE_DeltaM_Iqrk0_E_   = LHE->hepeup().PUP[i_iqrk_1][3];
		LHE_DeltaM_Iqrk0_pt_  = fs_Iqrk0.Pt();
		LHE_DeltaM_Iqrk0_eta_ = fs_Iqrk0.Eta();
		LHE_DeltaM_Iqrk0_phi_ = fs_Iqrk0.Phi();
		LHE_DeltaM_Iqrk0_M_   = fs_Iqrk0.M();
		LHE_DeltaM_Iqrk0_Mt_  = fs_Iqrk0.Mt();


        TLorentzVector fs_Iqrk1
        (
         LHE->hepeup().PUP[i_iqrk_2][0], //PG px
         LHE->hepeup().PUP[i_iqrk_2][1], //PG py
         LHE->hepeup().PUP[i_iqrk_2][2], //PG pz
         LHE->hepeup().PUP[i_iqrk_2][3] //PG E
         ) ;
	 	LHE_DeltaM_Iqrk1_E_   = LHE->hepeup().PUP[i_iqrk_2][3];
		LHE_DeltaM_Iqrk1_pt_  = fs_Iqrk1.Pt();
		LHE_DeltaM_Iqrk1_eta_ = fs_Iqrk1.Eta();
		LHE_DeltaM_Iqrk1_phi_ = fs_Iqrk1.Phi();
		LHE_DeltaM_Iqrk1_M_   = fs_Iqrk1.M();
		LHE_DeltaM_Iqrk1_Mt_  = fs_Iqrk1.Mt();

        TLorentzVector p4_WHad = fs_Wqrk0 + fs_Wqrk1;
        TLorentzVector p4_WLep = fs_lep0 + fs_lep1;        
        TLorentzVector p4_WW = p4_WHad + p4_WLep;
        
        double a_costheta1, a_costheta2, a_costhetastar, a_Phi, a_Phi1;
        computeAngles( p4_WW, p4_WLep, fs_lep0, fs_lep1, p4_WHad, fs_Wqrk0, fs_Wqrk1, 
                      a_costheta1, a_costheta2, a_Phi, a_costhetastar, a_Phi1);
        
        LHE_DeltaM_mWW_ = (float) p4_WW.M();
        LHE_DeltaM_mtWW_ = (float) p4_WW.Mt();
        LHE_DeltaM_mWLep_ = (float) p4_WLep.M();
        LHE_DeltaM_mtWLep_ = (float) p4_WLep.Mt();
        LHE_DeltaM_mWHad_ = (float) p4_WHad.M();        
        LHE_DeltaM_mtWHad_ = (float) p4_WHad.Mt();
        LHE_DeltaM_costheta1_ = (float) a_costheta1;                
        LHE_DeltaM_costheta2_ = (float) a_costheta2;
        LHE_DeltaM_phi_ = (float) a_Phi;
        LHE_DeltaM_costhetastar_ = (float) a_costhetastar;
        LHE_DeltaM_phi1_ = (float) a_Phi1;

        LHE_DeltaM_dEtajj_ = (float) fabs( fs_Iqrk0.Eta() - fs_Iqrk1.Eta() );
        LHE_DeltaM_dPhijj_ = (float) deltaPhi(fs_Iqrk0.Phi(),fs_Iqrk1.Phi());     
        LHE_DeltaM_mjj_ = (float) (fs_Iqrk0 + fs_Iqrk1).M();

	if (fabs(fs_Iqrk0.Eta()-fs_Iqrk1.Eta()) == 0.0)
		LHE_DeltaM_VBSCentrality_ = -999.0;
	else
		LHE_DeltaM_VBSCentrality_ = fabs(fs_Iqrk0.Eta()-(fs_Wqrk0.Rapidity()+fs_Wqrk1.Rapidity())-fs_Iqrk1.Eta())/fabs(fs_Iqrk0.Eta()-fs_Iqrk1.Eta());
	}
	}

        if (leptons.size() == 2){
	if ( WDaughter_MothInfo1 != 0 && WDaughter_MothInfo2 != 0 &&  FQuark_MothInfo1 != 0 && FQuark_MothInfo2 != 0){

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

	TLorentzVector MothInfo_Wqrk1
	(
	  LHE->hepeup().PUP[WDaughter_MothInfo1][0],
	  LHE->hepeup().PUP[WDaughter_MothInfo1][1],
	  LHE->hepeup().PUP[WDaughter_MothInfo1][2],
	  LHE->hepeup().PUP[WDaughter_MothInfo1][3]
	);
	 	LHE_MothInfo_Wqrk0_E_   = LHE->hepeup().PUP[WDaughter_MothInfo1][3];
		LHE_MothInfo_Wqrk0_pt_  = MothInfo_Wqrk1.Pt();
		LHE_MothInfo_Wqrk0_eta_ = MothInfo_Wqrk1.Eta();
		LHE_MothInfo_Wqrk0_phi_ = MothInfo_Wqrk1.Phi();
		LHE_MothInfo_Wqrk0_M_   = MothInfo_Wqrk1.M();
		LHE_MothInfo_Wqrk0_Mt_  = MothInfo_Wqrk1.Mt();
		
	TLorentzVector MothInfo_Wqrk2
	(
	  LHE->hepeup().PUP[WDaughter_MothInfo2][0],
	  LHE->hepeup().PUP[WDaughter_MothInfo2][1],
	  LHE->hepeup().PUP[WDaughter_MothInfo2][2],
	  LHE->hepeup().PUP[WDaughter_MothInfo2][3]
	);
	 	LHE_MothInfo_Wqrk1_E_   = LHE->hepeup().PUP[WDaughter_MothInfo2][3];
		LHE_MothInfo_Wqrk1_pt_  = MothInfo_Wqrk2.Pt();
		LHE_MothInfo_Wqrk1_eta_ = MothInfo_Wqrk2.Eta();
		LHE_MothInfo_Wqrk1_phi_ = MothInfo_Wqrk2.Phi();
		LHE_MothInfo_Wqrk1_M_   = MothInfo_Wqrk2.M();
		LHE_MothInfo_Wqrk1_Mt_  = MothInfo_Wqrk2.Mt();
	TLorentzVector MothInfo_Fqrk1
	(
	  LHE->hepeup().PUP[FQuark_MothInfo1][0],
	  LHE->hepeup().PUP[FQuark_MothInfo1][1],
	  LHE->hepeup().PUP[FQuark_MothInfo1][2],
	  LHE->hepeup().PUP[FQuark_MothInfo1][3]
	);
	 	LHE_MothInfo_Iqrk0_E_   = LHE->hepeup().PUP[FQuark_MothInfo1][3];
		LHE_MothInfo_Iqrk0_pt_  = MothInfo_Fqrk1.Pt();
		LHE_MothInfo_Iqrk0_eta_ = MothInfo_Fqrk1.Eta();
		LHE_MothInfo_Iqrk0_phi_ = MothInfo_Fqrk1.Phi();
		LHE_MothInfo_Iqrk0_M_   = MothInfo_Fqrk1.M();
		LHE_MothInfo_Iqrk0_Mt_  = MothInfo_Fqrk1.Mt();
	TLorentzVector MothInfo_Fqrk2
	(
	  LHE->hepeup().PUP[FQuark_MothInfo2][0],
	  LHE->hepeup().PUP[FQuark_MothInfo2][1],
	  LHE->hepeup().PUP[FQuark_MothInfo2][2],
	  LHE->hepeup().PUP[FQuark_MothInfo2][3]
	);
	 	LHE_MothInfo_Iqrk1_E_   = LHE->hepeup().PUP[FQuark_MothInfo2][3];
		LHE_MothInfo_Iqrk1_pt_  = MothInfo_Fqrk2.Pt();
		LHE_MothInfo_Iqrk1_eta_ = MothInfo_Fqrk2.Eta();
		LHE_MothInfo_Iqrk1_phi_ = MothInfo_Fqrk2.Phi();
		LHE_MothInfo_Iqrk1_M_   = MothInfo_Fqrk2.M();
		LHE_MothInfo_Iqrk1_Mt_  = MothInfo_Fqrk2.Mt();

		
        TLorentzVector p4_WHad_MothInfo = MothInfo_Wqrk1 + MothInfo_Wqrk2;
	TLorentzVector p4_WLep_MothInfo = fs_lep0 + fs_lep1;
        TLorentzVector p4_WLep = fs_lep0 + fs_lep1;        
        TLorentzVector p4_WW_MothInfo = p4_WHad_MothInfo + p4_WLep;

	if (MothInfo_Fqrk1.Pt()>30 && MothInfo_Fqrk2.Pt()>30 && fabs(MothInfo_Fqrk1.Eta())<5.0 && fabs( MothInfo_Fqrk2.Eta()) < 5.0 && p4_WHad_MothInfo.Pt() > 200 && fabs( p4_WHad_MothInfo.Eta()) < 2.0 &&  (p4_WHad_MothInfo.M()>60 && p4_WHad_MothInfo.M()<110) && (MothInfo_Fqrk1+MothInfo_Fqrk2).M()>800 && fabs(MothInfo_Fqrk1.Eta() - MothInfo_Fqrk2.Eta())>4.0 )
	{
	//cout<<p4_WLep_MothInfo.Px() << "\t" << p4_WLep_MothInfo.Py() << "\t" << p4_WLep_MothInfo.Pz() << "\t" << p4_WLep_MothInfo.E() << "\t" << p4_WHad_MothInfo.Px() << "\t" << p4_WHad_MothInfo.Py() << "\t" << p4_WHad_MothInfo.Pz() << "\t" << p4_WHad_MothInfo.E() << "\t" << MothInfo_Fqrk1.Px() << "\t" << MothInfo_Fqrk1.Py() << "\t" << MothInfo_Fqrk1.Pz() << "\t" << MothInfo_Fqrk1.E() << "\t" <<  MothInfo_Fqrk2.Px() << "\t" << MothInfo_Fqrk2.Py() << "\t" << MothInfo_Fqrk2.Pz() << "\t" << MothInfo_Fqrk2.E();
	file1<<p4_WLep_MothInfo.Px() << "\t" << p4_WLep_MothInfo.Py() << "\t" << p4_WLep_MothInfo.Pz() << "\t" << p4_WLep_MothInfo.E() << "\t" << p4_WHad_MothInfo.Px() << "\t" << p4_WHad_MothInfo.Py() << "\t" << p4_WHad_MothInfo.Pz() << "\t" << p4_WHad_MothInfo.E() << "\t" << MothInfo_Fqrk1.Px() << "\t" << MothInfo_Fqrk1.Py() << "\t" << MothInfo_Fqrk1.Pz() << "\t" << MothInfo_Fqrk1.E() << "\t" <<  MothInfo_Fqrk2.Px() << "\t" << MothInfo_Fqrk2.Py() << "\t" << MothInfo_Fqrk2.Pz() << "\t" << MothInfo_Fqrk2.E()<<endl;
	}
        
        double a_costheta1, a_costheta2, a_costhetastar, a_Phi, a_Phi1;
        computeAngles( p4_WW_MothInfo, p4_WLep, fs_lep0, fs_lep1, p4_WHad_MothInfo, MothInfo_Wqrk1, MothInfo_Wqrk2, 
                      a_costheta1, a_costheta2, a_Phi, a_costhetastar, a_Phi1);
        
        LHE_MothInfo_mWW_ = (float) p4_WW_MothInfo.M();
        LHE_MothInfo_mtWW_ = (float) p4_WW_MothInfo.Mt();
        LHE_MothInfo_mWLep_ = (float) p4_WLep.M();
        LHE_MothInfo_mtWLep_ = (float) p4_WLep.Mt();
        LHE_MothInfo_mWHad_ = (float) p4_WHad_MothInfo.M();        
        LHE_MothInfo_mtWHad_ = (float) p4_WHad_MothInfo.Mt();
        LHE_MothInfo_costheta1_ = (float) a_costheta1;                
        LHE_MothInfo_costheta2_ = (float) a_costheta2;
        LHE_MothInfo_phi_ = (float) a_Phi;
        LHE_MothInfo_costhetastar_ = (float) a_costhetastar;
        LHE_MothInfo_phi1_ = (float) a_Phi1;

        LHE_MothInfo_dEtajj_ = (float) fabs( MothInfo_Fqrk1.Eta() - MothInfo_Fqrk2.Eta() );
        LHE_MothInfo_dPhijj_ = (float) deltaPhi(MothInfo_Fqrk1.Phi(),MothInfo_Fqrk2.Phi());     
        LHE_MothInfo_mjj_ = (float) (MothInfo_Fqrk1 + MothInfo_Fqrk2).M();

	if (fabs(MothInfo_Fqrk1.Eta()-MothInfo_Fqrk2.Eta()) == 0.0)
		LHE_MothInfo_VBSCentrality_ = -999.0;
	else
		LHE_MothInfo_VBSCentrality_ = fabs(MothInfo_Fqrk1.Eta()-(MothInfo_Wqrk1.Rapidity()+MothInfo_Wqrk2.Rapidity())-MothInfo_Fqrk2.Eta())/fabs(MothInfo_Fqrk1.Eta()-MothInfo_Fqrk2.Eta());

        //isSignal = signalFlag;

	//initialQuarks_.clear();
		
	}
	}
  }

  int nGenParticle=0;

  
  double l_pt=0., l_eta=0., l_phi=0., l_mass=0., l_mother=0.;
  int l_pdgId=0, l_status=0, l_gmother=0;

  double nu_pt=0., nu_eta=0., nu_phi=0., nu_mass=0.;
  int nu_pdgId=0, nu_status=0, nu_mother=0, nu_gmother=0;

  TLorentzVector ELE, MU, TAU;
  TLorentzVector NU;
  TLorentzVector JET, VBFquarks, Wquarks;
  TLorentzVector VBF1, VBF2, VBFTOT;
  TLorentzVector Wjet1, Wjet2, WjetTOT;

  std::vector<TLorentzVector> tightLep, vNU, vJET, wJET;
  
  for(vector<reco::GenParticle>::const_iterator genps_it = genps_coll->begin(); genps_it != genps_coll->end(); genps_it++) 
  {
  	nGenParticle++;
	if((abs(genps_it->pdgId())==11 || abs(genps_it->pdgId())==13 || abs(genps_it->pdgId())==15) && (abs(genps_it->mother()->pdgId()) == 24) && genps_it->isHardProcess() )
	{
		if (Verbose_)
		cout<<"Status of leptons = "<<genps_it->status()<<endl;
		if (Verbose_)
		cout<<"Status of leptons monther = "<<genps_it->mother()->status()<<endl;
		//if (genps_it->pt() < 45.0) continue;
		//if (fabs(genps_it->eta()) > 2.1) continue;

		l_pdgId = genps_it->pdgId();
		l_status = genps_it->status();
		l_mother = genps_it->mother()->pdgId();
		l_gmother = genps_it->mother()->mother()->pdgId();

		ELE.SetPtEtaPhiE(genps_it->pt(), genps_it->eta(), genps_it->phi(), genps_it->energy());
		tightLep.push_back(ELE);
	}
	if((abs(genps_it->pdgId())==12 || abs(genps_it->pdgId())==14 || abs(genps_it->pdgId())==16) && (abs(genps_it->mother()->pdgId()) == 24)  && genps_it->isHardProcess() )
	{
		if (Verbose_)
		cout<<"Status of Neutrino = "<<genps_it->status()<<endl;

		nu_pdgId   = genps_it->pdgId();
		nu_status  = genps_it->status();
		nu_mother  = genps_it->mother()->pdgId();
		nu_gmother = genps_it->mother()->mother()->pdgId();

		NU.SetPtEtaPhiE(genps_it->pt(), genps_it->eta(), genps_it->phi(), genps_it->energy());
		vNU.push_back(NU);
	}
	if((abs(genps_it->pdgId())==1 || abs(genps_it->pdgId())==2 || abs(genps_it->pdgId())==3 || abs(genps_it->pdgId())==4 || abs(genps_it->pdgId())==5 || abs(genps_it->pdgId())==6) && (abs(genps_it->mother()->pdgId()) == 24) && (genps_it->status()==23) && genps_it->isHardProcess() )
	{
		Wquarks.SetPtEtaPhiE(genps_it->pt(), genps_it->eta(), genps_it->phi(), genps_it->energy());
		wJET.push_back(Wquarks);
	}
	if((abs(genps_it->pdgId())==1 || abs(genps_it->pdgId())==2 || abs(genps_it->pdgId())==3 || abs(genps_it->pdgId())==4 || abs(genps_it->pdgId())==5 || abs(genps_it->pdgId())==6) && (abs(genps_it->mother()->pdgId()) != 24) && (genps_it->status()==23) )
	{
		VBFquarks.SetPtEtaPhiE(genps_it->pt(), genps_it->eta(), genps_it->phi(), genps_it->energy());
		vJET.push_back(VBFquarks);
	}
  }

  if (tightLep.size()==1 && vNU.size()==1 && wJET.size()==2 && vJET.size()==2)
  {
		if (Verbose_)
  	cout<<"Event selected"<<nEVENT<<endl;
	nEVENT++;

	gen_LeptPt_ 	= ELE.Pt() ;
	gen_LeptEta_ 	= ELE.Eta();
	gen_LeptPhi_ 	= ELE.Phi();
	gen_LeptM_ 	= ELE.M() ;
	gen_LeptId_ 	= l_pdgId ;
	gen_LeptStatus_ 	= l_status ;
	gen_LeptMother_ 	= l_mother ;
	gen_LeptGrandMother_ = l_gmother;

	gen_NuPt_ 	= NU.Pt(); 
	gen_NuEta_ 	= NU.Eta();
	gen_NuPhi_ 	= NU.Phi();
	gen_NuM_ 	= NU.M();
	//gen_NuQ_ 	= ;
	gen_Nustatus_ 	= nu_status;
	gen_NuMother_ 	= nu_mother;
	gen_NuGrandMother_ = nu_gmother;
	gen_NuPdgId_ 	= nu_pdgId;	

	gen_VBFjet1_Pt_ = vJET[0].Pt();
	gen_VBFjet1_Eta_= vJET[0].Eta();
	gen_VBFjet1_Phi_= vJET[0].Phi();
	gen_VBFjet1_E_  = vJET[0].E();

	gen_VBFjet2_Pt_ = vJET[1].Pt();
	gen_VBFjet2_Eta_= vJET[1].Eta();
	gen_VBFjet2_Phi_= vJET[1].Phi();
	gen_VBFjet2_E_	 = vJET[1].E();

	gen_VBFjet1jet2_Pt_ = (vJET[0]+vJET[1]).Pt();
	gen_VBFjet1jet2_Eta_= (vJET[0]+vJET[1]).Eta();
	gen_VBFjet1jet2_Phi_= (vJET[0]+vJET[1]).Phi();
	gen_VBFjet1jet2_M_  = (vJET[0]+vJET[1]).M();
	gen_vbfjet_deltaR_= deltaR(vJET[0].Eta(),vJET[0].Phi(),vJET[1].Eta(),vJET[1].Phi());

	gen_WJet1_Pt_	= wJET[0].Pt();
	gen_WJet1_Eta_	= wJET[0].Eta();
	gen_WJet1_Phi_	= wJET[0].Phi();
	gen_WJet1_E_	= wJET[0].E();

	gen_WJet2_Pt_	= wJET[1].Pt();
	gen_WJet2_Eta_	= wJET[1].Eta();
	gen_WJet2_Phi_	= wJET[1].Phi();
	gen_WJet2_E_	= wJET[1].E();

	gen_WHad_Pt_	= (wJET[0]+wJET[1]).Pt();
	gen_WHad_M_= (wJET[0]+wJET[1]).M();
	gen_WHad_Mt_= (wJET[0]+wJET[1]).Mt();
	gen_WHad_deltaeta_= fabs(wJET[0].Eta()-wJET[1].Eta());
	gen_WHad_deltaphi_= fabs(wJET[0].Phi()-wJET[1].Phi());
	gen_WHad_deltar_= deltaR(wJET[0].Eta(),wJET[0].Phi(),wJET[1].Eta(),wJET[1].Phi());
	gen_deltaR_LepWHad_= deltaR((wJET[0]+wJET[1]).Eta(),(wJET[0]+wJET[1]).Phi(),l_eta, l_phi);
	gen_deltaphi_NuWHad_= fabs((wJET[0]+wJET[1]).Phi() - nu_phi);
	gen_deltaphi_WlepWHad_= fabs((wJET[0]+wJET[1]).Phi() - (ELE + NU).Phi());

		
        TLorentzVector p4_WHad_GEN = wJET[0] + wJET[1];
        TLorentzVector p4_WLep_GEN = ELE + NU;
        TLorentzVector p4_WW_GEN = p4_WHad_GEN + p4_WLep_GEN;
        
        double a_costheta1, a_costheta2, a_costhetastar, a_Phi, a_Phi1;
        computeAngles( p4_WW_GEN, p4_WLep_GEN, ELE, NU, p4_WHad_GEN, wJET[0], wJET[1], 
                      a_costheta1, a_costheta2, a_Phi, a_costhetastar, a_Phi1);
        
	gen_mWW_  = (float) p4_WW_GEN.M();
	gen_mtWW_ = (float) p4_WW_GEN.Mt();
        gen_mWLep_ = (float) p4_WLep_GEN.M();
        gen_mtWLep_ = (float) p4_WLep_GEN.Mt();
        gen_mWHad_ = (float) p4_WHad_GEN.M();        
        gen_mtWHad_ = (float) p4_WHad_GEN.Mt();
        gen_costheta1_ = (float) a_costheta1;                
        gen_costheta2_ = (float) a_costheta2;
        gen_phi_ = (float) a_Phi;
        gen_costhetastar_ = (float) a_costhetastar;
        gen_phi1_ = (float) a_Phi1;

        gen_dEtajj_ = (float) fabs( vJET[0].Eta() - vJET[1].Eta() );
        gen_dPhijj_ = (float) deltaPhi(vJET[0].Phi(),vJET[1].Phi());     
        gen_mjj_ = (float) (vJET[0] + vJET[1]).M();

	if (fabs(vJET[0].Eta()-vJET[1].Eta()) == 0.0)
		gen_VBSCentrality_ = -999.0;
	else
		gen_VBSCentrality_ = fabs(vJET[0].Eta()-(wJET[0].Rapidity()+wJET[1].Rapidity())-vJET[1].Eta())/fabs(vJET[0].Eta()-vJET[1].Eta());

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
    file1.open("out_TEMP_NAME.txt");

    SetBranches();
}

// ------------ method called once each job just after ending the event loop  ------------
void 
GenAnalyzer::endJob() 
{
    outputFile_->Write();
    outputFile_->Close();
    file1.close();
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
