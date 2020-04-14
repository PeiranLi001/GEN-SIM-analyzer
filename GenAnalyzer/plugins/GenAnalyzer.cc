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
	std::cout<<"size of LHEWeightIDS:\t"<<LHEWeightIDs_.size()<<std::endl;
	std::cout<<"size of LHEWeight: \t"<<LHEWeights_.size()<<std::endl;
	//std::cout<< " ID = " << LHEWeightIDs_[645] << "\t Weight = " << LHEWeights_[645] << std::endl;
	
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
            //if (LHE->hepeup().ISTUP.at (iPart) == 1){
            //    //PG leptons
            //    if (abs (LHE->hepeup().IDUP.at (iPart)) == 11 ||   //PG electron
            //        abs (LHE->hepeup().IDUP.at (iPart)) == 13 ||   //PG muon
            //        abs (LHE->hepeup().IDUP.at (iPart)) == 15 ||   //PG tau
            //        abs (LHE->hepeup().IDUP.at (iPart)) == 12 ||   //PG neutrino
            //        abs (LHE->hepeup().IDUP.at (iPart)) == 14 ||   //PG neutrino
            //        abs (LHE->hepeup().IDUP.at (iPart)) == 16)     //PG neutrino                    
            //        {
            //        leptons.push_back (iPart) ;
            //        } //PG leptons
            //    else
            //        {
            //        finalQuarks.push_back (iPart) ;
            //        }
            //    
            //} 
            
            //PG intermediates
            if (LHE->hepeup().ISTUP.at(iPart) == 2){
                intermediates.push_back (iPart) ;
            }
            
        } //PG loop over particles in the event
    
  }

  int nGenParticle=0;

  
  double l_pt=0., l_eta=0., l_phi=0., l_mass=0., l_mother=0.;
  int l_pdgId=0, l_status=0, l_gmother=0;

  double nu_pt=0., nu_eta=0., nu_phi=0., nu_mass=0.;
  int nu_pdgId=0, nu_status=0, nu_mother=0, nu_gmother=0;

  TLorentzVector PHO;
  TLorentzVector Wquarks;

  std::vector<TLorentzVector> photons;
  std::vector<TLorentzVector> wpJET;
  std::vector<TLorentzVector> wmJET;
  
  for(vector<reco::GenParticle>::const_iterator genps_it = genps_coll->begin(); genps_it != genps_coll->end(); genps_it++) 
  {
  	nGenParticle++;
	/*	Photon selection */
	if ( abs(genps_it->pdgId())==22 && (abs(genps_it->mother()->pdgId())==25) && genps_it->isHardProcess() && (genps_it->status()==23 || genps_it->status()==1) ) {
	  if (Verbose_) std::cout << "inphoton if condition" << std::endl;
	  //l_pdgId = genps_it->pdgId();
	  //l_status = genps_it->status();
	  //l_mother = genps_it->mother()->pdgId();
	  //l_gmother = genps_it->mother()->mother()->pdgId();

	  PHO.SetPtEtaPhiE(genps_it->pt(), genps_it->eta(), genps_it->phi(), genps_it->energy());
	  photons.push_back(PHO);
	} /* End if conditon for photon selection */

	/*	Quarks from W-boson	*/
	if ( abs(genps_it->pdgId())<7 && (genps_it->mother()->pdgId()==24) && genps_it->status()==23 && genps_it->isHardProcess() ) {
	  if (Verbose_) std::cout << "In quark loop... " << std::endl;
	  Wquarks.SetPtEtaPhiE(genps_it->pt(), genps_it->eta(), genps_it->phi(), genps_it->energy());
	  wpJET.push_back(Wquarks);
	} /* END if condition of wjETs */

	if ( abs(genps_it->pdgId())<7 && (genps_it->mother()->pdgId()==-24) && genps_it->status()==23 && genps_it->isHardProcess() ) {
	  if (Verbose_) std::cout << "In quark loop... " << std::endl;
	  Wquarks.SetPtEtaPhiE(genps_it->pt(), genps_it->eta(), genps_it->phi(), genps_it->energy());
	  wmJET.push_back(Wquarks);
	} /* END if condition of wjETs */
  }

  if (photons.size()==2 &&  wpJET.size()==2 && wmJET.size()==2)
  {
		if (Verbose_)
  	cout<<"Event selected"<<nEVENT<<endl;
	nEVENT++;

	gen_photon0_Pt_ 	= photons[0].Pt() ;
	gen_photon0_Eta_ 	= photons[0].Eta();
	gen_photon0_Phi_ 	= photons[0].Phi();
	gen_photon0_M_ 	= photons[0].M() ;

	gen_photon1_Pt_ 	= photons[1].Pt() ;
	gen_photon1_Eta_ 	= photons[1].Eta();
	gen_photon1_Phi_ 	= photons[1].Phi();
	gen_photon1_M_ 	= photons[1].M() ;

	gen_VBFjet1_Pt_ = wpJET[0].Pt();
	gen_VBFjet1_Eta_= wpJET[0].Eta();
	gen_VBFjet1_Phi_= wpJET[0].Phi();
	gen_VBFjet1_E_  = wpJET[0].E();

	gen_VBFjet2_Pt_ = wpJET[1].Pt();
	gen_VBFjet2_Eta_= wpJET[1].Eta();
	gen_VBFjet2_Phi_= wpJET[1].Phi();
	gen_VBFjet2_E_  = wpJET[1].E();

	//gen_vbfjet_deltaR_= deltaR(vJET[0].Eta(),vJET[0].Phi(),vJET[1].Eta(),vJET[1].Phi());
        //gen_dPhijj_ = (float) deltaPhi(vJET[0].Phi(),vJET[1].Phi());     
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
