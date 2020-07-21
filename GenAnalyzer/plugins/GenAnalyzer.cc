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
#include <typeinfo>

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
  
  // edm::Handle<reco::GenJetCollection> genAK4jetHandle;
  // iEvent.getByToken(genAK4jetToken, genAK4jetHandle);
  
  // edm::Handle<reco::GenJetCollection> genAK8jetHandle;
  // iEvent.getByToken(genAK8jetToken, genAK8jetHandle);
  
  const vector<reco::GenParticle>* genps_coll = genpsHandle.product();
  // const vector<reco::GenJet>* genAK4_coll = genAK4jetHandle.product();
  // const vector<reco::GenJet>* genAK8_coll = genAK8jetHandle.product();
  
  if (Verbose_)  std::cout<<"===================\n\n"<<endl;
  
  int nGenParticle=0;
  
  TLorentzVector PHO;
  TLorentzVector Wpquarks;
  TLorentzVector Wmquarks;
  TLorentzVector Wboson1;
  TLorentzVector Wboson2;
  TLorentzVector Higgs1;
  TLorentzVector Higgs2;
  
  std::vector<TLorentzVector> Vec_Photons;
  std::vector<TLorentzVector> Vec_wpJET;
  std::vector<TLorentzVector> Vec_wmJET;
  std::vector<TLorentzVector> Vec_Wboson;
  std::vector<TLorentzVector> Vec_Higgs;
    std::cout << "Size = " << genps_coll->size() << std::endl;
  
  for(vector<reco::GenParticle>::const_iterator genps_it = genps_coll->begin(); genps_it != genps_coll->end(); genps_it++)
  {
    nGenParticle++;
    
    /*	Photon selection */
    if ( (abs(genps_it->pdgId())==11 || abs(genps_it->pdgId())==13 || abs(genps_it->pdgId())==15) && (abs(genps_it->mother()->pdgId())== 23)) {
    //if ( (abs(genps_it->pdgId())==11 || abs(genps_it->pdgId())==13 || abs(genps_it->pdgId())==15)  && (abs(genps_it->mother()->pdgId())==25) && genps_it->isHardProcess() ) {
      if (Verbose_) std::cout << "inphoton if condition" << std::endl;

      std::cout<< "==> " << genps_it->pdgId() << "\t" << genps_it->mother()->pdgId() << std::endl;
      
      PHO.SetPtEtaPhiE(genps_it->pt(), genps_it->eta(), genps_it->phi(), genps_it->energy());
      Vec_Photons.push_back(PHO);
    } /* End if conditon for photon selection */
    
    }
    std::cout<< "size of vec = " << Vec_Photons.size() << endl;
  tree->Fill();
  
  Clear();
}

// ------------ method called once each job just before starting event loop  ------------
void 
GenAnalyzer::beginJob()
{
  std::cout<<"Inside beginJob()"<<std::endl;
  outputFile_ = new TFile(OutPutFileName_,"RECREATE");
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
