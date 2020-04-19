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
  
  edm::Handle<reco::GenJetCollection> genAK4jetHandle;
  iEvent.getByToken(genAK4jetToken, genAK4jetHandle);
  
  edm::Handle<reco::GenJetCollection> genAK8jetHandle;
  iEvent.getByToken(genAK8jetToken, genAK8jetHandle);
  
  const vector<reco::GenParticle>* genps_coll = genpsHandle.product();
  const vector<reco::GenJet>* genAK4_coll = genAK4jetHandle.product();
  const vector<reco::GenJet>* genAK8_coll = genAK8jetHandle.product();
  
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
  
  //  std::cout << "Size = " << genps_coll->size() << std::endl;
  
  for(vector<reco::GenParticle>::const_iterator genps_it = genps_coll->begin(); genps_it != genps_coll->end(); genps_it++)
  {
    nGenParticle++;
    
    /*	Photon selection */
    if ( abs(genps_it->pdgId())==22 && (abs(genps_it->mother()->pdgId())==25) && genps_it->isHardProcess() && (genps_it->status()==23 || genps_it->status()==1) ) {
      if (Verbose_) std::cout << "inphoton if condition" << std::endl;
      
      PHO.SetPtEtaPhiE(genps_it->pt(), genps_it->eta(), genps_it->phi(), genps_it->energy());
      Vec_Photons.push_back(PHO);
    } /* End if conditon for photon selection */
    
    /*	Quarks from W-boson	*/
    if ( abs(genps_it->pdgId())<7 && (genps_it->mother()->pdgId()==24) && genps_it->status()==23 && genps_it->isHardProcess() ) {
      if (Verbose_) std::cout << "In W+ loop... " << std::endl;
      Wpquarks.SetPtEtaPhiE(genps_it->pt(), genps_it->eta(), genps_it->phi(), genps_it->energy());
      Vec_wpJET.push_back(Wpquarks);
    } /* END if condition of wjETs */
    
    /*  Quarks from W- boson  */
    if ( abs(genps_it->pdgId())<7 && (genps_it->mother()->pdgId()==-24) && genps_it->status()==23 && genps_it->isHardProcess() ) {
      if (Verbose_) std::cout << "In W- loop... " << std::endl;
      Wmquarks.SetPtEtaPhiE(genps_it->pt(), genps_it->eta(), genps_it->phi(), genps_it->energy());
      Vec_wmJET.push_back(Wmquarks);
    } /* END if condition of wjETs */
  }
  
  if (Vec_Photons.size()==2 &&  Vec_wpJET.size()==2 && Vec_wmJET.size()==2)
  {
    if (Verbose_) cout<<"Event selected"<<nEVENT<<endl;
    nEVENT++;
    
    /*  Form the W-bosons 4-Vectors   */
    Vec_Wboson.push_back(Vec_wmJET[0]+Vec_wmJET[1]);
    Vec_Wboson.push_back(Vec_wpJET[0]+Vec_wpJET[1]);
    
    /*  Form the Higgs 4-Vectors  */
    Vec_Higgs.push_back(Vec_Wboson[0]+Vec_Wboson[1]);
    Vec_Higgs.push_back(Vec_Photons[0]+Vec_Photons[1]);
    
    /*  Sort each TLorentzVectors */
    std::sort(Vec_Photons.begin(), Vec_Photons.end(), GenAnalyzer::reorder);
    std::sort(Vec_wpJET.begin(), Vec_wpJET.end(), GenAnalyzer::reorder);
    std::sort(Vec_wmJET.begin(), Vec_wmJET.end(), GenAnalyzer::reorder);
    std::sort(Vec_Wboson.begin(), Vec_Wboson.end(), GenAnalyzer::reorder);
    /*	We did not sort the Higgs boson as I need to know which Higgs is reconstructed from W-boson and which one from photons	*/
    //std::sort(Vec_Higgs.begin(), Vec_Higgs.end(), GenAnalyzer::reorder);
    
    /************************************************************************
     **
     **         Get 4-moment for all final state GenParticles
     **
     *************************************************************************/
    gen_leading_photon_Pt_ 	= Vec_Photons[0].Pt() ;
    gen_leading_photon_Eta_ 	= Vec_Photons[0].Eta();
    gen_leading_photon_Phi_ 	= Vec_Photons[0].Phi();
    gen_leading_photon_M_ 	= Vec_Photons[0].M() ;
    
    gen_Subleading_photon_Pt_  = Vec_Photons[1].Pt() ;
    gen_Subleading_photon_Eta_     = Vec_Photons[1].Eta();
    gen_Subleading_photon_Phi_     = Vec_Photons[1].Phi();
    gen_Subleading_photon_M_   = Vec_Photons[1].M() ;
    
    gen_leading_WpJets_Pt_  = Vec_wpJET[0].Pt();
    gen_leading_WpJets_Eta_  = Vec_wpJET[0].Eta();
    gen_leading_WpJets_Phi_  = Vec_wpJET[0].Phi();
    gen_leading_WpJets_M_  = Vec_wpJET[0].M();
    
    gen_Subleading_WpJets_Pt_  = Vec_wpJET[1].Pt();
    gen_Subleading_WpJets_Eta_  = Vec_wpJET[1].Eta();
    gen_Subleading_WpJets_Phi_  = Vec_wpJET[1].Phi();
    gen_Subleading_WpJets_M_  = Vec_wpJET[1].M();
    
    gen_leading_WmJets_Pt_  = Vec_wmJET[0].Pt();
    gen_leading_WmJets_Eta_  = Vec_wmJET[0].Eta();
    gen_leading_WmJets_Phi_  = Vec_wmJET[0].Phi();
    gen_leading_WmJets_M_  = Vec_wmJET[0].M();
    
    gen_Subleading_WmJets_Pt_  = Vec_wmJET[1].Pt();
    gen_Subleading_WmJets_Eta_  = Vec_wmJET[1].Eta();
    gen_Subleading_WmJets_Phi_  = Vec_wmJET[1].Phi();
    gen_Subleading_WmJets_M_  = Vec_wmJET[1].M();
    
    gen_leading_WBoson_Pt_  = Vec_Wboson[0].Pt();
    gen_leading_WBoson_Eta_  = Vec_Wboson[0].Eta();
    gen_leading_WBoson_Phi_  = Vec_Wboson[0].Phi();
    gen_leading_WBoson_M_  = Vec_Wboson[0].M();
    
    gen_Subleading_WBoson_Pt_  = Vec_Wboson[1].Pt();
    gen_Subleading_WBoson_Eta_  = Vec_Wboson[1].Eta();
    gen_Subleading_WBoson_Phi_  = Vec_Wboson[1].Phi();
    gen_Subleading_WBoson_M_  = Vec_Wboson[1].M();
    
    gen_leading_Higgs_Pt_  = Vec_Higgs[0].Pt();
    gen_leading_Higgs_Eta_  = Vec_Higgs[0].Eta();
    gen_leading_Higgs_Phi_  = Vec_Higgs[0].Phi();
    gen_leading_Higgs_M_  = Vec_Higgs[0].M();
    
    gen_Subleading_Higgs_Pt_  = Vec_Higgs[1].Pt();
    gen_Subleading_Higgs_Eta_  = Vec_Higgs[1].Eta();
    gen_Subleading_Higgs_Phi_  = Vec_Higgs[1].Phi();
    gen_Subleading_Higgs_M_  = Vec_Higgs[1].M();
    
    
    /************************************************************************
     **
     **         Get Delta-R between each stage particles
     **
     *************************************************************************/
    gen_deltaR_Photon0_Photon1_ = deltaR(Vec_Photons[0].Eta(),Vec_Photons[0].Phi(), Vec_Photons[1].Eta(),Vec_Photons[1].Phi());
    
    gen_deltaR_Photon0_WmJ0_ = deltaR(Vec_Photons[0].Eta(),Vec_Photons[0].Phi(), Vec_wmJET[0].Eta(),Vec_wmJET[0].Phi());
    gen_deltaR_Photon0_WmJ1_ = deltaR(Vec_Photons[0].Eta(),Vec_Photons[0].Phi(), Vec_wmJET[1].Eta(),Vec_wmJET[1].Phi());
    gen_deltaR_Photon0_WpJ0_ = deltaR(Vec_Photons[0].Eta(),Vec_Photons[0].Phi(), Vec_wpJET[1].Eta(),Vec_wpJET[1].Phi());
    gen_deltaR_Photon0_WpJ1_ = deltaR(Vec_Photons[0].Eta(),Vec_Photons[0].Phi(), Vec_wpJET[1].Eta(),Vec_wpJET[1].Phi());
    
    gen_deltaR_Photon1_WmJ0_ = deltaR(Vec_Photons[1].Eta(),Vec_Photons[1].Phi(), Vec_wmJET[0].Eta(),Vec_wmJET[0].Phi());
    gen_deltaR_Photon1_WmJ1_ = deltaR(Vec_Photons[1].Eta(),Vec_Photons[1].Phi(), Vec_wmJET[1].Eta(),Vec_wmJET[1].Phi());
    gen_deltaR_Photon1_WpJ0_ = deltaR(Vec_Photons[1].Eta(),Vec_Photons[1].Phi(), Vec_wpJET[1].Eta(),Vec_wpJET[1].Phi());
    gen_deltaR_Photon1_WpJ1_ = deltaR(Vec_Photons[1].Eta(),Vec_Photons[1].Phi(), Vec_wpJET[1].Eta(),Vec_wpJET[1].Phi());
    
    gen_deltaR_WpJ0_WpJ1_    = deltaR(Vec_wpJET[0].Eta(),Vec_wpJET[0].Phi(), Vec_wpJET[1].Eta(),Vec_wpJET[1].Phi());
    gen_deltaR_WpJ0_WmJ0_    = deltaR(Vec_wpJET[0].Eta(),Vec_wpJET[0].Phi(), Vec_wmJET[0].Eta(),Vec_wmJET[0].Phi());
    gen_deltaR_WpJ0_WmJ1_    = deltaR(Vec_wpJET[0].Eta(),Vec_wpJET[0].Phi(), Vec_wmJET[1].Eta(),Vec_wmJET[1].Phi());
    gen_deltaR_WpJ1_WmJ0_    = deltaR(Vec_wpJET[1].Eta(),Vec_wpJET[1].Phi(), Vec_wmJET[0].Eta(),Vec_wmJET[0].Phi());
    gen_deltaR_WpJ1_WmJ1_    = deltaR(Vec_wpJET[1].Eta(),Vec_wpJET[1].Phi(), Vec_wmJET[1].Eta(),Vec_wmJET[1].Phi());
    gen_deltaR_WmJ0_WmJ1_    = deltaR(Vec_wmJET[0].Eta(),Vec_wmJET[0].Phi(), Vec_wmJET[1].Eta(),Vec_wmJET[1].Phi());
    
    gen_deltaR_Wp_Wm_    = deltaR(Vec_Wboson[0].Eta(), Vec_Wboson[0].Phi(), Vec_Wboson[1].Eta(), Vec_Wboson[1].Phi());
    gen_deltaR_H1_H2_    = deltaR(Vec_Higgs[0].Eta(), Vec_Higgs[0].Phi(), Vec_Higgs[1].Eta(), Vec_Higgs[1].Phi());
    //gen_dPhijj_ = (float) deltaPhi(vJET[0].Phi(),vJET[1].Phi());
    
  } else {
    std::cout << "WARNING:: the condition 2 photons and 2 jets from each w's are not satisfied...." << std::endl;
  }
  
  //genAK4_coll
  int nAK4jets = 0;
  double temp_AK4_jet_pt = 0.0;
  TLorentzVector genJetAK4;
  std::vector<TLorentzVector> Vec_genJetAK4;
  //std::vector<ROOT::Math::XYZTVector> Vec_genJetAK4;
  //std::vector<typename ROOT::Math::LorentzVector> Vec_genJetAK4;
  //std::vector<const TLorentzVector &> Vec_genJetAK4;
  for(vector<reco::GenJet>::const_iterator genjet = genAK4_coll->begin(); genjet != genAK4_coll->end(); genjet++) {
    if (deltaR(genjet->eta(),genjet->phi(), Vec_Photons[0].Eta(),Vec_Photons[0].Phi())>0.4 && deltaR(genjet->eta(),genjet->phi(), Vec_Photons[1].Eta(),Vec_Photons[1].Phi())>0.4 && genjet->pt()>10) {
      nAK4jets++;
      //std::cout << genjet->pt() << "\t"
      //<< Vec_wpJET[0].Pt() << "\t" << Vec_wpJET[1].Pt() << "\t"<< Vec_wmJET[0].Pt() << "\t" << Vec_wmJET[1].Pt() << "\t" << deltaR(genjet->eta(),genjet->phi(), Vec_wpJET[0].Eta(),Vec_wpJET[0].Phi()) << "\t"
      //<< deltaR(genjet->eta(),genjet->phi(), Vec_wpJET[1].Eta(),Vec_wpJET[1].Phi()) << "\t"
      //<< deltaR(genjet->eta(),genjet->phi(), Vec_wmJET[0].Eta(),Vec_wmJET[0].Phi()) << "\t"
      //<< deltaR(genjet->eta(),genjet->phi(), Vec_wmJET[1].Eta(),Vec_wmJET[1].Phi()) << std::endl;
      //std::cout<<"====> " << typeid(genjet->p4()).name() << std::endl;
      //std::cout<<"====> " << typeid(Vec_wmJET).name() << std::endl;
      if (genjet->pt()>temp_AK4_jet_pt) {
        genJetAK4.SetPtEtaPhiE(genjet->pt(), genjet->eta(), genjet->phi(), genjet->energy());
        Vec_genJetAK4.push_back(genJetAK4);
        //Vec_genJetAK4.push_back(genjet->p4());
        temp_AK4_jet_pt = genjet->pt();
      }
    }
    //    std::cout << genjet->pt() << std::endl;
  }
  /* Sort AK4 TLorentzVector */
  //std::sort(Vec_genJetAK4.begin(), Vec_genJetAK4.end(), GenAnalyzer::reorder);
  //std::cout << "Number of AK4 candidates = " << nAK4jets << std::endl;
  
  genJetAK4_njets_ = nAK4jets;
  
  //genJetAK4_leading_Pt_ = Vec_genJetAK4[0].Pt();
  //genJetAK4_leading_Eta_ = Vec_genJetAK4[0].Eta();
  //genJetAK4_leading_Phi_ = Vec_genJetAK4[0].Phi();
  //genJetAK4_leading_M_ = Vec_genJetAK4[0].M();
  //genJetAK4_leading_Energy_ = Vec_genJetAK4[0].Energy();
  
  //genJetAK4_Subleading_Pt_ = Vec_genJetAK4[0].Pt();
  //genJetAK4_Subleading_Eta_ = Vec_genJetAK4[0].Eta();
  //genJetAK4_Subleading_Phi_ = Vec_genJetAK4[0].Phi();
  //genJetAK4_Subleading_M_ = Vec_genJetAK4[0].M();
  //genJetAK4_Subleading_Energy_ = Vec_genJetAK4[0].Energy();
  
  int nAK8jets = 0;
  double temp_AK8jet_pt = -999.0;
  double temp_AK8jet_deltaM = 9999.0;
  TLorentzVector genJetAK8;
  TLorentzVector genJetAK8_minDMass;
  for(vector<reco::GenJet>::const_iterator genjet = genAK8_coll->begin(); genjet != genAK8_coll->end(); genjet++) {
    if (deltaR(genjet->eta(),genjet->phi(), Vec_Photons[0].Eta(),Vec_Photons[0].Phi())>0.4 && deltaR(genjet->eta(),genjet->phi(), Vec_Photons[1].Eta(),Vec_Photons[1].Phi())>0.4) {
      if ( jetCleaning(&(*genjet), genAK4_coll ) )
      {
        nAK8jets++;
        //std::cout << genjet->pt() << "\t"
        //<< Vec_wpJET[0].Pt() << "\t" << Vec_wpJET[1].Pt() << "\t" << Vec_wmJET[0].Pt() << "\t" << Vec_wmJET[1].Pt() << "\t" << deltaR(genjet->eta(),genjet->phi(), Vec_wpJET[0].Eta(),Vec_wpJET[0].Phi()) << "\t"
        //<< deltaR(genjet->eta(),genjet->phi(), Vec_wpJET[1].Eta(),Vec_wpJET[1].Phi()) << "\t"
        //<< deltaR(genjet->eta(),genjet->phi(), Vec_wmJET[0].Eta(),Vec_wmJET[0].Phi()) << "\t"
        //<< deltaR(genjet->eta(),genjet->phi(), Vec_wmJET[1].Eta(),Vec_wmJET[1].Phi()) << std::endl;
        
        if (genjet->pt()>temp_AK8jet_pt) {
          genJetAK8.SetPtEtaPhiE(genjet->pt(), genjet->eta(), genjet->phi(), genjet->energy());
          temp_AK8jet_pt = genjet->pt();
        }
        if ( abs(genjet->mass() - 125.0) < temp_AK8jet_deltaM)
        {
          genJetAK8_minDMass.SetPtEtaPhiE(genjet->pt(), genjet->eta(), genjet->phi(), genjet->energy());
          temp_AK8jet_deltaM = abs(genjet->mass() - 125.0);
        }
      }
    }
  }
  //std::cout << "Number of AK8 candidates = " << nAK8jets << std::endl;
  
  genJetAK8_njets_ = nAK8jets;
  genJetAK8_MaxPt_Pt_ = genJetAK8.Pt();
  genJetAK8_MaxPt_Eta_ = genJetAK8.Eta();
  genJetAK8_MaxPt_Phi_ = genJetAK8.Phi();
  genJetAK8_MaxPt_M_ = genJetAK8.M();
  genJetAK8_MaxPt_deltaR_H1_ = deltaR(genJetAK8.Eta(),genJetAK8.Phi(), Vec_Higgs[0].Eta(), Vec_Higgs[0].Phi());
  genJetAK8_MaxPt_deltaR_H2_ = deltaR(genJetAK8.Eta(),genJetAK8.Phi(), Vec_Higgs[1].Eta(), Vec_Higgs[1].Phi());
  
  genJetAK8_minDMass_Pt_ = genJetAK8_minDMass.Pt();
  genJetAK8_minDMass_Eta_ = genJetAK8_minDMass.Eta();
  genJetAK8_minDMass_Phi_ = genJetAK8_minDMass.Phi();
  genJetAK8_minDMass_M_ = genJetAK8_minDMass.M();
  genJetAK8_minDMass_deltaR_H1_ = deltaR(genJetAK8_minDMass.Eta(),genJetAK8_minDMass.Phi(), Vec_Higgs[0].Eta(), Vec_Higgs[0].Phi());
  genJetAK8_minDMass_deltaR_H2_ = deltaR(genJetAK8_minDMass.Eta(),genJetAK8_minDMass.Phi(), Vec_Higgs[1].Eta(), Vec_Higgs[1].Phi());
  
  
  tree->Fill();
  
  Vec_Wboson.clear();
  Vec_Higgs.clear();
  Vec_wmJET.clear();
  Vec_wpJET.clear();
  Vec_Photons.clear();
  
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
