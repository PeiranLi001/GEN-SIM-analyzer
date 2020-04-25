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
  
  /************************************************************************
   *
   *  Good AK4 jet selection
   *  Selection Procedure: First two AK4 jets selected having min delta
   *                        mass wrt the onshell higgs boson (i.e. 80 GeV)
   *                        Then select another pair of AK4 jets that when
   *                        combined with selected onshell AK4 jets and gives
   *                        the minimum delta mass wrt the Higgs boson.
   *
   *************************************************************************/  
  std::cout << "===================================================" << std::endl;
  int nAK4jets = 0;
  TLorentzVector temp_genJetAK4;
  std::vector<TLorentzVector> local_Vec_genJetAK4;
  for(vector<reco::GenJet>::const_iterator genjet = genAK4_coll->begin(); genjet != genAK4_coll->end(); genjet++) {
    if (deltaR(genjet->eta(),genjet->phi(), Vec_Photons[0].Eta(),Vec_Photons[0].Phi())>0.4 && deltaR(genjet->eta(),genjet->phi(), Vec_Photons[1].Eta(),Vec_Photons[1].Phi())>0.4) {
      if ( jetCleaning(&(*genjet), genAK8_coll ) ) {
        nAK4jets++;
        temp_genJetAK4.SetPtEtaPhiE(genjet->pt(), genjet->eta(), genjet->phi(), genjet->energy());
        local_Vec_genJetAK4.push_back(temp_genJetAK4);
        std::cout << "==> " << genjet->pt() << "\t" << genjet->eta() << "\t" << genjet->phi() << "\t" << genjet->energy()<< std::endl;
      }
    }
  }

  // std::vector<TLorentzVector> local_Vec_genJetAK41;
  // SortedCleanedJetVector( genAK4_coll,  genAK8_coll, Vec_Photons, local_Vec_genJetAK41);

  /* Sort AK4 TLorentzVector */
  std::sort(local_Vec_genJetAK4.begin(), local_Vec_genJetAK4.end(), GenAnalyzer::reorder);
  genJetAK4_njets_ = nAK4jets;

  const std::vector<TLorentzVector> &Vec_genJetAK4 = local_Vec_genJetAK4;
  local_Vec_genJetAK4.clear();

  if (Vec_genJetAK4.size()>4) {
    // std::cout << "total number of AK4 jets is " << Vec_genJetAK4.size() << std::endl;

  /**
   * Loop over all AK4 jets and selects a pair of jets having 
   * minimum delta mass wrt w-boson mass.
   */
  double tempMass1 = 9999.0;
  int onshell_WBoson_index1 = -1;
  int onshell_WBoson_index2 = -1;
  for ( int ak4_jet1 = 0; ak4_jet1 < int(Vec_genJetAK4.size())-1; ++ak4_jet1)
  {
    for ( int ak4_jet2 = ak4_jet1+1; ak4_jet2 < int(Vec_genJetAK4.size()); ++ak4_jet2)
    {
      double mass = (Vec_genJetAK4[ak4_jet1] + Vec_genJetAK4[ak4_jet2]).M();      
      if (abs(mass - 80.0) < tempMass1)
      {
        tempMass1 = abs(mass - 80.0);
        onshell_WBoson_index1 = ak4_jet1;
        onshell_WBoson_index2 = ak4_jet2;
      }
    }
  }
  
  /**
   *  Select another pair of AK4 jets that when combined with
   *  the selected onshell AK4 jets gives min delta M wrt the 
   *  Higgs boson.
   */
  double tempMass2 = 99999.0;
  int offshell_WBoson_index1 = -1;
  int offshell_WBoson_index2 = -1;
  for ( int ak4_jet1 = 0; ak4_jet1 < int(Vec_genJetAK4.size())-1; ++ak4_jet1)
  {
    if (ak4_jet1 == onshell_WBoson_index1) continue;
    if (ak4_jet1 == onshell_WBoson_index2) continue;
    // std::cout << "ak4_jet1 : " << ak4_jet1 << std::endl;
    for ( int ak4_jet2 = ak4_jet1+1; ak4_jet2 < int(Vec_genJetAK4.size()); ++ak4_jet2)
    {
      if (ak4_jet2 == onshell_WBoson_index1) continue;
      if (ak4_jet2 == onshell_WBoson_index2) continue;
    // std::cout << "ak4_jet2 : " << ak4_jet2 << std::endl;
      double mass = (Vec_genJetAK4[ak4_jet1] + Vec_genJetAK4[ak4_jet2] + 
                    Vec_genJetAK4[onshell_WBoson_index1] + Vec_genJetAK4[onshell_WBoson_index2]).M();      
      // std::cout << "ak4_jet1: " << ak4_jet1 << ", ak4_jet2: " << ak4_jet2 << ", mass = " << mass << std::endl;
      if (abs(mass - 125.0) < tempMass2)
      {
        tempMass2 = abs(mass - 125.0);
        offshell_WBoson_index1 = ak4_jet1;
        offshell_WBoson_index2 = ak4_jet2;
      }
    }
  }
  // std::cout << "onshell: " << onshell_WBoson_index1 << "\t" << onshell_WBoson_index2 << std::endl;
  // std::cout << "offshell: " << offshell_WBoson_index1 << "\t" << offshell_WBoson_index2 << std::endl;
  AK4GEN_AllResolved_onShellJet1_Pt_ = Vec_genJetAK4[onshell_WBoson_index1].Pt();
  AK4GEN_AllResolved_onShellJet1_Eta_ = Vec_genJetAK4[onshell_WBoson_index1].Eta();
  AK4GEN_AllResolved_onShellJet1_Phi_ = Vec_genJetAK4[onshell_WBoson_index1].Phi();
  AK4GEN_AllResolved_onShellJet1_M_ = Vec_genJetAK4[onshell_WBoson_index1].M();
  AK4GEN_AllResolved_onShellJet1_dR_q1_ = deltaR(Vec_genJetAK4[onshell_WBoson_index1].Eta(),Vec_genJetAK4[onshell_WBoson_index1].Phi(), Vec_wmJET[0].Eta(), Vec_wmJET[0].Phi());
  AK4GEN_AllResolved_onShellJet1_dR_q2_ = deltaR(Vec_genJetAK4[onshell_WBoson_index1].Eta(),Vec_genJetAK4[onshell_WBoson_index1].Phi(), Vec_wmJET[1].Eta(), Vec_wmJET[1].Phi());
  AK4GEN_AllResolved_onShellJet1_dR_q3_ = deltaR(Vec_genJetAK4[onshell_WBoson_index1].Eta(),Vec_genJetAK4[onshell_WBoson_index1].Phi(), Vec_wpJET[0].Eta(), Vec_wpJET[0].Phi());
  AK4GEN_AllResolved_onShellJet1_dR_q4_ = deltaR(Vec_genJetAK4[onshell_WBoson_index1].Eta(),Vec_genJetAK4[onshell_WBoson_index1].Phi(), Vec_wpJET[1].Eta(), Vec_wpJET[1].Phi());
  AK4GEN_AllResolved_onShellJet1_dR_g1_ = deltaR(Vec_genJetAK4[onshell_WBoson_index1].Eta(),Vec_genJetAK4[onshell_WBoson_index1].Phi(), Vec_Photons[0].Eta(), Vec_Photons[0].Phi());
  AK4GEN_AllResolved_onShellJet1_dR_g2_ = deltaR(Vec_genJetAK4[onshell_WBoson_index1].Eta(),Vec_genJetAK4[onshell_WBoson_index1].Phi(), Vec_Photons[1].Eta(), Vec_Photons[1].Phi());

  AK4GEN_AllResolved_onShellJet2_Pt_ = Vec_genJetAK4[onshell_WBoson_index2].Pt();
  AK4GEN_AllResolved_onShellJet2_Eta_ = Vec_genJetAK4[onshell_WBoson_index2].Eta();
  AK4GEN_AllResolved_onShellJet2_Phi_ = Vec_genJetAK4[onshell_WBoson_index2].Phi();
  AK4GEN_AllResolved_onShellJet2_M_ = Vec_genJetAK4[onshell_WBoson_index2].M();
  AK4GEN_AllResolved_onShellJets_dR_ = deltaR(Vec_genJetAK4[onshell_WBoson_index1].Eta(), Vec_genJetAK4[onshell_WBoson_index1].Phi(), Vec_genJetAK4[onshell_WBoson_index2].Eta(), Vec_genJetAK4[onshell_WBoson_index2].Phi());
  AK4GEN_AllResolved_onShellJet2_dR_q1_ = deltaR(Vec_genJetAK4[onshell_WBoson_index2].Eta(),Vec_genJetAK4[onshell_WBoson_index2].Phi(), Vec_wmJET[0].Eta(), Vec_wmJET[0].Phi());
  AK4GEN_AllResolved_onShellJet2_dR_q2_ = deltaR(Vec_genJetAK4[onshell_WBoson_index2].Eta(),Vec_genJetAK4[onshell_WBoson_index2].Phi(), Vec_wmJET[1].Eta(), Vec_wmJET[1].Phi());
  AK4GEN_AllResolved_onShellJet2_dR_q3_ = deltaR(Vec_genJetAK4[onshell_WBoson_index2].Eta(),Vec_genJetAK4[onshell_WBoson_index2].Phi(), Vec_wpJET[0].Eta(), Vec_wpJET[0].Phi());
  AK4GEN_AllResolved_onShellJet2_dR_q4_ = deltaR(Vec_genJetAK4[onshell_WBoson_index2].Eta(),Vec_genJetAK4[onshell_WBoson_index2].Phi(), Vec_wpJET[1].Eta(), Vec_wpJET[1].Phi());
  AK4GEN_AllResolved_onShellJet2_dR_g1_ = deltaR(Vec_genJetAK4[onshell_WBoson_index2].Eta(),Vec_genJetAK4[onshell_WBoson_index2].Phi(), Vec_Photons[0].Eta(), Vec_Photons[0].Phi());
  AK4GEN_AllResolved_onShellJet2_dR_g2_ = deltaR(Vec_genJetAK4[onshell_WBoson_index2].Eta(),Vec_genJetAK4[onshell_WBoson_index2].Phi(), Vec_Photons[1].Eta(), Vec_Photons[1].Phi());

  AK4GEN_AllResolved_onShellWboson_Pt_ = (Vec_genJetAK4[onshell_WBoson_index1]+Vec_genJetAK4[onshell_WBoson_index2]).Pt();
  AK4GEN_AllResolved_onShellWboson_Eta_ = (Vec_genJetAK4[onshell_WBoson_index1]+Vec_genJetAK4[onshell_WBoson_index2]).Eta();
  AK4GEN_AllResolved_onShellWboson_Phi_ = (Vec_genJetAK4[onshell_WBoson_index1]+Vec_genJetAK4[onshell_WBoson_index2]).Phi();
  AK4GEN_AllResolved_onShellWboson_M_ = (Vec_genJetAK4[onshell_WBoson_index1]+Vec_genJetAK4[onshell_WBoson_index2]).M();
  double tempDeltaR = deltaR ((Vec_genJetAK4[onshell_WBoson_index1]+Vec_genJetAK4[onshell_WBoson_index2]).Eta(),
                              (Vec_genJetAK4[onshell_WBoson_index1]+Vec_genJetAK4[onshell_WBoson_index2]).Phi(),
                              Vec_Wboson[0].Eta(), Vec_Wboson[0].Phi());
  AK4GEN_AllResolved_onShellWboson_dR_W0PID_ = tempDeltaR;
  tempDeltaR = deltaR ((Vec_genJetAK4[onshell_WBoson_index1]+Vec_genJetAK4[onshell_WBoson_index2]).Eta(),
                              (Vec_genJetAK4[onshell_WBoson_index1]+Vec_genJetAK4[onshell_WBoson_index2]).Phi(),
                              Vec_Wboson[1].Eta(), Vec_Wboson[1].Phi());  
  AK4GEN_AllResolved_onShellWboson_dR_W1PID_ = tempDeltaR;

  AK4GEN_AllResolved_offShellJet1_Pt_ = Vec_genJetAK4[offshell_WBoson_index1].Pt();
  AK4GEN_AllResolved_offShellJet1_Eta_ = Vec_genJetAK4[offshell_WBoson_index1].Eta();
  AK4GEN_AllResolved_offShellJet1_Phi_ = Vec_genJetAK4[offshell_WBoson_index1].Phi();
  AK4GEN_AllResolved_offShellJet1_M_ = Vec_genJetAK4[offshell_WBoson_index1].M();
  AK4GEN_AllResolved_offShellJet1_dR_q1_ = deltaR(Vec_genJetAK4[offshell_WBoson_index1].Eta(),Vec_genJetAK4[offshell_WBoson_index1].Phi(), Vec_wmJET[0].Eta(), Vec_wmJET[0].Phi());
  AK4GEN_AllResolved_offShellJet1_dR_q2_ = deltaR(Vec_genJetAK4[offshell_WBoson_index1].Eta(),Vec_genJetAK4[offshell_WBoson_index1].Phi(), Vec_wmJET[1].Eta(), Vec_wmJET[1].Phi());
  AK4GEN_AllResolved_offShellJet1_dR_q3_ = deltaR(Vec_genJetAK4[offshell_WBoson_index1].Eta(),Vec_genJetAK4[offshell_WBoson_index1].Phi(), Vec_wpJET[0].Eta(), Vec_wpJET[0].Phi());
  AK4GEN_AllResolved_offShellJet1_dR_q4_ = deltaR(Vec_genJetAK4[offshell_WBoson_index1].Eta(),Vec_genJetAK4[offshell_WBoson_index1].Phi(), Vec_wpJET[1].Eta(), Vec_wpJET[1].Phi());
  AK4GEN_AllResolved_offShellJet1_dR_g1_ = deltaR(Vec_genJetAK4[offshell_WBoson_index1].Eta(),Vec_genJetAK4[offshell_WBoson_index1].Phi(), Vec_Photons[0].Eta(), Vec_Photons[0].Phi());
  AK4GEN_AllResolved_offShellJet1_dR_g2_ = deltaR(Vec_genJetAK4[offshell_WBoson_index1].Eta(),Vec_genJetAK4[offshell_WBoson_index1].Phi(), Vec_Photons[1].Eta(), Vec_Photons[1].Phi());

  AK4GEN_AllResolved_offShellJet2_Pt_ = Vec_genJetAK4[offshell_WBoson_index2].Pt();
  AK4GEN_AllResolved_offShellJet2_Eta_ = Vec_genJetAK4[offshell_WBoson_index2].Eta();
  AK4GEN_AllResolved_offShellJet2_Phi_ = Vec_genJetAK4[offshell_WBoson_index2].Phi();
  AK4GEN_AllResolved_offShellJet2_M_ = Vec_genJetAK4[offshell_WBoson_index2].M();
  AK4GEN_AllResolved_offShellJets_dR_ = deltaR(Vec_genJetAK4[offshell_WBoson_index1].Eta(), Vec_genJetAK4[offshell_WBoson_index1].Phi(), Vec_genJetAK4[offshell_WBoson_index2].Eta(), Vec_genJetAK4[offshell_WBoson_index2].Phi());
  double deltaRonoff = ((Vec_genJetAK4[onshell_WBoson_index1]+Vec_genJetAK4[onshell_WBoson_index2]).Pt(),
                        (Vec_genJetAK4[onshell_WBoson_index1]+Vec_genJetAK4[onshell_WBoson_index2]).Eta(),
                        (Vec_genJetAK4[offshell_WBoson_index1]+Vec_genJetAK4[offshell_WBoson_index2]).Pt(),
                        (Vec_genJetAK4[offshell_WBoson_index1]+Vec_genJetAK4[offshell_WBoson_index2]).Eta());
  AK4GEN_AllResolved_onShelloffShellJets_dR_ = deltaRonoff;
  AK4GEN_AllResolved_offShellJet2_dR_q1_ = deltaR(Vec_genJetAK4[offshell_WBoson_index2].Eta(),Vec_genJetAK4[offshell_WBoson_index2].Phi(), Vec_wmJET[0].Eta(), Vec_wmJET[0].Phi());
  AK4GEN_AllResolved_offShellJet2_dR_q2_ = deltaR(Vec_genJetAK4[offshell_WBoson_index2].Eta(),Vec_genJetAK4[offshell_WBoson_index2].Phi(), Vec_wmJET[1].Eta(), Vec_wmJET[1].Phi());
  AK4GEN_AllResolved_offShellJet2_dR_q3_ = deltaR(Vec_genJetAK4[offshell_WBoson_index2].Eta(),Vec_genJetAK4[offshell_WBoson_index2].Phi(), Vec_wpJET[0].Eta(), Vec_wpJET[0].Phi());
  AK4GEN_AllResolved_offShellJet2_dR_q4_ = deltaR(Vec_genJetAK4[offshell_WBoson_index2].Eta(),Vec_genJetAK4[offshell_WBoson_index2].Phi(), Vec_wpJET[1].Eta(), Vec_wpJET[1].Phi());
  AK4GEN_AllResolved_offShellJet2_dR_g1_ = deltaR(Vec_genJetAK4[offshell_WBoson_index2].Eta(),Vec_genJetAK4[offshell_WBoson_index2].Phi(), Vec_Photons[0].Eta(), Vec_Photons[0].Phi());
  AK4GEN_AllResolved_offShellJet2_dR_g2_ = deltaR(Vec_genJetAK4[offshell_WBoson_index2].Eta(),Vec_genJetAK4[offshell_WBoson_index2].Phi(), Vec_Photons[1].Eta(), Vec_Photons[1].Phi());
  

  AK4GEN_AllResolved_offShellWboson_Pt_ = (Vec_genJetAK4[offshell_WBoson_index1]+Vec_genJetAK4[offshell_WBoson_index2]).Pt();
  AK4GEN_AllResolved_offShellWboson_Eta_ = (Vec_genJetAK4[offshell_WBoson_index1]+Vec_genJetAK4[offshell_WBoson_index2]).Eta();
  AK4GEN_AllResolved_offShellWboson_Phi_ = (Vec_genJetAK4[offshell_WBoson_index1]+Vec_genJetAK4[offshell_WBoson_index2]).Phi();
  AK4GEN_AllResolved_offShellWboson_M_ = (Vec_genJetAK4[offshell_WBoson_index1]+Vec_genJetAK4[offshell_WBoson_index2]).M();
  tempDeltaR = deltaR ((Vec_genJetAK4[offshell_WBoson_index1]+Vec_genJetAK4[offshell_WBoson_index2]).Eta(),
                              (Vec_genJetAK4[offshell_WBoson_index1]+Vec_genJetAK4[offshell_WBoson_index2]).Phi(),
                              Vec_Wboson[0].Eta(), Vec_Wboson[0].Phi());
  AK4GEN_AllResolved_offShellWboson_dR_W0PID_ = tempDeltaR;
  tempDeltaR = deltaR ((Vec_genJetAK4[offshell_WBoson_index1]+Vec_genJetAK4[offshell_WBoson_index2]).Eta(),
                              (Vec_genJetAK4[offshell_WBoson_index1]+Vec_genJetAK4[offshell_WBoson_index2]).Phi(),
                              Vec_Wboson[1].Eta(), Vec_Wboson[1].Phi());  
  AK4GEN_AllResolved_offShellWboson_dR_W1PID_ = tempDeltaR;

  AK4GEN_AllResolved_Higgs_Pt_ = (Vec_genJetAK4[offshell_WBoson_index1]+Vec_genJetAK4[offshell_WBoson_index2]+
                                  Vec_genJetAK4[onshell_WBoson_index1]+Vec_genJetAK4[onshell_WBoson_index2]).Pt();
  AK4GEN_AllResolved_Higgs_Eta_ = (Vec_genJetAK4[offshell_WBoson_index1]+Vec_genJetAK4[offshell_WBoson_index2]+
                                  Vec_genJetAK4[onshell_WBoson_index1]+Vec_genJetAK4[onshell_WBoson_index2]).Eta();
  AK4GEN_AllResolved_Higgs_Phi_ = (Vec_genJetAK4[offshell_WBoson_index1]+Vec_genJetAK4[offshell_WBoson_index2]+
                                  Vec_genJetAK4[onshell_WBoson_index1]+Vec_genJetAK4[onshell_WBoson_index2]).Phi();
  AK4GEN_AllResolved_Higgs_M_ = (Vec_genJetAK4[offshell_WBoson_index1]+Vec_genJetAK4[offshell_WBoson_index2]+
                                  Vec_genJetAK4[onshell_WBoson_index1]+Vec_genJetAK4[onshell_WBoson_index2]).M();
  AK4GEN_AllResolved_Higgs_DR_Higgs0PID_ =((Vec_genJetAK4[offshell_WBoson_index1]+Vec_genJetAK4[offshell_WBoson_index2]+
                                           Vec_genJetAK4[onshell_WBoson_index1]+Vec_genJetAK4[onshell_WBoson_index2]).Eta(),
                                          (Vec_genJetAK4[offshell_WBoson_index1]+Vec_genJetAK4[offshell_WBoson_index2]+
                                           Vec_genJetAK4[onshell_WBoson_index1]+Vec_genJetAK4[onshell_WBoson_index2]).Phi(),
                                          Vec_Higgs[0].Eta(), Vec_Higgs[0].Phi());
  AK4GEN_AllResolved_Higgs_DR_Higgs1PID_ =((Vec_genJetAK4[offshell_WBoson_index1]+Vec_genJetAK4[offshell_WBoson_index2]+
                                           Vec_genJetAK4[onshell_WBoson_index1]+Vec_genJetAK4[onshell_WBoson_index2]).Eta(),
                                          (Vec_genJetAK4[offshell_WBoson_index1]+Vec_genJetAK4[offshell_WBoson_index2]+
                                           Vec_genJetAK4[onshell_WBoson_index1]+Vec_genJetAK4[onshell_WBoson_index2]).Phi(),
                                          Vec_Higgs[1].Eta(), Vec_Higgs[1].Phi());
  }

  /**********************************************************
   * 
   *              ak8 jet selection
   * 
   **********************************************************/
  int nAK8jets = 0;
  double temp_AK8jet_pt = -999.0;
  double temp_AK8jet_deltaM = 9999.0;
  TLorentzVector AK8Gen_HiggsJet_MaxPt;
  TLorentzVector AK8Gen_HiggsJet_minDMass;
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
          AK8Gen_HiggsJet_MaxPt.SetPtEtaPhiE(genjet->pt(), genjet->eta(), genjet->phi(), genjet->energy());
          temp_AK8jet_pt = genjet->pt();
        }
        if ( abs(genjet->mass() - 125.0) < temp_AK8jet_deltaM)
        {
          AK8Gen_HiggsJet_minDMass.SetPtEtaPhiE(genjet->pt(), genjet->eta(), genjet->phi(), genjet->energy());
          temp_AK8jet_deltaM = abs(genjet->mass() - 125.0);
        }
      }
    }
  }
  
  //std::cout << "Number of AK8 candidates = " << nAK8jets << std::endl;
  
  AK8Gen_HiggsJet_njets_ = nAK8jets;
  AK8Gen_HiggsJet_MaxPt_Pt_ = AK8Gen_HiggsJet_MaxPt.Pt();
  AK8Gen_HiggsJet_MaxPt_Eta_ = AK8Gen_HiggsJet_MaxPt.Eta();
  AK8Gen_HiggsJet_MaxPt_Phi_ = AK8Gen_HiggsJet_MaxPt.Phi();
  AK8Gen_HiggsJet_MaxPt_M_ = AK8Gen_HiggsJet_MaxPt.M();
  AK8Gen_HiggsJet_MaxPt_deltaR_H1_ = deltaR(AK8Gen_HiggsJet_MaxPt.Eta(),AK8Gen_HiggsJet_MaxPt.Phi(), Vec_Higgs[0].Eta(), Vec_Higgs[0].Phi());
  AK8Gen_HiggsJet_MaxPt_deltaR_H2_ = deltaR(AK8Gen_HiggsJet_MaxPt.Eta(),AK8Gen_HiggsJet_MaxPt.Phi(), Vec_Higgs[1].Eta(), Vec_Higgs[1].Phi());
  
  AK8Gen_HiggsJet_minDMass_Pt_ = AK8Gen_HiggsJet_minDMass.Pt();
  AK8Gen_HiggsJet_minDMass_Eta_ = AK8Gen_HiggsJet_minDMass.Eta();
  AK8Gen_HiggsJet_minDMass_Phi_ = AK8Gen_HiggsJet_minDMass.Phi();
  AK8Gen_HiggsJet_minDMass_M_ = AK8Gen_HiggsJet_minDMass.M();
  AK8Gen_HiggsJet_minDMass_deltaR_H1_ = deltaR(AK8Gen_HiggsJet_minDMass.Eta(),AK8Gen_HiggsJet_minDMass.Phi(), Vec_Higgs[0].Eta(), Vec_Higgs[0].Phi());
  AK8Gen_HiggsJet_minDMass_deltaR_H2_ = deltaR(AK8Gen_HiggsJet_minDMass.Eta(),AK8Gen_HiggsJet_minDMass.Phi(), Vec_Higgs[1].Eta(), Vec_Higgs[1].Phi());
  
  
  /************************************************************************
   *
   *  Two AK8 jet selection
   *  The case where Higgs will be selected as two FAT jets
   *  CASE: 1: Two leading AK8 jets are selected as W-boson jets
   *
   *************************************************************************/
  TLorentzVector genJetResolvedAK8;
  std::vector<TLorentzVector> Vec_AK8Gen_MergedWjets;
  for(vector<reco::GenJet>::const_iterator genjet = genAK8_coll->begin(); genjet != genAK8_coll->end(); genjet++) {
    if (deltaR(genjet->eta(),genjet->phi(), Vec_Photons[0].Eta(),Vec_Photons[0].Phi())>0.4 && deltaR(genjet->eta(),genjet->phi(), Vec_Photons[1].Eta(),Vec_Photons[1].Phi())>0.4) {
      if ( jetCleaning(&(*genjet), genAK4_coll ) )
      {
        genJetResolvedAK8.SetPtEtaPhiE(genjet->pt(), genjet->eta(), genjet->phi(), genjet->energy());
        Vec_AK8Gen_MergedWjets.push_back(genJetResolvedAK8);
      }
    }
  }
  std::sort(Vec_AK8Gen_MergedWjets.begin(), Vec_AK8Gen_MergedWjets.end(), GenAnalyzer::reorder);
  
  if (Vec_AK8Gen_MergedWjets.size()>=2)
  {
    AK8Gen_MergedWjets_MaxPt_Leading_Pt_ = Vec_AK8Gen_MergedWjets[0].Pt();
    AK8Gen_MergedWjets_MaxPt_Leading_Eta_ = Vec_AK8Gen_MergedWjets[0].Eta();
    AK8Gen_MergedWjets_MaxPt_Leading_Phi_ = Vec_AK8Gen_MergedWjets[0].Phi();
    AK8Gen_MergedWjets_MaxPt_Leading_M_ = Vec_AK8Gen_MergedWjets[0].M();
    AK8Gen_MergedWjets_MaxPt_Leading_deltaR_H1_ = deltaR(Vec_AK8Gen_MergedWjets[0].Eta(),Vec_AK8Gen_MergedWjets[0].Phi(), Vec_Higgs[0].Eta(), Vec_Higgs[0].Phi());
    AK8Gen_MergedWjets_MaxPt_Leading_deltaR_H2_ = deltaR(Vec_AK8Gen_MergedWjets[0].Eta(),Vec_AK8Gen_MergedWjets[0].Phi(), Vec_Higgs[1].Eta(), Vec_Higgs[1].Phi());
    
    AK8Gen_MergedWjets_MaxPt_SubLeading_Pt_ = Vec_AK8Gen_MergedWjets[1].Pt();
    AK8Gen_MergedWjets_MaxPt_SubLeading_Eta_ = Vec_AK8Gen_MergedWjets[1].Eta();
    AK8Gen_MergedWjets_MaxPt_SubLeading_Phi_ = Vec_AK8Gen_MergedWjets[1].Phi();
    AK8Gen_MergedWjets_MaxPt_SubLeading_M_ = Vec_AK8Gen_MergedWjets[1].M();
    AK8Gen_MergedWjets_MaxPt_SubLeading_deltaR_H1_ = deltaR(Vec_AK8Gen_MergedWjets[1].Eta(),Vec_AK8Gen_MergedWjets[1].Phi(), Vec_Higgs[0].Eta(), Vec_Higgs[0].Phi());
    AK8Gen_MergedWjets_MaxPt_SubLeading_deltaR_H2_ = deltaR(Vec_AK8Gen_MergedWjets[1].Eta(),Vec_AK8Gen_MergedWjets[1].Phi(), Vec_Higgs[1].Eta(), Vec_Higgs[1].Phi());
    AK8Gen_MergedWjets_MaxPt_LeadingSubLeading_DR_ = deltaR(Vec_AK8Gen_MergedWjets[1].Eta(),Vec_AK8Gen_MergedWjets[1].Phi(),Vec_AK8Gen_MergedWjets[0].Eta(),Vec_AK8Gen_MergedWjets[0].Phi());

    AK8Gen_MergedWjets_MaxPt_Higgs_Pt_ = (Vec_AK8Gen_MergedWjets[0]+Vec_AK8Gen_MergedWjets[1]).Pt();
    AK8Gen_MergedWjets_MaxPt_Higgs_Eta_ = (Vec_AK8Gen_MergedWjets[0]+Vec_AK8Gen_MergedWjets[1]).Eta();
    AK8Gen_MergedWjets_MaxPt_Higgs_Phi_ = (Vec_AK8Gen_MergedWjets[0]+Vec_AK8Gen_MergedWjets[1]).Phi();
    AK8Gen_MergedWjets_MaxPt_Higgs_M_ = (Vec_AK8Gen_MergedWjets[0]+Vec_AK8Gen_MergedWjets[1]).M();
  
    AK8Gen_MergedWjets_MaxPt_Higgs_deltaR_H1_ = deltaR((Vec_AK8Gen_MergedWjets[0]+Vec_AK8Gen_MergedWjets[1]).Eta(),(Vec_AK8Gen_MergedWjets[0]+Vec_AK8Gen_MergedWjets[1]).Phi(), Vec_Higgs[0].Eta(), Vec_Higgs[0].Phi());
    AK8Gen_MergedWjets_MaxPt_Higgs_deltaR_H2_ = deltaR((Vec_AK8Gen_MergedWjets[0]+Vec_AK8Gen_MergedWjets[1]).Eta(),(Vec_AK8Gen_MergedWjets[0]+Vec_AK8Gen_MergedWjets[1]).Phi(), Vec_Higgs[1].Eta(), Vec_Higgs[1].Phi());

  }
  
  /************************************************************************
   *
   *  Two AK8 jet selection
   *  The case where Higgs will be selected as two FAT jets
   *  CASE: 2: Two jets with min delta mass wrt Higgs are selected.
   *
   *************************************************************************/  
  double temp_ResolvedJetAK8_deltaMass = 99999.0;
  TLorentzVector temp_ReolvedJetAK8_LV1, temp_ReolvedJetAK8_LV2;
  TLorentzVector LV_leadingWjet, LV_SubleadingWjet;
  int foundjets = 0;
  for(vector<reco::GenJet>::const_iterator genjet1 = genAK8_coll->begin(); genjet1 != genAK8_coll->end()-1; genjet1++) {
      for(vector<reco::GenJet>::const_iterator genjet2 = genjet1+1; genjet2 != genAK8_coll->end(); genjet2++) {
        temp_ReolvedJetAK8_LV1.SetPtEtaPhiE(genjet1->pt(), genjet1->eta(), genjet1->phi(), genjet1->energy());
        temp_ReolvedJetAK8_LV2.SetPtEtaPhiE(genjet2->pt(), genjet2->eta(), genjet2->phi(), genjet2->energy());

        if ( (abs((temp_ReolvedJetAK8_LV1+temp_ReolvedJetAK8_LV2).M() - 125.0) < temp_ResolvedJetAK8_deltaMass) ) {
          temp_ResolvedJetAK8_deltaMass = abs((temp_ReolvedJetAK8_LV1+temp_ReolvedJetAK8_LV2).M() - 125.0);
          foundjets = 1;
          if (genjet1->pt()>genjet2->pt())
            LV_leadingWjet = temp_ReolvedJetAK8_LV1;
          else
            LV_SubleadingWjet = temp_ReolvedJetAK8_LV2;
        }
      }
  }

  if (foundjets) {
  AK8Gen_MergedWjets_minDMass_Leading_Pt_ = LV_leadingWjet.Pt();
  AK8Gen_MergedWjets_minDMass_Leading_Eta_ = LV_leadingWjet.Eta();
  AK8Gen_MergedWjets_minDMass_Leading_Phi_ = LV_leadingWjet.Phi();
  AK8Gen_MergedWjets_minDMass_Leading_M_ = LV_leadingWjet.M();
  AK8Gen_MergedWjets_minDMass_Leading_deltaR_H1_ = deltaR(LV_leadingWjet.Eta(),LV_leadingWjet.Phi(), Vec_Higgs[0].Eta(), Vec_Higgs[0].Phi());
  AK8Gen_MergedWjets_minDMass_Leading_deltaR_H2_ = deltaR(LV_leadingWjet.Eta(),LV_leadingWjet.Phi(), Vec_Higgs[1].Eta(), Vec_Higgs[1].Phi());
  
  AK8Gen_MergedWjets_minDMass_SubLeading_Pt_ = LV_SubleadingWjet.Pt();
  AK8Gen_MergedWjets_minDMass_SubLeading_Eta_ = LV_SubleadingWjet.Eta();
  AK8Gen_MergedWjets_minDMass_SubLeading_Phi_ = LV_SubleadingWjet.Phi();
  AK8Gen_MergedWjets_minDMass_SubLeading_M_ = LV_SubleadingWjet.M();
  AK8Gen_MergedWjets_minDMass_SubLeading_deltaR_H1_ = deltaR(LV_SubleadingWjet.Eta(),LV_SubleadingWjet.Phi(), Vec_Higgs[0].Eta(), Vec_Higgs[0].Phi());
  AK8Gen_MergedWjets_minDMass_SubLeading_deltaR_H2_ = deltaR(LV_SubleadingWjet.Eta(),LV_SubleadingWjet.Phi(), Vec_Higgs[1].Eta(), Vec_Higgs[1].Phi());
  AK8Gen_MergedWjets_minDMass_LeadingSubLeading_DR_ = deltaR(LV_SubleadingWjet.Eta(),LV_SubleadingWjet.Phi(),LV_leadingWjet.Eta(),LV_leadingWjet.Phi());

  AK8Gen_MergedWjets_minDMass_Higgs_Pt_ = (LV_leadingWjet+LV_SubleadingWjet).Pt();
  AK8Gen_MergedWjets_minDMass_Higgs_Eta_ = (LV_leadingWjet+LV_SubleadingWjet).Eta();
  AK8Gen_MergedWjets_minDMass_Higgs_Phi_ = (LV_leadingWjet+LV_SubleadingWjet).Phi();
  AK8Gen_MergedWjets_minDMass_Higgs_M_ = (LV_leadingWjet+LV_SubleadingWjet).M();
  
  AK8Gen_MergedWjets_minDMass_Higgs_deltaR_H1_ = deltaR((LV_leadingWjet+LV_SubleadingWjet).Eta(),(LV_leadingWjet+LV_SubleadingWjet).Phi(), Vec_Higgs[0].Eta(), Vec_Higgs[0].Phi());
  AK8Gen_MergedWjets_minDMass_Higgs_deltaR_H2_ = deltaR((LV_leadingWjet+LV_SubleadingWjet).Eta(),(LV_leadingWjet+LV_SubleadingWjet).Phi(), Vec_Higgs[1].Eta(), Vec_Higgs[1].Phi());
  }

  /************************************************************************
   *
   *  Two AK8 jet selection
   *  The case where Higgs will be selected as two FAT jets
   *  CASE: 3:  One of jets is selected having min delta M wrt W-boson mass
   *            out of remaning another one is selected that forms min Higgs
   *            mass with the one selected for w-boson.
   *
   *************************************************************************/  
  double temp_minWBosonMass = 9999.0;
  double temp_minHMass = 9999.0;
  TLorentzVector temp_ReolvedJetAK8_LV1_a, temp_ReolvedJetAK8_LV2_a;
  TLorentzVector LV_OnShell_WBoson;
  TLorentzVector LV_OffShell_WBoson;  
  int counter_of_selected_wjet=-1;
  int temp_counter_of_selected_wjet=0;
  int foundjets1 = 0;
  int foundjets2 = 0;
  /* Iterator over genAK8 jets for W-jet selection  */
  for(vector<reco::GenJet>::const_iterator genjet1 = genAK8_coll->begin(); genjet1 != genAK8_coll->end(); genjet1++) {
    temp_counter_of_selected_wjet++;
    temp_ReolvedJetAK8_LV1_a.SetPtEtaPhiE(genjet1->pt(), genjet1->eta(), genjet1->phi(), genjet1->energy());
    if ( abs(temp_ReolvedJetAK8_LV1_a.M() - 80.0) < temp_minWBosonMass) {
      foundjets1 = 1;
      LV_OnShell_WBoson = temp_ReolvedJetAK8_LV1_a;
      counter_of_selected_wjet = temp_counter_of_selected_wjet;
    }
  }
  /* Iterator for off-shell w-jet selection after on-shell w-jet selection  */
  temp_counter_of_selected_wjet = 0;
  for(vector<reco::GenJet>::const_iterator genjet1 = genAK8_coll->begin(); genjet1 != genAK8_coll->end(); genjet1++) {
    temp_counter_of_selected_wjet++;
    if (temp_counter_of_selected_wjet == counter_of_selected_wjet) continue;
    temp_ReolvedJetAK8_LV2_a.SetPtEtaPhiE(genjet1->pt(), genjet1->eta(), genjet1->phi(), genjet1->energy());
    if ( abs((temp_ReolvedJetAK8_LV2_a+LV_OnShell_WBoson).M() - 125.0) < temp_minHMass){
      foundjets2 = 1;
      LV_OffShell_WBoson = temp_ReolvedJetAK8_LV2_a;
    }
  }

  if (LV_OnShell_WBoson.Pt() > LV_OffShell_WBoson.Pt()) {
    LV_leadingWjet = LV_OnShell_WBoson;
    LV_SubleadingWjet = LV_OffShell_WBoson;
  } else {
    LV_leadingWjet = LV_OffShell_WBoson;
    LV_SubleadingWjet = LV_OnShell_WBoson;
  }

  if (foundjets1 == 1 && foundjets2 == 1)
  {
  AK8Gen_MergedWjets_minWminHmass_Leading_Pt_ = LV_leadingWjet.Pt();
  AK8Gen_MergedWjets_minWminHmass_Leading_Eta_ = LV_leadingWjet.Eta();
  AK8Gen_MergedWjets_minWminHmass_Leading_Phi_ = LV_leadingWjet.Phi();
  AK8Gen_MergedWjets_minWminHmass_Leading_M_ = LV_leadingWjet.M();
  AK8Gen_MergedWjets_minWminHmass_Leading_deltaR_H1_ = deltaR(LV_leadingWjet.Eta(),LV_leadingWjet.Phi(), Vec_Higgs[0].Eta(), Vec_Higgs[0].Phi());
  AK8Gen_MergedWjets_minWminHmass_Leading_deltaR_H2_ = deltaR(LV_leadingWjet.Eta(),LV_leadingWjet.Phi(), Vec_Higgs[1].Eta(), Vec_Higgs[1].Phi());
  
  AK8Gen_MergedWjets_minWminHmass_SubLeading_Pt_ = LV_SubleadingWjet.Pt();
  AK8Gen_MergedWjets_minWminHmass_SubLeading_Eta_ = LV_SubleadingWjet.Eta();
  AK8Gen_MergedWjets_minWminHmass_SubLeading_Phi_ = LV_SubleadingWjet.Phi();
  AK8Gen_MergedWjets_minWminHmass_SubLeading_M_ = LV_SubleadingWjet.M();
  AK8Gen_MergedWjets_minWminHmass_SubLeading_deltaR_H1_ = deltaR(LV_SubleadingWjet.Eta(),LV_SubleadingWjet.Phi(), Vec_Higgs[0].Eta(), Vec_Higgs[0].Phi());
  AK8Gen_MergedWjets_minWminHmass_SubLeading_deltaR_H2_ = deltaR(LV_SubleadingWjet.Eta(),LV_SubleadingWjet.Phi(), Vec_Higgs[1].Eta(), Vec_Higgs[1].Phi());
  AK8Gen_MergedWjets_minWminHmass_LeadingSubLeading_DR_ = deltaR(LV_SubleadingWjet.Eta(),LV_SubleadingWjet.Phi(),LV_leadingWjet.Eta(),LV_leadingWjet.Phi());

  AK8Gen_MergedWjets_minWminHmass_Higgs_Pt_ = (LV_leadingWjet+LV_SubleadingWjet).Pt();
  AK8Gen_MergedWjets_minWminHmass_Higgs_Eta_ = (LV_leadingWjet+LV_SubleadingWjet).Eta();
  AK8Gen_MergedWjets_minWminHmass_Higgs_Phi_ = (LV_leadingWjet+LV_SubleadingWjet).Phi();
  AK8Gen_MergedWjets_minWminHmass_Higgs_M_ = (LV_leadingWjet+LV_SubleadingWjet).M();
  
  AK8Gen_MergedWjets_minWminHmass_Higgs_deltaR_H1_ = deltaR((LV_leadingWjet+LV_SubleadingWjet).Eta(),(LV_leadingWjet+LV_SubleadingWjet).Phi(), Vec_Higgs[0].Eta(), Vec_Higgs[0].Phi());
  AK8Gen_MergedWjets_minWminHmass_Higgs_deltaR_H2_ = deltaR((LV_leadingWjet+LV_SubleadingWjet).Eta(),(LV_leadingWjet+LV_SubleadingWjet).Phi(), Vec_Higgs[1].Eta(), Vec_Higgs[1].Phi());
  }




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
