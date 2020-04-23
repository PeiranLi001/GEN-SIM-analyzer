/*
 * =====================================================================================
 *
 *       Filename:  GenAnalyzer.h
 *
 *    Description:
 *
 *        Version:  1.0
 *        Created:  12.10.2016 09:03:18
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Ramkrishna Sharma (Ram), ramkrishna.sharma71@gmail.com
 *   Organization:  University Of Delhi, Delhi, India
 *
 * =====================================================================================
 */

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include "TTree.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "Math/GenVector/VectorUtil.h"
#include "DataFormats/GeometryVector/interface/VectorUtil.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/METReco/interface/GenMETCollection.h"
#include "DataFormats/METReco/interface/GenMET.h"

#include "TLorentzVector.h"

using namespace edm;
using namespace std;
using namespace reco;
using namespace ROOT::Math::VectorUtil;

class GenAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {                                                                  
public:
  explicit GenAnalyzer(const edm::ParameterSet&);
  ~GenAnalyzer();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  static bool reorder(const TLorentzVector &a, const TLorentzVector &b);
  //https://stackoverflow.com/a/26230635/2302094
  static bool jetCleaning(const reco::GenJet  * genAK8jet,const vector<reco::GenJet>* genAK4_coll, const double r_seperation=0.8);

  void SetBranches();
  void Clear();
  
private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  
  // ----------member data ---------------------------
  virtual void AddBranch(std::vector<std::string>*, std::string name);
  virtual void AddBranch(std::vector<double>*, std::string name);
  virtual void AddBranch(std::vector<int>*, std::string name);
  virtual void AddBranch(int* vec, std::string name);
  virtual void AddBranch(double* vec, std::string name);
  virtual void computeAngles(TLorentzVector thep4H, TLorentzVector thep4Z1, TLorentzVector thep4M11, TLorentzVector thep4M12, TLorentzVector thep4Z2, TLorentzVector thep4M21, TLorentzVector thep4M22, double& costheta1, double& costheta2, double& Phi, double& costhetastar, double& Phi1);
  
  edm::EDGetTokenT<reco::GenParticleCollection> genParticlesToken;
  edm::EDGetTokenT<LHEEventProduct> LHEEventToken;
  edm::EDGetTokenT<reco::GenJetCollection> genAK4jetToken;
  edm::EDGetTokenT<reco::GenJetCollection> genAK8jetToken;
  edm::EDGetTokenT<reco::GenMETCollection> genMetCaloToken;
  
  bool	Verbose_;
  TString OutPutFileName_;
  
  edm::Service<TFileService> fs;
  TFile * outputFile_;
  TTree* tree;
  std::ofstream file1;
  
  //std::vector<int> pdgID_;
  
  //  std::vector<std::string> LHEWeightIDs_;
  //  std::vector<double> LHEWeights_;
  
  int nEVENT=-999;
  
  //  int		isMuMinus_ = -999;
  //  double 	LHELeptPt_ = -999.0;
  //  double 	LHELeptEta_ = -999.0;
  //  double 	LHELeptPhi_ = -999.0;
  //  double 	LHELeptM_ = -999.0;
  //  double 	LHELeptE_ = -999.0;
  
  double gen_leading_photon_Pt_   = -999.0;
  double gen_leading_photon_Eta_  = -999.0;
  double gen_leading_photon_Phi_  = -999.0;
  double gen_leading_photon_M_    = -999.0;
  double gen_Subleading_photon_Pt_    = -999.0;
  double gen_Subleading_photon_Eta_   = -999.0;
  double gen_Subleading_photon_Phi_   = -999.0;
  double gen_Subleading_photon_M_     = -999.0;
  double gen_leading_WpJets_Pt_   = -999.0;
  double gen_leading_WpJets_Eta_  = -999.0;
  double gen_leading_WpJets_Phi_  = -999.0;
  double gen_leading_WpJets_M_    = -999.0;
  double gen_Subleading_WpJets_Pt_    = -999.0;
  double gen_Subleading_WpJets_Eta_   = -999.0;
  double gen_Subleading_WpJets_Phi_   = -999.0;
  double gen_Subleading_WpJets_M_     = -999.0;
  double gen_leading_WmJets_Pt_   = -999.0;
  double gen_leading_WmJets_Eta_  = -999.0;
  double gen_leading_WmJets_Phi_  = -999.0;
  double gen_leading_WmJets_M_    = -999.0;
  double gen_Subleading_WmJets_Pt_    = -999.0;
  double gen_Subleading_WmJets_Eta_   = -999.0;
  double gen_Subleading_WmJets_Phi_   = -999.0;
  double gen_Subleading_WmJets_M_     = -999.0;
  double gen_leading_WBoson_Pt_   = -999.0;
  double gen_leading_WBoson_Eta_  = -999.0;
  double gen_leading_WBoson_Phi_  = -999.0;
  double gen_leading_WBoson_M_    = -999.0;
  double gen_Subleading_WBoson_Pt_    = -999.0;
  double gen_Subleading_WBoson_Eta_   = -999.0;
  double gen_Subleading_WBoson_Phi_   = -999.0;
  double gen_Subleading_WBoson_M_     = -999.0;
  double gen_leading_Higgs_Pt_    = -999.0;
  double gen_leading_Higgs_Eta_   = -999.0;
  double gen_leading_Higgs_Phi_   = -999.0;
  double gen_leading_Higgs_M_     = -999.0;
  double gen_Subleading_Higgs_Pt_     = -999.0;
  double gen_Subleading_Higgs_Eta_    = -999.0;
  double gen_Subleading_Higgs_Phi_    = -999.0;
  double gen_Subleading_Higgs_M_  = -999.0;
  double gen_deltaR_Photon0_Photon1_  = -999.0;
  double gen_deltaR_Photon0_WmJ0_     = -999.0;
  double gen_deltaR_Photon0_WmJ1_     = -999.0;
  double gen_deltaR_Photon0_WpJ0_     = -999.0;
  double gen_deltaR_Photon0_WpJ1_     = -999.0;
  double gen_deltaR_Photon1_WmJ0_     = -999.0;
  double gen_deltaR_Photon1_WmJ1_     = -999.0;
  double gen_deltaR_Photon1_WpJ0_     = -999.0;
  double gen_deltaR_Photon1_WpJ1_     = -999.0;
  double gen_deltaR_WpJ0_WpJ1_     = -999.0;
  double gen_deltaR_WpJ0_WmJ0_     = -999.0;
  double gen_deltaR_WpJ0_WmJ1_     = -999.0;
  double gen_deltaR_WpJ1_WmJ0_     = -999.0;
  double gen_deltaR_WpJ1_WmJ1_     = -999.0;
  double gen_deltaR_WmJ0_WmJ1_     = -999.0;
  double gen_deltaR_Wp_Wm_     = -999.0;
  double gen_deltaR_H1_H2_     = -999.0;
  
  double genJetAK4_njets_ = -999.0;
  double genJetAK4_leading_Pt_ = -999.0;
  double genJetAK4_leading_Eta_ = -999.0;
  double genJetAK4_leading_Phi_ = -999.0;
  double genJetAK4_leading_M_ = -999.0;
  double genJetAK4_leading_Energy_ = -999.0;
  double genJetAK4_Subleading_Pt_ = -999.0;
  double genJetAK4_Subleading_Eta_ = -999.0;
  double genJetAK4_Subleading_Phi_ = -999.0;
  double genJetAK4_Subleading_M_ = -999.0;
  double genJetAK4_Subleading_Energy_ = -999.0;
  
  double AK8Gen_HiggsJet_njets_ = -999.0;
  double AK8Gen_HiggsJet_MaxPt_Pt_ = -999.0;
  double AK8Gen_HiggsJet_MaxPt_Eta_ = -999.0;
  double AK8Gen_HiggsJet_MaxPt_Phi_ = -999.0;
  double AK8Gen_HiggsJet_MaxPt_M_   = -999.0;
  double AK8Gen_HiggsJet_MaxPt_deltaR_H1_ = -999.0;
  double AK8Gen_HiggsJet_MaxPt_deltaR_H2_ = -999.0;
  
  double AK8Gen_HiggsJet_minDMass_Pt_ = -999.0;
  double AK8Gen_HiggsJet_minDMass_Eta_ = -999.0;
  double AK8Gen_HiggsJet_minDMass_Phi_ = -999.0;
  double AK8Gen_HiggsJet_minDMass_M_   = -999.0;
  double AK8Gen_HiggsJet_minDMass_deltaR_H1_ = -999.0;
  double AK8Gen_HiggsJet_minDMass_deltaR_H2_ = -999.0;
  
  double AK8Gen_MergedWjets_MaxPt_Leading_Pt_ = -999.0;
  double AK8Gen_MergedWjets_MaxPt_Leading_Eta_ = -999.0;
  double AK8Gen_MergedWjets_MaxPt_Leading_Phi_ = -999.0;
  double AK8Gen_MergedWjets_MaxPt_Leading_M_ = -999.0;
  double AK8Gen_MergedWjets_MaxPt_Leading_deltaR_H1_ = -999.0;
  double AK8Gen_MergedWjets_MaxPt_Leading_deltaR_H2_ = -999.0;
  double AK8Gen_MergedWjets_MaxPt_SubLeading_Pt_ = -999.0;
  double AK8Gen_MergedWjets_MaxPt_SubLeading_Eta_ = -999.0;
  double AK8Gen_MergedWjets_MaxPt_SubLeading_Phi_ = -999.0;
  double AK8Gen_MergedWjets_MaxPt_SubLeading_M_ = -999.0;
  double AK8Gen_MergedWjets_MaxPt_SubLeading_deltaR_H1_ = -999.0;
  double AK8Gen_MergedWjets_MaxPt_SubLeading_deltaR_H2_ = -999.0;
  double AK8Gen_MergedWjets_MaxPt_LeadingSubLeading_DR_ = -999.0;

  double AK8Gen_MergedWjets_minDMass_Leading_Pt_ = -999.0;
  double AK8Gen_MergedWjets_minDMass_Leading_Eta_ = -999.0;
  double AK8Gen_MergedWjets_minDMass_Leading_Phi_ = -999.0;
  double AK8Gen_MergedWjets_minDMass_Leading_M_ = -999.0;
  double AK8Gen_MergedWjets_minDMass_Leading_deltaR_H1_ = -999.0;
  double AK8Gen_MergedWjets_minDMass_Leading_deltaR_H2_ = -999.0;
  double AK8Gen_MergedWjets_minDMass_SubLeading_Pt_ = -999.0;
  double AK8Gen_MergedWjets_minDMass_SubLeading_Eta_ = -999.0;
  double AK8Gen_MergedWjets_minDMass_SubLeading_Phi_ = -999.0;
  double AK8Gen_MergedWjets_minDMass_SubLeading_M_ = -999.0;
  double AK8Gen_MergedWjets_minDMass_SubLeading_deltaR_H1_ = -999.0;
  double AK8Gen_MergedWjets_minDMass_SubLeading_deltaR_H2_ = -999.0;
  double AK8Gen_MergedWjets_minDMass_LeadingSubLeading_DR_   = -999.0;

  double AK8Gen_MergedWjets_minWminHmass_Leading_Pt_ = -999.0;
  double AK8Gen_MergedWjets_minWminHmass_Leading_Eta_ = -999.0;
  double AK8Gen_MergedWjets_minWminHmass_Leading_Phi_ = -999.0;
  double AK8Gen_MergedWjets_minWminHmass_Leading_M_ = -999.0;
  double AK8Gen_MergedWjets_minWminHmass_Leading_deltaR_H1_ = -999.0;
  double AK8Gen_MergedWjets_minWminHmass_Leading_deltaR_H2_ = -999.0;
  double AK8Gen_MergedWjets_minWminHmass_SubLeading_Pt_ = -999.0;
  double AK8Gen_MergedWjets_minWminHmass_SubLeading_Eta_ = -999.0;
  double AK8Gen_MergedWjets_minWminHmass_SubLeading_Phi_ = -999.0;
  double AK8Gen_MergedWjets_minWminHmass_SubLeading_M_ = -999.0;
  double AK8Gen_MergedWjets_minWminHmass_SubLeading_deltaR_H1_ = -999.0;
  double AK8Gen_MergedWjets_minWminHmass_SubLeading_deltaR_H2_ = -999.0;
  double AK8Gen_MergedWjets_minWminHmass_LeadingSubLeading_DR_   = -999.0;

  double AK4GEN_AllResolved_onShellJet1_Pt_ = -999.0;
  double AK4GEN_AllResolved_onShellJet1_Eta_  = -999.0;
  double AK4GEN_AllResolved_onShellJet1_Phi_  = -999.0;
  double AK4GEN_AllResolved_onShellJet1_M_  = -999.0;
  double AK4GEN_AllResolved_onShellJet1_dR_q1_  = -999.0;
  double AK4GEN_AllResolved_onShellJet1_dR_q2_  = -999.0;
  double AK4GEN_AllResolved_onShellJet1_dR_q3_  = -999.0;
  double AK4GEN_AllResolved_onShellJet1_dR_q4_  = -999.0;
  double AK4GEN_AllResolved_onShellJet1_dR_g1_  = -999.0;
  double AK4GEN_AllResolved_onShellJet1_dR_g2_  = -999.0;
  double AK4GEN_AllResolved_onShellJet2_Pt_ = -999.0;
  double AK4GEN_AllResolved_onShellJet2_Eta_  = -999.0;
  double AK4GEN_AllResolved_onShellJet2_Phi_  = -999.0;
  double AK4GEN_AllResolved_onShellJet2_M_  = -999.0;
  double AK4GEN_AllResolved_onShellJets_dR_ = -999.0;
  double AK4GEN_AllResolved_onShellJet2_dR_q1_  = -999.0;
  double AK4GEN_AllResolved_onShellJet2_dR_q2_  = -999.0;
  double AK4GEN_AllResolved_onShellJet2_dR_q3_  = -999.0;
  double AK4GEN_AllResolved_onShellJet2_dR_q4_  = -999.0;
  double AK4GEN_AllResolved_onShellJet2_dR_g1_  = -999.0;
  double AK4GEN_AllResolved_onShellJet2_dR_g2_  = -999.0;
  double AK4GEN_AllResolved_offShellJet1_Pt_  = -999.0;
  double AK4GEN_AllResolved_offShellJet1_Eta_ = -999.0;
  double AK4GEN_AllResolved_offShellJet1_Phi_ = -999.0;
  double AK4GEN_AllResolved_offShellJet1_M_ = -999.0;
  double AK4GEN_AllResolved_offShellJet1_dR_q1_ = -999.0;
  double AK4GEN_AllResolved_offShellJet1_dR_q2_ = -999.0;
  double AK4GEN_AllResolved_offShellJet1_dR_q3_ = -999.0;
  double AK4GEN_AllResolved_offShellJet1_dR_q4_ = -999.0;
  double AK4GEN_AllResolved_offShellJet1_dR_g1_ = -999.0;
  double AK4GEN_AllResolved_offShellJet1_dR_g2_ = -999.0;
  double AK4GEN_AllResolved_offShellJet2_Pt_  = -999.0;
  double AK4GEN_AllResolved_offShellJet2_Eta_ = -999.0;
  double AK4GEN_AllResolved_offShellJet2_Phi_ = -999.0;
  double AK4GEN_AllResolved_offShellJet2_M_ = -999.0;
  double AK4GEN_AllResolved_offShellJets_dR_  = -999.0;
  double AK4GEN_AllResolved_onShelloffShellJets_dR_ = -999.0;
  double AK4GEN_AllResolved_offShellJet2_dR_q1_ = -999.0;
  double AK4GEN_AllResolved_offShellJet2_dR_q2_ = -999.0;
  double AK4GEN_AllResolved_offShellJet2_dR_q3_ = -999.0;
  double AK4GEN_AllResolved_offShellJet2_dR_q4_ = -999.0;
  double AK4GEN_AllResolved_offShellJet2_dR_g1_ = -999.0;
  double AK4GEN_AllResolved_offShellJet2_dR_g2_ = -999.0;

  double AK4GEN_AllResolved_onShellWboson_Pt_ = -999.0;
  double AK4GEN_AllResolved_onShellWboson_Eta_  = -999.0;
  double AK4GEN_AllResolved_onShellWboson_Phi_  = -999.0;
  double AK4GEN_AllResolved_onShellWboson_M_  = -999.0;
  double AK4GEN_AllResolved_onShellWboson_dR_W0PID_ = -999.0;
  double AK4GEN_AllResolved_onShellWboson_dR_W1PID_ = -999.0;
  double AK4GEN_AllResolved_offShellWboson_Pt_  = -999.0;
  double AK4GEN_AllResolved_offShellWboson_Eta_ = -999.0;
  double AK4GEN_AllResolved_offShellWboson_Phi_ = -999.0;
  double AK4GEN_AllResolved_offShellWboson_M_ = -999.0;
  double AK4GEN_AllResolved_offShellWboson_dR_W0PID_  = -999.0;
  double AK4GEN_AllResolved_offShellWboson_dR_W1PID_  = -999.0;
  double AK4GEN_AllResolved_Higgs_Pt_ = -999.0;
  double AK4GEN_AllResolved_Higgs_Eta_  = -999.0;
  double AK4GEN_AllResolved_Higgs_Phi_  = -999.0;
  double AK4GEN_AllResolved_Higgs_M_  = -999.0;
  double AK4GEN_AllResolved_Higgs_DR_Higgs0PID_ = -999.0;
  double AK4GEN_AllResolved_Higgs_DR_Higgs1PID_ = -999.0;

};

//
//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
GenAnalyzer::GenAnalyzer(const edm::ParameterSet& iConfig)
{
  //now do what ever initialization is needed
  usesResource("TFileService");
  LHEEventToken = consumes<LHEEventProduct>(iConfig.getParameter<edm::InputTag>("LHEEventInputTag"));
  genParticlesToken = consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genParticlesInputTag"));
  genAK4jetToken = consumes<reco::GenJetCollection>(iConfig.getParameter<edm::InputTag>("ak4GenJetsInputTag"));
  genAK8jetToken = consumes<reco::GenJetCollection>(iConfig.getParameter<edm::InputTag>("ak8GenJetsInputTag"));
  Verbose_ = iConfig.getParameter<bool>("Verbose");
  OutPutFileName_ = iConfig.getParameter<std::string>("OutPutFileName");
}


GenAnalyzer::~GenAnalyzer()
{
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  delete outputFile_;
}

void GenAnalyzer::AddBranch(std::vector<std::string>* vec, std::string name){
  tree->Branch(name.c_str(),vec);
}

void GenAnalyzer::AddBranch(std::vector<double>* vec, std::string name){
  tree->Branch(name.c_str(),vec);
}

void GenAnalyzer::AddBranch(std::vector<int>* vec, std::string name){
  tree->Branch(name.c_str(),vec);
}

void GenAnalyzer::AddBranch(int* vec, std::string name){
  tree->Branch(name.c_str(),vec,(name+"/I").c_str());
}

void GenAnalyzer::AddBranch(double* vec, std::string name){
  tree->Branch(name.c_str(),vec,(name+"/D").c_str());
}


//////////////////////////////////
//// P A P E R   4 - V E C T O R   D E F I N I T I O N   O F   P H I   A N D   P H I 1
//////////////////////////////////
void GenAnalyzer::computeAngles(TLorentzVector thep4H, TLorentzVector thep4Z1, TLorentzVector thep4M11, TLorentzVector thep4M12, TLorentzVector thep4Z2, TLorentzVector thep4M21, TLorentzVector thep4M22, double& costheta1, double& costheta2, double& Phi, double& costhetastar, double& Phi1){
  ///////////////////////////////////////////////
  // check for z1/z2 convention, redefine all 4 vectors with convention
  ///////////////////////////////////////////////
  TLorentzVector p4H, p4Z1, p4M11, p4M12, p4Z2, p4M21, p4M22;
  p4H = thep4H;
  
  p4Z1 = thep4Z1; p4M11 = thep4M11; p4M12 = thep4M12;
  p4Z2 = thep4Z2; p4M21 = thep4M21; p4M22 = thep4M22;
  //// costhetastar
  TVector3 boostX = -(thep4H.BoostVector());
  TLorentzVector thep4Z1inXFrame( p4Z1 );
  TLorentzVector thep4Z2inXFrame( p4Z2 );
  thep4Z1inXFrame.Boost( boostX );
  thep4Z2inXFrame.Boost( boostX );
  TVector3 theZ1X_p3 = TVector3( thep4Z1inXFrame.X(), thep4Z1inXFrame.Y(), thep4Z1inXFrame.Z() );
  TVector3 theZ2X_p3 = TVector3( thep4Z2inXFrame.X(), thep4Z2inXFrame.Y(), thep4Z2inXFrame.Z() );
  costhetastar = theZ1X_p3.CosTheta();
  
  //// --------------------------- costheta1
  TVector3 boostV1 = -(thep4Z1.BoostVector());
  TLorentzVector p4M11_BV1( p4M11 );
  TLorentzVector p4M12_BV1( p4M12 );
  TLorentzVector p4M21_BV1( p4M21 );
  TLorentzVector p4M22_BV1( p4M22 );
  p4M11_BV1.Boost( boostV1 );
  p4M12_BV1.Boost( boostV1 );
  p4M21_BV1.Boost( boostV1 );
  p4M22_BV1.Boost( boostV1 );
  
  TLorentzVector p4V2_BV1 = p4M21_BV1 + p4M22_BV1;
  //// costheta1
  costheta1 = -p4V2_BV1.Vect().Dot( p4M11_BV1.Vect() )/p4V2_BV1.Vect().Mag()/p4M11_BV1.Vect().Mag();
  
  //// --------------------------- costheta2
  TVector3 boostV2 = -(thep4Z2.BoostVector());
  TLorentzVector p4M11_BV2( p4M11 );
  TLorentzVector p4M12_BV2( p4M12 );
  TLorentzVector p4M21_BV2( p4M21 );
  TLorentzVector p4M22_BV2( p4M22 );
  p4M11_BV2.Boost( boostV2 );
  p4M12_BV2.Boost( boostV2 );
  p4M21_BV2.Boost( boostV2 );
  p4M22_BV2.Boost( boostV2 );
  
  TLorentzVector p4V1_BV2 = p4M11_BV2 + p4M12_BV2;
  //// costheta2
  costheta2 = -p4V1_BV2.Vect().Dot( p4M21_BV2.Vect() )/p4V1_BV2.Vect().Mag()/p4M21_BV2.Vect().Mag();
  
  //// --------------------------- Phi and Phi1
  //    TVector3 boostX = -(thep4H.BoostVector());
  TLorentzVector p4M11_BX( p4M11 );
  TLorentzVector p4M12_BX( p4M12 );
  TLorentzVector p4M21_BX( p4M21 );
  TLorentzVector p4M22_BX( p4M22 );
  
  p4M11_BX.Boost( boostX );
  p4M12_BX.Boost( boostX );
  p4M21_BX.Boost( boostX );
  p4M22_BX.Boost( boostX );
  
  TVector3 tmp1 = p4M11_BX.Vect().Cross( p4M12_BX.Vect() );
  TVector3 tmp2 = p4M21_BX.Vect().Cross( p4M22_BX.Vect() );
  
  TVector3 normal1_BX( tmp1.X()/tmp1.Mag(), tmp1.Y()/tmp1.Mag(), tmp1.Z()/tmp1.Mag() );
  TVector3 normal2_BX( tmp2.X()/tmp2.Mag(), tmp2.Y()/tmp2.Mag(), tmp2.Z()/tmp2.Mag() );
  
  //// Phi
  TLorentzVector p4Z1_BX = p4M11_BX + p4M12_BX;
  double tmpSgnPhi = p4Z1_BX.Vect().Dot( normal1_BX.Cross( normal2_BX) );
  double sgnPhi = tmpSgnPhi/fabs(tmpSgnPhi);
  Phi = sgnPhi * acos( -1.*normal1_BX.Dot( normal2_BX) );
  
  //////////////
  TVector3 beamAxis(0,0,1);
  TVector3 tmp3 = (p4M11_BX + p4M12_BX).Vect();
  
  TVector3 p3V1_BX( tmp3.X()/tmp3.Mag(), tmp3.Y()/tmp3.Mag(), tmp3.Z()/tmp3.Mag() );
  TVector3 tmp4 = beamAxis.Cross( p3V1_BX );
  TVector3 normalSC_BX( tmp4.X()/tmp4.Mag(), tmp4.Y()/tmp4.Mag(), tmp4.Z()/tmp4.Mag() );
  
  //// Phi1
  double tmpSgnPhi1 = p4Z1_BX.Vect().Dot( normal1_BX.Cross( normalSC_BX) );
  double sgnPhi1 = tmpSgnPhi1/fabs(tmpSgnPhi1);
  Phi1 = sgnPhi1 * acos( normal1_BX.Dot( normalSC_BX) );
  
  //    std::cout << "extractAngles: " << std::endl;
  //    std::cout << "costhetastar = " << costhetastar << ", costheta1 = " << costheta1 << ", costheta2 = " << costheta2 << std::endl;
  //    std::cout << "Phi = " << Phi << ", Phi1 = " << Phi1 << std::endl;
}


void GenAnalyzer::SetBranches(){
  //AddBranch(&pdgID_,	"pdgID");
  //  AddBranch(&isMuMinus_ , "isMuMinus");
  //  AddBranch(&LHELeptPt_ ,	"LHELeptPt");
  //  AddBranch(&LHELeptEta_ ,	"LHELeptEta");
  //  AddBranch(&LHELeptPhi_ ,	"LHELeptPhi");
  //  AddBranch(&LHELeptM_ ,	"LHELeptM");
  //  AddBranch(&LHELeptE_ ,	"LHELeptE");
  //
  //  AddBranch(&LHEWeightIDs_, "LHEWeightIDs");
  //  AddBranch(&LHEWeights_, "LHEWeights");
  
  AddBranch(&gen_leading_photon_Pt_, "gen_leading_photon_Pt");
  AddBranch(&gen_leading_photon_Eta_, "gen_leading_photon_Eta");
  AddBranch(&gen_leading_photon_Phi_, "gen_leading_photon_Phi");
  AddBranch(&gen_leading_photon_M_, "gen_leading_photon_M");
  AddBranch(&gen_Subleading_photon_Pt_, "gen_Subleading_photon_Pt");
  AddBranch(&gen_Subleading_photon_Eta_, "gen_Subleading_photon_Eta");
  AddBranch(&gen_Subleading_photon_Phi_, "gen_Subleading_photon_Phi");
  AddBranch(&gen_Subleading_photon_M_, "gen_Subleading_photon_M");
  AddBranch(&gen_leading_WpJets_Pt_, "gen_leading_WpJets_Pt");
  AddBranch(&gen_leading_WpJets_Eta_, "gen_leading_WpJets_Eta");
  AddBranch(&gen_leading_WpJets_Phi_, "gen_leading_WpJets_Phi");
  AddBranch(&gen_leading_WpJets_M_, "gen_leading_WpJets_M");
  AddBranch(&gen_Subleading_WpJets_Pt_, "gen_Subleading_WpJets_Pt");
  AddBranch(&gen_Subleading_WpJets_Eta_, "gen_Subleading_WpJets_Eta");
  AddBranch(&gen_Subleading_WpJets_Phi_, "gen_Subleading_WpJets_Phi");
  AddBranch(&gen_Subleading_WpJets_M_, "gen_Subleading_WpJets_M");
  AddBranch(&gen_leading_WmJets_Pt_, "gen_leading_WmJets_Pt");
  AddBranch(&gen_leading_WmJets_Eta_, "gen_leading_WmJets_Eta");
  AddBranch(&gen_leading_WmJets_Phi_, "gen_leading_WmJets_Phi");
  AddBranch(&gen_leading_WmJets_M_, "gen_leading_WmJets_M");
  AddBranch(&gen_Subleading_WmJets_Pt_, "gen_Subleading_WmJets_Pt");
  AddBranch(&gen_Subleading_WmJets_Eta_, "gen_Subleading_WmJets_Eta");
  AddBranch(&gen_Subleading_WmJets_Phi_, "gen_Subleading_WmJets_Phi");
  AddBranch(&gen_Subleading_WmJets_M_, "gen_Subleading_WmJets_M");
  AddBranch(&gen_leading_WBoson_Pt_, "gen_leading_WBoson_Pt");
  AddBranch(&gen_leading_WBoson_Eta_, "gen_leading_WBoson_Eta");
  AddBranch(&gen_leading_WBoson_Phi_, "gen_leading_WBoson_Phi");
  AddBranch(&gen_leading_WBoson_M_, "gen_leading_WBoson_M");
  AddBranch(&gen_Subleading_WBoson_Pt_, "gen_Subleading_WBoson_Pt");
  AddBranch(&gen_Subleading_WBoson_Eta_, "gen_Subleading_WBoson_Eta");
  AddBranch(&gen_Subleading_WBoson_Phi_, "gen_Subleading_WBoson_Phi");
  AddBranch(&gen_Subleading_WBoson_M_, "gen_Subleading_WBoson_M");
  AddBranch(&gen_leading_Higgs_Pt_, "gen_leading_Higgs_Pt");
  AddBranch(&gen_leading_Higgs_Eta_, "gen_leading_Higgs_Eta");
  AddBranch(&gen_leading_Higgs_Phi_, "gen_leading_Higgs_Phi");
  AddBranch(&gen_leading_Higgs_M_, "gen_leading_Higgs_M");
  AddBranch(&gen_Subleading_Higgs_Pt_, "gen_Subleading_Higgs_Pt");
  AddBranch(&gen_Subleading_Higgs_Eta_, "gen_Subleading_Higgs_Eta");
  AddBranch(&gen_Subleading_Higgs_Phi_, "gen_Subleading_Higgs_Phi");
  AddBranch(&gen_Subleading_Higgs_M_, "gen_Subleading_Higgs_M");
  AddBranch(&gen_deltaR_Photon0_Photon1_, "gen_deltaR_Photon0_Photon1");
  AddBranch(&gen_deltaR_Photon0_WmJ0_, "gen_deltaR_Photon0_WmJ0");
  AddBranch(&gen_deltaR_Photon0_WmJ1_, "gen_deltaR_Photon0_WmJ1");
  AddBranch(&gen_deltaR_Photon0_WpJ0_, "gen_deltaR_Photon0_WpJ0");
  AddBranch(&gen_deltaR_Photon0_WpJ1_, "gen_deltaR_Photon0_WpJ1");
  AddBranch(&gen_deltaR_Photon1_WmJ0_, "gen_deltaR_Photon1_WmJ0");
  AddBranch(&gen_deltaR_Photon1_WmJ1_, "gen_deltaR_Photon1_WmJ1");
  AddBranch(&gen_deltaR_Photon1_WpJ0_, "gen_deltaR_Photon1_WpJ0");
  AddBranch(&gen_deltaR_Photon1_WpJ1_, "gen_deltaR_Photon1_WpJ1");
  AddBranch(&gen_deltaR_WpJ0_WpJ1_, "gen_deltaR_WpJ0_WpJ1");
  AddBranch(&gen_deltaR_WpJ0_WmJ0_, "gen_deltaR_WpJ0_WmJ0");
  AddBranch(&gen_deltaR_WpJ0_WmJ1_, "gen_deltaR_WpJ0_WmJ1");
  AddBranch(&gen_deltaR_WpJ1_WmJ0_, "gen_deltaR_WpJ1_WmJ0");
  AddBranch(&gen_deltaR_WpJ1_WmJ1_, "gen_deltaR_WpJ1_WmJ1");
  AddBranch(&gen_deltaR_WmJ0_WmJ1_, "gen_deltaR_WmJ0_WmJ1");
  AddBranch(&gen_deltaR_Wp_Wm_, "gen_deltaR_Wp_Wm");
  AddBranch(&gen_deltaR_H1_H2_  , "gen_deltaR_H1_H2");
  
  AddBranch(&genJetAK4_njets_, "genJetAK4_njets");
  AddBranch(&genJetAK4_leading_Pt_, "genJetAK4_leading_Pt");
  AddBranch(&genJetAK4_leading_Eta_, "genJetAK4_leading_Eta");
  AddBranch(&genJetAK4_leading_Phi_, "genJetAK4_leading_Phi");
  AddBranch(&genJetAK4_leading_M_, "genJetAK4_leading_M");
  AddBranch(&genJetAK4_leading_Energy_, "genJetAK4_leading_Energy");
  AddBranch(&genJetAK4_Subleading_Pt_, "genJetAK4_Subleading_Pt");
  AddBranch(&genJetAK4_Subleading_Eta_, "genJetAK4_Subleading_Eta");
  AddBranch(&genJetAK4_Subleading_Phi_, "genJetAK4_Subleading_Phi");
  AddBranch(&genJetAK4_Subleading_M_, "genJetAK4_Subleading_M");
  AddBranch(&genJetAK4_Subleading_Energy_, "genJetAK4_Subleading_Energy");
  
  AddBranch(&AK8Gen_HiggsJet_njets_, "genJetAK8_njets");
  AddBranch(&AK8Gen_HiggsJet_MaxPt_Pt_, "AK8Gen_HiggsJet_MaxPt_Pt");
  AddBranch(&AK8Gen_HiggsJet_MaxPt_Eta_, "AK8Gen_HiggsJet_MaxPt_Eta");
  AddBranch(&AK8Gen_HiggsJet_MaxPt_Phi_, "AK8Gen_HiggsJet_MaxPt_Phi");
  AddBranch(&AK8Gen_HiggsJet_MaxPt_M_  , "AK8Gen_HiggsJet_MaxPt_M");
  
  AddBranch(&AK8Gen_HiggsJet_MaxPt_deltaR_H1_, "AK8Gen_HiggsJet_MaxPt_deltaR_H1");
  AddBranch(&AK8Gen_HiggsJet_MaxPt_deltaR_H2_, "AK8Gen_HiggsJet_MaxPt_deltaR_H2");
  
  AddBranch(&AK8Gen_HiggsJet_minDMass_Pt_, "AK8Gen_HiggsJet_minDMass_Pt");
  AddBranch(&AK8Gen_HiggsJet_minDMass_Eta_, "AK8Gen_HiggsJet_minDMass_Eta");
  AddBranch(&AK8Gen_HiggsJet_minDMass_Phi_, "AK8Gen_HiggsJet_minDMass_Phi");
  AddBranch(&AK8Gen_HiggsJet_minDMass_M_  , "AK8Gen_HiggsJet_minDMass_M");
  
  AddBranch(&AK8Gen_HiggsJet_minDMass_deltaR_H1_, "AK8Gen_HiggsJet_minDMass_deltaR_H1");
  AddBranch(&AK8Gen_HiggsJet_minDMass_deltaR_H2_, "AK8Gen_HiggsJet_minDMass_deltaR_H2");
  
  AddBranch(&AK8Gen_MergedWjets_MaxPt_Leading_Pt_,"AK8Gen_MergedWjets_MaxPt_Leading_Pt");
  AddBranch(&AK8Gen_MergedWjets_MaxPt_Leading_Eta_,"AK8Gen_MergedWjets_MaxPt_Leading_Eta");
  AddBranch(&AK8Gen_MergedWjets_MaxPt_Leading_Phi_,"AK8Gen_MergedWjets_MaxPt_Leading_Phi");
  AddBranch(&AK8Gen_MergedWjets_MaxPt_Leading_M_,"AK8Gen_MergedWjets_MaxPt_Leading_M");
  AddBranch(&AK8Gen_MergedWjets_MaxPt_Leading_deltaR_H1_,"AK8Gen_MergedWjets_MaxPt_Leading_deltaR_H1");
  AddBranch(&AK8Gen_MergedWjets_MaxPt_Leading_deltaR_H2_,"AK8Gen_MergedWjets_MaxPt_Leading_deltaR_H2");
  AddBranch(&AK8Gen_MergedWjets_MaxPt_SubLeading_Pt_,"AK8Gen_MergedWjets_MaxPt_SubLeading_Pt");
  AddBranch(&AK8Gen_MergedWjets_MaxPt_SubLeading_Eta_,"AK8Gen_MergedWjets_MaxPt_SubLeading_Eta");
  AddBranch(&AK8Gen_MergedWjets_MaxPt_SubLeading_Phi_,"AK8Gen_MergedWjets_MaxPt_SubLeading_Phi");
  AddBranch(&AK8Gen_MergedWjets_MaxPt_SubLeading_M_,"AK8Gen_MergedWjets_MaxPt_SubLeading_M");
  AddBranch(&AK8Gen_MergedWjets_MaxPt_SubLeading_deltaR_H1_,"AK8Gen_MergedWjets_MaxPt_SubLeading_deltaR_H1");
  AddBranch(&AK8Gen_MergedWjets_MaxPt_SubLeading_deltaR_H2_,"AK8Gen_MergedWjets_MaxPt_SubLeading_deltaR_H2");
  AddBranch(&AK8Gen_MergedWjets_MaxPt_LeadingSubLeading_DR_,"AK8Gen_MergedWjets_MaxPt_LeadingSubLeading_DR");

  AddBranch(&AK8Gen_MergedWjets_minDMass_Leading_Pt_,"AK8Gen_MergedWjets_minDMass_Leading_Pt");
  AddBranch(&AK8Gen_MergedWjets_minDMass_Leading_Eta_,"AK8Gen_MergedWjets_minDMass_Leading_Eta");
  AddBranch(&AK8Gen_MergedWjets_minDMass_Leading_Phi_,"AK8Gen_MergedWjets_minDMass_Leading_Phi");
  AddBranch(&AK8Gen_MergedWjets_minDMass_Leading_M_,"AK8Gen_MergedWjets_minDMass_Leading_M");
  AddBranch(&AK8Gen_MergedWjets_minDMass_Leading_deltaR_H1_,"AK8Gen_MergedWjets_minDMass_Leading_deltaR_H1");
  AddBranch(&AK8Gen_MergedWjets_minDMass_Leading_deltaR_H2_,"AK8Gen_MergedWjets_minDMass_Leading_deltaR_H2");
  AddBranch(&AK8Gen_MergedWjets_minDMass_SubLeading_Pt_,"AK8Gen_MergedWjets_minDMass_SubLeading_Pt");
  AddBranch(&AK8Gen_MergedWjets_minDMass_SubLeading_Eta_,"AK8Gen_MergedWjets_minDMass_SubLeading_Eta");
  AddBranch(&AK8Gen_MergedWjets_minDMass_SubLeading_Phi_,"AK8Gen_MergedWjets_minDMass_SubLeading_Phi");
  AddBranch(&AK8Gen_MergedWjets_minDMass_SubLeading_M_,"AK8Gen_MergedWjets_minDMass_SubLeading_M");
  AddBranch(&AK8Gen_MergedWjets_minDMass_SubLeading_deltaR_H1_,"AK8Gen_MergedWjets_minDMass_SubLeading_deltaR_H1");
  AddBranch(&AK8Gen_MergedWjets_minDMass_SubLeading_deltaR_H2_,"AK8Gen_MergedWjets_minDMass_SubLeading_deltaR_H2");
  AddBranch(&AK8Gen_MergedWjets_minDMass_LeadingSubLeading_DR_,"AK8Gen_MergedWjets_minDMass_LeadingSubLeading_DR");

  AddBranch(&AK8Gen_MergedWjets_minWminHmass_Leading_Pt_,"AK8Gen_MergedWjets_minWminHmass_Leading_Pt");
  AddBranch(&AK8Gen_MergedWjets_minWminHmass_Leading_Eta_,"AK8Gen_MergedWjets_minWminHmass_Leading_Eta");
  AddBranch(&AK8Gen_MergedWjets_minWminHmass_Leading_Phi_,"AK8Gen_MergedWjets_minWminHmass_Leading_Phi");
  AddBranch(&AK8Gen_MergedWjets_minWminHmass_Leading_M_,"AK8Gen_MergedWjets_minWminHmass_Leading_M");
  AddBranch(&AK8Gen_MergedWjets_minWminHmass_Leading_deltaR_H1_,"AK8Gen_MergedWjets_minWminHmass_Leading_deltaR_H1");
  AddBranch(&AK8Gen_MergedWjets_minWminHmass_Leading_deltaR_H2_,"AK8Gen_MergedWjets_minWminHmass_Leading_deltaR_H2");
  AddBranch(&AK8Gen_MergedWjets_minWminHmass_SubLeading_Pt_,"AK8Gen_MergedWjets_minWminHmass_SubLeading_Pt");
  AddBranch(&AK8Gen_MergedWjets_minWminHmass_SubLeading_Eta_,"AK8Gen_MergedWjets_minWminHmass_SubLeading_Eta");
  AddBranch(&AK8Gen_MergedWjets_minWminHmass_SubLeading_Phi_,"AK8Gen_MergedWjets_minWminHmass_SubLeading_Phi");
  AddBranch(&AK8Gen_MergedWjets_minWminHmass_SubLeading_M_,"AK8Gen_MergedWjets_minWminHmass_SubLeading_M");
  AddBranch(&AK8Gen_MergedWjets_minWminHmass_SubLeading_deltaR_H1_,"AK8Gen_MergedWjets_minWminHmass_SubLeading_deltaR_H1");
  AddBranch(&AK8Gen_MergedWjets_minWminHmass_SubLeading_deltaR_H2_,"AK8Gen_MergedWjets_minWminHmass_SubLeading_deltaR_H2");
  AddBranch(&AK8Gen_MergedWjets_minWminHmass_LeadingSubLeading_DR_,"AK8Gen_MergedWjets_minWminHmass_LeadingSubLeading_DR");

  AddBranch(&AK4GEN_AllResolved_onShellJet1_Pt_,"AK4GEN_AllResolved_onShellJet1_Pt");
  AddBranch(&AK4GEN_AllResolved_onShellJet1_Eta_,"AK4GEN_AllResolved_onShellJet1_Eta");
  AddBranch(&AK4GEN_AllResolved_onShellJet1_Phi_,"AK4GEN_AllResolved_onShellJet1_Phi");
  AddBranch(&AK4GEN_AllResolved_onShellJet1_M_,"AK4GEN_AllResolved_onShellJet1_M");
  AddBranch(&AK4GEN_AllResolved_onShellJet1_dR_q1_,"AK4GEN_AllResolved_onShellJet1_dR_q1");
  AddBranch(&AK4GEN_AllResolved_onShellJet1_dR_q2_,"AK4GEN_AllResolved_onShellJet1_dR_q2");
  AddBranch(&AK4GEN_AllResolved_onShellJet1_dR_q3_,"AK4GEN_AllResolved_onShellJet1_dR_q3");
  AddBranch(&AK4GEN_AllResolved_onShellJet1_dR_q4_,"AK4GEN_AllResolved_onShellJet1_dR_q4");
  AddBranch(&AK4GEN_AllResolved_onShellJet1_dR_g1_,"AK4GEN_AllResolved_onShellJet1_dR_g1");
  AddBranch(&AK4GEN_AllResolved_onShellJet1_dR_g2_,"AK4GEN_AllResolved_onShellJet1_dR_g2");
  AddBranch(&AK4GEN_AllResolved_onShellJet2_Pt_,"AK4GEN_AllResolved_onShellJet2_Pt");
  AddBranch(&AK4GEN_AllResolved_onShellJet2_Eta_,"AK4GEN_AllResolved_onShellJet2_Eta");
  AddBranch(&AK4GEN_AllResolved_onShellJet2_Phi_,"AK4GEN_AllResolved_onShellJet2_Phi");
  AddBranch(&AK4GEN_AllResolved_onShellJet2_M_,"AK4GEN_AllResolved_onShellJet2_M");
  AddBranch(&AK4GEN_AllResolved_onShellJets_dR_,"AK4GEN_AllResolved_onShellJets_dR");
  AddBranch(&AK4GEN_AllResolved_onShellJet2_dR_q1_,"AK4GEN_AllResolved_onShellJet2_dR_q1");
  AddBranch(&AK4GEN_AllResolved_onShellJet2_dR_q2_,"AK4GEN_AllResolved_onShellJet2_dR_q2");
  AddBranch(&AK4GEN_AllResolved_onShellJet2_dR_q3_,"AK4GEN_AllResolved_onShellJet2_dR_q3");
  AddBranch(&AK4GEN_AllResolved_onShellJet2_dR_q4_,"AK4GEN_AllResolved_onShellJet2_dR_q4");
  AddBranch(&AK4GEN_AllResolved_onShellJet2_dR_g1_,"AK4GEN_AllResolved_onShellJet2_dR_g1");
  AddBranch(&AK4GEN_AllResolved_onShellJet2_dR_g2_,"AK4GEN_AllResolved_onShellJet2_dR_g2");
  AddBranch(&AK4GEN_AllResolved_offShellJet1_Pt_,"AK4GEN_AllResolved_offShellJet1_Pt");
  AddBranch(&AK4GEN_AllResolved_offShellJet1_Eta_,"AK4GEN_AllResolved_offShellJet1_Eta");
  AddBranch(&AK4GEN_AllResolved_offShellJet1_Phi_,"AK4GEN_AllResolved_offShellJet1_Phi");
  AddBranch(&AK4GEN_AllResolved_offShellJet1_M_,"AK4GEN_AllResolved_offShellJet1_M");
  AddBranch(&AK4GEN_AllResolved_offShellJet1_dR_q1_,"AK4GEN_AllResolved_offShellJet1_dR_q1");
  AddBranch(&AK4GEN_AllResolved_offShellJet1_dR_q2_,"AK4GEN_AllResolved_offShellJet1_dR_q2");
  AddBranch(&AK4GEN_AllResolved_offShellJet1_dR_q3_,"AK4GEN_AllResolved_offShellJet1_dR_q3");
  AddBranch(&AK4GEN_AllResolved_offShellJet1_dR_q4_,"AK4GEN_AllResolved_offShellJet1_dR_q4");
  AddBranch(&AK4GEN_AllResolved_offShellJet1_dR_g1_,"AK4GEN_AllResolved_offShellJet1_dR_g1");
  AddBranch(&AK4GEN_AllResolved_offShellJet1_dR_g2_,"AK4GEN_AllResolved_offShellJet1_dR_g2");
  AddBranch(&AK4GEN_AllResolved_offShellJet2_Pt_,"AK4GEN_AllResolved_offShellJet2_Pt");
  AddBranch(&AK4GEN_AllResolved_offShellJet2_Eta_,"AK4GEN_AllResolved_offShellJet2_Eta");
  AddBranch(&AK4GEN_AllResolved_offShellJet2_Phi_,"AK4GEN_AllResolved_offShellJet2_Phi");
  AddBranch(&AK4GEN_AllResolved_offShellJet2_M_,"AK4GEN_AllResolved_offShellJet2_M");
  AddBranch(&AK4GEN_AllResolved_offShellJets_dR_,"AK4GEN_AllResolved_offShellJets_dR");
  AddBranch(&AK4GEN_AllResolved_onShelloffShellJets_dR_,"AK4GEN_AllResolved_onShelloffShellJets_dR");
  AddBranch(&AK4GEN_AllResolved_offShellJet2_dR_q1_,"AK4GEN_AllResolved_offShellJet2_dR_q1");
  AddBranch(&AK4GEN_AllResolved_offShellJet2_dR_q2_,"AK4GEN_AllResolved_offShellJet2_dR_q2");
  AddBranch(&AK4GEN_AllResolved_offShellJet2_dR_q3_,"AK4GEN_AllResolved_offShellJet2_dR_q3");
  AddBranch(&AK4GEN_AllResolved_offShellJet2_dR_q4_,"AK4GEN_AllResolved_offShellJet2_dR_q4");
  AddBranch(&AK4GEN_AllResolved_offShellJet2_dR_g1_,"AK4GEN_AllResolved_offShellJet2_dR_g1");
  AddBranch(&AK4GEN_AllResolved_offShellJet2_dR_g2_,"AK4GEN_AllResolved_offShellJet2_dR_g2");

  AddBranch(&AK4GEN_AllResolved_onShellWboson_Pt_,"AK4GEN_AllResolved_onShellWboson_Pt");
  AddBranch(&AK4GEN_AllResolved_onShellWboson_Eta_,"AK4GEN_AllResolved_onShellWboson_Eta");
  AddBranch(&AK4GEN_AllResolved_onShellWboson_Phi_,"AK4GEN_AllResolved_onShellWboson_Phi");
  AddBranch(&AK4GEN_AllResolved_onShellWboson_M_,"AK4GEN_AllResolved_onShellWboson_M");
  AddBranch(&AK4GEN_AllResolved_onShellWboson_dR_W0PID_,"AK4GEN_AllResolved_onShellWboson_dR_W0PID");
  AddBranch(&AK4GEN_AllResolved_onShellWboson_dR_W1PID_,"AK4GEN_AllResolved_onShellWboson_dR_W1PID");
  AddBranch(&AK4GEN_AllResolved_offShellWboson_Pt_,"AK4GEN_AllResolved_offShellWboson_Pt");
  AddBranch(&AK4GEN_AllResolved_offShellWboson_Eta_,"AK4GEN_AllResolved_offShellWboson_Eta");
  AddBranch(&AK4GEN_AllResolved_offShellWboson_Phi_,"AK4GEN_AllResolved_offShellWboson_Phi");
  AddBranch(&AK4GEN_AllResolved_offShellWboson_M_,"AK4GEN_AllResolved_offShellWboson_M");
  AddBranch(&AK4GEN_AllResolved_offShellWboson_dR_W0PID_,"AK4GEN_AllResolved_offShellWboson_dR_W0PID");
  AddBranch(&AK4GEN_AllResolved_offShellWboson_dR_W1PID_,"AK4GEN_AllResolved_offShellWboson_dR_W1PID");
  AddBranch(&AK4GEN_AllResolved_Higgs_Pt_,"AK4GEN_AllResolved_Higgs_Pt");
  AddBranch(&AK4GEN_AllResolved_Higgs_Eta_,"AK4GEN_AllResolved_Higgs_Eta");
  AddBranch(&AK4GEN_AllResolved_Higgs_Phi_,"AK4GEN_AllResolved_Higgs_Phi");
  AddBranch(&AK4GEN_AllResolved_Higgs_M_,"AK4GEN_AllResolved_Higgs_M");
  AddBranch(&AK4GEN_AllResolved_Higgs_DR_Higgs0PID_,"AK4GEN_AllResolved_Higgs_DR_Higgs0PID");
  AddBranch(&AK4GEN_AllResolved_Higgs_DR_Higgs1PID_,"AK4GEN_AllResolved_Higgs_DR_Higgs1PID");

}

void GenAnalyzer::Clear(){
  //pdgID_.clear();
  //  LHEWeightIDs_.clear();
  //  LHEWeights_.clear();
  
  gen_leading_photon_Pt_ = -999.0;
  gen_leading_photon_Eta_ = -999.0;
  gen_leading_photon_Phi_ = -999.0;
  gen_leading_photon_M_ = -999.0;
  gen_Subleading_photon_Pt_ = -999.0;
  gen_Subleading_photon_Eta_ = -999.0;
  gen_Subleading_photon_Phi_ = -999.0;
  gen_Subleading_photon_M_ = -999.0;
  gen_leading_WpJets_Pt_ = -999.0;
  gen_leading_WpJets_Eta_ = -999.0;
  gen_leading_WpJets_Phi_ = -999.0;
  gen_leading_WpJets_M_ = -999.0;
  gen_Subleading_WpJets_Pt_ = -999.0;
  gen_Subleading_WpJets_Eta_ = -999.0;
  gen_Subleading_WpJets_Phi_ = -999.0;
  gen_Subleading_WpJets_M_ = -999.0;
  gen_leading_WmJets_Pt_ = -999.0;
  gen_leading_WmJets_Eta_ = -999.0;
  gen_leading_WmJets_Phi_ = -999.0;
  gen_leading_WmJets_M_ = -999.0;
  gen_Subleading_WmJets_Pt_ = -999.0;
  gen_Subleading_WmJets_Eta_ = -999.0;
  gen_Subleading_WmJets_Phi_ = -999.0;
  gen_Subleading_WmJets_M_ = -999.0;
  gen_leading_WBoson_Pt_ = -999.0;
  gen_leading_WBoson_Eta_ = -999.0;
  gen_leading_WBoson_Phi_ = -999.0;
  gen_leading_WBoson_M_ = -999.0;
  gen_Subleading_WBoson_Pt_ = -999.0;
  gen_Subleading_WBoson_Eta_ = -999.0;
  gen_Subleading_WBoson_Phi_ = -999.0;
  gen_Subleading_WBoson_M_ = -999.0;
  gen_leading_Higgs_Pt_ = -999.0;
  gen_leading_Higgs_Eta_ = -999.0;
  gen_leading_Higgs_Phi_ = -999.0;
  gen_leading_Higgs_M_ = -999.0;
  gen_Subleading_Higgs_Pt_ = -999.0;
  gen_Subleading_Higgs_Eta_ = -999.0;
  gen_Subleading_Higgs_Phi_ = -999.0;
  gen_Subleading_Higgs_M_ = -999.0;
  gen_deltaR_Photon0_Photon1_ = -999.0;
  gen_deltaR_Photon0_WmJ0_ = -999.0;
  gen_deltaR_Photon0_WmJ1_ = -999.0;
  gen_deltaR_Photon0_WpJ0_ = -999.0;
  gen_deltaR_Photon0_WpJ1_ = -999.0;
  gen_deltaR_Photon1_WmJ0_ = -999.0;
  gen_deltaR_Photon1_WmJ1_ = -999.0;
  gen_deltaR_Photon1_WpJ0_ = -999.0;
  gen_deltaR_Photon1_WpJ1_ = -999.0;
  gen_deltaR_WpJ0_WpJ1_ = -999.0;
  gen_deltaR_WpJ0_WmJ0_ = -999.0;
  gen_deltaR_WpJ0_WmJ1_ = -999.0;
  gen_deltaR_WpJ1_WmJ0_ = -999.0;
  gen_deltaR_WpJ1_WmJ1_ = -999.0;
  gen_deltaR_WmJ0_WmJ1_ = -999.0;
  gen_deltaR_Wp_Wm_ = -999.0;
  gen_deltaR_H1_H2_ = -999.0;
  
  genJetAK4_njets_ = -999.0;
  genJetAK4_leading_Pt_ = -999.0;
  genJetAK4_leading_Eta_ = -999.0;
  genJetAK4_leading_Phi_ = -999.0;
  genJetAK4_leading_M_ = -999.0;
  genJetAK4_leading_Energy_ = -999.0;
  genJetAK4_Subleading_Pt_ = -999.0;
  genJetAK4_Subleading_Eta_ = -999.0;
  genJetAK4_Subleading_Phi_ = -999.0;
  genJetAK4_Subleading_M_ = -999.0;
  genJetAK4_Subleading_Energy_ = -999.0;
  
  AK8Gen_HiggsJet_njets_ = -999.0;
  AK8Gen_HiggsJet_MaxPt_Pt_ = -999.0;
  AK8Gen_HiggsJet_MaxPt_Eta_ = -999.0;
  AK8Gen_HiggsJet_MaxPt_Phi_ = -999.0;
  AK8Gen_HiggsJet_MaxPt_M_ = -999.0;
  
  AK8Gen_HiggsJet_MaxPt_deltaR_H1_ = -999.0;
  AK8Gen_HiggsJet_MaxPt_deltaR_H2_ = -999.0;
  
  AK8Gen_HiggsJet_minDMass_Pt_ = -999.0;
  AK8Gen_HiggsJet_minDMass_Eta_ = -999.0;
  AK8Gen_HiggsJet_minDMass_Phi_ = -999.0;
  AK8Gen_HiggsJet_minDMass_M_ = -999.0;
  
  AK8Gen_HiggsJet_minDMass_deltaR_H1_ = -999.0;
  AK8Gen_HiggsJet_minDMass_deltaR_H2_ = -999.0;

  AK8Gen_MergedWjets_MaxPt_Leading_Pt_ = -999.0;
  AK8Gen_MergedWjets_MaxPt_Leading_Eta_ = -999.0;
  AK8Gen_MergedWjets_MaxPt_Leading_Phi_ = -999.0;
  AK8Gen_MergedWjets_MaxPt_Leading_M_ = -999.0;
  AK8Gen_MergedWjets_MaxPt_Leading_deltaR_H1_ = -999.0;
  AK8Gen_MergedWjets_MaxPt_Leading_deltaR_H2_ = -999.0;
  AK8Gen_MergedWjets_MaxPt_SubLeading_Pt_ = -999.0;
  AK8Gen_MergedWjets_MaxPt_SubLeading_Eta_ = -999.0;
  AK8Gen_MergedWjets_MaxPt_SubLeading_Phi_ = -999.0;
  AK8Gen_MergedWjets_MaxPt_SubLeading_M_ = -999.0;
  AK8Gen_MergedWjets_MaxPt_SubLeading_deltaR_H1_ = -999.0;
  AK8Gen_MergedWjets_MaxPt_SubLeading_deltaR_H2_   = -999.0;
  AK8Gen_MergedWjets_MaxPt_LeadingSubLeading_DR_ = -999.0;

  AK8Gen_MergedWjets_minDMass_Leading_Pt_ = -999.0;
  AK8Gen_MergedWjets_minDMass_Leading_Eta_ = -999.0;
  AK8Gen_MergedWjets_minDMass_Leading_Phi_ = -999.0;
  AK8Gen_MergedWjets_minDMass_Leading_M_ = -999.0;
  AK8Gen_MergedWjets_minDMass_Leading_deltaR_H1_ = -999.0;
  AK8Gen_MergedWjets_minDMass_Leading_deltaR_H2_ = -999.0;
  AK8Gen_MergedWjets_minDMass_SubLeading_Pt_ = -999.0;
  AK8Gen_MergedWjets_minDMass_SubLeading_Eta_ = -999.0;
  AK8Gen_MergedWjets_minDMass_SubLeading_Phi_ = -999.0;
  AK8Gen_MergedWjets_minDMass_SubLeading_M_ = -999.0;
  AK8Gen_MergedWjets_minDMass_SubLeading_deltaR_H1_ = -999.0;
  AK8Gen_MergedWjets_minDMass_SubLeading_deltaR_H2_ = -999.0;
  AK8Gen_MergedWjets_minDMass_LeadingSubLeading_DR_  = -999.0;

  AK8Gen_MergedWjets_minWminHmass_Leading_Pt_ = -999.0;
  AK8Gen_MergedWjets_minWminHmass_Leading_Eta_ = -999.0;
  AK8Gen_MergedWjets_minWminHmass_Leading_Phi_ = -999.0;
  AK8Gen_MergedWjets_minWminHmass_Leading_M_ = -999.0;
  AK8Gen_MergedWjets_minWminHmass_Leading_deltaR_H1_ = -999.0;
  AK8Gen_MergedWjets_minWminHmass_Leading_deltaR_H2_ = -999.0;
  AK8Gen_MergedWjets_minWminHmass_SubLeading_Pt_ = -999.0;
  AK8Gen_MergedWjets_minWminHmass_SubLeading_Eta_ = -999.0;
  AK8Gen_MergedWjets_minWminHmass_SubLeading_Phi_ = -999.0;
  AK8Gen_MergedWjets_minWminHmass_SubLeading_M_ = -999.0;
  AK8Gen_MergedWjets_minWminHmass_SubLeading_deltaR_H1_ = -999.0;
  AK8Gen_MergedWjets_minWminHmass_SubLeading_deltaR_H2_ = -999.0;
  AK8Gen_MergedWjets_minWminHmass_LeadingSubLeading_DR_   = -999.0;

  AK4GEN_AllResolved_onShellJet1_Pt_  = -999.0;
  AK4GEN_AllResolved_onShellJet1_Eta_ = -999.0;
  AK4GEN_AllResolved_onShellJet1_Phi_ = -999.0;
  AK4GEN_AllResolved_onShellJet1_M_ = -999.0;
  AK4GEN_AllResolved_onShellJet1_dR_q1_ = -999.0;
  AK4GEN_AllResolved_onShellJet1_dR_q2_ = -999.0;
  AK4GEN_AllResolved_onShellJet1_dR_q3_ = -999.0;
  AK4GEN_AllResolved_onShellJet1_dR_q4_ = -999.0;
  AK4GEN_AllResolved_onShellJet1_dR_g1_ = -999.0;
  AK4GEN_AllResolved_onShellJet1_dR_g2_ = -999.0;
  AK4GEN_AllResolved_onShellJet2_Pt_  = -999.0;
  AK4GEN_AllResolved_onShellJet2_Eta_ = -999.0;
  AK4GEN_AllResolved_onShellJet2_Phi_ = -999.0;
  AK4GEN_AllResolved_onShellJet2_M_ = -999.0;
  AK4GEN_AllResolved_onShellJets_dR_  = -999.0;
  AK4GEN_AllResolved_onShellJet2_dR_q1_ = -999.0;
  AK4GEN_AllResolved_onShellJet2_dR_q2_ = -999.0;
  AK4GEN_AllResolved_onShellJet2_dR_q3_ = -999.0;
  AK4GEN_AllResolved_onShellJet2_dR_q4_ = -999.0;
  AK4GEN_AllResolved_onShellJet2_dR_g1_ = -999.0;
  AK4GEN_AllResolved_onShellJet2_dR_g2_ = -999.0;
  AK4GEN_AllResolved_offShellJet1_Pt_ = -999.0;
  AK4GEN_AllResolved_offShellJet1_Eta_  = -999.0;
  AK4GEN_AllResolved_offShellJet1_Phi_  = -999.0;
  AK4GEN_AllResolved_offShellJet1_M_  = -999.0;
  AK4GEN_AllResolved_offShellJet1_dR_q1_  = -999.0;
  AK4GEN_AllResolved_offShellJet1_dR_q2_  = -999.0;
  AK4GEN_AllResolved_offShellJet1_dR_q3_  = -999.0;
  AK4GEN_AllResolved_offShellJet1_dR_q4_  = -999.0;
  AK4GEN_AllResolved_offShellJet1_dR_g1_  = -999.0;
  AK4GEN_AllResolved_offShellJet1_dR_g2_  = -999.0;
  AK4GEN_AllResolved_offShellJet2_Pt_ = -999.0;
  AK4GEN_AllResolved_offShellJet2_Eta_  = -999.0;
  AK4GEN_AllResolved_offShellJet2_Phi_  = -999.0;
  AK4GEN_AllResolved_offShellJet2_M_  = -999.0;
  AK4GEN_AllResolved_offShellJets_dR_ = -999.0;
  AK4GEN_AllResolved_onShelloffShellJets_dR_  = -999.0;
  AK4GEN_AllResolved_offShellJet2_dR_q1_  = -999.0;
  AK4GEN_AllResolved_offShellJet2_dR_q2_  = -999.0;
  AK4GEN_AllResolved_offShellJet2_dR_q3_  = -999.0;
  AK4GEN_AllResolved_offShellJet2_dR_q4_  = -999.0;
  AK4GEN_AllResolved_offShellJet2_dR_g1_  = -999.0;
  AK4GEN_AllResolved_offShellJet2_dR_g2_  = -999.0;

  AK4GEN_AllResolved_onShellWboson_Pt_  = -999.0;
  AK4GEN_AllResolved_onShellWboson_Eta_ = -999.0;
  AK4GEN_AllResolved_onShellWboson_Phi_ = -999.0;
  AK4GEN_AllResolved_onShellWboson_M_ = -999.0;
  AK4GEN_AllResolved_onShellWboson_dR_W0PID_  = -999.0;
  AK4GEN_AllResolved_onShellWboson_dR_W1PID_  = -999.0;
  AK4GEN_AllResolved_offShellWboson_Pt_ = -999.0;
  AK4GEN_AllResolved_offShellWboson_Eta_  = -999.0;
  AK4GEN_AllResolved_offShellWboson_Phi_  = -999.0;
  AK4GEN_AllResolved_offShellWboson_M_  = -999.0;
  AK4GEN_AllResolved_offShellWboson_dR_W0PID_ = -999.0;
  AK4GEN_AllResolved_offShellWboson_dR_W1PID_ = -999.0;
  AK4GEN_AllResolved_Higgs_Pt_  = -999.0;
  AK4GEN_AllResolved_Higgs_Eta_ = -999.0;
  AK4GEN_AllResolved_Higgs_Phi_ = -999.0;
  AK4GEN_AllResolved_Higgs_M_ = -999.0;
  AK4GEN_AllResolved_Higgs_DR_Higgs0PID_  = -999.0;
  AK4GEN_AllResolved_Higgs_DR_Higgs1PID_  = -999.0;
}


bool GenAnalyzer::reorder(const TLorentzVector &a, const TLorentzVector &b)
{
  return a.Pt() > b.Pt();
}


bool GenAnalyzer::jetCleaning(const reco::GenJet  * genAK8jet,const vector<reco::GenJet>* genAK4_coll, const double r_seperation)
{
  for(vector<reco::GenJet>::const_iterator genjet = genAK4_coll->begin(); genjet != genAK4_coll->end(); genjet++) {
    if (deltaR(genAK8jet->pt(), genAK8jet->eta(), genjet->pt(), genjet->eta()) < r_seperation) return false;
  }
  return true;
}
