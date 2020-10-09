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
  static bool jetCleaning(const reco::GenJet  * genAK8jet,const vector<reco::GenJet>* firstJetCollection, const double r_seperation=0.8);
  static void SortedCleanedJetVector(const std::vector<reco::GenJet>* jetCollection1, const std::vector<reco::GenJet>* jetCollection2, const std::vector<TLorentzVector> &photons, std::vector<TLorentzVector> &outputLorentzVector, const double r_seperation=0.8);
  static void SortedCleanedJetVectorAK8(const std::vector<reco::GenJet>* jetCollection1, const std::vector<reco::GenJet>* jetCollection2, const std::vector<TLorentzVector> &photons, std::vector<TLorentzVector> &outputLorentzVector, const double r_seperation=0.8);
  static void indexOfSelectedJet(const std::vector<TLorentzVector> &inputLorentzVector, int &index1, int &index2);
  static void indexOfSelectedJet(const std::vector<TLorentzVector> &inputLorentzVector, double massComp, int &index1, int &index2, int in_index1=-1, int in_index2=-1);
  static TLorentzVector maxPtLorentzVector(const std::vector<TLorentzVector> &inputLorentzVector);
  static TLorentzVector minMassLorentzVector(const std::vector<TLorentzVector> &inputLorentzVector, const double mass);
  static TLorentzVector minMassLorentzVector(const std::vector<TLorentzVector> &inputLorentzVector, const double mass, int &position, bool skip);
  static void minMassLorentzVector(const std::vector<TLorentzVector> &inputLorentzVector, const double mass, TLorentzVector &leadingJet, TLorentzVector &subleadingJet);
  static std::vector<TLorentzVector> minMassLorentzVector(const std::vector<TLorentzVector> &input_AK4LorentzVector, const std::vector<TLorentzVector> &input_AK8LorentzVector);

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

  std::vector<int> pdgID_;

  //  std::vector<std::string> LHEWeightIDs_;
  //  std::vector<double> LHEWeights_;

  int nEVENT=-999;

  //  int		isMuMinus_ = -999;
  //  double 	LHELeptPt_ = -999.0;
  //  double 	LHELeptEta_ = -999.0;
  //  double 	LHELeptPhi_ = -999.0;
  //  double 	LHELeptM_ = -999.0;
  //  double 	LHELeptE_ = -999.0;

  //--------------------------------------------------------------
  //  Photon variables
  double gen_leading_photon_Pt_   = -999.0;
  double gen_leading_photon_Eta_  = -999.0;
  double gen_leading_photon_Phi_  = -999.0;
  double gen_leading_photon_M_    = -999.0;
  double gen_Subleading_photon_Pt_    = -999.0;
  double gen_Subleading_photon_Eta_   = -999.0;
  double gen_Subleading_photon_Phi_   = -999.0;
  double gen_Subleading_photon_M_     = -999.0;
  double gen_HiggsGG_Pt_  = -999.0;
  double gen_HiggsGG_Eta_  = -999.0;
  double gen_HiggsGG_Phi_  = -999.0;
  double gen_HiggsGG_M_  = -999.0;

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
  double gen_HH_Pt_   = -999.0;
  double gen_HH_Eta_  = -999.0;
  double gen_HH_Phi_  = -999.0;
  double gen_HH_M_    = -999.0;

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

  double gen_deltaPhi_Photon0_Photon1_  = -999.0;
  double gen_deltaPhi_Wp_Wm_  = -999.0;
  double gen_deltaPhi_H1_H2_  = -999.0;
  double gen_deltaPhi_WpJ0_WpJ1_  = -999.0;
  double gen_deltaPhi_WmJ0_WmJ1_  = -999.0;
  double gen_deltaEta_Photon0_Photon1_  = -999.0;
  double gen_deltaEta_Wp_Wm_  = -999.0;
  double gen_deltaEta_H1_H2_  = -999.0;
  double gen_deltaEta_WpJ0_WpJ1_  = -999.0;
  double gen_deltaEta_WmJ0_WmJ1_    = -999.0;

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

  double AK8Gen_MergedWjets_MaxPt_Higgs_Pt_ = -999.0;
  double AK8Gen_MergedWjets_MaxPt_Higgs_Eta_ = -999.0;
  double AK8Gen_MergedWjets_MaxPt_Higgs_Phi_  = -999.0;
  double AK8Gen_MergedWjets_MaxPt_Higgs_M_  = -999.0;
  double AK8Gen_MergedWjets_MaxPt_Higgs_deltaR_H1_  = -999.0;
  double AK8Gen_MergedWjets_MaxPt_Higgs_deltaR_H2_  = -999.0;

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

  double AK8Gen_MergedWjets_minDMass_Higgs_Pt_  = -999.0;
  double AK8Gen_MergedWjets_minDMass_Higgs_Eta_ = -999.0;
  double AK8Gen_MergedWjets_minDMass_Higgs_Phi_ = -999.0;
  double AK8Gen_MergedWjets_minDMass_Higgs_M_ = -999.0;
  double AK8Gen_MergedWjets_minDMass_Higgs_deltaR_H1_ = -999.0;
  double AK8Gen_MergedWjets_minDMass_Higgs_deltaR_H2_ = -999.0;

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

  double AK8Gen_MergedWjets_minWminHmass_Higgs_Pt_  = -999.0;
  double AK8Gen_MergedWjets_minWminHmass_Higgs_Eta_ = -999.0;
  double AK8Gen_MergedWjets_minWminHmass_Higgs_Phi_ = -999.0;
  double AK8Gen_MergedWjets_minWminHmass_Higgs_M_ = -999.0;
  double AK8Gen_MergedWjets_minWminHmass_Higgs_deltaR_H1_ = -999.0;
  double AK8Gen_MergedWjets_minWminHmass_Higgs_deltaR_H2_ = -999.0;

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

  double AK4GEN_AllResolved_HH_Pt_  = -999.;
  double AK4GEN_AllResolved_HH_Eta_ = -999.;
  double AK4GEN_AllResolved_HH_Phi_ = -999.;
  double AK4GEN_AllResolved_HH_M_ = -999.;

  double OneAK8TwoAK4_pTMax_AK8_Pt_ = -999.0;
  double OneAK8TwoAK4_pTMax_AK8_Eta_  = -999.0;
  double OneAK8TwoAK4_pTMax_AK8_Phi_  = -999.0;
  double OneAK8TwoAK4_pTMax_AK8_M_  = -999.0;
  double OneAK8TwoAK4_pTMax_AK8_dR_W1_  = -999.0;
  double OneAK8TwoAK4_pTMax_AK8_dR_W2_  = -999.0;
  double OneAK8TwoAK4_pTMax_AK8_dR_H1_  = -999.0;
  double OneAK8TwoAK4_pTMax_AK8_dR_H2_  = -999.0;
  double OneAK8TwoAK4_pTMax_leadingAK4_Pt_  = -999.0;
  double OneAK8TwoAK4_pTMax_leadingAK4_Eta_ = -999.0;
  double OneAK8TwoAK4_pTMax_leadingAK4_Phi_ = -999.0;
  double OneAK8TwoAK4_pTMax_leadingAK4_M_ = -999.0;
  double OneAK8TwoAK4_pTMax_leadingAK4_dR_W1_ = -999.0;
  double OneAK8TwoAK4_pTMax_leadingAK4_dR_W2_ = -999.0;
  double OneAK8TwoAK4_pTMax_leadingAK4_dR_H1_ = -999.0;
  double OneAK8TwoAK4_pTMax_leadingAK4_dR_H2_ = -999.0;
  double OneAK8TwoAK4_pTMax_subleadingAK4_Pt_ = -999.0;
  double OneAK8TwoAK4_pTMax_subleadingAK4_Eta_  = -999.0;
  double OneAK8TwoAK4_pTMax_subleadingAK4_Phi_  = -999.0;
  double OneAK8TwoAK4_pTMax_subleadingAK4_M_  = -999.0;
  double OneAK8TwoAK4_pTMax_subleadingAK4_dR_W1_  = -999.0;
  double OneAK8TwoAK4_pTMax_subleadingAK4_dR_W2_  = -999.0;
  double OneAK8TwoAK4_pTMax_subleadingAK4_dR_H1_  = -999.0;
  double OneAK8TwoAK4_pTMax_subleadingAK4_dR_H2_  = -999.0;
  double OneAK8TwoAK4_pTMax_ReconsW_AK4_Pt_ = -999.0;
  double OneAK8TwoAK4_pTMax_ReconsW_AK4_Eta_  = -999.0;
  double OneAK8TwoAK4_pTMax_ReconsW_AK4_Phi_  = -999.0;
  double OneAK8TwoAK4_pTMax_ReconsW_AK4_M_  = -999.0;
  double OneAK8TwoAK4_pTMax_ReconsW_AK4_dR_W1_  = -999.0;
  double OneAK8TwoAK4_pTMax_ReconsW_AK4_dR_W2_  = -999.0;
  double OneAK8TwoAK4_pTMax_ReconsW_AK4_dR_H1_  = -999.0;
  double OneAK8TwoAK4_pTMax_ReconsW_AK4_dR_H2_  = -999.0;
  double OneAK8TwoAK4_pTMax_ReconsH_Pt_ = -999.0;
  double OneAK8TwoAK4_pTMax_ReconsH_Eta_  = -999.0;
  double OneAK8TwoAK4_pTMax_ReconsH_Phi_  = -999.0;
  double OneAK8TwoAK4_pTMax_ReconsH_M_  = -999.0;
  double OneAK8TwoAK4_pTMax_ReconsH_dR_W1_  = -999.0;
  double OneAK8TwoAK4_pTMax_ReconsH_dR_W2_  = -999.0;
  double OneAK8TwoAK4_pTMax_ReconsH_dR_H1_  = -999.0;
  double OneAK8TwoAK4_pTMax_ReconsH_dR_H2_  = -999.0;

  double OneAK8TwoAK4_minMass_AK8_Pt_ = -999.0;
  double OneAK8TwoAK4_minMass_AK8_Eta_  = -999.0;
  double OneAK8TwoAK4_minMass_AK8_Phi_  = -999.0;
  double OneAK8TwoAK4_minMass_AK8_M_  = -999.0;
  double OneAK8TwoAK4_minMass_AK8_dR_W1_  = -999.0;
  double OneAK8TwoAK4_minMass_AK8_dR_W2_  = -999.0;
  double OneAK8TwoAK4_minMass_AK8_dR_H1_  = -999.0;
  double OneAK8TwoAK4_minMass_AK8_dR_H2_  = -999.0;
  double OneAK8TwoAK4_minMass_leadingAK4_Pt_  = -999.0;
  double OneAK8TwoAK4_minMass_leadingAK4_Eta_ = -999.0;
  double OneAK8TwoAK4_minMass_leadingAK4_Phi_ = -999.0;
  double OneAK8TwoAK4_minMass_leadingAK4_M_ = -999.0;
  double OneAK8TwoAK4_minMass_leadingAK4_dR_W1_ = -999.0;
  double OneAK8TwoAK4_minMass_leadingAK4_dR_W2_ = -999.0;
  double OneAK8TwoAK4_minMass_leadingAK4_dR_H1_ = -999.0;
  double OneAK8TwoAK4_minMass_leadingAK4_dR_H2_ = -999.0;
  double OneAK8TwoAK4_minMass_subleadingAK4_Pt_ = -999.0;
  double OneAK8TwoAK4_minMass_subleadingAK4_Eta_  = -999.0;
  double OneAK8TwoAK4_minMass_subleadingAK4_Phi_  = -999.0;
  double OneAK8TwoAK4_minMass_subleadingAK4_M_  = -999.0;
  double OneAK8TwoAK4_minMass_subleadingAK4_dR_W1_  = -999.0;
  double OneAK8TwoAK4_minMass_subleadingAK4_dR_W2_  = -999.0;
  double OneAK8TwoAK4_minMass_subleadingAK4_dR_H1_  = -999.0;
  double OneAK8TwoAK4_minMass_subleadingAK4_dR_H2_  = -999.0;
  double OneAK8TwoAK4_minMass_ReconsW_AK4_Pt_ = -999.0;
  double OneAK8TwoAK4_minMass_ReconsW_AK4_Eta_  = -999.0;
  double OneAK8TwoAK4_minMass_ReconsW_AK4_Phi_  = -999.0;
  double OneAK8TwoAK4_minMass_ReconsW_AK4_M_  = -999.0;
  double OneAK8TwoAK4_minMass_ReconsW_AK4_dR_W1_  = -999.0;
  double OneAK8TwoAK4_minMass_ReconsW_AK4_dR_W2_  = -999.0;
  double OneAK8TwoAK4_minMass_ReconsW_AK4_dR_H1_  = -999.0;
  double OneAK8TwoAK4_minMass_ReconsW_AK4_dR_H2_  = -999.0;
  double OneAK8TwoAK4_minMass_ReconsH_Pt_ = -999.0;
  double OneAK8TwoAK4_minMass_ReconsH_Eta_  = -999.0;
  double OneAK8TwoAK4_minMass_ReconsH_Phi_  = -999.0;
  double OneAK8TwoAK4_minMass_ReconsH_M_  = -999.0;
  double OneAK8TwoAK4_minMass_ReconsH_dR_W1_  = -999.0;
  double OneAK8TwoAK4_minMass_ReconsH_dR_W2_  = -999.0;
  double OneAK8TwoAK4_minMass_ReconsH_dR_H1_  = -999.0;
  double OneAK8TwoAK4_minMass_ReconsH_dR_H2_  = -999.0;

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
  AddBranch(&pdgID_,	"pdgID");
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
  AddBranch(&gen_HiggsGG_Pt_,"gen_HiggsGG_Pt");
  AddBranch(&gen_HiggsGG_Eta_,"gen_HiggsGG_Eta");
  AddBranch(&gen_HiggsGG_Phi_,"gen_HiggsGG_Phi");
  AddBranch(&gen_HiggsGG_M_,"gen_HiggsGG_M");

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
  AddBranch(&gen_HH_Pt_, "gen_HH_Pt");
  AddBranch(&gen_HH_Eta_, "gen_HH_Eta");
  AddBranch(&gen_HH_Phi_, "gen_HH_Phi");
  AddBranch(&gen_HH_M_  , "gen_HH_M");
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

  AddBranch(&gen_deltaPhi_Photon0_Photon1_, "gen_deltaPhi_Photon0_Photon1");
  AddBranch(&gen_deltaPhi_Wp_Wm_, "gen_deltaPhi_Wp_Wm");
  AddBranch(&gen_deltaPhi_H1_H2_, "gen_deltaPhi_H1_H2");
  AddBranch(&gen_deltaPhi_WpJ0_WpJ1_, "gen_deltaPhi_WpJ0_WpJ1");
  AddBranch(&gen_deltaPhi_WmJ0_WmJ1_, "gen_deltaPhi_WmJ0_WmJ1");
  AddBranch(&gen_deltaEta_Photon0_Photon1_, "gen_deltaEta_Photon0_Photon1");
  AddBranch(&gen_deltaEta_Wp_Wm_, "gen_deltaEta_Wp_Wm");
  AddBranch(&gen_deltaEta_H1_H2_, "gen_deltaEta_H1_H2");
  AddBranch(&gen_deltaEta_WpJ0_WpJ1_, "gen_deltaEta_WpJ0_WpJ1");
  AddBranch(&gen_deltaEta_WmJ0_WmJ1_  , "gen_deltaEta_WmJ0_WmJ1");

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

  AddBranch(&AK8Gen_MergedWjets_MaxPt_Higgs_Pt_,"AK8Gen_MergedWjets_MaxPt_Higgs_Pt");
  AddBranch(&AK8Gen_MergedWjets_MaxPt_Higgs_Eta_,"AK8Gen_MergedWjets_MaxPt_Higgs_Eta");
  AddBranch(&AK8Gen_MergedWjets_MaxPt_Higgs_Phi_,"AK8Gen_MergedWjets_MaxPt_Higgs_Phi");
  AddBranch(&AK8Gen_MergedWjets_MaxPt_Higgs_M_,"AK8Gen_MergedWjets_MaxPt_Higgs_M");
  AddBranch(&AK8Gen_MergedWjets_MaxPt_Higgs_deltaR_H1_,"AK8Gen_MergedWjets_MaxPt_Higgs_deltaR_H1");
  AddBranch(&AK8Gen_MergedWjets_MaxPt_Higgs_deltaR_H2_,"AK8Gen_MergedWjets_MaxPt_Higgs_deltaR_H2");

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

  AddBranch(&AK8Gen_MergedWjets_minDMass_Higgs_Pt_,"AK8Gen_MergedWjets_minDMass_Higgs_Pt");
  AddBranch(&AK8Gen_MergedWjets_minDMass_Higgs_Eta_,"AK8Gen_MergedWjets_minDMass_Higgs_Eta");
  AddBranch(&AK8Gen_MergedWjets_minDMass_Higgs_Phi_,"AK8Gen_MergedWjets_minDMass_Higgs_Phi");
  AddBranch(&AK8Gen_MergedWjets_minDMass_Higgs_M_,"AK8Gen_MergedWjets_minDMass_Higgs_M");
  AddBranch(&AK8Gen_MergedWjets_minDMass_Higgs_deltaR_H1_,"AK8Gen_MergedWjets_minDMass_Higgs_deltaR_H1");
  AddBranch(&AK8Gen_MergedWjets_minDMass_Higgs_deltaR_H2_,"AK8Gen_MergedWjets_minDMass_Higgs_deltaR_H2");

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

  AddBranch(&AK8Gen_MergedWjets_minWminHmass_Higgs_Pt_,"AK8Gen_MergedWjets_minWminHmass_Higgs_Pt");
  AddBranch(&AK8Gen_MergedWjets_minWminHmass_Higgs_Eta_,"AK8Gen_MergedWjets_minWminHmass_Higgs_Eta");
  AddBranch(&AK8Gen_MergedWjets_minWminHmass_Higgs_Phi_,"AK8Gen_MergedWjets_minWminHmass_Higgs_Phi");
  AddBranch(&AK8Gen_MergedWjets_minWminHmass_Higgs_M_,"AK8Gen_MergedWjets_minWminHmass_Higgs_M");
  AddBranch(&AK8Gen_MergedWjets_minWminHmass_Higgs_deltaR_H1_,"AK8Gen_MergedWjets_minWminHmass_Higgs_deltaR_H1");
  AddBranch(&AK8Gen_MergedWjets_minWminHmass_Higgs_deltaR_H2_,"AK8Gen_MergedWjets_minWminHmass_Higgs_deltaR_H2");

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
  AddBranch(&AK4GEN_AllResolved_HH_Pt_,"AK4GEN_AllResolved_HH_Pt");
  AddBranch(&AK4GEN_AllResolved_HH_Eta_,"AK4GEN_AllResolved_HH_Eta");
  AddBranch(&AK4GEN_AllResolved_HH_Phi_,"AK4GEN_AllResolved_HH_Phi");
  AddBranch(&AK4GEN_AllResolved_HH_M_,"AK4GEN_AllResolved_HH_M");

  AddBranch(&OneAK8TwoAK4_pTMax_AK8_Pt_,"OneAK8TwoAK4_pTMax_AK8_Pt");
  AddBranch(&OneAK8TwoAK4_pTMax_AK8_Eta_,"OneAK8TwoAK4_pTMax_AK8_Eta");
  AddBranch(&OneAK8TwoAK4_pTMax_AK8_Phi_,"OneAK8TwoAK4_pTMax_AK8_Phi");
  AddBranch(&OneAK8TwoAK4_pTMax_AK8_M_,"OneAK8TwoAK4_pTMax_AK8_M");
  AddBranch(&OneAK8TwoAK4_pTMax_AK8_dR_W1_,"OneAK8TwoAK4_pTMax_AK8_dR_W1");
  AddBranch(&OneAK8TwoAK4_pTMax_AK8_dR_W2_,"OneAK8TwoAK4_pTMax_AK8_dR_W2");
  AddBranch(&OneAK8TwoAK4_pTMax_AK8_dR_H1_,"OneAK8TwoAK4_pTMax_AK8_dR_H1");
  AddBranch(&OneAK8TwoAK4_pTMax_AK8_dR_H2_,"OneAK8TwoAK4_pTMax_AK8_dR_H2");
  AddBranch(&OneAK8TwoAK4_pTMax_leadingAK4_Pt_,"OneAK8TwoAK4_pTMax_leadingAK4_Pt");
  AddBranch(&OneAK8TwoAK4_pTMax_leadingAK4_Eta_,"OneAK8TwoAK4_pTMax_leadingAK4_Eta");
  AddBranch(&OneAK8TwoAK4_pTMax_leadingAK4_Phi_,"OneAK8TwoAK4_pTMax_leadingAK4_Phi");
  AddBranch(&OneAK8TwoAK4_pTMax_leadingAK4_M_,"OneAK8TwoAK4_pTMax_leadingAK4_M");
  AddBranch(&OneAK8TwoAK4_pTMax_leadingAK4_dR_W1_,"OneAK8TwoAK4_pTMax_leadingAK4_dR_W1");
  AddBranch(&OneAK8TwoAK4_pTMax_leadingAK4_dR_W2_,"OneAK8TwoAK4_pTMax_leadingAK4_dR_W2");
  AddBranch(&OneAK8TwoAK4_pTMax_leadingAK4_dR_H1_,"OneAK8TwoAK4_pTMax_leadingAK4_dR_H1");
  AddBranch(&OneAK8TwoAK4_pTMax_leadingAK4_dR_H2_,"OneAK8TwoAK4_pTMax_leadingAK4_dR_H2");
  AddBranch(&OneAK8TwoAK4_pTMax_subleadingAK4_Pt_,"OneAK8TwoAK4_pTMax_subleadingAK4_Pt");
  AddBranch(&OneAK8TwoAK4_pTMax_subleadingAK4_Eta_,"OneAK8TwoAK4_pTMax_subleadingAK4_Eta");
  AddBranch(&OneAK8TwoAK4_pTMax_subleadingAK4_Phi_,"OneAK8TwoAK4_pTMax_subleadingAK4_Phi");
  AddBranch(&OneAK8TwoAK4_pTMax_subleadingAK4_M_,"OneAK8TwoAK4_pTMax_subleadingAK4_M");
  AddBranch(&OneAK8TwoAK4_pTMax_subleadingAK4_dR_W1_,"OneAK8TwoAK4_pTMax_subleadingAK4_dR_W1");
  AddBranch(&OneAK8TwoAK4_pTMax_subleadingAK4_dR_W2_,"OneAK8TwoAK4_pTMax_subleadingAK4_dR_W2");
  AddBranch(&OneAK8TwoAK4_pTMax_subleadingAK4_dR_H1_,"OneAK8TwoAK4_pTMax_subleadingAK4_dR_H1");
  AddBranch(&OneAK8TwoAK4_pTMax_subleadingAK4_dR_H2_,"OneAK8TwoAK4_pTMax_subleadingAK4_dR_H2");
  AddBranch(&OneAK8TwoAK4_pTMax_ReconsW_AK4_Pt_,"OneAK8TwoAK4_pTMax_ReconsW_AK4_Pt");
  AddBranch(&OneAK8TwoAK4_pTMax_ReconsW_AK4_Eta_,"OneAK8TwoAK4_pTMax_ReconsW_AK4_Eta");
  AddBranch(&OneAK8TwoAK4_pTMax_ReconsW_AK4_Phi_,"OneAK8TwoAK4_pTMax_ReconsW_AK4_Phi");
  AddBranch(&OneAK8TwoAK4_pTMax_ReconsW_AK4_M_,"OneAK8TwoAK4_pTMax_ReconsW_AK4_M");
  AddBranch(&OneAK8TwoAK4_pTMax_ReconsW_AK4_dR_W1_,"OneAK8TwoAK4_pTMax_ReconsW_AK4_dR_W1");
  AddBranch(&OneAK8TwoAK4_pTMax_ReconsW_AK4_dR_W2_,"OneAK8TwoAK4_pTMax_ReconsW_AK4_dR_W2");
  AddBranch(&OneAK8TwoAK4_pTMax_ReconsW_AK4_dR_H1_,"OneAK8TwoAK4_pTMax_ReconsW_AK4_dR_H1");
  AddBranch(&OneAK8TwoAK4_pTMax_ReconsW_AK4_dR_H2_,"OneAK8TwoAK4_pTMax_ReconsW_AK4_dR_H2");
  AddBranch(&OneAK8TwoAK4_pTMax_ReconsH_Pt_,"OneAK8TwoAK4_pTMax_ReconsH_Pt");
  AddBranch(&OneAK8TwoAK4_pTMax_ReconsH_Eta_,"OneAK8TwoAK4_pTMax_ReconsH_Eta");
  AddBranch(&OneAK8TwoAK4_pTMax_ReconsH_Phi_,"OneAK8TwoAK4_pTMax_ReconsH_Phi");
  AddBranch(&OneAK8TwoAK4_pTMax_ReconsH_M_,"OneAK8TwoAK4_pTMax_ReconsH_M");
  AddBranch(&OneAK8TwoAK4_pTMax_ReconsH_dR_W1_,"OneAK8TwoAK4_pTMax_ReconsH_dR_W1");
  AddBranch(&OneAK8TwoAK4_pTMax_ReconsH_dR_W2_,"OneAK8TwoAK4_pTMax_ReconsH_dR_W2");
  AddBranch(&OneAK8TwoAK4_pTMax_ReconsH_dR_H1_,"OneAK8TwoAK4_pTMax_ReconsH_dR_H1");
  AddBranch(&OneAK8TwoAK4_pTMax_ReconsH_dR_H2_,"OneAK8TwoAK4_pTMax_ReconsH_dR_H2");

  AddBranch(&OneAK8TwoAK4_minMass_AK8_Pt_,"OneAK8TwoAK4_minMass_AK8_Pt");
  AddBranch(&OneAK8TwoAK4_minMass_AK8_Eta_,"OneAK8TwoAK4_minMass_AK8_Eta");
  AddBranch(&OneAK8TwoAK4_minMass_AK8_Phi_,"OneAK8TwoAK4_minMass_AK8_Phi");
  AddBranch(&OneAK8TwoAK4_minMass_AK8_M_,"OneAK8TwoAK4_minMass_AK8_M");
  AddBranch(&OneAK8TwoAK4_minMass_AK8_dR_W1_,"OneAK8TwoAK4_minMass_AK8_dR_W1");
  AddBranch(&OneAK8TwoAK4_minMass_AK8_dR_W2_,"OneAK8TwoAK4_minMass_AK8_dR_W2");
  AddBranch(&OneAK8TwoAK4_minMass_AK8_dR_H1_,"OneAK8TwoAK4_minMass_AK8_dR_H1");
  AddBranch(&OneAK8TwoAK4_minMass_AK8_dR_H2_,"OneAK8TwoAK4_minMass_AK8_dR_H2");
  AddBranch(&OneAK8TwoAK4_minMass_leadingAK4_Pt_,"OneAK8TwoAK4_minMass_leadingAK4_Pt");
  AddBranch(&OneAK8TwoAK4_minMass_leadingAK4_Eta_,"OneAK8TwoAK4_minMass_leadingAK4_Eta");
  AddBranch(&OneAK8TwoAK4_minMass_leadingAK4_Phi_,"OneAK8TwoAK4_minMass_leadingAK4_Phi");
  AddBranch(&OneAK8TwoAK4_minMass_leadingAK4_M_,"OneAK8TwoAK4_minMass_leadingAK4_M");
  AddBranch(&OneAK8TwoAK4_minMass_leadingAK4_dR_W1_,"OneAK8TwoAK4_minMass_leadingAK4_dR_W1");
  AddBranch(&OneAK8TwoAK4_minMass_leadingAK4_dR_W2_,"OneAK8TwoAK4_minMass_leadingAK4_dR_W2");
  AddBranch(&OneAK8TwoAK4_minMass_leadingAK4_dR_H1_,"OneAK8TwoAK4_minMass_leadingAK4_dR_H1");
  AddBranch(&OneAK8TwoAK4_minMass_leadingAK4_dR_H2_,"OneAK8TwoAK4_minMass_leadingAK4_dR_H2");
  AddBranch(&OneAK8TwoAK4_minMass_subleadingAK4_Pt_,"OneAK8TwoAK4_minMass_subleadingAK4_Pt");
  AddBranch(&OneAK8TwoAK4_minMass_subleadingAK4_Eta_,"OneAK8TwoAK4_minMass_subleadingAK4_Eta");
  AddBranch(&OneAK8TwoAK4_minMass_subleadingAK4_Phi_,"OneAK8TwoAK4_minMass_subleadingAK4_Phi");
  AddBranch(&OneAK8TwoAK4_minMass_subleadingAK4_M_,"OneAK8TwoAK4_minMass_subleadingAK4_M");
  AddBranch(&OneAK8TwoAK4_minMass_subleadingAK4_dR_W1_,"OneAK8TwoAK4_minMass_subleadingAK4_dR_W1");
  AddBranch(&OneAK8TwoAK4_minMass_subleadingAK4_dR_W2_,"OneAK8TwoAK4_minMass_subleadingAK4_dR_W2");
  AddBranch(&OneAK8TwoAK4_minMass_subleadingAK4_dR_H1_,"OneAK8TwoAK4_minMass_subleadingAK4_dR_H1");
  AddBranch(&OneAK8TwoAK4_minMass_subleadingAK4_dR_H2_,"OneAK8TwoAK4_minMass_subleadingAK4_dR_H2");
  AddBranch(&OneAK8TwoAK4_minMass_ReconsW_AK4_Pt_,"OneAK8TwoAK4_minMass_ReconsW_AK4_Pt");
  AddBranch(&OneAK8TwoAK4_minMass_ReconsW_AK4_Eta_,"OneAK8TwoAK4_minMass_ReconsW_AK4_Eta");
  AddBranch(&OneAK8TwoAK4_minMass_ReconsW_AK4_Phi_,"OneAK8TwoAK4_minMass_ReconsW_AK4_Phi");
  AddBranch(&OneAK8TwoAK4_minMass_ReconsW_AK4_M_,"OneAK8TwoAK4_minMass_ReconsW_AK4_M");
  AddBranch(&OneAK8TwoAK4_minMass_ReconsW_AK4_dR_W1_,"OneAK8TwoAK4_minMass_ReconsW_AK4_dR_W1");
  AddBranch(&OneAK8TwoAK4_minMass_ReconsW_AK4_dR_W2_,"OneAK8TwoAK4_minMass_ReconsW_AK4_dR_W2");
  AddBranch(&OneAK8TwoAK4_minMass_ReconsW_AK4_dR_H1_,"OneAK8TwoAK4_minMass_ReconsW_AK4_dR_H1");
  AddBranch(&OneAK8TwoAK4_minMass_ReconsW_AK4_dR_H2_,"OneAK8TwoAK4_minMass_ReconsW_AK4_dR_H2");
  AddBranch(&OneAK8TwoAK4_minMass_ReconsH_Pt_,"OneAK8TwoAK4_minMass_ReconsH_Pt");
  AddBranch(&OneAK8TwoAK4_minMass_ReconsH_Eta_,"OneAK8TwoAK4_minMass_ReconsH_Eta");
  AddBranch(&OneAK8TwoAK4_minMass_ReconsH_Phi_,"OneAK8TwoAK4_minMass_ReconsH_Phi");
  AddBranch(&OneAK8TwoAK4_minMass_ReconsH_M_,"OneAK8TwoAK4_minMass_ReconsH_M");
  AddBranch(&OneAK8TwoAK4_minMass_ReconsH_dR_W1_,"OneAK8TwoAK4_minMass_ReconsH_dR_W1");
  AddBranch(&OneAK8TwoAK4_minMass_ReconsH_dR_W2_,"OneAK8TwoAK4_minMass_ReconsH_dR_W2");
  AddBranch(&OneAK8TwoAK4_minMass_ReconsH_dR_H1_,"OneAK8TwoAK4_minMass_ReconsH_dR_H1");
  AddBranch(&OneAK8TwoAK4_minMass_ReconsH_dR_H2_,"OneAK8TwoAK4_minMass_ReconsH_dR_H2");

}

void GenAnalyzer::Clear(){
  pdgID_.clear();
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
  gen_HiggsGG_Pt_ = -999.0;
  gen_HiggsGG_Eta_ = -999.0;
  gen_HiggsGG_Phi_ = -999.0;
  gen_HiggsGG_M_ = -999.0;

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
  gen_HH_Pt_ = -999.0;
  gen_HH_Eta_ = -999.0;
  gen_HH_Phi_ = -999.0;
  gen_HH_M_   = -999.0;
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

  gen_deltaPhi_Photon0_Photon1_ = -999.0;
  gen_deltaPhi_Wp_Wm_ = -999.0;
  gen_deltaPhi_H1_H2_ = -999.0;
  gen_deltaPhi_WpJ0_WpJ1_ = -999.0;
  gen_deltaPhi_WmJ0_WmJ1_ = -999.0;
  gen_deltaEta_Photon0_Photon1_ = -999.0;
  gen_deltaEta_Wp_Wm_ = -999.0;
  gen_deltaEta_H1_H2_ = -999.0;
  gen_deltaEta_WpJ0_WpJ1_ = -999.0;
  gen_deltaEta_WmJ0_WmJ1_   = -999.0;

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

  AK8Gen_MergedWjets_MaxPt_Higgs_Pt_  = -999.0;
  AK8Gen_MergedWjets_MaxPt_Higgs_Eta_ = -999.0;
  AK8Gen_MergedWjets_MaxPt_Higgs_Phi_ = -999.0;
  AK8Gen_MergedWjets_MaxPt_Higgs_M_ = -999.0;
  AK8Gen_MergedWjets_MaxPt_Higgs_deltaR_H1_ = -999.0;
  AK8Gen_MergedWjets_MaxPt_Higgs_deltaR_H2_ = -999.0;

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

  AK8Gen_MergedWjets_minDMass_Higgs_Pt_ = -999.0;
  AK8Gen_MergedWjets_minDMass_Higgs_Eta_  = -999.0;
  AK8Gen_MergedWjets_minDMass_Higgs_Phi_  = -999.0;
  AK8Gen_MergedWjets_minDMass_Higgs_M_  = -999.0;
  AK8Gen_MergedWjets_minDMass_Higgs_deltaR_H1_  = -999.0;
  AK8Gen_MergedWjets_minDMass_Higgs_deltaR_H2_  = -999.0;

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

  AK8Gen_MergedWjets_minWminHmass_Higgs_Pt_ = -999.0;
  AK8Gen_MergedWjets_minWminHmass_Higgs_Eta_  = -999.0;
  AK8Gen_MergedWjets_minWminHmass_Higgs_Phi_  = -999.0;
  AK8Gen_MergedWjets_minWminHmass_Higgs_M_  = -999.0;
  AK8Gen_MergedWjets_minWminHmass_Higgs_deltaR_H1_  = -999.0;
  AK8Gen_MergedWjets_minWminHmass_Higgs_deltaR_H2_  = -999.0;

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
  AK4GEN_AllResolved_HH_Pt_ = -999.0;
  AK4GEN_AllResolved_HH_Eta_ = -999.0;
  AK4GEN_AllResolved_HH_Phi_ = -999.0;
  AK4GEN_AllResolved_HH_M_ = -999.0;

  OneAK8TwoAK4_pTMax_AK8_Pt_  = -999.0;
  OneAK8TwoAK4_pTMax_AK8_Eta_ = -999.0;
  OneAK8TwoAK4_pTMax_AK8_Phi_ = -999.0;
  OneAK8TwoAK4_pTMax_AK8_M_ = -999.0;
  OneAK8TwoAK4_pTMax_AK8_dR_W1_ = -999.0;
  OneAK8TwoAK4_pTMax_AK8_dR_W2_ = -999.0;
  OneAK8TwoAK4_pTMax_AK8_dR_H1_ = -999.0;
  OneAK8TwoAK4_pTMax_AK8_dR_H2_ = -999.0;
  OneAK8TwoAK4_pTMax_leadingAK4_Pt_ = -999.0;
  OneAK8TwoAK4_pTMax_leadingAK4_Eta_  = -999.0;
  OneAK8TwoAK4_pTMax_leadingAK4_Phi_  = -999.0;
  OneAK8TwoAK4_pTMax_leadingAK4_M_  = -999.0;
  OneAK8TwoAK4_pTMax_leadingAK4_dR_W1_  = -999.0;
  OneAK8TwoAK4_pTMax_leadingAK4_dR_W2_  = -999.0;
  OneAK8TwoAK4_pTMax_leadingAK4_dR_H1_  = -999.0;
  OneAK8TwoAK4_pTMax_leadingAK4_dR_H2_  = -999.0;
  OneAK8TwoAK4_pTMax_subleadingAK4_Pt_  = -999.0;
  OneAK8TwoAK4_pTMax_subleadingAK4_Eta_ = -999.0;
  OneAK8TwoAK4_pTMax_subleadingAK4_Phi_ = -999.0;
  OneAK8TwoAK4_pTMax_subleadingAK4_M_ = -999.0;
  OneAK8TwoAK4_pTMax_subleadingAK4_dR_W1_ = -999.0;
  OneAK8TwoAK4_pTMax_subleadingAK4_dR_W2_ = -999.0;
  OneAK8TwoAK4_pTMax_subleadingAK4_dR_H1_ = -999.0;
  OneAK8TwoAK4_pTMax_subleadingAK4_dR_H2_ = -999.0;
  OneAK8TwoAK4_pTMax_ReconsW_AK4_Pt_  = -999.0;
  OneAK8TwoAK4_pTMax_ReconsW_AK4_Eta_ = -999.0;
  OneAK8TwoAK4_pTMax_ReconsW_AK4_Phi_ = -999.0;
  OneAK8TwoAK4_pTMax_ReconsW_AK4_M_ = -999.0;
  OneAK8TwoAK4_pTMax_ReconsW_AK4_dR_W1_ = -999.0;
  OneAK8TwoAK4_pTMax_ReconsW_AK4_dR_W2_ = -999.0;
  OneAK8TwoAK4_pTMax_ReconsW_AK4_dR_H1_ = -999.0;
  OneAK8TwoAK4_pTMax_ReconsW_AK4_dR_H2_ = -999.0;
  OneAK8TwoAK4_pTMax_ReconsH_Pt_  = -999.0;
  OneAK8TwoAK4_pTMax_ReconsH_Eta_ = -999.0;
  OneAK8TwoAK4_pTMax_ReconsH_Phi_ = -999.0;
  OneAK8TwoAK4_pTMax_ReconsH_M_ = -999.0;
  OneAK8TwoAK4_pTMax_ReconsH_dR_W1_ = -999.0;
  OneAK8TwoAK4_pTMax_ReconsH_dR_W2_ = -999.0;
  OneAK8TwoAK4_pTMax_ReconsH_dR_H1_ = -999.0;
  OneAK8TwoAK4_pTMax_ReconsH_dR_H2_   = -999.0;

  OneAK8TwoAK4_minMass_AK8_Pt_  = -999.0;
  OneAK8TwoAK4_minMass_AK8_Eta_ = -999.0;
  OneAK8TwoAK4_minMass_AK8_Phi_ = -999.0;
  OneAK8TwoAK4_minMass_AK8_M_ = -999.0;
  OneAK8TwoAK4_minMass_AK8_dR_W1_ = -999.0;
  OneAK8TwoAK4_minMass_AK8_dR_W2_ = -999.0;
  OneAK8TwoAK4_minMass_AK8_dR_H1_ = -999.0;
  OneAK8TwoAK4_minMass_AK8_dR_H2_ = -999.0;
  OneAK8TwoAK4_minMass_leadingAK4_Pt_ = -999.0;
  OneAK8TwoAK4_minMass_leadingAK4_Eta_  = -999.0;
  OneAK8TwoAK4_minMass_leadingAK4_Phi_  = -999.0;
  OneAK8TwoAK4_minMass_leadingAK4_M_  = -999.0;
  OneAK8TwoAK4_minMass_leadingAK4_dR_W1_  = -999.0;
  OneAK8TwoAK4_minMass_leadingAK4_dR_W2_  = -999.0;
  OneAK8TwoAK4_minMass_leadingAK4_dR_H1_  = -999.0;
  OneAK8TwoAK4_minMass_leadingAK4_dR_H2_  = -999.0;
  OneAK8TwoAK4_minMass_subleadingAK4_Pt_  = -999.0;
  OneAK8TwoAK4_minMass_subleadingAK4_Eta_ = -999.0;
  OneAK8TwoAK4_minMass_subleadingAK4_Phi_ = -999.0;
  OneAK8TwoAK4_minMass_subleadingAK4_M_ = -999.0;
  OneAK8TwoAK4_minMass_subleadingAK4_dR_W1_ = -999.0;
  OneAK8TwoAK4_minMass_subleadingAK4_dR_W2_ = -999.0;
  OneAK8TwoAK4_minMass_subleadingAK4_dR_H1_ = -999.0;
  OneAK8TwoAK4_minMass_subleadingAK4_dR_H2_ = -999.0;
  OneAK8TwoAK4_minMass_ReconsW_AK4_Pt_  = -999.0;
  OneAK8TwoAK4_minMass_ReconsW_AK4_Eta_ = -999.0;
  OneAK8TwoAK4_minMass_ReconsW_AK4_Phi_ = -999.0;
  OneAK8TwoAK4_minMass_ReconsW_AK4_M_ = -999.0;
  OneAK8TwoAK4_minMass_ReconsW_AK4_dR_W1_ = -999.0;
  OneAK8TwoAK4_minMass_ReconsW_AK4_dR_W2_ = -999.0;
  OneAK8TwoAK4_minMass_ReconsW_AK4_dR_H1_ = -999.0;
  OneAK8TwoAK4_minMass_ReconsW_AK4_dR_H2_ = -999.0;
  OneAK8TwoAK4_minMass_ReconsH_Pt_  = -999.0;
  OneAK8TwoAK4_minMass_ReconsH_Eta_ = -999.0;
  OneAK8TwoAK4_minMass_ReconsH_Phi_ = -999.0;
  OneAK8TwoAK4_minMass_ReconsH_M_ = -999.0;
  OneAK8TwoAK4_minMass_ReconsH_dR_W1_ = -999.0;
  OneAK8TwoAK4_minMass_ReconsH_dR_W2_ = -999.0;
  OneAK8TwoAK4_minMass_ReconsH_dR_H1_ = -999.0;
  OneAK8TwoAK4_minMass_ReconsH_dR_H2_ = -999.0;

}

/**
 * This helps to identify the higher pT lorentzVector.
 * @param  a first lorentz vector
 * @param  b second lorentz vector
 * @return   [bool] True if first LV has higher pT than second otherwise False.
 */
bool GenAnalyzer::reorder(const TLorentzVector &a, const TLorentzVector &b)
{
  return a.Pt() > b.Pt();
}

/**
 * This takes one AK4 (AK8) jet and checks from AK8 (AK4) jet collection if its cleaned or not.
 * @param  genAK8jet          [reco::GenJet] one gen jet
 * @param  firstJetCollection [const vector<reco::GenJet>] vector of genJet collection
 * @param  r_seperation       [double] this is the value of delta R
 * @return                    [bool] returns true if jet is cleaned else false.
 */
bool GenAnalyzer::jetCleaning(const reco::GenJet  * genAK8jet,const vector<reco::GenJet>* firstJetCollection, const double r_seperation)
{
  for(vector<reco::GenJet>::const_iterator genjet = firstJetCollection->begin(); genjet != firstJetCollection->end(); genjet++) {
    if (deltaR(genAK8jet->pt(), genAK8jet->eta(), genjet->pt(), genjet->eta()) < r_seperation) return false;
  }
  return true;
}

/**
 * This function returns the sorted TLorentzVector that contains information from first
 * passed jetCollection.
 * @param firstJetCollection         first genjetCollection whose information will be passed to
 *                            TLorentzVector local_Vec_genJetAK5
 * @param secondJetCollection         Another genjetCollection which from which the first jet
 *                            collection is going to be cleaned.
 * @param Vec_Photons         Sorted TLorentzVector that contains the photon collection.
 * @param local_Vec_genJetAK4 vector of TLorentzVector that we need
 * @param r_seperation        This is the delta-R seperation distance beetween the two
 *                            jets.
 */
void GenAnalyzer::SortedCleanedJetVector(const std::vector<reco::GenJet>* firstJetCollection, const std::vector<reco::GenJet>* secondJetCollection, const std::vector<TLorentzVector> &Vec_Photons, std::vector<TLorentzVector> &local_Vec_genJetAK4, const double r_seperation)
{
  // int nAK4jets = 0;
  TLorentzVector temp_genJetAK4;
  for(vector<reco::GenJet>::const_iterator genjet = firstJetCollection->begin(); genjet != firstJetCollection->end(); genjet++)
  {
    if (genjet->pt()<15) continue;
    if (deltaR(genjet->eta(),genjet->phi(), Vec_Photons[0].Eta(),Vec_Photons[0].Phi())>0.4 && deltaR(genjet->eta(),genjet->phi(), Vec_Photons[1].Eta(),Vec_Photons[1].Phi())>0.4) 
    {
      // if ( GenAnalyzer::jetCleaning(&(*genjet), secondJetCollection, r_seperation ) )
      // {
        temp_genJetAK4.SetPtEtaPhiE(genjet->pt(), genjet->eta(), genjet->phi(), genjet->energy());
        local_Vec_genJetAK4.push_back(temp_genJetAK4);
        // std::cout << "**> " << genjet->pt() << "\t" << genjet->eta() << "\t" << genjet->phi() << "\t" << genjet->energy()<< std::endl;
      // }
    }
  }
  std::sort(local_Vec_genJetAK4.begin(), local_Vec_genJetAK4.end(), GenAnalyzer::reorder);
}

void GenAnalyzer::SortedCleanedJetVector(const std::vector<reco::GenJet>* firstJetCollection, const std::vector<reco::GenJet>* secondJetCollection, const std::vector<TLorentzVector> &Vec_Photons, std::vector<TLorentzVector> &local_Vec_genJetAK4, const double r_seperation)
{
  // int nAK4jets = 0;
  TLorentzVector temp_genJetAK8;
  for(vector<reco::GenJet>::const_iterator genjet = firstJetCollection->begin(); genjet != firstJetCollection->end(); genjet++)
  {
    if (genjet->pt()<30) continue;
    if (deltaR(genjet->eta(),genjet->phi(), Vec_Photons[0].Eta(),Vec_Photons[0].Phi())>0.8 && deltaR(genjet->eta(),genjet->phi(), Vec_Photons[1].Eta(),Vec_Photons[1].Phi())>0.8) 
    {
      // if ( GenAnalyzer::jetCleaning(&(*genjet), secondJetCollection, r_seperation ) )
      // {
        temp_genJetAK8.SetPtEtaPhiE(genjet->pt(), genjet->eta(), genjet->phi(), genjet->energy());
        local_Vec_genJetAK8.push_back(temp_genJetAK8);
        // std::cout << "**> " << genjet->pt() << "\t" << genjet->eta() << "\t" << genjet->phi() << "\t" << genjet->energy()<< std::endl;
      // }
    }
  }
  std::sort(local_Vec_genJetAK8.begin(), local_Vec_genJetAK8.end(), GenAnalyzer::reorder);
}

/**
 * This takes a vector of TLorentzVector and give a TLorentzVector having highest pT
 * @param  inputLorentzVector Vector of TLorentzVector
 * @return                    TLorentzVector having highest pT
 */
TLorentzVector GenAnalyzer::maxPtLorentzVector(const std::vector<TLorentzVector> &inputLorentzVector) {
  double temp_pt = -999.0;
  TLorentzVector AK8Gen_HiggsJet_MaxPt;
  for (std::vector<TLorentzVector>::const_iterator jet = inputLorentzVector.begin(); jet != inputLorentzVector.end(); ++jet)
  {
    if (jet->Pt()>temp_pt)
    {
          AK8Gen_HiggsJet_MaxPt.SetPtEtaPhiE(jet->Pt(), jet->Eta(), jet->Phi(), jet->Energy());
          temp_pt = jet->Pt();
    }
  }
  return AK8Gen_HiggsJet_MaxPt;
}

/**
 * This takes a vector of TLorentzVector and gives a TLorentzVector having minimum mass difference from the specified mass.
 * @param  inputLorentzVector [TLorentzVector] Input vector of TLorentzVector
 * @param  mass               [double] mass from which it will compare.
 * @return                    TLorentzVector having minimum mass difference w.r.t. specified mass.
 */
TLorentzVector GenAnalyzer::minMassLorentzVector(const std::vector<TLorentzVector> &inputLorentzVector, const double mass) {
  double temp_AK8jet_deltaM = 9999.0;
  TLorentzVector AK8Gen_HiggsJet_MaxPt;
  for (std::vector<TLorentzVector>::const_iterator jet = inputLorentzVector.begin(); jet != inputLorentzVector.end(); ++jet)
  {
    if (abs(jet->M()-mass)<temp_AK8jet_deltaM)
    {
          AK8Gen_HiggsJet_MaxPt.SetPtEtaPhiE(jet->Pt(), jet->Eta(), jet->Phi(), jet->Energy());
          temp_AK8jet_deltaM = abs(jet->M()-mass);
    }
  }
  return AK8Gen_HiggsJet_MaxPt;
}

/**
 * This takes a vector of TLorentzVector and returns a TLorentzVector having minimum delta mass w.r.t. the
 * specified mass. Also, it takes a input positional arguments. It just skip that while checking.
 * @param  inputLorentzVector Input vector of TLorentzVector.
 * @param  mass               mass from which it will compare
 * @param  position           position of the vector that need to be skip from calculation
 * @param  skip               [description]
 * @return                    TLorentzVector having minimum mass difference w.r.t. specified mass.
 */
TLorentzVector GenAnalyzer::minMassLorentzVector(const std::vector<TLorentzVector> &inputLorentzVector, const double mass, int &position, bool skip) {
  double temp_AK8jet_deltaM = 9999.0;
  TLorentzVector AK8Gen_HiggsJet_MaxPt;
  int counter = -1;
  // std::cout << "===================" << position << std::endl;
  for (std::vector<TLorentzVector>::const_iterator jet = inputLorentzVector.begin(); jet != inputLorentzVector.end(); ++jet)
  {
    counter++;
    // std::cout << "==> " << counter << "\t" << position << std::endl;
    if (skip)
    {
      if (counter == position) continue;
    }
    // std::cout << counter << "\t" << position << "\t" << jet->Pt() <<  "\t" << abs(jet->M()-mass) << "\t" << temp_AK8jet_deltaM<< std::endl;
    if (abs(jet->M()-mass)<temp_AK8jet_deltaM)
    {
      AK8Gen_HiggsJet_MaxPt.SetPtEtaPhiE(jet->Pt(), jet->Eta(), jet->Phi(), jet->Energy());
      temp_AK8jet_deltaM = abs(jet->M()-mass);
      if (!skip)
      {
        position = counter;
        // std::cout << "position : " << position << std::endl;
      }
    }
  }
  return AK8Gen_HiggsJet_MaxPt;
}

/**
 * [GenAnalyzer::minMassLorentzVector description]
 * @param inputLorentzVector [description]
 * @param mass               [description]
 * @param leadingJet         [description]
 * @param subleadingJet      [description]
 */
void GenAnalyzer::minMassLorentzVector(const std::vector<TLorentzVector> &inputLorentzVector, const double mass, TLorentzVector &leadingJet, TLorentzVector &subleadingJet) {
  double temp_AK8jet_deltaM = 9999.0;
  for (std::vector<TLorentzVector>::const_iterator jet = inputLorentzVector.begin(); jet != inputLorentzVector.end()-1; ++jet)
    for (std::vector<TLorentzVector>::const_iterator jet1 = jet+1; jet1 != inputLorentzVector.end(); ++jet1)
    {
      if (abs((*jet+*jet1).M()-mass)<temp_AK8jet_deltaM)
      {
        temp_AK8jet_deltaM = abs((*jet+*jet1).M() - mass);
        if (jet->Pt()>jet1->Pt()) {
          leadingJet.SetPtEtaPhiE((jet)->Pt(), (jet)->Eta(), (jet)->Phi(), (jet)->Energy());
          subleadingJet.SetPtEtaPhiE((jet1)->Pt(), (jet1)->Eta(), (jet1)->Phi(), (jet1)->Energy());
        } else {
          leadingJet.SetPtEtaPhiE((jet1)->Pt(), (jet1)->Eta(), (jet1)->Phi(), (jet1)->Energy());
          subleadingJet.SetPtEtaPhiE((jet)->Pt(), (jet)->Eta(), (jet)->Phi(), (jet)->Energy());
        }
      }
    }
}


/**
 * Returns the two index from jet collection that satisfies condition abs(mass-80) minimum.
 * @param inputLorentzVector Input TLorentzVector that contains AK8 jet information
 * @param index1             Index of first on-shell selected jet
 * @param index2             Index of second on-shell selected jet
 */
void GenAnalyzer::indexOfSelectedJet(const std::vector<TLorentzVector> &inputLorentzVector, int &index1, int &index2) {
  double tempMass1 = 9999.0;
  for ( int ak4_jet1 = 0; ak4_jet1 < int(inputLorentzVector.size())-1; ++ak4_jet1)
  {
    for ( int ak4_jet2 = ak4_jet1+1; ak4_jet2 < int(inputLorentzVector.size()); ++ak4_jet2)
    {
      double mass = (inputLorentzVector[ak4_jet1] + inputLorentzVector[ak4_jet2]).M();
      if (abs(mass - 80.0) < tempMass1)
      {
        tempMass1 = abs(mass - 80.0);
        index1 = ak4_jet1;
        index2 = ak4_jet2;
      }
    }
  }
}

/**
 * Return the two index from the jet collection vector of type TLorentzVector.
 * @param inputLorentzVector Input TLorentzVector that contains jet information
 * @param massComp           Value of mass for which the minimum delta mass is going to be calculated.
 * @param index1             Index of first selected jet.
 * @param index2             Index of second selected jet.
 * @param in_index1          Index of fist jet to be ignored.
 * @param in_index2          Index of second jet to be ignored.
 */
void GenAnalyzer::indexOfSelectedJet(const std::vector<TLorentzVector> &inputLorentzVector, double massComp, int &index1, int &index2, int in_index1, int in_index2) {
  double tempMass1 = 9999.0;
  for ( int ak4_jet1 = 0; ak4_jet1 < int(inputLorentzVector.size())-1; ++ak4_jet1)
  {
    double mass = 0.0;
    if (ak4_jet1 == in_index1) continue;
    if (ak4_jet1 == in_index2) continue;
    for ( int ak4_jet2 = ak4_jet1+1; ak4_jet2 < int(inputLorentzVector.size()); ++ak4_jet2)
    {
      if (ak4_jet2 == in_index1) continue;
      if (ak4_jet2 == in_index2) continue;
      if (in_index1 == -1 && in_index2 == -1) {
          mass = (inputLorentzVector[ak4_jet1] + inputLorentzVector[ak4_jet2]).M();
      } else {
          mass = (inputLorentzVector[ak4_jet1] + inputLorentzVector[ak4_jet2] +
          inputLorentzVector[in_index1] + inputLorentzVector[in_index2]).M();
      }
      if (abs(mass - massComp) < tempMass1)
      {
        tempMass1 = abs(mass - massComp);
        index1 = ak4_jet1;
        index2 = ak4_jet2;
      }
    }
  }
}


/**
 * @brief      Selection of two AK4 jet and one AK8 jet. This selection is done
 *             by searching the combination of two AK4 and one AK8 which gives
 *             us the minimum delta mass w.r.t. the Higgs mass.
 *
 * @param      input_AK4LorentzVector  The input vector of AK4 lorentz vector
 * @param      input_AK8LorentzVector  The input vector of AK8 lorentz vector
 *
 * @return     Returns the vector of TLorentzVector having size 3. Whose first two elements
 *             contains information about the AK4 jets and the third element contains
 *             information of AK8 jet.
 */
std::vector<TLorentzVector> GenAnalyzer::minMassLorentzVector(const std::vector<TLorentzVector> &input_AK4LorentzVector, const std::vector<TLorentzVector> &input_AK8LorentzVector)
{
  std::vector<TLorentzVector> outputLorentzVector;
  /**
   * Fill the outputLorentzVector with three blank TLorentzVector
   * So, that if it did not get any element it won't crash.
   */
  outputLorentzVector.push_back(TLorentzVector(0,0,0,0));
  outputLorentzVector.push_back(TLorentzVector(0,0,0,0));
  outputLorentzVector.push_back(TLorentzVector(0,0,0,0));

  double tempMass1 = 9999.0;

  for (std::vector<TLorentzVector>::const_iterator i = input_AK4LorentzVector.begin(); i != input_AK4LorentzVector.end()-1; ++i)
  {
    for (std::vector<TLorentzVector>::const_iterator j = i+1; j != input_AK4LorentzVector.end(); ++j)
    {
      for (std::vector<TLorentzVector>::const_iterator k = input_AK8LorentzVector.begin(); k != input_AK8LorentzVector.end(); ++k)
      {
        double mass = (*i + *j + *k).M();
        if ((*i).Pt()<15 || (*j).Pt()<15) continue;
        if ((*k).Pt()<30) continue;
        if (abs(mass - 125.0) < tempMass1)
        {
          outputLorentzVector.clear();
          outputLorentzVector.push_back(*i);
          outputLorentzVector.push_back(*j);
          outputLorentzVector.push_back(*k);
          tempMass1 = abs(mass - 125.0);
        }
      }
    }
  }
  if (outputLorentzVector.size()!=3)
  {
    std::cout << "size of output vector seems not equal to three... please check the code. " << outputLorentzVector.size() << std::endl;
    exit(EXIT_FAILURE);
  }
  return outputLorentzVector;
}
