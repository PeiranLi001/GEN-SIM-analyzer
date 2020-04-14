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
  edm::EDGetTokenT<reco::GenMETCollection> genMetCaloToken;
  
  bool	Verbose_;
  
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
  Verbose_ = iConfig.getParameter<bool>("Verbose");
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
}


bool GenAnalyzer::reorder(const TLorentzVector &a, const TLorentzVector &b)
{
  return a.Pt() > b.Pt();
}
