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

  std::vector<std::string> LHEWeightIDs_;
  std::vector<double> LHEWeights_;

int nEVENT=-999;

  int		isMuMinus_ = -999;
  double 	LHELeptPt_ = -999.0;
  double 	LHELeptEta_ = -999.0;
  double 	LHELeptPhi_ = -999.0;
  double 	LHELeptM_ = -999.0;
  double 	LHELeptE_ = -999.0;

 double LHENuPt_ = -999.0;
 double LHENuEta_ = -999.0;
 double LHENuPhi_ = -999.0;
 double LHENuM_ = -999.0;
 double LHENuE_ = -999.0;

 double LHE_DeltaM_Wqrk0_pt_ = -999.0;
 double LHE_DeltaM_Wqrk0_eta_ = -999.0;
 double LHE_DeltaM_Wqrk0_phi_ = -999.0;
 double LHE_DeltaM_Wqrk0_M_ = -999.0;
 double LHE_DeltaM_Wqrk0_E_ = -999.0;
 double LHE_DeltaM_Wqrk0_Mt_ = -999.0;

 double LHE_DeltaM_Wqrk1_pt_ = -999.0;
 double LHE_DeltaM_Wqrk1_eta_ = -999.0;
 double LHE_DeltaM_Wqrk1_phi_ = -999.0;
 double LHE_DeltaM_Wqrk1_M_ = -999.0;
 double LHE_DeltaM_Wqrk1_E_ = -999.0;
 double LHE_DeltaM_Wqrk1_Mt_ = -999.0;

 double LHE_DeltaM_Iqrk0_pt_ = -999.0; 
 double LHE_DeltaM_Iqrk0_eta_ = -999.0;
 double LHE_DeltaM_Iqrk0_phi_ = -999.0;
 double LHE_DeltaM_Iqrk0_E_ = -999.0;
 double LHE_DeltaM_Iqrk0_M_ = -999.0;
 double LHE_DeltaM_Iqrk0_Mt_ = -999.0;

 double LHE_DeltaM_Iqrk1_pt_ = -999.0; 
 double LHE_DeltaM_Iqrk1_eta_ = -999.0;
 double LHE_DeltaM_Iqrk1_phi_ = -999.0;
 double LHE_DeltaM_Iqrk1_E_ = -999.0;
 double LHE_DeltaM_Iqrk1_M_ = -999.0;
 double LHE_DeltaM_Iqrk1_Mt_ = -999.0;

 double  LHE_DeltaM_mWW_ = -999.0;
 double  LHE_DeltaM_mtWW_ = -999.0;
 double  LHE_DeltaM_mWLep_ = -999.0;
 double  LHE_DeltaM_mtWLep_ = -999.0;
 double  LHE_DeltaM_mWHad_ = -999.0;
 double  LHE_DeltaM_mtWHad_ = -999.0;
 double  LHE_DeltaM_costheta1_ = -999.0;
 double  LHE_DeltaM_costheta2_ = -999.0;
 double  LHE_DeltaM_phi_ = -999.0;
 double  LHE_DeltaM_costhetastar_ = -999.0;
 double  LHE_DeltaM_phi1_ = -999.0;
 double  LHE_DeltaM_dEtajj_ = -999.0;
 double  LHE_DeltaM_dPhijj_ = -999.0;
 double  LHE_DeltaM_mjj_ = -999.0;
 double  LHE_DeltaM_VBSCentrality_ = -999.0;


 double LHE_MothInfo_Wqrk0_pt_ = -999.0;
 double LHE_MothInfo_Wqrk0_eta_ = -999.0;
 double LHE_MothInfo_Wqrk0_phi_ = -999.0;
 double LHE_MothInfo_Wqrk0_M_ = -999.0;
 double LHE_MothInfo_Wqrk0_E_ = -999.0;
 double LHE_MothInfo_Wqrk0_Mt_ = -999.0;

 double LHE_MothInfo_Wqrk1_pt_ = -999.0;
 double LHE_MothInfo_Wqrk1_eta_ = -999.0;
 double LHE_MothInfo_Wqrk1_phi_ = -999.0;
 double LHE_MothInfo_Wqrk1_M_ = -999.0;
 double LHE_MothInfo_Wqrk1_E_ = -999.0;
 double LHE_MothInfo_Wqrk1_Mt_ = -999.0;

 double LHE_MothInfo_Iqrk0_pt_ = -999.0; 
 double LHE_MothInfo_Iqrk0_eta_ = -999.0;
 double LHE_MothInfo_Iqrk0_phi_ = -999.0;
 double LHE_MothInfo_Iqrk0_E_ = -999.0;
 double LHE_MothInfo_Iqrk0_M_ = -999.0;
 double LHE_MothInfo_Iqrk0_Mt_ = -999.0;

 double LHE_MothInfo_Iqrk1_pt_ = -999.0; 
 double LHE_MothInfo_Iqrk1_eta_ = -999.0;
 double LHE_MothInfo_Iqrk1_phi_ = -999.0;
 double LHE_MothInfo_Iqrk1_E_ = -999.0;
 double LHE_MothInfo_Iqrk1_M_ = -999.0;
 double LHE_MothInfo_Iqrk1_Mt_ = -999.0;

 double  LHE_MothInfo_mWW_ = -999.0;
 double  LHE_MothInfo_mtWW_ = -999.0;
 double  LHE_MothInfo_mWLep_ = -999.0;
 double  LHE_MothInfo_mtWLep_ = -999.0;
 double  LHE_MothInfo_mWHad_ = -999.0;
 double  LHE_MothInfo_mtWHad_ = -999.0;
 double  LHE_MothInfo_costheta1_ = -999.0;
 double  LHE_MothInfo_costheta2_ = -999.0;
 double  LHE_MothInfo_phi_ = -999.0;
 double  LHE_MothInfo_costhetastar_ = -999.0;
 double  LHE_MothInfo_phi1_ = -999.0;
 double  LHE_MothInfo_dEtajj_ = -999.0;
 double  LHE_MothInfo_dPhijj_ = -999.0;
 double  LHE_MothInfo_mjj_ = -999.0;
 double  LHE_MothInfo_VBSCentrality_ = -999.0;

  int ngen_Lept_;
  double gen_LeptPt_;
  double gen_LeptEta_;
  double gen_LeptPhi_;
  double gen_LeptM_;
  int gen_LeptId_;
  int gen_LeptStatus_;
  double gen_LeptMother_;
  int gen_LeptGrandMother_;
    
 int ngen_Nu_;
 double gen_NuPt_;
 double gen_NuEta_;
 double gen_NuPhi_;
 double gen_NuM_;
 double gen_NuQ_;
 int gen_Nustatus_;
 int gen_NuMother_;
 int gen_NuGrandMother_;
 int gen_NuPdgId_;

 int ngen_WJet1__;
 double gen_WJet1_Pt_;
 double gen_WJet1_Eta_;
 double gen_WJet1_Phi_;
 double gen_WJet1_M_;
 double gen_WJet1_E_;
 double gen_WJet1_Q_;
 int gen_WJet1_status_;
 int gen_WJet1_Mother_;
 int gen_WJet1_GrandMother_;
 int gen_WJet1_PdgId_;

 int ngen_WJet2__;
 double gen_WJet2_Pt_;
 double gen_WJet2_Eta_;
 double gen_WJet2_Phi_;
 double gen_WJet2_M_;
 double gen_WJet2_E_;
 double gen_WJet2_Q_;
 int gen_WJet2_status_;
 int gen_WJet2_Mother_;
 int gen_WJet2_GrandMother_;
 int gen_WJet2_PdgId_;

 int ngen_VBFjet1__;
 double gen_VBFjet1_Pt_;
 double gen_VBFjet1_Eta_;
 double gen_VBFjet1_Phi_;
 double gen_VBFjet1_M_;
 double gen_VBFjet1_E_;
 double gen_VBFjet1_Q_;
 int gen_VBFjet1_status_;
 int gen_VBFjet1_Mother_;
 int gen_VBFjet1_GrandMother_;
 int gen_VBFjet1_PdgId_;

 int ngen_VBFjet2__;
 double gen_VBFjet2_Pt_;
 double gen_VBFjet2_Eta_;
 double gen_VBFjet2_Phi_;
 double gen_VBFjet2_M_;
 double gen_VBFjet2_E_;
 double gen_VBFjet2_Q_;
 int gen_VBFjet2_status_;
 int gen_VBFjet2_Mother_;
 int gen_VBFjet2_GrandMother_;
 int gen_VBFjet2_PdgId_;

double gen_VBFjet1jet2_Pt_;
double gen_VBFjet1jet2_Eta_;
double gen_VBFjet1jet2_Phi_;
double gen_VBFjet1jet2_M_;
double gen_vbfjet_deltaR_;

double gen_WHad_Pt_;
double gen_WHad_M_;
double gen_WHad_Mt_;
double gen_WHad_deltaeta_;
double gen_WHad_deltaphi_;
double gen_WHad_deltar_;

double gen_mWW_;
double gen_mtWW_;
double gen_mWLep_;
double gen_mtWLep_;
double gen_mWHad_;
double gen_mtWHad_;
double gen_costheta1_;
double gen_costheta2_;
double gen_phi_;
double gen_costhetastar_;
double gen_phi1_;
double gen_dEtajj_;
double gen_dPhijj_;
double gen_mjj_;
double gen_VBSCentrality_;
//double 

double gen_deltaR_LepWHad_;
double gen_deltaphi_NuWHad_;
double gen_deltaphi_WlepWHad_;


  int ngenJet_;
  int nVBFJet_;


  std::vector<double> genQuarkStatus_;
    
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

class PtGreater {
         public:
         template <typename T> bool operator () (const T& i, const T& j) {
         return (i->pt() > j->pt());
         }
};


 
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
AddBranch(&isMuMinus_ , "isMuMinus");
AddBranch(&LHELeptPt_ ,	"LHELeptPt");
AddBranch(&LHELeptEta_ ,	"LHELeptEta");
AddBranch(&LHELeptPhi_ ,	"LHELeptPhi");
AddBranch(&LHELeptM_ ,	"LHELeptM");
AddBranch(&LHELeptE_ ,	"LHELeptE");
AddBranch(&LHENuPt_ ,	"LHENuPt");
AddBranch(&LHENuEta_ ,	"LHENuEta");
AddBranch(&LHENuPhi_ ,	"LHENuPhi");
AddBranch(&LHENuM_ ,	"LHENuM");
AddBranch(&LHENuE_ ,	"LHENuE");
AddBranch(&LHE_DeltaM_Wqrk0_pt_ ,	"LHE_DeltaM_Wqrk0_pt");
AddBranch(&LHE_DeltaM_Wqrk0_eta_ ,	"LHE_DeltaM_Wqrk0_eta");
AddBranch(&LHE_DeltaM_Wqrk0_phi_ ,	"LHE_DeltaM_Wqrk0_phi");
AddBranch(&LHE_DeltaM_Wqrk0_M_ ,	"LHE_DeltaM_Wqrk0_M");
AddBranch(&LHE_DeltaM_Wqrk0_E_ ,	"LHE_DeltaM_Wqrk0_E");
AddBranch(&LHE_DeltaM_Wqrk0_Mt_ ,	"LHE_DeltaM_Wqrk0_Mt");
AddBranch(&LHE_DeltaM_Wqrk1_pt_ ,	"LHE_DeltaM_Wqrk1_pt");
AddBranch(&LHE_DeltaM_Wqrk1_eta_ ,	"LHE_DeltaM_Wqrk1_eta");
AddBranch(&LHE_DeltaM_Wqrk1_phi_ ,	"LHE_DeltaM_Wqrk1_phi");
AddBranch(&LHE_DeltaM_Wqrk1_M_ ,	"LHE_DeltaM_Wqrk1_M");
AddBranch(&LHE_DeltaM_Wqrk1_E_ ,	"LHE_DeltaM_Wqrk1_E");
AddBranch(&LHE_DeltaM_Wqrk1_Mt_ ,	"LHE_DeltaM_Wqrk1_Mt");
AddBranch(&LHE_DeltaM_Iqrk0_pt_ ,	"LHE_DeltaM_Iqrk0_pt");
AddBranch(&LHE_DeltaM_Iqrk0_eta_ ,	"LHE_DeltaM_Iqrk0_eta");
AddBranch(&LHE_DeltaM_Iqrk0_phi_ ,	"LHE_DeltaM_Iqrk0_phi");
AddBranch(&LHE_DeltaM_Iqrk0_E_ ,	"LHE_DeltaM_Iqrk0_E");
AddBranch(&LHE_DeltaM_Iqrk0_M_ ,	"LHE_DeltaM_Iqrk0_M");
AddBranch(&LHE_DeltaM_Iqrk0_Mt_ ,	"LHE_DeltaM_Iqrk0_Mt");
AddBranch(&LHE_DeltaM_Iqrk1_pt_ ,	"LHE_DeltaM_Iqrk1_pt");
AddBranch(&LHE_DeltaM_Iqrk1_eta_ ,	"LHE_DeltaM_Iqrk1_eta");
AddBranch(&LHE_DeltaM_Iqrk1_phi_ ,	"LHE_DeltaM_Iqrk1_phi");
AddBranch(&LHE_DeltaM_Iqrk1_E_ ,	"LHE_DeltaM_Iqrk1_E");
AddBranch(&LHE_DeltaM_Iqrk1_M_ ,	"LHE_DeltaM_Iqrk1_M");
AddBranch(&LHE_DeltaM_Iqrk1_Mt_ ,	"LHE_DeltaM_Iqrk1_Mt");
AddBranch(&LHE_DeltaM_mWW_ ,	"LHE_DeltaM_mWW");
AddBranch(&LHE_DeltaM_mtWW_ ,	"LHE_DeltaM_mtWW");
AddBranch(&LHE_DeltaM_mWLep_ ,	"LHE_DeltaM_mWLep");
AddBranch(&LHE_DeltaM_mtWLep_ ,	"LHE_DeltaM_mtWLep");
AddBranch(&LHE_DeltaM_mWHad_ ,	"LHE_DeltaM_mWHad");
AddBranch(&LHE_DeltaM_mtWHad_ ,	"LHE_DeltaM_mtWHad");
AddBranch(&LHE_DeltaM_costheta1_ ,	"LHE_DeltaM_costheta1");
AddBranch(&LHE_DeltaM_costheta2_ ,	"LHE_DeltaM_costheta2");
AddBranch(&LHE_DeltaM_phi_ ,	"LHE_DeltaM_phi");
AddBranch(&LHE_DeltaM_costhetastar_ ,	"LHE_DeltaM_costhetastar");
AddBranch(&LHE_DeltaM_phi1_ ,	"LHE_DeltaM_phi1");
AddBranch(&LHE_DeltaM_dEtajj_ ,	"LHE_DeltaM_dEtajj");
AddBranch(&LHE_DeltaM_dPhijj_ ,	"LHE_DeltaM_dPhijj");
AddBranch(&LHE_DeltaM_mjj_ ,	"LHE_DeltaM_mjj");
AddBranch(&LHE_DeltaM_VBSCentrality_ ,	"LHE_DeltaM_VBSCentrality");


AddBranch(&LHE_MothInfo_Wqrk0_pt_ ,	"LHE_MothInfo_Wqrk0_pt");
AddBranch(&LHE_MothInfo_Wqrk0_eta_ ,	"LHE_MothInfo_Wqrk0_eta");
AddBranch(&LHE_MothInfo_Wqrk0_phi_ ,	"LHE_MothInfo_Wqrk0_phi");
AddBranch(&LHE_MothInfo_Wqrk0_M_ ,	"LHE_MothInfo_Wqrk0_M");
AddBranch(&LHE_MothInfo_Wqrk0_E_ ,	"LHE_MothInfo_Wqrk0_E");
AddBranch(&LHE_MothInfo_Wqrk0_Mt_ ,	"LHE_MothInfo_Wqrk0_Mt");
AddBranch(&LHE_MothInfo_Wqrk1_pt_ ,	"LHE_MothInfo_Wqrk1_pt");
AddBranch(&LHE_MothInfo_Wqrk1_eta_ ,	"LHE_MothInfo_Wqrk1_eta");
AddBranch(&LHE_MothInfo_Wqrk1_phi_ ,	"LHE_MothInfo_Wqrk1_phi");
AddBranch(&LHE_MothInfo_Wqrk1_M_ ,	"LHE_MothInfo_Wqrk1_M");
AddBranch(&LHE_MothInfo_Wqrk1_E_ ,	"LHE_MothInfo_Wqrk1_E");
AddBranch(&LHE_MothInfo_Wqrk1_Mt_ ,	"LHE_MothInfo_Wqrk1_Mt");
AddBranch(&LHE_MothInfo_Iqrk0_pt_ ,	"LHE_MothInfo_Iqrk0_pt");
AddBranch(&LHE_MothInfo_Iqrk0_eta_ ,	"LHE_MothInfo_Iqrk0_eta");
AddBranch(&LHE_MothInfo_Iqrk0_phi_ ,	"LHE_MothInfo_Iqrk0_phi");
AddBranch(&LHE_MothInfo_Iqrk0_E_ ,	"LHE_MothInfo_Iqrk0_E");
AddBranch(&LHE_MothInfo_Iqrk0_M_ ,	"LHE_MothInfo_Iqrk0_M");
AddBranch(&LHE_MothInfo_Iqrk0_Mt_ ,	"LHE_MothInfo_Iqrk0_Mt");
AddBranch(&LHE_MothInfo_Iqrk1_pt_ ,	"LHE_MothInfo_Iqrk1_pt");
AddBranch(&LHE_MothInfo_Iqrk1_eta_ ,	"LHE_MothInfo_Iqrk1_eta");
AddBranch(&LHE_MothInfo_Iqrk1_phi_ ,	"LHE_MothInfo_Iqrk1_phi");
AddBranch(&LHE_MothInfo_Iqrk1_E_ ,	"LHE_MothInfo_Iqrk1_E");
AddBranch(&LHE_MothInfo_Iqrk1_M_ ,	"LHE_MothInfo_Iqrk1_M");
AddBranch(&LHE_MothInfo_Iqrk1_Mt_ ,	"LHE_MothInfo_Iqrk1_Mt");
AddBranch(&LHE_MothInfo_mWW_ ,	"LHE_MothInfo_mWW");
AddBranch(&LHE_MothInfo_mtWW_ ,	"LHE_MothInfo_mtWW");
AddBranch(&LHE_MothInfo_mWLep_ ,	"LHE_MothInfo_mWLep");
AddBranch(&LHE_MothInfo_mtWLep_ ,	"LHE_MothInfo_mtWLep");
AddBranch(&LHE_MothInfo_mWHad_ ,	"LHE_MothInfo_mWHad");
AddBranch(&LHE_MothInfo_mtWHad_ ,	"LHE_MothInfo_mtWHad");
AddBranch(&LHE_MothInfo_costheta1_ ,	"LHE_MothInfo_costheta1");
AddBranch(&LHE_MothInfo_costheta2_ ,	"LHE_MothInfo_costheta2");
AddBranch(&LHE_MothInfo_phi_ ,	"LHE_MothInfo_phi");
AddBranch(&LHE_MothInfo_costhetastar_ ,	"LHE_MothInfo_costhetastar");
AddBranch(&LHE_MothInfo_phi1_ ,	"LHE_MothInfo_phi1");
AddBranch(&LHE_MothInfo_dEtajj_ ,	"LHE_MothInfo_dEtajj");
AddBranch(&LHE_MothInfo_dPhijj_ ,	"LHE_MothInfo_dPhijj");
AddBranch(&LHE_MothInfo_mjj_ ,	"LHE_MothInfo_mjj");
AddBranch(&LHE_MothInfo_VBSCentrality_ ,	"LHE_MothInfo_VBSCentrality");
	
  	AddBranch(&LHEWeightIDs_, "LHEWeightIDs");
  	AddBranch(&LHEWeights_, "LHEWeights");

	AddBranch(&ngen_Lept_, "ngen_Lept");
	AddBranch(&gen_LeptPt_, "gen_LeptPt");
	AddBranch(&gen_LeptEta_,"gen_LeptEta");
	AddBranch(&gen_LeptPhi_,"gen_LeptPhi");
	AddBranch(&gen_LeptM_,"gen_LeptM");
	AddBranch(&gen_LeptStatus_,"gen_LeptStatus");
	AddBranch(&gen_LeptId_,"gen_LeptId");
	AddBranch(&gen_LeptMother_,"gen_LeptMother");
	AddBranch(&gen_LeptGrandMother_,"gen_LeptGrandMother");

  	AddBranch(&gen_NuPdgId_,"gen_NuPdgId");
  	AddBranch(&ngen_Nu_,"ngen_Nu");
  	AddBranch(&gen_NuPt_, "gen_NuPt");
  	AddBranch(&gen_NuEta_,"gen_NuEta");
  	AddBranch(&gen_NuPhi_,"gen_NuPhi");
  	AddBranch(&gen_NuM_,"gen_NuM");
  	AddBranch(&gen_NuQ_,"gen_NuQ");
  	AddBranch(&gen_Nustatus_,"gen_Nustatus");
  	AddBranch(&gen_NuMother_,"gen_NuMother");
  	AddBranch(&gen_NuGrandMother_,"gen_NuGrandMother");

  	AddBranch(&gen_WJet1_PdgId_,"gen_WJet1_PdgId");
  	AddBranch(&ngen_WJet1__,"ngen_WJet1_");
  	AddBranch(&gen_WJet1_Pt_, "gen_WJet1_Pt");
  	AddBranch(&gen_WJet1_Eta_,"gen_WJet1_Eta");
  	AddBranch(&gen_WJet1_Phi_,"gen_WJet1_Phi");
  	AddBranch(&gen_WJet1_M_,"gen_WJet1_M");
	AddBranch(&gen_WJet1_E_, "gen_WJet1_E");
  	AddBranch(&gen_WJet1_Q_,"gen_WJet1_Q");
  	AddBranch(&gen_WJet1_status_,"gen_WJet1_status");
  	AddBranch(&gen_WJet1_Mother_,"gen_WJet1_Mother");
  	AddBranch(&gen_WJet1_GrandMother_,"gen_WJet1_GrandMother");

  	AddBranch(&gen_WJet2_PdgId_,"gen_WJet2_PdgId");
  	AddBranch(&ngen_WJet2__,"ngen_WJet2_");
  	AddBranch(&gen_WJet2_Pt_, "gen_WJet2_Pt");
  	AddBranch(&gen_WJet2_Eta_,"gen_WJet2_Eta");
  	AddBranch(&gen_WJet2_Phi_,"gen_WJet2_Phi");
  	AddBranch(&gen_WJet2_M_,"gen_WJet2_M");
	AddBranch(&gen_WJet2_E_, "gen_WJet2_E");
  	AddBranch(&gen_WJet2_Q_,"gen_WJet2_Q");
  	AddBranch(&gen_WJet2_status_,"gen_WJet2_status");
  	AddBranch(&gen_WJet2_Mother_,"gen_WJet2_Mother");
  	AddBranch(&gen_WJet2_GrandMother_,"gen_WJet2_GrandMother");

  	AddBranch(&gen_VBFjet1_PdgId_,"gen_VBFjet1_PdgId");
  	AddBranch(&ngen_VBFjet1__,"ngen_VBFjet1_");
  	AddBranch(&gen_VBFjet1_Pt_, "gen_VBFjet1_Pt");
  	AddBranch(&gen_VBFjet1_Eta_,"gen_VBFjet1_Eta");
  	AddBranch(&gen_VBFjet1_Phi_,"gen_VBFjet1_Phi");
  	AddBranch(&gen_VBFjet1_M_,"gen_VBFjet1_M");
  	AddBranch(&gen_VBFjet1_Q_,"gen_VBFjet1_Q");
  	AddBranch(&gen_VBFjet1_status_,"gen_VBFjet1_status");
  	AddBranch(&gen_VBFjet1_Mother_,"gen_VBFjet1_Mother");
  	AddBranch(&gen_VBFjet1_GrandMother_,"gen_VBFjet1_GrandMother");

  	AddBranch(&gen_VBFjet2_PdgId_,"gen_VBFjet2_PdgId");
  	AddBranch(&ngen_VBFjet2__,"ngen_VBFjet2_");
  	AddBranch(&gen_VBFjet2_Pt_, "gen_VBFjet2_Pt");
  	AddBranch(&gen_VBFjet2_Eta_,"gen_VBFjet2_Eta");
  	AddBranch(&gen_VBFjet2_Phi_,"gen_VBFjet2_Phi");
  	AddBranch(&gen_VBFjet2_M_,"gen_VBFjet2_M");
  	AddBranch(&gen_VBFjet2_Q_,"gen_VBFjet2_Q");
  	AddBranch(&gen_VBFjet2_status_,"gen_VBFjet2_status");
  	AddBranch(&gen_VBFjet2_Mother_,"gen_VBFjet2_Mother");
  	AddBranch(&gen_VBFjet2_GrandMother_,"gen_VBFjet2_GrandMother");

AddBranch(&gen_VBFjet1jet2_Pt_, "gen_VBFjet1jet2_Pt");
AddBranch(&gen_VBFjet1jet2_Eta_, "gen_VBFjet1jet2_Eta");
AddBranch(&gen_VBFjet1jet2_Phi_, "gen_VBFjet1jet2_Phi");
AddBranch(&gen_VBFjet1jet2_M_, "gen_VBFjet1jet2_M");
AddBranch(&gen_vbfjet_deltaR_, "gen_vbfjet_deltaR");
AddBranch(&gen_WHad_Pt_, "gen_WHad_Pt");
AddBranch(&gen_WHad_M_, "gen_WHad_M");
AddBranch(&gen_WHad_Mt_, "gen_WHad_Mt");
AddBranch(&gen_WHad_deltaeta_, "gen_WHad_deltaeta");
AddBranch(&gen_WHad_deltaphi_, "gen_WHad_deltaphi");
AddBranch(&gen_WHad_deltar_, "gen_WHad_deltar");
AddBranch(&gen_deltaR_LepWHad_, "gen_deltaR_LepWHad");
AddBranch(&gen_deltaphi_NuWHad_, "gen_deltaphi_NuWHad");
AddBranch(&gen_deltaphi_WlepWHad_, "gen_deltaphi_WlepWHad");

AddBranch(&gen_mWW_, "gen_mWW");
AddBranch(&gen_mtWW_, "gen_mtWW");
AddBranch(&gen_mWLep_, "gen_mWLep");
AddBranch(&gen_mtWLep_, "gen_mtWLep");
AddBranch(&gen_mWHad_, "gen_mWHad");
AddBranch(&gen_mtWHad_, "gen_mtWHad");
AddBranch(&gen_costheta1_, "gen_costheta1");
AddBranch(&gen_costheta2_, "gen_costheta2");
AddBranch(&gen_phi_, "gen_phi");
AddBranch(&gen_costhetastar_, "gen_costhetastar");
AddBranch(&gen_phi1_, "gen_phi1");
AddBranch(&gen_dEtajj_, "gen_dEtajj");
AddBranch(&gen_dPhijj_, "gen_dPhijj");
AddBranch(&gen_mjj_, "gen_mjj");
AddBranch(&gen_VBSCentrality_, "gen_VBSCentrality");

  AddBranch(&genQuarkStatus_, "genQuarkStatus");
  AddBranch(&ngenJet_, "ngenJet");
  AddBranch(&nVBFJet_, "nVBFJet");



}

void GenAnalyzer::Clear(){
	//pdgID_.clear();
LHEWeightIDs_.clear();
LHEWeights_.clear();
isMuMinus_	= -999.0;
LHELeptPt_	= -999.0;
LHELeptEta_	= -999.0;
LHELeptPhi_	= -999.0;
LHELeptM_	= -999.0;
LHELeptE_	= -999.0;
LHENuPt_	= -999.0;
LHENuEta_	= -999.0;
LHENuPhi_	= -999.0;
LHENuM_	= -999.0;
LHENuE_	= -999.0;
LHE_DeltaM_Wqrk0_pt_	= -999.0;
LHE_DeltaM_Wqrk0_eta_	= -999.0;
LHE_DeltaM_Wqrk0_phi_	= -999.0;
LHE_DeltaM_Wqrk0_M_	= -999.0;
LHE_DeltaM_Wqrk0_E_	= -999.0;
LHE_DeltaM_Wqrk0_Mt_	= -999.0;
LHE_DeltaM_Wqrk1_pt_	= -999.0;
LHE_DeltaM_Wqrk1_eta_	= -999.0;
LHE_DeltaM_Wqrk1_phi_	= -999.0;
LHE_DeltaM_Wqrk1_M_	= -999.0;
LHE_DeltaM_Wqrk1_E_	= -999.0;
LHE_DeltaM_Wqrk1_Mt_	= -999.0;
LHE_DeltaM_Iqrk0_pt_	= -999.0;
LHE_DeltaM_Iqrk0_eta_	= -999.0;
LHE_DeltaM_Iqrk0_phi_	= -999.0;
LHE_DeltaM_Iqrk0_E_	= -999.0;
LHE_DeltaM_Iqrk0_M_	= -999.0;
LHE_DeltaM_Iqrk0_Mt_	= -999.0;
LHE_DeltaM_Iqrk1_pt_	= -999.0;
LHE_DeltaM_Iqrk1_eta_	= -999.0;
LHE_DeltaM_Iqrk1_phi_	= -999.0;
LHE_DeltaM_Iqrk1_E_	= -999.0;
LHE_DeltaM_Iqrk1_M_	= -999.0;
LHE_DeltaM_Iqrk1_Mt_	= -999.0;
LHE_DeltaM_mWW_	= -999.0;
LHE_DeltaM_mtWW_	= -999.0;
LHE_DeltaM_mWLep_	= -999.0;
LHE_DeltaM_mtWLep_	= -999.0;
LHE_DeltaM_mWHad_	= -999.0;
LHE_DeltaM_mtWHad_	= -999.0;
LHE_DeltaM_costheta1_	= -999.0;
LHE_DeltaM_costheta2_	= -999.0;
LHE_DeltaM_phi_	= -999.0;
LHE_DeltaM_costhetastar_	= -999.0;
LHE_DeltaM_phi1_	= -999.0;
LHE_DeltaM_dEtajj_	= -999.0;
LHE_DeltaM_dPhijj_	= -999.0;
LHE_DeltaM_mjj_	= -999.0;
LHE_DeltaM_VBSCentrality_	= -999.0;


LHE_MothInfo_Wqrk0_pt_	= -999.0;
LHE_MothInfo_Wqrk0_eta_	= -999.0;
LHE_MothInfo_Wqrk0_phi_	= -999.0;
LHE_MothInfo_Wqrk0_M_	= -999.0;
LHE_MothInfo_Wqrk0_E_	= -999.0;
LHE_MothInfo_Wqrk0_Mt_	= -999.0;
LHE_MothInfo_Wqrk1_pt_	= -999.0;
LHE_MothInfo_Wqrk1_eta_	= -999.0;
LHE_MothInfo_Wqrk1_phi_	= -999.0;
LHE_MothInfo_Wqrk1_M_	= -999.0;
LHE_MothInfo_Wqrk1_E_	= -999.0;
LHE_MothInfo_Wqrk1_Mt_	= -999.0;
LHE_MothInfo_Iqrk0_pt_	= -999.0;
LHE_MothInfo_Iqrk0_eta_	= -999.0;
LHE_MothInfo_Iqrk0_phi_	= -999.0;
LHE_MothInfo_Iqrk0_E_	= -999.0;
LHE_MothInfo_Iqrk0_M_	= -999.0;
LHE_MothInfo_Iqrk0_Mt_	= -999.0;
LHE_MothInfo_Iqrk1_pt_	= -999.0;
LHE_MothInfo_Iqrk1_eta_	= -999.0;
LHE_MothInfo_Iqrk1_phi_	= -999.0;
LHE_MothInfo_Iqrk1_E_	= -999.0;
LHE_MothInfo_Iqrk1_M_	= -999.0;
LHE_MothInfo_Iqrk1_Mt_	= -999.0;
LHE_MothInfo_mWW_	= -999.0;
LHE_MothInfo_mtWW_	= -999.0;
LHE_MothInfo_mWLep_	= -999.0;
LHE_MothInfo_mtWLep_	= -999.0;
LHE_MothInfo_mWHad_	= -999.0;
LHE_MothInfo_mtWHad_	= -999.0;
LHE_MothInfo_costheta1_	= -999.0;
LHE_MothInfo_costheta2_	= -999.0;
LHE_MothInfo_phi_	= -999.0;
LHE_MothInfo_costhetastar_	= -999.0;
LHE_MothInfo_phi1_	= -999.0;
LHE_MothInfo_dEtajj_	= -999.0;
LHE_MothInfo_dPhijj_	= -999.0;
LHE_MothInfo_mjj_	= -999.0;
LHE_MothInfo_VBSCentrality_	= -999.0;



	ngen_Lept_ = -999;
	gen_LeptPt_ = -999.0;
	gen_LeptEta_ = -999.0;
	gen_LeptPhi_ = -999.0;
	gen_LeptStatus_ = -999;
	gen_LeptMother_ = -999.0;
	gen_LeptGrandMother_ = -999;
	gen_LeptId_ = -999;
	gen_LeptM_ = -999.0;

  ngen_Nu_ = -999;
  gen_NuPt_ = -999.0;
  gen_NuEta_ = -999.0;
  gen_NuPhi_ = -999.0;
  gen_NuM_ = -999.0;
  gen_NuQ_ = -999.0;
  gen_Nustatus_ = -999;
  gen_NuMother_ = -999;
  gen_NuGrandMother_ = -999;
  gen_NuPdgId_ = -999;

  ngen_WJet1__ = -999;
  gen_WJet1_Pt_ = -999.0;
  gen_WJet1_Eta_ = -999.0;
  gen_WJet1_Phi_ = -999.0;
  gen_WJet1_M_ = -999.0;
  gen_WJet1_E_ = -999.0;
  gen_WJet1_Q_ = -999.0;
  gen_WJet1_status_ = -999;
  gen_WJet1_Mother_ = -999;
  gen_WJet1_GrandMother_ = -999;
  gen_WJet1_PdgId_ = -999;

  ngen_WJet2__ = -999;
  gen_WJet2_Pt_ = -999.0;
  gen_WJet2_Eta_ = -999.0;
  gen_WJet2_Phi_ = -999.0;
  gen_WJet2_M_ = -999.0;
  gen_WJet2_E_ = -999.0;
  gen_WJet2_Q_ = -999.0;
  gen_WJet2_status_ = -999;
  gen_WJet2_Mother_ = -999;
  gen_WJet2_GrandMother_ = -999;
  gen_WJet2_PdgId_ = -999;

  ngen_VBFjet1__ = -999;
  gen_VBFjet1_Pt_ = -999.0;
  gen_VBFjet1_Eta_ = -999.0;
  gen_VBFjet1_Phi_ = -999.0;
  gen_VBFjet1_M_ = -999.0;
  gen_VBFjet1_Q_ = -999.0;
  gen_VBFjet1_status_ = -999;
  gen_VBFjet1_Mother_ = -999;
  gen_VBFjet1_GrandMother_ = -999;
  gen_VBFjet1_PdgId_ = -999;

  ngen_VBFjet2__ = -999;
  gen_VBFjet2_Pt_ = -999.0;
  gen_VBFjet2_Eta_ = -999.0;
  gen_VBFjet2_Phi_ = -999.0;
  gen_VBFjet2_M_ = -999.0;
  gen_VBFjet2_Q_ = -999.0;
  gen_VBFjet2_status_ = -999;
  gen_VBFjet2_Mother_ = -999;
  gen_VBFjet2_GrandMother_ = -999;
  gen_VBFjet2_PdgId_ = -999;


gen_VBFjet1jet2_Pt_ = -999.0;
gen_VBFjet1jet2_Eta_ = -999.0;
gen_VBFjet1jet2_Phi_ = -999.0;
gen_VBFjet1jet2_M_ = -999.0;
gen_vbfjet_deltaR_ = -999.0;
gen_WHad_Pt_ = -999.0;
gen_WHad_M_ = -999.0;
gen_WHad_Mt_ = -999.0;
gen_WHad_deltaeta_ = -999.0;
gen_WHad_deltaphi_ = -999.0;
gen_WHad_deltar_ = -999.0;
gen_deltaR_LepWHad_ = -999.0;
gen_deltaphi_NuWHad_ = -999.0;
gen_deltaphi_WlepWHad_ = -999.0;


gen_mWLep_ = -999.0;
gen_mtWLep_ = -999.0;
gen_mWHad_ = -999.0;
gen_mtWHad_ = -999.0;
gen_costheta1_ = -999.0;
gen_costheta2_ = -999.0;
gen_phi_ = -999.0;
gen_costhetastar_ = -999.0;
gen_phi1_ = -999.0;
gen_dEtajj_ = -999.0;
gen_dPhijj_ = -999.0;
gen_mjj_ = -999.0;
gen_VBSCentrality_ = -999.0;

nVBFJet_ = -999;
ngenJet_ = -999;
LHEWeightIDs_.clear();
LHEWeights_.clear();
genQuarkStatus_.clear();


}

