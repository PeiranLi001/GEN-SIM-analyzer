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

 double LHE_Wqrk0_pt_ = -999.0;
 double LHE_Wqrk0_eta_ = -999.0;
 double LHE_Wqrk0_phi_ = -999.0;
 double LHE_Wqrk0_M_ = -999.0;
 double LHE_Wqrk0_E_ = -999.0;
 double LHE_Wqrk0_Mt_ = -999.0;

 double LHE_Wqrk1_pt_ = -999.0;
 double LHE_Wqrk1_eta_ = -999.0;
 double LHE_Wqrk1_phi_ = -999.0;
 double LHE_Wqrk1_M_ = -999.0;
 double LHE_Wqrk1_E_ = -999.0;
 double LHE_Wqrk1_Mt_ = -999.0;

 double LHE_Iqrk0_pt_ = -999.0; 
 double LHE_Iqrk0_eta_ = -999.0;
 double LHE_Iqrk0_phi_ = -999.0;
 double LHE_Iqrk0_E_ = -999.0;
 double LHE_Iqrk0_M_ = -999.0;
 double LHE_Iqrk0_Mt_ = -999.0;

 double LHE_Iqrk1_pt_ = -999.0; 
 double LHE_Iqrk1_eta_ = -999.0;
 double LHE_Iqrk1_phi_ = -999.0;
 double LHE_Iqrk1_E_ = -999.0;
 double LHE_Iqrk1_M_ = -999.0;
 double LHE_Iqrk1_Mt_ = -999.0;

 double  LHE_mWW_ = -999.0;
 double  LHE_mtWW_ = -999.0;
 double  LHE_mWLep_ = -999.0;
 double  LHE_mtWLep_ = -999.0;
 double  LHE_mWHad_ = -999.0;
 double  LHE_mtWHad_ = -999.0;
 double  LHE_costheta1_ = -999.0;
 double  LHE_costheta2_ = -999.0;
 double  LHE_phi_ = -999.0;
 double  LHE_costhetastar_ = -999.0;
 double  LHE_phi1_ = -999.0;
 double  LHE_dEtajj_ = -999.0;
 double  LHE_dPhijj_ = -999.0;
 double  LHE_mjj_ = -999.0;
 double  LHE_VBSCentrality_ = -999.0;
    
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
	//genParticles_(iConfig.getParameter<edm::InputTag>( "genParticle"))
	LHEEventToken = consumes<LHEEventProduct>(iConfig.getParameter<edm::InputTag>("LHEEventInputTag"));

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
AddBranch(&LHE_Wqrk0_pt_ ,	"LHE_Wqrk0_pt");
AddBranch(&LHE_Wqrk0_eta_ ,	"LHE_Wqrk0_eta");
AddBranch(&LHE_Wqrk0_phi_ ,	"LHE_Wqrk0_phi");
AddBranch(&LHE_Wqrk0_M_ ,	"LHE_Wqrk0_M");
AddBranch(&LHE_Wqrk0_E_ ,	"LHE_Wqrk0_E");
AddBranch(&LHE_Wqrk0_Mt_ ,	"LHE_Wqrk0_Mt");
AddBranch(&LHE_Wqrk1_pt_ ,	"LHE_Wqrk1_pt");
AddBranch(&LHE_Wqrk1_eta_ ,	"LHE_Wqrk1_eta");
AddBranch(&LHE_Wqrk1_phi_ ,	"LHE_Wqrk1_phi");
AddBranch(&LHE_Wqrk1_M_ ,	"LHE_Wqrk1_M");
AddBranch(&LHE_Wqrk1_E_ ,	"LHE_Wqrk1_E");
AddBranch(&LHE_Wqrk1_Mt_ ,	"LHE_Wqrk1_Mt");
AddBranch(&LHE_Iqrk0_pt_ ,	"LHE_Iqrk0_pt");
AddBranch(&LHE_Iqrk0_eta_ ,	"LHE_Iqrk0_eta");
AddBranch(&LHE_Iqrk0_phi_ ,	"LHE_Iqrk0_phi");
AddBranch(&LHE_Iqrk0_E_ ,	"LHE_Iqrk0_E");
AddBranch(&LHE_Iqrk0_M_ ,	"LHE_Iqrk0_M");
AddBranch(&LHE_Iqrk0_Mt_ ,	"LHE_Iqrk0_Mt");
AddBranch(&LHE_Iqrk1_pt_ ,	"LHE_Iqrk1_pt");
AddBranch(&LHE_Iqrk1_eta_ ,	"LHE_Iqrk1_eta");
AddBranch(&LHE_Iqrk1_phi_ ,	"LHE_Iqrk1_phi");
AddBranch(&LHE_Iqrk1_E_ ,	"LHE_Iqrk1_E");
AddBranch(&LHE_Iqrk1_M_ ,	"LHE_Iqrk1_M");
AddBranch(&LHE_Iqrk1_Mt_ ,	"LHE_Iqrk1_Mt");
AddBranch(&LHE_mWW_ ,	"LHE_mWW");
AddBranch(&LHE_mtWW_ ,	"LHE_mtWW");
AddBranch(&LHE_mWLep_ ,	"LHE_mWLep");
AddBranch(&LHE_mtWLep_ ,	"LHE_mtWLep");
AddBranch(&LHE_mWHad_ ,	"LHE_mWHad");
AddBranch(&LHE_mtWHad_ ,	"LHE_mtWHad");
AddBranch(&LHE_costheta1_ ,	"LHE_costheta1");
AddBranch(&LHE_costheta2_ ,	"LHE_costheta2");
AddBranch(&LHE_phi_ ,	"LHE_phi");
AddBranch(&LHE_costhetastar_ ,	"LHE_costhetastar");
AddBranch(&LHE_phi1_ ,	"LHE_phi1");
AddBranch(&LHE_dEtajj_ ,	"LHE_dEtajj");
AddBranch(&LHE_dPhijj_ ,	"LHE_dPhijj");
AddBranch(&LHE_mjj_ ,	"LHE_mjj");
AddBranch(&LHE_VBSCentrality_ ,	"LHE_VBSCentrality");
	
  	AddBranch(&LHEWeightIDs_, "LHEWeightIDs");
  	AddBranch(&LHEWeights_, "LHEWeights");

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
LHE_Wqrk0_pt_	= -999.0;
LHE_Wqrk0_eta_	= -999.0;
LHE_Wqrk0_phi_	= -999.0;
LHE_Wqrk0_M_	= -999.0;
LHE_Wqrk0_E_	= -999.0;
LHE_Wqrk0_Mt_	= -999.0;
LHE_Wqrk1_pt_	= -999.0;
LHE_Wqrk1_eta_	= -999.0;
LHE_Wqrk1_phi_	= -999.0;
LHE_Wqrk1_M_	= -999.0;
LHE_Wqrk1_E_	= -999.0;
LHE_Wqrk1_Mt_	= -999.0;
LHE_Iqrk0_pt_	= -999.0;
LHE_Iqrk0_eta_	= -999.0;
LHE_Iqrk0_phi_	= -999.0;
LHE_Iqrk0_E_	= -999.0;
LHE_Iqrk0_M_	= -999.0;
LHE_Iqrk0_Mt_	= -999.0;
LHE_Iqrk1_pt_	= -999.0;
LHE_Iqrk1_eta_	= -999.0;
LHE_Iqrk1_phi_	= -999.0;
LHE_Iqrk1_E_	= -999.0;
LHE_Iqrk1_M_	= -999.0;
LHE_Iqrk1_Mt_	= -999.0;
LHE_mWW_	= -999.0;
LHE_mtWW_	= -999.0;
LHE_mWLep_	= -999.0;
LHE_mtWLep_	= -999.0;
LHE_mWHad_	= -999.0;
LHE_mtWHad_	= -999.0;
LHE_costheta1_	= -999.0;
LHE_costheta2_	= -999.0;
LHE_phi_	= -999.0;
LHE_costhetastar_	= -999.0;
LHE_phi1_	= -999.0;
LHE_dEtajj_	= -999.0;
LHE_dPhijj_	= -999.0;
LHE_mjj_	= -999.0;
LHE_VBSCentrality_	= -999.0;

}

