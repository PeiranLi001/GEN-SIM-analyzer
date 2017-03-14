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

int nEVENT=0;
  int ngenLept_;
  double genLeptPt_;
  double genLeptEta_;
  double genLeptPhi_;
  double genLeptM_;
  int genLeptId_;
  int genLeptStatus_;
  double genLeptMother_;
  int genLeptGrandMother_;
    
 int ngenNu_;
 double genNuPt_;
 double genNuEta_;
 double genNuPhi_;
 double genNuM_;
 double genNuQ_;
 int genNustatus_;
 int genNuMother_;
 int genNuGrandMother_;
 int genNuPdgId_;

 int ngenWJet1__;
 double genWJet1_Pt_;
 double genWJet1_Eta_;
 double genWJet1_Phi_;
 double genWJet1_M_;
 double genWJet1_Q_;
 int genWJet1_status_;
 int genWJet1_Mother_;
 int genWJet1_GrandMother_;
 int genWJet1_PdgId_;

 int ngenWJet2__;
 double genWJet2_Pt_;
 double genWJet2_Eta_;
 double genWJet2_Phi_;
 double genWJet2_M_;
 double genWJet2_Q_;
 int genWJet2_status_;
 int genWJet2_Mother_;
 int genWJet2_GrandMother_;
 int genWJet2_PdgId_;

 int ngenVBFjet1__;
 double genVBFjet1_Pt_;
 double genVBFjet1_Eta_;
 double genVBFjet1_Phi_;
 double genVBFjet1_M_;
 double genVBFjet1_Q_;
 int genVBFjet1_status_;
 int genVBFjet1_Mother_;
 int genVBFjet1_GrandMother_;
 int genVBFjet1_PdgId_;

 int ngenVBFjet2__;
 double genVBFjet2_Pt_;
 double genVBFjet2_Eta_;
 double genVBFjet2_Phi_;
 double genVBFjet2_M_;
 double genVBFjet2_Q_;
 int genVBFjet2_status_;
 int genVBFjet2_Mother_;
 int genVBFjet2_GrandMother_;
 int genVBFjet2_PdgId_;

  int ngenJet_;
  int nVBFJet_;
  double vbf_maxpt_j1_pt_;
  double vbf_maxpt_j1_eta_;
  double vbf_maxpt_j1_phi_;
  double vbf_maxpt_j1_e_;
  double vbf_maxpt_j1_bDiscriminatorCSV_;
  double vbf_maxpt_j2_pt_;
  double vbf_maxpt_j2_eta_;
  double vbf_maxpt_j2_phi_;
  double vbf_maxpt_j2_e_;
  double vbf_maxpt_j2_bDiscriminatorCSV_;
  double vbf_maxpt_jj_pt_;
  double vbf_maxpt_jj_eta_;
  double vbf_maxpt_jj_phi_;
  double vbf_maxpt_jj_m_;
  double vbf_maxpt_deltaR_;
  double AK4_jet1_pt_;
  double AK4_jet1_eta_;
  double AK4_jet1_phi_;
  double AK4_jet1_e_;
  double AK4_jet1_bDiscriminatorCSV_;
  double AK4_jet2_pt_;
  double AK4_jet2_eta_;
  double AK4_jet2_phi_;
  double AK4_jet2_e_;
  double AK4_jet2_bDiscriminatorCSV_;
  double AK4_jetjet_pt_;
  double AK4_jetjet_mass_;
  double AK4_jetjet_deltaeta_;
  double AK4_jetjet_deltaphi_;
  double AK4_jetjet_deltar_;
  double deltaR_lak4jetjet_;
  double deltaphi_METak4jetjet_;
  double deltaphi_Vak4jetjet_;
  double mass_lvjj_run2_AK4_;

  std::vector<double> genQuarkStatus_;
  std::vector<double> genJetPt_;
  std::vector<double> genJetEta_;
  std::vector<double> genJetPhi_;
  std::vector<double> genJetMass_;
  std::vector<double> genCaloMET_;
  std::vector<double> genCaloMETPhi_;
  std::vector<double> genTrueMET_;
  std::vector<double> genTrueMETPhi_;

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
	genParticlesToken = consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genParticlesInputTag"));
	LHEEventToken = consumes<LHEEventProduct>(iConfig.getParameter<edm::InputTag>("LHEEventInputTag"));
	genAK4jetToken = consumes<reco::GenJetCollection>(iConfig.getParameter<edm::InputTag>("genJetsAK4jetsInputTag"));
	genMetCaloToken= consumes<reco::GenMETCollection>(iConfig.getParameter<edm::InputTag>("genMetCaloInputTag"));
	//genMetCaloToken= consumes<reco::GenMETCollection>(iConfig.getParameter<edm::InputTag>("genMetCaloInputTag"));

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

void GenAnalyzer::SetBranches(){
	//AddBranch(&pdgID_,	"pdgID");
	
  	AddBranch(&LHEWeightIDs_, "LHEWeightIDs");
  	AddBranch(&LHEWeights_, "LHEWeights");

	AddBranch(&ngenLept_, "ngenLept");
	AddBranch(&genLeptPt_, "genLeptPt");
	AddBranch(&genLeptEta_,"genLeptEta");
	AddBranch(&genLeptPhi_,"genLeptPhi");
	AddBranch(&genLeptM_,"genLeptM");
	AddBranch(&genLeptStatus_,"genLeptStatus");
	AddBranch(&genLeptId_,"genLeptId");
	AddBranch(&genLeptMother_,"genLeptMother");
	AddBranch(&genLeptGrandMother_,"genLeptGrandMother");

  	AddBranch(&genNuPdgId_,"genNuPdgId");
  	AddBranch(&ngenNu_,"ngenNu");
  	AddBranch(&genNuPt_, "genNuPt");
  	AddBranch(&genNuEta_,"genNuEta");
  	AddBranch(&genNuPhi_,"genNuPhi");
  	AddBranch(&genNuM_,"genNuM");
  	AddBranch(&genNuQ_,"genNuQ");
  	AddBranch(&genNustatus_,"genNustatus");
  	AddBranch(&genNuMother_,"genNuMother");
  	AddBranch(&genNuGrandMother_,"genNuGrandMother");

  	AddBranch(&genWJet1_PdgId_,"genWJet1_PdgId");
  	AddBranch(&ngenWJet1__,"ngenWJet1_");
  	AddBranch(&genWJet1_Pt_, "genWJet1_Pt");
  	AddBranch(&genWJet1_Eta_,"genWJet1_Eta");
  	AddBranch(&genWJet1_Phi_,"genWJet1_Phi");
  	AddBranch(&genWJet1_M_,"genWJet1_M");
  	AddBranch(&genWJet1_Q_,"genWJet1_Q");
  	AddBranch(&genWJet1_status_,"genWJet1_status");
  	AddBranch(&genWJet1_Mother_,"genWJet1_Mother");
  	AddBranch(&genWJet1_GrandMother_,"genWJet1_GrandMother");

  	AddBranch(&genWJet2_PdgId_,"genWJet2_PdgId");
  	AddBranch(&ngenWJet2__,"ngenWJet2_");
  	AddBranch(&genWJet2_Pt_, "genWJet2_Pt");
  	AddBranch(&genWJet2_Eta_,"genWJet2_Eta");
  	AddBranch(&genWJet2_Phi_,"genWJet2_Phi");
  	AddBranch(&genWJet2_M_,"genWJet2_M");
  	AddBranch(&genWJet2_Q_,"genWJet2_Q");
  	AddBranch(&genWJet2_status_,"genWJet2_status");
  	AddBranch(&genWJet2_Mother_,"genWJet2_Mother");
  	AddBranch(&genWJet2_GrandMother_,"genWJet2_GrandMother");

  	AddBranch(&genVBFjet1_PdgId_,"genVBFjet1_PdgId");
  	AddBranch(&ngenVBFjet1__,"ngenVBFjet1_");
  	AddBranch(&genVBFjet1_Pt_, "genVBFjet1_Pt");
  	AddBranch(&genVBFjet1_Eta_,"genVBFjet1_Eta");
  	AddBranch(&genVBFjet1_Phi_,"genVBFjet1_Phi");
  	AddBranch(&genVBFjet1_M_,"genVBFjet1_M");
  	AddBranch(&genVBFjet1_Q_,"genVBFjet1_Q");
  	AddBranch(&genVBFjet1_status_,"genVBFjet1_status");
  	AddBranch(&genVBFjet1_Mother_,"genVBFjet1_Mother");
  	AddBranch(&genVBFjet1_GrandMother_,"genVBFjet1_GrandMother");

  	AddBranch(&genVBFjet2_PdgId_,"genVBFjet2_PdgId");
  	AddBranch(&ngenVBFjet2__,"ngenVBFjet2_");
  	AddBranch(&genVBFjet2_Pt_, "genVBFjet2_Pt");
  	AddBranch(&genVBFjet2_Eta_,"genVBFjet2_Eta");
  	AddBranch(&genVBFjet2_Phi_,"genVBFjet2_Phi");
  	AddBranch(&genVBFjet2_M_,"genVBFjet2_M");
  	AddBranch(&genVBFjet2_Q_,"genVBFjet2_Q");
  	AddBranch(&genVBFjet2_status_,"genVBFjet2_status");
  	AddBranch(&genVBFjet2_Mother_,"genVBFjet2_Mother");
  	AddBranch(&genVBFjet2_GrandMother_,"genVBFjet2_GrandMother");

  AddBranch(&genQuarkStatus_, "genQuarkStatus");
  AddBranch(&genJetPt_, "genJetPt");
  AddBranch(&genJetEta_, "genJetEta");
  AddBranch(&genJetPhi_, "genJetPhi");
  AddBranch(&genJetMass_, "genJetMass");
  AddBranch(&ngenJet_, "ngenJet");
  AddBranch(&nVBFJet_, "nVBFJet");

  AddBranch(&genCaloMET_, "genCaloMET");
  AddBranch(&genCaloMETPhi_, "genCaloMETPhi");
  AddBranch(&genTrueMET_, "genTrueMET");
  AddBranch(&genTrueMETPhi_, "genTrueMETPhi");

AddBranch(&vbf_maxpt_j1_pt_, "vbf_maxpt_j1_pt");
AddBranch(&vbf_maxpt_j1_eta_, "vbf_maxpt_j1_eta");
AddBranch(&vbf_maxpt_j1_phi_, "vbf_maxpt_j1_phi");
AddBranch(&vbf_maxpt_j1_e_, "vbf_maxpt_j1_e");
AddBranch(&vbf_maxpt_j1_bDiscriminatorCSV_, "vbf_maxpt_j1_bDiscriminatorCSV");
AddBranch(&vbf_maxpt_j2_pt_, "vbf_maxpt_j2_pt");
AddBranch(&vbf_maxpt_j2_eta_, "vbf_maxpt_j2_eta");
AddBranch(&vbf_maxpt_j2_phi_, "vbf_maxpt_j2_phi");
AddBranch(&vbf_maxpt_j2_e_, "vbf_maxpt_j2_e");
AddBranch(&vbf_maxpt_j2_bDiscriminatorCSV_, "vbf_maxpt_j2_bDiscriminatorCSV");
AddBranch(&vbf_maxpt_jj_pt_, "vbf_maxpt_jj_pt");
AddBranch(&vbf_maxpt_jj_eta_, "vbf_maxpt_jj_eta");
AddBranch(&vbf_maxpt_jj_phi_, "vbf_maxpt_jj_phi");
AddBranch(&vbf_maxpt_jj_m_, "vbf_maxpt_jj_m");
AddBranch(&vbf_maxpt_deltaR_, "vbf_maxpt_deltaR");
AddBranch(&AK4_jet1_pt_, "AK4_jet1_pt");
AddBranch(&AK4_jet1_eta_, "AK4_jet1_eta");
AddBranch(&AK4_jet1_phi_, "AK4_jet1_phi");
AddBranch(&AK4_jet1_e_, "AK4_jet1_e");
AddBranch(&AK4_jet1_bDiscriminatorCSV_, "AK4_jet1_bDiscriminatorCSV");
AddBranch(&AK4_jet2_pt_, "AK4_jet2_pt");
AddBranch(&AK4_jet2_eta_, "AK4_jet2_eta");
AddBranch(&AK4_jet2_phi_, "AK4_jet2_phi");
AddBranch(&AK4_jet2_e_, "AK4_jet2_e");
AddBranch(&AK4_jet2_bDiscriminatorCSV_, "AK4_jet2_bDiscriminatorCSV");
AddBranch(&AK4_jetjet_pt_, "AK4_jetjet_pt");
AddBranch(&AK4_jetjet_mass_, "AK4_jetjet_mass");
AddBranch(&AK4_jetjet_deltaeta_, "AK4_jetjet_deltaeta");
AddBranch(&AK4_jetjet_deltaphi_, "AK4_jetjet_deltaphi");
AddBranch(&AK4_jetjet_deltar_, "AK4_jetjet_deltar");
AddBranch(&deltaR_lak4jetjet_, "deltaR_lak4jetjet");
AddBranch(&deltaphi_METak4jetjet_, "deltaphi_METak4jetjet");
AddBranch(&deltaphi_Vak4jetjet_, "deltaphi_Vak4jetjet");
AddBranch(&mass_lvjj_run2_AK4_, "mass_lvjj_run2_AK4");

}

void GenAnalyzer::Clear(){
	//pdgID_.clear();

	ngenLept_ = -999;
	genLeptPt_ = -999.0;
	genLeptEta_ = -999.0;
	genLeptPhi_ = -999.0;
	genLeptStatus_ = -999;
	genLeptMother_ = -999.0;
	genLeptGrandMother_ = -999;
	genLeptId_ = -999;
	genLeptM_ = -999.0;

  ngenNu_ = -999;
  genNuPt_ = -999.0;
  genNuEta_ = -999.0;
  genNuPhi_ = -999.0;
  genNuM_ = -999.0;
  genNuQ_ = -999.0;
  genNustatus_ = -999;
  genNuMother_ = -999;
  genNuGrandMother_ = -999;
  genNuPdgId_ = -999;

  ngenWJet1__ = -999;
  genWJet1_Pt_ = -999.0;
  genWJet1_Eta_ = -999.0;
  genWJet1_Phi_ = -999.0;
  genWJet1_M_ = -999.0;
  genWJet1_Q_ = -999.0;
  genWJet1_status_ = -999;
  genWJet1_Mother_ = -999;
  genWJet1_GrandMother_ = -999;
  genWJet1_PdgId_ = -999;

  ngenWJet2__ = -999;
  genWJet2_Pt_ = -999.0;
  genWJet2_Eta_ = -999.0;
  genWJet2_Phi_ = -999.0;
  genWJet2_M_ = -999.0;
  genWJet2_Q_ = -999.0;
  genWJet2_status_ = -999;
  genWJet2_Mother_ = -999;
  genWJet2_GrandMother_ = -999;
  genWJet2_PdgId_ = -999;

  ngenVBFjet1__ = -999;
  genVBFjet1_Pt_ = -999.0;
  genVBFjet1_Eta_ = -999.0;
  genVBFjet1_Phi_ = -999.0;
  genVBFjet1_M_ = -999.0;
  genVBFjet1_Q_ = -999.0;
  genVBFjet1_status_ = -999;
  genVBFjet1_Mother_ = -999;
  genVBFjet1_GrandMother_ = -999;
  genVBFjet1_PdgId_ = -999;

  ngenVBFjet2__ = -999;
  genVBFjet2_Pt_ = -999.0;
  genVBFjet2_Eta_ = -999.0;
  genVBFjet2_Phi_ = -999.0;
  genVBFjet2_M_ = -999.0;
  genVBFjet2_Q_ = -999.0;
  genVBFjet2_status_ = -999;
  genVBFjet2_Mother_ = -999;
  genVBFjet2_GrandMother_ = -999;
  genVBFjet2_PdgId_ = -999;

nVBFJet_ = -999;
ngenJet_ = -999;
LHEWeightIDs_.clear();
LHEWeights_.clear();
genQuarkStatus_.clear();
genJetPt_.clear();
genJetEta_.clear();
genJetPhi_.clear();
genJetMass_.clear();

genCaloMET_.clear();
genCaloMETPhi_.clear();
genTrueMET_.clear();
genTrueMETPhi_.clear();


vbf_maxpt_j1_pt_ = -999.0;
vbf_maxpt_j1_eta_ = -999.0;
vbf_maxpt_j1_phi_ = -999.0;
vbf_maxpt_j1_e_ = -999.0;
vbf_maxpt_j1_bDiscriminatorCSV_ = -999.0;
vbf_maxpt_j2_pt_ = -999.0;
vbf_maxpt_j2_eta_ = -999.0;
vbf_maxpt_j2_phi_ = -999.0;
vbf_maxpt_j2_e_ = -999.0;
vbf_maxpt_j2_bDiscriminatorCSV_ = -999.0;
vbf_maxpt_jj_pt_ = -999.0;
vbf_maxpt_jj_pt_ = -999.0;
vbf_maxpt_jj_eta_ = -999.0;
vbf_maxpt_jj_phi_ = -999.0;
vbf_maxpt_jj_m_ = -999.0;
vbf_maxpt_deltaR_ = -999.0;
AK4_jet1_pt_ = -999.0;
AK4_jet1_eta_ = -999.0;
AK4_jet1_phi_ = -999.0;
AK4_jet1_e_ = -999.0;
AK4_jet1_bDiscriminatorCSV_ = -999.0;
AK4_jet2_pt_ = -999.0;
AK4_jet2_eta_ = -999.0;
AK4_jet2_phi_ = -999.0;
AK4_jet2_e_ = -999.0;
AK4_jet2_bDiscriminatorCSV_ = -999.0;
AK4_jetjet_pt_ = -999.0;
AK4_jetjet_mass_ = -999.0;
AK4_jetjet_deltaeta_ = -999.0;
AK4_jetjet_deltaphi_ = -999.0;
AK4_jetjet_deltar_ = -999.0;
deltaR_lak4jetjet_ = -999.0;
deltaphi_METak4jetjet_ = -999.0;
deltaphi_Vak4jetjet_ = -999.0;
mass_lvjj_run2_AK4_ = -999.0;
}

