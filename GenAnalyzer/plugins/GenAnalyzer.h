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
      edm::EDGetTokenT<reco::GenJetCollection> genAK4jetToken;
      edm::EDGetTokenT<reco::GenMETCollection> genMetCaloToken;
      
      bool	Verbose_;
      
      edm::Service<TFileService> fs;
      TFile * outputFile_;
      TTree* tree;


    std::vector<int> pdgID_;

  int ngenLept_;
  std::vector<double> genLeptPt_;
  std::vector<double> genLeptEta_;
  std::vector<double> genLeptPhi_;
  std::vector<double> genLeptM_;
  std::vector<int> genLeptId_;
  std::vector<int> genLeptStatus_;
  std::vector<double> genLeptMother_;
  std::vector<int> genLeptGrandMother_;
    
 int ngenNu_;
 std::vector<double> genNuPt_;
 std::vector<double> genNuEta_;
 std::vector<double> genNuPhi_;
 std::vector<double> genNuQ_;
 std::vector<int> genNustatus_;
 std::vector<int> genNuMother_;
 std::vector<int> genNuGrandMother_;
 std::vector<int> genNuPdgId_;

  int ngenJet_;
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
