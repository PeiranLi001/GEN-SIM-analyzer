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

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "GenAnalyzer.h"

#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/METReco/interface/GenMETCollection.h"
#include "DataFormats/METReco/interface/GenMET.h"

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

  int nNu = 0;
  int nLept = 0;
  int njets = 0;

  edm::Handle<reco::GenParticleCollection> genpsHandle;
  iEvent.getByToken(genParticlesToken, genpsHandle);
  
  edm::Handle<reco::GenJetCollection> genAK4JetsHandle;
  iEvent.getByToken(genAK4jetToken, genAK4JetsHandle);

  edm::Handle<reco::GenMETCollection> genMetCaloHandle;
  iEvent.getByToken(genMetCaloToken , genMetCaloHandle);

  //  edm::Handle<reco::GenMETCollection> genMetTrueHandle;
  //  iEvent.getByToken(genMetTrueToken , genMetTrueHandle);
	
  const vector<reco::GenParticle>* genps_coll = genpsHandle.product();
  //const vector<reco::GenJet>* genJetColl= genAK4JetsHandle.product();
  //const vector<reco::GenMET>* genMETCaloColl= genMetCaloHandle.product();
  const reco::GenJetCollection* genJetColl= &(*genAK4JetsHandle);
  const reco::GenMETCollection* genMETColl = &(*genMetCaloHandle);


  int nGenParticle=0;
  
  for(vector<reco::GenParticle>::const_iterator genps_it = genps_coll->begin(); genps_it != genps_coll->end(); genps_it++) 
  {
  	nGenParticle++;

	int id = genps_it->pdgId();
	int st = genps_it->status();
	const Candidate * mother = genps_it->mother();

	if(fabs(id)==11 && st==23 && mother->pdgId()==24){
	//cout << "mass =   " << mother->mass() << endl;  
	}

	if (Verbose_)
		cout<<"id = "<<id<<endl;
	//pdgID_.push_back(id);

	if((abs(genps_it->pdgId())==11 || abs(genps_it->pdgId())==13 || abs(genps_it->pdgId())==15) && (abs(genps_it->mother()->pdgId()) == 24 || (abs(genps_it->mother()->pdgId()) == 15 && abs(genps_it->mother()->mother()->pdgId()) == 24))){

	genLeptPt_.push_back(genps_it->pt());
	genLeptEta_.push_back(genps_it->eta());
	genLeptPhi_.push_back(genps_it->phi());
	genLeptM_.push_back(genps_it->mass());
	genLeptId_.push_back(genps_it->pdgId());
	genLeptStatus_.push_back(genps_it->status());
	genLeptMother_.push_back(genps_it->mother()->pdgId());
	genLeptGrandMother_.push_back(genps_it->mother()->mother()->pdgId());

	nLept++;
	}
	ngenLept_ = nLept;

	if((abs(genps_it->pdgId())==12 || abs(genps_it->pdgId())==14 || abs(genps_it->pdgId())==16) && (genps_it->status() == 1))
	{
	
	       genNuPt_.push_back(genps_it->pt());
	       genNuEta_.push_back(genps_it->eta());
	       genNuPhi_.push_back(genps_it->phi());
	       genNuQ_.push_back(genps_it->charge());
	       genNustatus_.push_back(genps_it->status());
	       genNuMother_.push_back(genps_it->mother()->pdgId());
	       genNuGrandMother_.push_back(genps_it->mother()->mother()->pdgId());
	       genNuPdgId_.push_back(genps_it->pdgId());
	nNu++;
	}
	ngenNu_ = nNu;
   }

	
//  for(vector<reco::GenJet>::const_iterator genps_it = genJetColl->begin(); genps_it != genJetColl->end(); genps_it++) 
//  {
// 	nNu++;
//}
//  for(vector<reco::GenMET>::const_iterator genps_it = genMETCaloColl->begin(); genps_it != genMETCaloColl->end(); genps_it++) 
//  {
// 	nNu++;
//}

   std::vector<const reco::GenJet *> sortedjets;
   sortedjets.reserve(genJetColl->size());
   
   for (const reco::GenJet &g : *genJetColl) { sortedjets.push_back(&g); }
   
   std::sort(sortedjets.begin(),sortedjets.end(),PtGreater());
   
   for (auto const & genPtrjets : sortedjets) 
   {
   	auto const & gjets = *genPtrjets;
	//cout << "Jets PT = " << gjets.pt() << "   Jets Eta = " << gjets.eta() << endl;
	if(gjets.pt() > 30. && fabs(gjets.eta()) < 5.0) {
		genJetPt_.push_back(gjets.pt());
		genJetEta_.push_back(gjets.eta());
		genJetPhi_.push_back(gjets.phi());
		genJetMass_.push_back(gjets.mass());
		njets++;
	}
	//cout << gjets.mass() << endl;
   }
   ngenJet_ = njets;
   
   
   reco::GenMETCollection::const_iterator i;
   
   for(i=genMETColl->begin(); i!=genMETColl->end(); i++)
   {
   	genCaloMET_.push_back(i->pt());
	genCaloMETPhi_.push_back(i->phi());
   }

   //cout<<"nGenParticle = "<<nGenParticle<<endl;
   tree->Fill();
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

void GenAnalyzer::SetBranches(){
	AddBranch(&pdgID_,	"pdgID");

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
  AddBranch(&genNuQ_,"genNuQ");
  AddBranch(&genNustatus_,"genNustatus");
  AddBranch(&genNuMother_,"genNuMother");
  AddBranch(&genNuGrandMother_,"genNuGrandMother");

  AddBranch(&genJetPt_, "genJetPt");
  AddBranch(&genJetEta_, "genJetEta");
  AddBranch(&genJetPhi_, "genJetPhi");
  AddBranch(&genJetMass_, "genJetMass");
  AddBranch(&ngenJet_, "ngenJet");

  AddBranch(&genCaloMET_, "genCaloMET");
  AddBranch(&genCaloMETPhi_, "genCaloMETPhi");
  AddBranch(&genTrueMET_, "genTrueMET");
  AddBranch(&genTrueMETPhi_, "genTrueMETPhi");


}

void GenAnalyzer::Clear(){
	pdgID_.clear();

	ngenLept_ = -999;
	genLeptPt_.clear();
	genLeptEta_.clear();
	genLeptPhi_.clear();
	genLeptStatus_.clear();
	genLeptMother_.clear();
	genLeptGrandMother_.clear();
	genLeptId_.clear();
	genLeptM_.clear();

  ngenNu_ = -999;
  genNuPt_.clear();
  genNuEta_.clear();
  genNuPhi_.clear();
  genNuQ_.clear();
  genNustatus_.clear();
  genNuMother_.clear();
  genNuGrandMother_.clear();
  genNuPdgId_.clear();

ngenJet_ = -999;
genJetPt_.clear();
genJetEta_.clear();
genJetPhi_.clear();
genJetMass_.clear();

genCaloMET_.clear();
genCaloMETPhi_.clear();
genTrueMET_.clear();
genTrueMETPhi_.clear();
}

// ------------ method called once each job just before starting event loop  ------------
void 
GenAnalyzer::beginJob()
{
    std::cout<<"Inside beginJob()"<<std::endl;
    outputFile_ = new TFile("Test.root","RECREATE"); 
    outputFile_->SetCompressionLevel(2);
    tree = new TTree("genParticles","GenParticles Basic Info"); 

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
