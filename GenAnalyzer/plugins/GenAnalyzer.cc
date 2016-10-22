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
  //const reco::GenMETCollection* genMETColl = &(*genMetCaloHandle);


  int nGenParticle=0;

  
  double l_pt=0., l_eta=0., l_phi=0., l_mass=0., l_mother=0.;
  int l_pdgId=0, l_status=0, l_gmother=0;

  double nu_pt=0., nu_eta=0., nu_phi=0., nu_mass=0.;
  int nu_pdgId=0, nu_status=0, nu_mother=0, nu_gmother=0;

  TLorentzVector ELE, MU, TAU;
  TLorentzVector NU;
  TLorentzVector JET;
  TLorentzVector VBF1, VBF2, VBFTOT;
  TLorentzVector Wjet1, Wjet2, WjetTOT;

  std::vector<TLorentzVector> tightEle, vNU, vJET;
  
  float tempPt=0., tempWm=100.0;

  //genLeptPt_.clear();
  for(vector<reco::GenParticle>::const_iterator genps_it = genps_coll->begin(); genps_it != genps_coll->end(); genps_it++) 
  {
  	nGenParticle++;

	int id = genps_it->pdgId();

	if (Verbose_)
		cout<<"id = "<<id<<endl;

	if((abs(genps_it->pdgId())==11 || abs(genps_it->pdgId())==13 || abs(genps_it->pdgId())==15) && (abs(genps_it->mother()->pdgId()) == 24 || (abs(genps_it->mother()->pdgId()) == 15 && abs(genps_it->mother()->mother()->pdgId()) == 24)))
	{
		if (genps_it->pt() < 30.0) continue;
		if (fabs(genps_it->eta()) > 2.1) continue;
		if (genps_it->pt() < tempPt ) continue;

		ELE.SetPtEtaPhiE(genps_it->pt(), genps_it->eta(), genps_it->phi(), genps_it->energy());
		tightEle.push_back(ELE);


		l_pt      = genps_it->pt();
		l_eta     = genps_it->eta();
		l_phi     = genps_it->phi();
		l_mass    = genps_it->mass();
		l_pdgId   = genps_it->pdgId();
		l_status  = genps_it->status();
		l_mother  = genps_it->mother()->pdgId();
		l_gmother = genps_it->mother()->mother()->pdgId();
		
		tempPt = genps_it->pt();

		nLept++;
		
		if (Verbose_)
			cout<<"n lpe = "<<nLept<<"\tpt = "<<genps_it->pt()<<endl;
	}
	ngenLept_ = nLept;
	//if((abs(genps_it->pdgId())==12 || abs(genps_it->pdgId())==14 || abs(genps_it->pdgId())==16) && (genps_it->status() == 1))
	if((abs(genps_it->pdgId())==12 || abs(genps_it->pdgId())==14 || abs(genps_it->pdgId())==16))
	{
		if (tempWm < fabs(80.385-(l_mass+genps_it->mass()))) continue;

		NU.SetPtEtaPhiE(genps_it->pt(), genps_it->eta(), genps_it->phi(), genps_it->energy());
		vNU.push_back(ELE);

		tempWm = fabs(80.385-(l_mass+genps_it->mass()));


		nu_pt      = genps_it->pt();
		nu_eta     = genps_it->eta();
		nu_phi     = genps_it->phi();
		nu_mass    = genps_it->mass();
		nu_pdgId   = genps_it->pdgId();
		nu_status  = genps_it->status();
		nu_mother  = genps_it->mother()->pdgId();
		nu_gmother = genps_it->mother()->mother()->pdgId();
		//genNuQ_.push_back(genps_it->charge());
	nNu++;
	}
	ngenNu_ = nNu;
   }

   std::vector<const reco::GenJet *> sortedjets;
   sortedjets.reserve(genJetColl->size());
   
   for (const reco::GenJet &g : *genJetColl) { sortedjets.push_back(&g); }
   
   std::sort(sortedjets.begin(),sortedjets.end(),PtGreater());
   
   for (auto const & genPtrjets : sortedjets) 
   {
   	auto const & gjets = *genPtrjets;
	//cout << "Jets PT = " << gjets.pt() << "   Jets Eta = " << gjets.eta() << endl;
	if (nLept != 1) continue;
	if (gjets.pt() < 30.) continue;
	if (fabs(gjets.eta()) > 2.4 ) continue;
	if (deltaR(l_eta,l_phi,gjets.eta(),gjets.phi()) < 0.3) continue;

	JET.SetPtEtaPhiE(gjets.pt(), gjets.eta(), gjets.phi(), gjets.energy());
	vJET.push_back(JET);
	njets++;
   }
   ngenJet_ = njets;


		if (Verbose_)
  if (njets >=4)
	cout << "4 jets*******************************************" << endl;
   
  int nVBFjets = 0;
  float DeltaEta = 0.;
  unsigned int nVBF1=0,nVBF2=0;
  int nGoodAK4Wjets = 0;
  if (njets >=4)
  {
   for(unsigned int i = 0; i<vJET.size()-1; i++)
   {
   	for(unsigned int j = i+1; j<vJET.size(); j++)
	{
		if (Verbose_)
			cout<<DeltaEta<<"\tDelta eta = "<<abs(vJET[i].Eta()-vJET[j].Eta())<<"\t opp hemi = "<<vJET[i].Eta()*vJET[j].Eta()<<"\tinv Mass = "<< (vJET[i]+vJET[j]).M()<<endl;
		if (DeltaEta > abs(vJET[i].Eta()-vJET[j].Eta()) || vJET[i].Eta()*vJET[j].Eta()>0 || (vJET[i]+vJET[j]).M()<300) continue;
		if (abs(vJET[i].Eta()-vJET[j].Eta()) < 3.0) continue;
		VBF1 = vJET[i];
		VBF2 = vJET[j];
		VBFTOT = VBF1 + VBF2 ;
		DeltaEta = abs(vJET[i].Eta()-vJET[j].Eta());
		nVBFjets++;
		nVBF1 = i;
		nVBF2 = j;
	}
   }
	if (Verbose_)
   if (nVBFjets != 0)
   cout<<"jet size = "<<vJET.size()<<"\tnVBF1 = "<<nVBF1<<"\tnVBF2 = "<<nVBF2<<endl;

   if (nVBFjets != 0)
   {
   for(unsigned int i = 0; i<vJET.size(); i++)
   {
	if (i == nVBF1) continue;
	if (i == nVBF2) continue;
	nGoodAK4Wjets++;
	if  (nGoodAK4Wjets == 1)	
		Wjet1 = vJET[i];
	if  (nGoodAK4Wjets == 2)	
		Wjet2 = vJET[i];
   }
   WjetTOT = Wjet1 + Wjet2 ;
   }   
   }
   nVBFJet_ = nVBFjets;

   if (nVBFjets != 0 && nLept==1 && nNu == 1 && nGoodAK4Wjets >= 2)
   {
	if (Verbose_)
	cout<<"Number of leptons = "<<nLept<<endl;

	genLeptPt_ 	= l_pt ;
	genLeptEta_ 	= l_eta ;
	genLeptPhi_ 	= l_phi ;
	genLeptM_ 	= l_mass ;
	genLeptId_ 	= l_pdgId ;
	genLeptStatus_ 	= l_status ;
	genLeptMother_ 	= l_mother ;
	genLeptGrandMother_ = l_gmother;

	genNuPt_ 	= nu_pt;
	genNuEta_ 	= nu_eta;
	genNuPhi_ 	= nu_phi;
	genNuM_ 	= nu_mass;
	//genNuQ_ 	= ;
	genNustatus_ 	= nu_status;
	genNuMother_ 	= nu_mother;
	genNuGrandMother_ = nu_gmother;
	genNuPdgId_ 	= nu_pdgId;	

	vbf_maxpt_j1_pt_ = VBF1.Pt();
	vbf_maxpt_j1_eta_= VBF1.Eta();
	vbf_maxpt_j1_phi_= VBF1.Phi();
	vbf_maxpt_j1_e_  = VBF1.E();

	vbf_maxpt_j2_pt_ = VBF1.Pt();
	vbf_maxpt_j2_eta_= VBF1.Eta();
	vbf_maxpt_j2_phi_= VBF1.Phi();
	vbf_maxpt_j2_e_	 = VBF1.E();

	vbf_maxpt_jj_pt_ = VBFTOT.Pt();
	vbf_maxpt_jj_eta_= VBFTOT.Eta();
	vbf_maxpt_jj_phi_= VBFTOT.Phi();
	vbf_maxpt_jj_m_  = VBFTOT.M();
	vbf_maxpt_deltaR_= deltaR(VBF1.Eta(),VBF1.Phi(),VBF2.Eta(),VBF2.Phi());

	AK4_jet1_pt_	= Wjet1.Pt();
	AK4_jet1_eta_	= Wjet1.Eta();
	AK4_jet1_phi_	= Wjet1.Phi();
	AK4_jet1_e_	= Wjet1.E();

	AK4_jet2_pt_	= Wjet2.Pt();
	AK4_jet2_eta_	= Wjet2.Eta();
	AK4_jet2_phi_	= Wjet2.Phi();
	AK4_jet2_e_	= Wjet2.E();

	AK4_jetjet_pt_	= WjetTOT.Pt();
	AK4_jetjet_mass_= WjetTOT.Eta();
	AK4_jetjet_deltaeta_= fabs(Wjet1.Eta()-Wjet2.Eta());
	AK4_jetjet_deltaphi_= fabs(Wjet1.Phi()-Wjet2.Phi());
	AK4_jetjet_deltar_= deltaR(Wjet1.Eta(),Wjet1.Phi(),Wjet2.Eta(),Wjet2.Phi());
	deltaR_lak4jetjet_= deltaR(WjetTOT.Eta(),WjetTOT.Phi(),l_eta, l_phi);
	deltaphi_METak4jetjet_= fabs(WjetTOT.Phi() - nu_phi);
	deltaphi_Vak4jetjet_= fabs(WjetTOT.Phi() - l_phi);
	mass_lvjj_run2_AK4_= (ELE + NU + WjetTOT + VBFTOT ).M();


   tree->Fill();
   }
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
	//AddBranch(&pdgID_,	"pdgID");

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

nVBFJet_ = -999;
ngenJet_ = -999;
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
