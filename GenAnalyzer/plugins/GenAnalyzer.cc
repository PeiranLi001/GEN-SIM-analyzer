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
  
  const vector<reco::GenParticle>* genps_coll = genpsHandle.product();

  edm::Handle<LHEEventProduct> LHEEventHandle;
  iEvent.getByToken(LHEEventToken, LHEEventHandle);
  const LHEEventProduct* LHE = 0;

  if(LHEEventHandle.isValid()){
  	LHE = LHEEventHandle.product();

	for(const auto& weight : LHE->weights()) {
		LHEWeightIDs_.push_back(weight.id);
		LHEWeights_.push_back(weight.wgt);
	}
  }


  int nGenParticle=0;

  
  double l_pt=0., l_eta=0., l_phi=0., l_mass=0., l_mother=0.;
  int l_pdgId=0, l_status=0, l_gmother=0;

  double nu_pt=0., nu_eta=0., nu_phi=0., nu_mass=0.;
  int nu_pdgId=0, nu_status=0, nu_mother=0, nu_gmother=0;

  TLorentzVector ELE, MU, TAU;
  TLorentzVector NU;
  TLorentzVector JET, VBFquarks, Wquarks;
  TLorentzVector VBF1, VBF2, VBFTOT;
  TLorentzVector Wjet1, Wjet2, WjetTOT;

  std::vector<TLorentzVector> tightLep, vNU, vJET, wJET;
  
  for(vector<reco::GenParticle>::const_iterator genps_it = genps_coll->begin(); genps_it != genps_coll->end(); genps_it++) 
  {
  	nGenParticle++;
	if((abs(genps_it->pdgId())==11 || abs(genps_it->pdgId())==13 || abs(genps_it->pdgId())==15) && (abs(genps_it->mother()->pdgId()) == 24) )
	{
		if (Verbose_)
		cout<<"Status of leptons = "<<genps_it->status()<<endl;
		if (Verbose_)
		cout<<"Status of leptons monther = "<<genps_it->mother()->status()<<endl;
		if (genps_it->pt() < 45.0) continue;
		if (fabs(genps_it->eta()) > 2.1) continue;

		l_pdgId = genps_it->pdgId();
		l_status = genps_it->status();
		l_mother = genps_it->mother()->pdgId();
		l_gmother = genps_it->mother()->mother()->pdgId();

		ELE.SetPtEtaPhiE(genps_it->pt(), genps_it->eta(), genps_it->phi(), genps_it->energy());
		tightLep.push_back(ELE);
	}
	if((abs(genps_it->pdgId())==12 || abs(genps_it->pdgId())==14 || abs(genps_it->pdgId())==16) && (abs(genps_it->mother()->pdgId()) == 24) )
	{
		if (Verbose_)
		cout<<"Status of Neutrino = "<<genps_it->status()<<endl;

		nu_pdgId   = genps_it->pdgId();
		nu_status  = genps_it->status();
		nu_mother  = genps_it->mother()->pdgId();
		nu_gmother = genps_it->mother()->mother()->pdgId();

		NU.SetPtEtaPhiE(genps_it->pt(), genps_it->eta(), genps_it->phi(), genps_it->energy());
		vNU.push_back(NU);
	}
	if((abs(genps_it->pdgId())==1 || abs(genps_it->pdgId())==2 || abs(genps_it->pdgId())==3 || abs(genps_it->pdgId())==4 || abs(genps_it->pdgId())==5 || abs(genps_it->pdgId())==6) && (abs(genps_it->mother()->pdgId()) == 24) && (genps_it->status()==23) )
	{
		Wquarks.SetPtEtaPhiE(genps_it->pt(), genps_it->eta(), genps_it->phi(), genps_it->energy());
		wJET.push_back(Wquarks);
	}
	if((abs(genps_it->pdgId())==1 || abs(genps_it->pdgId())==2 || abs(genps_it->pdgId())==3 || abs(genps_it->pdgId())==4 || abs(genps_it->pdgId())==5 || abs(genps_it->pdgId())==6) && (abs(genps_it->mother()->pdgId()) != 24) && (genps_it->status()==23) )
	{
		VBFquarks.SetPtEtaPhiE(genps_it->pt(), genps_it->eta(), genps_it->phi(), genps_it->energy());
		vJET.push_back(Wquarks);
	}
  }

  if (tightLep.size()==1 && vNU.size()==1 && wJET.size()==2 && vJET.size()==2)
  {
		if (Verbose_)
  	cout<<"Event selected"<<nEVENT<<endl;
	nEVENT++;

	genLeptPt_ 	= ELE.Pt() ;
	genLeptEta_ 	= ELE.Eta();
	genLeptPhi_ 	= ELE.Phi();
	genLeptM_ 	= ELE.M() ;
	genLeptId_ 	= l_pdgId ;
	genLeptStatus_ 	= l_status ;
	genLeptMother_ 	= l_mother ;
	genLeptGrandMother_ = l_gmother;

	genNuPt_ 	= NU.Pt(); 
	genNuEta_ 	= NU.Eta();
	genNuPhi_ 	= NU.Phi();
	genNuM_ 	= NU.M();
	//genNuQ_ 	= ;
	genNustatus_ 	= nu_status;
	genNuMother_ 	= nu_mother;
	genNuGrandMother_ = nu_gmother;
	genNuPdgId_ 	= nu_pdgId;	

	vbf_maxpt_j1_pt_ = vJET[0].Pt();
	vbf_maxpt_j1_eta_= vJET[0].Eta();
	vbf_maxpt_j1_phi_= vJET[0].Phi();
	vbf_maxpt_j1_e_  = vJET[0].E();

	vbf_maxpt_j2_pt_ = vJET[1].Pt();
	vbf_maxpt_j2_eta_= vJET[1].Eta();
	vbf_maxpt_j2_phi_= vJET[1].Phi();
	vbf_maxpt_j2_e_	 = vJET[1].E();

	vbf_maxpt_jj_pt_ = (vJET[0]+vJET[1]).Pt();
	vbf_maxpt_jj_eta_= (vJET[0]+vJET[1]).Eta();
	vbf_maxpt_jj_phi_= (vJET[0]+vJET[1]).Phi();
	vbf_maxpt_jj_m_  = (vJET[0]+vJET[1]).M();
	vbf_maxpt_deltaR_= deltaR(vJET[0].Eta(),vJET[0].Phi(),vJET[1].Eta(),vJET[1].Phi());

	AK4_jet1_pt_	= wJET[0].Pt();
	AK4_jet1_eta_	= wJET[0].Eta();
	AK4_jet1_phi_	= wJET[0].Phi();
	AK4_jet1_e_	= wJET[0].E();

	AK4_jet2_pt_	= wJET[1].Pt();
	AK4_jet2_eta_	= wJET[1].Eta();
	AK4_jet2_phi_	= wJET[1].Phi();
	AK4_jet2_e_	= wJET[1].E();

	AK4_jetjet_pt_	= (wJET[0]+wJET[1]).Pt();
	AK4_jetjet_mass_= (wJET[0]+wJET[1]).Eta();
	AK4_jetjet_deltaeta_= fabs(wJET[0].Eta()-wJET[1].Eta());
	AK4_jetjet_deltaphi_= fabs(wJET[0].Phi()-wJET[1].Phi());
	AK4_jetjet_deltar_= deltaR(wJET[0].Eta(),wJET[0].Phi(),wJET[1].Eta(),wJET[1].Phi());
	deltaR_lak4jetjet_= deltaR((wJET[0]+wJET[1]).Eta(),(wJET[0]+wJET[1]).Phi(),l_eta, l_phi);
	deltaphi_METak4jetjet_= fabs((wJET[0]+wJET[1]).Phi() - nu_phi);
	deltaphi_Vak4jetjet_= fabs((wJET[0]+wJET[1]).Phi() - l_phi);
	mass_lvjj_run2_AK4_= (ELE + NU + (wJET[0]+wJET[1]) + VBFTOT ).M();


  }

tree->Fill();
Clear();   
}

// ------------ method called once each job just before starting event loop  ------------
void 
GenAnalyzer::beginJob()
{
    std::cout<<"Inside beginJob()"<<std::endl;
    outputFile_ = new TFile("Test.root","RECREATE"); 
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
